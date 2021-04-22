#include "short_read_mapper.h"

#include <unistd.h>

#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <random>
#include <string>

#include "utils.h"

#define READ_NOT_MAPPED 0b00
#define READ_MAPPED 0b01
#define READ_SATELLITE 0b10

using namespace std;

void ShortReadMapper::genSeedMask() {
    uint64_t seed_mask = 0;
    for (int i = 0; i < _seed_len; ++i) {
        seed_mask = (seed_mask << 2) + 3;
    }
    _seed_mask = seed_mask;
}

void ShortReadMapper::updateSeed(char& base, uint64_t& seed) {
    if (base == 'A' || base == 'a')
        seed = (seed << 2) + 0;
    else if (base == 'C' || base == 'c')
        seed = (seed << 2) + 1;
    else if (base == 'G' || base == 'g')
        seed = (seed << 2) + 2;
    else if (base == 'T' || base == 't')
        seed = (seed << 2) + 3;
    seed &= _seed_mask;
}

void ShortReadMapper::updateRefSeq(char& base, long base_cnt) {
    if (base == 'A' || base == 'a')
        _ref_seq[base_cnt] = 'A';
    else if (base == 'C' || base == 'c')
        _ref_seq[base_cnt] = 'C';
    else if (base == 'G' || base == 'g')
        _ref_seq[base_cnt] = 'G';
    else if (base == 'T' || base == 't')
        _ref_seq[base_cnt] = 'T';
}

string ShortReadMapper::getRefSeqFromLoc(long loc, int len) {
    // TODO: char* does not have terminating char now
    char* seq = new char[len];
    for (int i = 0; i < len; i++) {
        seq[i] = _ref_seq[loc + i];
    }
    string rv = string(seq);
    delete seq;
    return rv;
}

bool ShortReadMapper::isSatellite(int layer_id, int hit_cnt[]) {
    for (int i = 0; i < _bf_amount[layer_id]; i++) {
        if (hit_cnt[i] >= _hit_threshold) {
            _layer_hit_cnt[layer_id] += 1;
        }
    }
    return _layer_hit_cnt[layer_id] > _satellite_threshold;
}

void ShortReadMapper::initQuery() {
    _bml_sel->init();
    for (int i = 0; i < _layer_num; i++) {
        _layer_hit_cnt[i] = 0;
    }
}

int ShortReadMapper::queryLayer(string& read, int layer_id, long hier_offset,
                                long base_offset) {
    /*
    In each layer, query every seeds of the read and
    record the hit count. If hit count > threshold, recursively
    qurey the next layer.

    Return value:
    0th bit: read mapped
    1st bit: satellite
    */

    // Whether it is the last layer
    bool last_layer = layer_id == _layer_num - 1;

    // Initialize return value
    int rv = READ_NOT_MAPPED;

    // Build hit count array
    long bf_amount = _bf_amount[layer_id];
    int hit_cnt[bf_amount];
    for (int i = 0; i < bf_amount; i++) {
        hit_cnt[i] = 0;
    }

    // Generate seed and query the layer
    uint64_t seed = 0;
    for (int i = 0; i < read.length(); i++) {
        updateSeed(read[i], seed);
        seed &= _seed_mask;

        // When the seed is not long enough, keep reading base
        if (i < _seed_len - 1) continue;

        // Query the layer
        // If it's the last layer, OR the nearby Bloom filter
        if (last_layer)
            _layers[layer_id]->query(seed, hit_cnt, hier_offset, true);
        else
            _layers[layer_id]->query(seed, hit_cnt, hier_offset, false);
    }

    long hit_threshold;
    if (layer_id == 0)
        hit_threshold = min(meanPlusStdev(hit_cnt, 14, 1), _hit_threshold);
    else
        hit_threshold = _hit_threshold;

    if (layer_id == 1) {
        int max_hit_cnt = findMax(hit_cnt, bf_amount);
        if (countHitBF(hit_cnt, bf_amount, hit_threshold) == 0) {
            hit_threshold = max_hit_cnt;
        }
    }

    // printHitCnt(layer_id, hit_cnt, bf_amount);

    // If too many hits, return satelllite code
    if (layer_id != 0)
        if (isSatellite(layer_id, hit_cnt)) return READ_SATELLITE;

    // For each Bloom filter
    for (int i = 0; i < bf_amount; i++) {
        if (hit_cnt[i] >= hit_threshold) {
            // If it is the last layer,
            // use the BML selector to calculate the score.
            if (last_layer) {
                rv = READ_MAPPED;
                // Calculate the CML location
                long cml_loc = base_offset + i * _seed_range[layer_id];
                int seq_len = _seed_range[layer_id] * 2;
                string ref_seq = getRefSeqFromLoc(cml_loc, seq_len);

                // Send one CML to the BML engine
                _seeding_sw->pause();
                _seed_extraction_sw->start();
                _bml_sel->update(ref_seq, read, cml_loc);
                _seeding_sw->start();
                _seed_extraction_sw->pause();
            }
            // If not the last layer, query the next layer
            else {
                long hier_offset_next = hier_offset + i * _bf_size[layer_id];
                long base_offset_next = base_offset + i * _seed_range[layer_id];

                rv |= queryLayer(read, layer_id + 1, hier_offset_next,
                                 base_offset_next);
                // If we found the read is satellite at the child layer,
                // return immediately.
                if (rv & READ_SATELLITE) return rv;
            }
        }
    }

    return rv;
}

void ShortReadMapper::updateScoreboard(int& rv, long& golden_loc,
                                       long& mapped_loc, bool verbose) {
    /*
    Return value:
    0th bit: read mapped
    1st bit: satellite
    */

    if (rv & READ_SATELLITE) {
        // Satellite
        _satellite += 1;
        if (verbose) cout << "Satellite" << endl;
    }
    else if (rv & READ_MAPPED) {
        // Mapped
        if (abs(golden_loc - mapped_loc) <= _ans_margin) {
            _correctly_mapped += 1;
            if (verbose) cout << "Correctly mapped" << endl;
        }
        else {
            _wrongly_mapped += 1;
            if (verbose) cout << "Wrongly mapped" << endl;
        }
    }
    else {
        // Not mapped
        _not_mapped += 1;
        if (verbose) cout << "Not mapped" << endl;
    }
}

ShortReadMapper::ShortReadMapper(string& ref_path, string& read_path,
                                 long read_len, long seed_len,
                                 long query_shift_amt, long hit_threshold,
                                 long ans_margin, long satellite_threshold) {
    // File paths
    _ref_path = ref_path;
    _read_path = read_path;

    // Configuration
    _read_len = read_len;
    _seed_len = seed_len;
    genSeedMask();
    _query_skip_amt = query_shift_amt;
    _hit_threshold = hit_threshold;
    _ans_margin = ans_margin;
    _satellite_threshold = satellite_threshold;
    _layer_num = 3;

    // Bloom filters configuration
    // 16 MB for each Bloom filter in layer 1
    // 64 kB for each Bloom filter in layer 2
    // 256 Bytes for each Bloom filter in layer 3
    _bf_size = new long[_layer_num];
    _bf_size[0] = 16 * 1024 * 1024 * 8;
    _bf_size[1] = 64 * 1024 * 8;
    _bf_size[2] = 256 * 8;

    // 256 layer 1 Bloom filters
    // 256 layer 2 Bloom filters associated with one layer 1 BF
    // 256 layer 3 Bloom filters associated with one layer 2 BF
    _bf_amount = new long[_layer_num];
    _bf_amount[0] = 256;
    _bf_amount[1] = 256;
    _bf_amount[2] = 256;

    _bf_total = new long[_layer_num];
    _bf_total[0] = 256;
    _bf_total[1] = 256 * 256;
    _bf_total[2] = 256 * 256 * 256;

    // Save 256 * 256 * 256 seeds in a Bloom filter in layer 1
    // Save 256 * 256 seeds in a Bloom filter in layer 2
    // Save 256 seeds in a Bloom filter in layer 3
    _seed_range = new long[_layer_num];
    _seed_range[0] = 256 * 256 * 256;
    _seed_range[1] = 256 * 256;
    _seed_range[2] = 256;

    // Mapping configuration
    _test_num = 10000;
    _ref_size = 2948627755;

    // Generate hash factor
    int rand_seed = 666;
    srand(1);
    random_device rd;
    default_random_engine generator(rand_seed);
    uniform_int_distribution<uint64_t> distribution(0, 0xFFFFFFFFFFFFFFFF);
    uint64_t hash_factors[_layer_num];
    for (int i = 0; i < _layer_num; i++) {
        hash_factors[i] = distribution(generator);
    }

    // Instantiate N layers
    _layers = new Layer*[_layer_num];
    for (int i = 0; i < _layer_num; i++) {
        _layers[i] = new Layer(_bf_size[i], _bf_amount[i], _bf_total[i],
                               _seed_range[i], hash_factors[i]);
    }

    // Total hit count in each layer
    // If hit_cnt > _satellite_threshold, read is satellite
    _layer_hit_cnt = new int[_layer_num];

    // Initialize _ref_seq
    _ref_seq = new char[_ref_size];

    // Instantiate BML selector
    _bml_sel = new BMLSelector();

    // Scoreboard
    _correctly_mapped = 0;
    _wrongly_mapped = 0;
    _satellite = 0;
    _not_mapped = 0;

    // Stopwatch
    _training_sw = new Stopwatch();
    _seeding_sw = new Stopwatch();
    _seed_extraction_sw = new Stopwatch();
    _training_sw->reset();
    _seeding_sw->reset();
    _seed_extraction_sw->reset();
}

ShortReadMapper::~ShortReadMapper() {
    // Layer configuration
    delete[] _bf_size;
    delete[] _bf_amount;
    delete[] _bf_total;
    delete[] _seed_range;
    for (int i = 0; i < _layer_num; i++) {
        delete _layers[i];
    }
    delete[] _layers;
    delete[] _layer_hit_cnt;

    delete _ref_seq;
    delete _bml_sel;

    // Stopwatch
    delete _training_sw;
    delete _seeding_sw;
    delete _seed_extraction_sw;
}

void ShortReadMapper::trainBF(bool ignoreSatellite) {
    cout << "[trainBF] Start training the Bloom filter" << endl;
    if (ignoreSatellite) cout << "[trainBF] Ignore satellite DNA" << endl;

    // Start stopwatch
    _training_sw->start();

    // Seed
    uint64_t seed = 0;
    long base_cnt = 0;

    // Open ref file
    ifstream ref_seq_fs(_ref_path);
    if (!ref_seq_fs.is_open()) {
        cerr << "[trainBF] Cannot open the reference sequence file." << endl;
        exit(1);
    }

    string line;

    // If ignoreSatellite, start building the _seed_cnt map
    if (ignoreSatellite) {
        // First pass seed and base_cnt
        uint64_t seed_fp = 0;
        long base_cnt_fp = 0;

        while (ref_seq_fs >> line) {
            // If the line starts with '>', ignore it
            if (line[0] == '>') continue;

            // For each character, generate a seed
            for (int i = 0; i < line.size(); i++) {
                updateSeed(line[i], seed_fp);

                // If the seed variable contains more than seed_len seeds,
                // start updating seed count.
                if (base_cnt_fp >= _seed_len - 1) {
                    findAndIncrement(_seed_cnt, seed_fp);
                }

                base_cnt_fp += 1;
                if (base_cnt_fp % 10000000 == 0)
                    cout << "[trainBF] Counted for " << base_cnt_fp << " seeds"
                         << endl;
                if (base_cnt_fp == _ref_size) break;
            }
        }

        // Set file stream position to the beginning
        ref_seq_fs.clear();
        ref_seq_fs.seekg(0, ios::beg);
    }

    // Parse the ref file line by line
    while (ref_seq_fs >> line) {
        // If the line starts with '>', ignore it
        if (line[0] == '>') continue;

        // For each character, generate a seed
        for (int i = 0; i < line.size(); i++) {
            updateSeed(line[i], seed);
            updateRefSeq(line[i], base_cnt);

            // If the seed variable contains more than seed_len seeds,
            // start updating the Bloom filter.
            if (base_cnt >= _seed_len - 1) {
                if (ignoreSatellite) {
                    int cnt = getCount(_seed_cnt, seed);
                    if (cnt <= _satellite_threshold) {
                        _layers[0]->update(seed, base_cnt);
                        _layers[1]->update(seed, base_cnt);
                        _layers[2]->update(seed, base_cnt);
                    }
                }
                else {
                    // Normal mode
                    _layers[0]->update(seed, base_cnt);
                    _layers[1]->update(seed, base_cnt);
                    _layers[2]->update(seed, base_cnt);
                }
            }

            base_cnt += 1;
            if (base_cnt % 10000000 == 0)
                cout << "[trainBF] Processed " << base_cnt << " seeds" << endl;
            if (base_cnt == _ref_size) break;
        }
    }

    // printSeedCnt(_seed_cnt);

    // Pause stopwatch
    _training_sw->pause();
}

void ShortReadMapper::writeBF() {
    // Currently not used
    _layers[0]->write_bf_bin("layer1_bf.dat");
    _layers[1]->write_bf_bin("layer2_bf.dat");
    _layers[2]->write_bf_bin("layer3_bf.dat");
}

void ShortReadMapper::readBF() {
    // Currently not used
    _layers[0]->read_bf_bin("layer1_bf.dat");
    _layers[1]->read_bf_bin("layer2_bf.dat");
    _layers[2]->read_bf_bin("layer3_bf.dat");
}

void ShortReadMapper::mapRead() {
    cout << "[mapRead] Start mapping the reads" << endl;

    // Start stopwatch
    _seeding_sw->start();

    // Open read file
    ifstream read_seq_fs(_read_path);
    if (!read_seq_fs.is_open()) {
        cerr << "[mapRead] Cannot open the read sequence file." << endl;
        exit(1);
    }

    // Read header
    string line;
    while (getline(read_seq_fs, line)) {
        if (line == "##Header End") break;
    }

    // Read reads
    int read_cnt = 0;
    string token;
    string golden_loc_s;
    string fwd_rev;
    string read;

    while (read_cnt < _test_num) {
        /* Read format:
        >chr1  chr1-1536540  116446253  -
        <Original sequence>
        <Simulater generated sequence>
        */

        read_seq_fs >> token;         // Ignore reference sequence name
        read_seq_fs >> token;         // Ignore read name
        read_seq_fs >> golden_loc_s;  // Get golden location
        read_seq_fs >> fwd_rev;       // Get forward/reverse
        read_seq_fs >> token;         // Ignore original sequence
        read_seq_fs >> read;          // Get simulator generated read sequence

        long golden_loc = stol(golden_loc_s);

        // Only map the forward sequence, ignore the reverse sequence
        if (fwd_rev == "-") continue;

        long bf_idx1 = golden_loc / (256 * 256 * 256);
        long bf_idx2 = (golden_loc % (256 * 256 * 256)) / (256 * 256);
        long bf_idx3 = (golden_loc % (256 * 256)) / (256);
        // cout << "\n" << read << endl;
        // cout << "golden loc: " << golden_loc;
        // cout << " (" << bf_idx1 << ", " << bf_idx2 << ", " << bf_idx3 << ")"
        //      << endl;

        // Query the read in each layer recursively
        initQuery();
        int layer_id = 0;
        long hier_offset = 0;
        long base_offset = 0;
        int rv = queryLayer(read, layer_id, hier_offset, base_offset);

        // Get mapped location from the BML selector
        long mapped_loc = _bml_sel->getMapLoc();
        bool verbose = false;
        updateScoreboard(rv, golden_loc, mapped_loc, verbose);

        read_cnt += 1;
        if (read_cnt % 1000 == 0)
            cout << "[mapRead] Processed " << read_cnt << " reads" << endl;
    }

    // pause stopwatch
    _seeding_sw->pause();
}

void ShortReadMapper::displayResult() {
    int sum = _correctly_mapped + _wrongly_mapped + _satellite + _not_mapped;

    cout << "\n---- Mapping Result ----" << endl;
    cout << "Correctly mapped: " << setw(5) << _correctly_mapped << endl;
    cout << "Wrongly mapped:   " << setw(5) << _wrongly_mapped << endl;
    cout << "Satellite:        " << setw(5) << _satellite << endl;
    cout << "Not mapped:       " << setw(5) << _not_mapped << endl;
    cout << "Total:            " << setw(5) << sum << endl;

    cout << "\n---- Duration (sec) ----" << endl;
    cout << fixed << setprecision(2);
    cout << "Training:         " << setw(5) << _training_sw->getSec() << endl;
    cout << "Seeding:          " << setw(5) << _seeding_sw->getSec() << endl;
    cout << "Seed extraction:  " << setw(5) << _seed_extraction_sw->getSec()
         << endl;
}
