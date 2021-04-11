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
    for (int i = 0; i < _seed_len; ++i) seed_mask = (seed_mask << 2) + 3;
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

bool ShortReadMapper::isSatellite(int layer, int hit_cnt[]) {
    if (layer == 1) {
        for (int i = 0; i < _bf_amount1; i++) {
            if (hit_cnt[i] >= _hit_threshold) {
                _layer1_hit_cnt += 1;
            }
        }
        return _layer1_hit_cnt > _satellite_threshold;
    }
    else if (layer == 2) {
        for (int i = 0; i < _bf_amount2; i++) {
            if (hit_cnt[i] >= _hit_threshold) {
                _layer2_hit_cnt += 1;
            }
        }
        return _layer2_hit_cnt > _satellite_threshold;
    }
    else if (layer == 3) {
        for (int i = 0; i < _bf_amount3; i++) {
            if (hit_cnt[i] >= _hit_threshold) {
                _layer3_hit_cnt += 1;
            }
        }
        return _layer3_hit_cnt > _satellite_threshold;
    }
    return false;
}

void ShortReadMapper::initQuery() {
    _bml_sel->init();
    _layer1_hit_cnt = 0;
    _layer2_hit_cnt = 0;
    _layer3_hit_cnt = 0;
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

    int rv = READ_NOT_MAPPED;

    if (layer_id == 1) {
        // Build layer 1 hit count array
        int hit_cnt[_bf_amount1];
        for (int i = 0; i < _bf_amount1; i++) hit_cnt[i] = 0;

        // Generate seed
        uint64_t seed = 0;
        for (int i = 0; i < read.length(); i++) {
            updateSeed(read[i], seed);
            seed &= _seed_mask;

            // When the seed is not long enough, keep reading base
            if (i < _seed_len - 1) continue;

            // Query the layer 1
            _layer[0]->query(seed, hit_cnt, hier_offset, false);
        }

        // If too many hits, return satelllite code
        if (isSatellite(1, hit_cnt)) return READ_SATELLITE;

        // For each Bloom filter that has more than _hit_threshold hits,
        // query the next layer.
        for (int i = 0; i < _bf_amount1; i++) {
            if (hit_cnt[i] >= _hit_threshold) {
                long hier_offset_next = i * _bf_amount2 * _bf_size2;
                long base_offset_next = i * _seed_range1;

                rv |= queryLayer(read, 2, hier_offset_next, base_offset_next);
                // If we found the read is satellite at the child layer,
                // return immediately.
                if (rv & READ_SATELLITE) return rv;
            }
        }
    }
    else if (layer_id == 2) {
        // Build layer 2 hit count array
        int hit_cnt[_bf_amount2];
        for (int i = 0; i < _bf_amount2; i++) hit_cnt[i] = 0;

        // Generate seed
        uint64_t seed = 0;
        for (int i = 0; i < read.length(); i++) {
            updateSeed(read[i], seed);

            // When the seed is not long enough, keep reading base
            if (i < _seed_len - 1) continue;

            // Query the layer 2
            _layer[1]->query(seed, hit_cnt, hier_offset, false);
        }

        // If too many hits, return satelllite code
        if (isSatellite(2, hit_cnt)) return READ_SATELLITE;

        // For each Bloom filter that has more than _hit_threshold hits,
        // query the next layer.
        for (int i = 0; i < _bf_amount1; i++) {
            if (hit_cnt[i] >= _hit_threshold) {
                long hier_offset_next =
                    hier_offset + i * _bf_amount3 * _bf_size3;
                long base_offset_next = base_offset + i * _seed_range2;

                rv |= queryLayer(read, 3, hier_offset_next, base_offset_next);
                // If we found the read is satellite at the child layer,
                // return immediately.
                if (rv & READ_SATELLITE) return rv;
            }
        }
    }
    else if (layer_id == 3) {
        // Build layer 3 hit count array
        int hit_cnt_or[_bf_amount3];
        for (int i = 0; i < _bf_amount3; i++) hit_cnt_or[i] = 0;

        // Generate seed
        uint64_t seed = 0;
        for (int i = 0; i < read.length(); i++) {
            updateSeed(read[i], seed);

            // When the seed is not long enough, keep reading base
            if (i < _seed_len - 1) continue;

            // Query the layer 3
            _layer[2]->query(seed, hit_cnt_or, hier_offset, true);
        }

        // If too many hits, return satelllite code
        if (isSatellite(3, hit_cnt_or)) return READ_SATELLITE;

        // Because this is the last layer, we can now use the BML
        // selector to calculate the score.
        for (int i = 0; i < _bf_amount3 - 1; i++) {
            if (hit_cnt_or[i] >= _hit_threshold) {
                rv = READ_MAPPED;
                // Calculate the CML location
                long cml_loc = base_offset + i * _seed_range3;
                int seq_len = _seed_range3 * 2;
                string ref_seq = getRefSeqFromLoc(cml_loc, seq_len);

                // Send one CML to the BML engine
                _bml_sel->update(ref_seq, read, cml_loc);
            }
        }
    }
    else {
        cerr << "Layer ID " << layer_id << " is not recognized" << endl;
        exit(1);
    }
    return rv;
}

void ShortReadMapper::updateScoreboard(int& rv, long& golden_loc,
                                       long& mapped_loc) {
    /*
    Return value:
    0th bit: read mapped
    1st bit: satellite
    */

    if (rv & READ_SATELLITE) {
        // Satellite
        _satellite += 1;
    }
    else if (rv & READ_MAPPED) {
        // Mapped
        if (abs(golden_loc - mapped_loc) <= _ans_margin)
            _correctly_mapped += 1;
        else
            _wrongly_mapped += 1;
    }
    else {
        // Not mapped
        _not_mapped += 1;
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

    // Bloom filters configuration
    // 16 MB for each Bloom filter in layer 1
    // 64 kB for each Bloom filter in layer 2
    // 256 Bytes for each Bloom filter in layer 3
    _bf_size1 = 16 * 1024 * 1024 * 8;
    _bf_size2 = 64 * 1024 * 8;
    _bf_size3 = 256 * 8;

    // 256 layer 1 Bloom filters
    // 256 layer 2 Bloom filters associated with one layer 1 BF
    // 256 layer 3 Bloom filters associated with one layer 2 BF
    _bf_amount1 = 256;
    _bf_amount2 = 256;
    _bf_amount3 = 256;

    _bf_total1 = 256;
    _bf_total2 = 256 * 256;
    _bf_total3 = 256 * 256 * 256;

    // Save 256 * 256 * 256 seeds in a Bloom filter in layer 1
    // Save 256 * 256 seeds in a Bloom filter in layer 2
    // Save 256 seeds in a Bloom filter in layer 3
    _seed_range1 = 256 * 256 * 256;
    _seed_range2 = 256 * 256;
    _seed_range3 = 256;

    // Mapping configuration
    _test_num = 10000;
    _ref_size = 2948627755;

    // Generate hash factor
    int rand_seed = 666;
    srand(1);
    random_device rd;
    default_random_engine generator(rand_seed);
    uniform_int_distribution<uint64_t> distribution(0, 0xFFFFFFFFFFFFFFFF);
    uint64_t hash_factor1 = distribution(generator);
    uint64_t hash_factor2 = distribution(generator);
    uint64_t hash_factor3 = distribution(generator);

    // Instantiate 3 layers
    _layer[0] = new Layer(_bf_size1, _bf_amount1, _bf_total1, _seed_range1,
                          hash_factor1);
    _layer[1] = new Layer(_bf_size2, _bf_amount2, _bf_total2, _seed_range2,
                          hash_factor2);
    _layer[2] = new Layer(_bf_size3, _bf_amount3, _bf_total3, _seed_range3,
                          hash_factor3);

    // Initialize _ref_seq
    _ref_seq = new char[_ref_size];

    // Instantiate BML selector
    _bml_sel = new BMLSelector();

    // Scoreboard
    _correctly_mapped = 0;
    _wrongly_mapped = 0;
    _satellite = 0;
    _not_mapped = 0;
}

ShortReadMapper::~ShortReadMapper() {
    for (int i = 0; i < 3; i++) {
        delete _layer[i];
    }
    delete _ref_seq;
    delete _bml_sel;
}

void ShortReadMapper::trainBF(bool ignoreSatellite) {
    cout << "Start training the Bloom filter" << endl;
    if (ignoreSatellite) cout << "Ignore satellite DNA" << endl;

    // Seed
    uint64_t seed = 0;
    long base_cnt = 0;

    // Open ref file
    ifstream ref_seq_fs(_ref_path);
    if (!ref_seq_fs.is_open()) {
        cerr << "Cannot open the reference sequence file." << endl;
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
                    cout << "Counted for " << base_cnt_fp << " seeds" << endl;
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
                        _layer[0]->update(seed, base_cnt);
                        _layer[1]->update(seed, base_cnt);
                        _layer[2]->update(seed, base_cnt);
                    }
                }
                else {
                    // Normal mode
                    _layer[0]->update(seed, base_cnt);
                    _layer[1]->update(seed, base_cnt);
                    _layer[2]->update(seed, base_cnt);
                }
            }

            base_cnt += 1;
            if (base_cnt % 10000000 == 0)
                cout << "Processed " << base_cnt << " seeds" << endl;
            if (base_cnt == _ref_size) break;
        }
    }

    // printSeedCnt(_seed_cnt);
}

void ShortReadMapper::writeBF() {
    // Currently not used
    _layer[0]->write_bf_bin("layer1_bf.dat");
    _layer[1]->write_bf_bin("layer2_bf.dat");
    _layer[2]->write_bf_bin("layer3_bf.dat");
}

void ShortReadMapper::readBF() {
    // Currently not used
    _layer[0]->read_bf_bin("layer1_bf.dat");
    _layer[1]->read_bf_bin("layer2_bf.dat");
    _layer[2]->read_bf_bin("layer3_bf.dat");
}

void ShortReadMapper::mapRead() {
    cout << "Start mapping the reads" << endl;

    // Open read file
    ifstream read_seq_fs(_read_path);
    if (!read_seq_fs.is_open()) {
        cerr << "Cannot open the read sequence file." << endl;
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

        // Query the read in each layer recursively
        initQuery();
        int rv = queryLayer(read, 1, 0, 0);

        // Get mapped location from the BML selector
        long mapped_loc = _bml_sel->getMapLoc();
        updateScoreboard(rv, golden_loc, mapped_loc);

        read_cnt += 1;
        if (read_cnt % 1000 == 0)
            cout << "Processed " << read_cnt << " reads" << endl;
    }
}

void ShortReadMapper::displayResult() {
    int sum = _correctly_mapped + _wrongly_mapped + _satellite + _not_mapped;

    cout << "Correctly mapped: " << setw(5) << _correctly_mapped << endl;
    cout << "Wrongly mapped:   " << setw(5) << _wrongly_mapped << endl;
    cout << "Satellite:        " << setw(5) << _satellite << endl;
    cout << "Not mapped:       " << setw(5) << _not_mapped << endl;
    cout << "Total:            " << setw(5) << sum << endl;
}
