#include "short_read_mapper.h"

#include <unistd.h>

#include <fstream>
#include <iostream>
#include <random>
#include <string>

using namespace std;

void ShortReadMapper::genSeedMask() {
    uint64_t seed_mask = 0;
    for (int i = 0; i < _seed_len; ++i) seed_mask = (seed_mask << 2) + 3;
    _seed_mask = seed_mask;
}

void ShortReadMapper::readBase(uint64_t& out, char& base) {
    if (base == 'A' || base == 'a')
        out = (out << 2) + 0;
    else if (base == 'C' || base == 'c')
        out = (out << 2) + 1;
    else if (base == 'G' || base == 'g')
        out = (out << 2) + 2;
    else if (base == 'T' || base == 't')
        out = (out << 2) + 3;
}

void ShortReadMapper::queryLayer(string& read, int layer_id, long hier_offset,
                                 long base_offset) {
    /* In each layer, query every seeds of the read and
    record the hit count. If hit count > threshold, recursively
    qurey the next layer */

    if (layer_id == 1) {
        // Build layer 1 hit count array
        int hit_cnt[_bf_amount1];
        for (int i = 0; i < _bf_amount1; i++) hit_cnt[i] = 0;

        // Generate seed
        uint64_t seed = 0;
        for (int i = 0; i < read.length(); i++) {
            readBase(seed, read[i]);
            seed &= _seed_mask;

            // When the seed is not long enough, keep reading base
            if (i < _seed_len - 1) continue;

            // Query the layer 1
            _layer[0]->query(seed, hit_cnt, hier_offset);
        }

        // For each Bloom filter that has more than _hit_threshold hits,
        // query the next layer.
        for (int i = 0; i < _bf_amount1; i++) {
            if (hit_cnt[i] > _hit_threshold) {
                long hier_offset_next = i * _bf_amount2 * _bf_size2;
                long base_offset_next = i * _seed_range1;
                queryLayer(read, 2, hier_offset_next, base_offset_next);
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
            readBase(seed, read[i]);
            seed &= _seed_mask;

            // When the seed is not long enough, keep reading base
            if (i < _seed_len - 1) continue;

            // Query the layer 2
            _layer[1]->query(seed, hit_cnt, hier_offset);
        }

        // For each Bloom filter that has more than _hit_threshold hits,
        // query the next layer.
        for (int i = 0; i < _bf_amount1; i++) {
            if (hit_cnt[i] > _hit_threshold) {
                long hier_offset_next =
                    hier_offset + i * _bf_amount3 * _bf_size3;
                long base_offset_next = base_offset + i * _seed_range2;
                queryLayer(read, 3, hier_offset_next, base_offset_next);
            }
        }
    }
    else if (layer_id == 3) {
        // Build layer 3 hit count array
        int hit_cnt[_bf_amount3];
        for (int i = 0; i < _bf_amount3; i++) hit_cnt[i] = 0;

        // Generate seed
        uint64_t seed = 0;
        for (int i = 0; i < read.length(); i++) {
            readBase(seed, read[i]);
            seed &= _seed_mask;

            // When the seed is not long enough, keep reading base
            if (i < _seed_len - 1) continue;

            // Query the layer 3
            _layer[2]->query(seed, hit_cnt, hier_offset);
        }

        // Because this is the last layer, we can now use the BML
        // selector to calculate the score.
        for (int i = 0; i < _bf_amount3; i++) {
            if (hit_cnt[i] > _hit_threshold) {
                // Calculate the CML (base_offset)
                long cml = base_offset + i * _seed_range3;
            }
        }
    }
    else {
        cerr << "Layer ID " << layer_id << " is not recognized" << endl;
        exit(1);
    }
}

ShortReadMapper::ShortReadMapper(string ref_path, string read_path,
                                 long read_len, long seed_len,
                                 long query_shift_amt, long hit_threshold,
                                 long ans_margin, long satellate_threshold) {
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
    _satellate_threshold = satellate_threshold;

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
}

ShortReadMapper::~ShortReadMapper() {
    for (int i = 0; i < 3; i++) {
        delete _layer[i];
    }
}

void ShortReadMapper::trainBF() {
    cout << "Start training the Bloom filter" << endl;

    // Seed
    uint64_t seed = 0;
    int base_cnt = 0;

    // Open ref file
    ifstream ref_seq_fs(_ref_path);
    if (!ref_seq_fs.is_open()) {
        cerr << "Cannot open the reference sequence file." << endl;
        exit(1);
    }

    // Parse the ref file line by line
    string line;
    while (ref_seq_fs >> line) {
        // If the line starts with '>', ignore it
        if (line[0] == '>') continue;

        // For each character, generate a seed
        for (int i = 0; i < line.size(); i++) {
            if (base_cnt && base_cnt % 10000000 == 0)
                printf("Processed %d seeds\n", base_cnt);

            readBase(seed, line[i]);
            seed &= _seed_mask;

            // If the seed variable contains more than seed_len seeds,
            // start updating the Bloom filter.
            if (base_cnt >= _seed_len - 1) {
                _layer[0]->update(seed, base_cnt);
                _layer[1]->update(seed, base_cnt);
                _layer[2]->update(seed, base_cnt);
            }

            base_cnt += 1;
            if (base_cnt == _ref_size) break;
        }
    }
}

void ShortReadMapper::writeBF() {
    _layer[0]->write_bf_bin("layer1_bf.dat");
    _layer[1]->write_bf_bin("layer2_bf.dat");
    _layer[2]->write_bf_bin("layer3_bf.dat");
}

void ShortReadMapper::readBF() {
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
    int line_cnt = 0;
    string token;
    string golden_loc_str;
    string fwd_rev;
    string read;
    while (true) {
        /* Read format:
        >chr1   chr1-1536540    116446253       -
        <Original sequence>
        <Simulater generated sequence>
        */

        // Reference sequence name
        read_seq_fs >> token;
        // Read name
        read_seq_fs >> token;
        // Location
        read_seq_fs >> golden_loc_str;
        long golden_loc = stol(golden_loc_str);

        // Forward/reverse
        read_seq_fs >> fwd_rev;
        // Original sequence
        read_seq_fs >> token;
        // Simulator generated sequence
        read_seq_fs >> read;

        // Only map the forward sequence, ignore the reverse sequence
        if (fwd_rev == "-") continue;

        cout << read << endl;

        // Query the read in each layer recursively
        // queryLayer(read, layer_id, hier_offset, base_offset)
        queryLayer(read, 1, 0, 0);

        read_cnt += 1;
        line_cnt = (line_cnt + 1) % 3;
        if (read_cnt > _test_num) break;
    }
}
