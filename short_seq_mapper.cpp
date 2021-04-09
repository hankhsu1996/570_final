#include "short_seq_mapper.h"

#include <unistd.h>

#include <fstream>
#include <iostream>
#include <random>
#include <string>

using namespace std;

void ShortSeqMapper::genSeedMask() {
    uint64_t seed_mask = 0;
    for (int i = 0; i < _seed_len; ++i) seed_mask = (seed_mask << 2) + 3;
    _seed_mask = seed_mask;
}

ShortSeqMapper::ShortSeqMapper(string ref_path, string read_path, long read_len,
                               long seed_len, long query_shift_amt,
                               long hit_threshold, long ans_margin,
                               long satellate_threshold) {
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

    _layer[0] = new Layer(_bf_size1, _bf_amount1, _bf_total1, _seed_range1,
                          hash_factor1);
    _layer[1] = new Layer(_bf_size2, _bf_amount2, _bf_total2, _seed_range2,
                          hash_factor2);
    _layer[2] = new Layer(_bf_size3, _bf_amount3, _bf_total3, _seed_range3,
                          hash_factor3);
}

ShortSeqMapper::~ShortSeqMapper() {
    for (int i = 0; i < 3; i++) {
        delete _layer[i];
    }
}

void ShortSeqMapper::train() {
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

            if (line[i] == 'A' || line[i] == 'a')
                seed = (seed << 2) + 0;
            else if (line[i] == 'C' || line[i] == 'c')
                seed = (seed << 2) + 1;
            else if (line[i] == 'G' || line[i] == 'g')
                seed = (seed << 2) + 2;
            else if (line[i] == 'T' || line[i] == 't')
                seed = (seed << 2) + 3;
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

#ifdef WRITE_BF
    _layer[0]->write_bf("layer1_bf.hex");
    _layer[1]->write_bf("layer2_bf.hex");
    _layer[2]->write_bf("layer3_bf.hex");
#endif
}

void ShortSeqMapper::map() {
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
    while (getline(read_seq_fs, line)) {
        if (read_cnt > _test_num) break;
        usleep(100000);
    }
}
