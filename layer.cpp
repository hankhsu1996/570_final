
#include "layer.h"

#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>

void Layer::genBFMask() {
    _bf_bitwidth = log2(double(_bf_size));

    for (int i = 0; i < _bf_bitwidth; i++) {
        _bf_mask = (_bf_mask << 1) + 1;
    }
}

bool Layer::isHit(int mem_content, int mem_bit) {
    return (mem_content & (1 << (31 - mem_bit))) != 0;
}

Layer::Layer(long bf_size, long bf_amount, long bf_total, long seed_range,
             uint64_t& hash_factor) {
    _bf_size = bf_size;
    _bf_amount = bf_amount;
    _bf_total = bf_total;
    _seed_range = seed_range;
    genBFMask();
    _hash_factor = hash_factor;

    // Memory array
    _mem_size = (bf_size / 32) * bf_total;
    _memory = new int[_mem_size];

    // Memory arrangement
    _mem_arrangement = INTERLEAVED;
}

Layer::~Layer() { delete _memory; }

void Layer::update(uint64_t& seed, long base_cnt) {
    // hash function
    uint64_t hash_val = (seed ^ _hash_factor) & _bf_mask;

    // Hierarchical Bloom filter
    long last_layer_range = (long)_seed_range * (long)_bf_amount;
    long nth_last_layer = base_cnt / last_layer_range;
    long hier_offset = nth_last_layer * _bf_amount * _bf_size;

    if (_mem_arrangement == INORDERED) {
        /* Memory content:
        N = _bf_size
        bf[0] bit 0, 1, ..., N
        bf[1] bit 0, 1, ..., N
        ...
        bf[255] bit 0, 1, ..., N
        */
        long bit_offset = hash_val;
        long bf_offset =
            ((base_cnt % last_layer_range) / _seed_range) * _bf_size;
        long mem_idx = hier_offset + bit_offset + bf_offset;
        long mem_addr = mem_idx / 32;
        long mem_bit = mem_idx % 32;
        _memory[mem_addr] |= 1 << (31 - mem_bit);
    }
    else if (_mem_arrangement == INTERLEAVED) {
        /* Memory content:
        bit 0 of bf[0] bf[1] ... bf[255]
        bit 1 of bf[0] bf[1] ... bf[255]
        ...
        bit N of bf[0] bf[1] ... bf[255]
        */
        long bit_offset = hash_val * _bf_amount;
        long bf_offset = (base_cnt % last_layer_range) / _seed_range;
        long mem_idx = hier_offset + bit_offset + bf_offset;
        long mem_addr = mem_idx / 32;
        long mem_bit = mem_idx % 32;
        _memory[mem_addr] |= 1 << (31 - mem_bit);
    }
}

void Layer::query(uint64_t& seed, int hit_cnt[], long hier_offset,
                  bool or_next) {
    // hash_function
    uint64_t hash_val = (seed ^ _hash_factor) & _bf_mask;

    if (_mem_arrangement == INORDERED) {
        long bit_offset = hash_val;
        for (int i = 0; i < _bf_amount; i++) {
            long bf_offset = i * _bf_size;
            long mem_idx = hier_offset + bit_offset + bf_offset;
            long mem_addr = mem_idx / 32;
            int mem_bit = mem_idx % 32;

            bool hit = isHit(_memory[mem_addr], mem_bit);

            if (or_next && i < _bf_amount - 1) {
                mem_idx += _bf_size;
                mem_addr = mem_idx / 32;
                mem_bit = mem_idx % 32;
                hit |= isHit(_memory[mem_addr], mem_bit);
            }

            hit_cnt[i] += hit;
        }
    }
    else if (_mem_arrangement == INTERLEAVED) {
        long bit_offset = hash_val * _bf_amount;
        for (int i = 0; i < _bf_amount; i++) {
            long bf_offset = i;
            long mem_idx = hier_offset + bit_offset + bf_offset;
            long mem_addr = mem_idx / 32;
            int mem_bit = mem_idx % 32;

            bool hit = isHit(_memory[mem_addr], mem_bit);

            if (or_next && i < _bf_amount - 1) {
                mem_idx += 1;
                mem_addr = mem_idx / 32;
                mem_bit = mem_idx % 32;
                hit |= isHit(_memory[mem_addr], mem_bit);
            }

            hit_cnt[i] += hit;
        }
    }
}

void Layer::write_bf_hex(string path) {
    ofstream bf_os(path);

    if (!bf_os.is_open()) {
        cerr << "Cannot open " << path << endl;
        exit(1);
    }

    cout << "Write Bloom filter content to file " << path << endl;

    for (int i = 0; i < _mem_size; i++) {
        bf_os << hex << setw(8) << setfill('0') << _memory[i];
        if (i % 16 == 15) bf_os << endl;
    }
}

void Layer::write_bf_bin(string path) {
    ofstream bf_os(path, ios::out | ios::binary);

    if (!bf_os.is_open()) {
        cerr << "Cannot open " << path << endl;
        exit(1);
    }

    cout << "Write Bloom filter content to file " << path << endl;

    for (int i = 0; i < _mem_size; i++) {
        bf_os.write((char*)&_memory[i], sizeof(int));
    }
    bf_os.close();
}

void Layer::read_bf_bin(string path) {
    ifstream bf_is(path, ios::in | ios::binary);

    if (!bf_is.is_open()) {
        cerr << "Cannot open " << path << endl;
        exit(1);
    }

    cout << "Read Bloom filter content from file " << path << endl;

    for (int i = 0; i < _mem_size; i++) {
        bf_is.read((char*)&_memory[i], sizeof(int));
    }
    bf_is.close();
}
