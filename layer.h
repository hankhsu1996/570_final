
#include <string>
using namespace std;

#ifndef __LAYER__
#define __LAYER__

class Layer {
   private:
    long _bf_size;
    long _bf_bitwidth;
    uint64_t _bf_mask;
    long _bf_amount;
    long _bf_total;
    long _seed_range;
    long _mem_size;

    // hash function parameters
    uint64_t _hash_factor;

    // Bloom filter memory
    int* _memory;
    void genBFMask();

   public:
    Layer(long bf_size, long bf_amount, long bf_total, long seed_range,
          uint64_t hash_factor);
    ~Layer();
    void update(uint64_t seed, long base_cnt);
    void write_bf(string path);
};

#endif
