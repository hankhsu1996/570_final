
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
    bool isHit(int&, int&);

   public:
    Layer(long, long, long, long, uint64_t);
    ~Layer();
    void update(uint64_t, long);
    void query(uint64_t&, int[], long);
    void write_bf(string);
};

#endif
