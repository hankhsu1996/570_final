
#include <cstdlib>
#include <string>

#include "layer.h"

using namespace std;

#ifndef __SHORT_SEQ_MAPPER__
#define __SHORT_SEQ_MAPPER__

class ShortReadMapper {
   private:
    // File paths
    string _ref_path;
    string _read_path;

    // Basic configuration
    long _read_len;
    long _seed_len;
    uint64_t _seed_mask;
    long _query_skip_amt;
    long _hit_threshold;
    long _ans_margin;
    long _satellate_threshold;

    // Bloom filters configuration
    // 16 MB for each Bloom filter in layer 1
    // 64 kB for each Bloom filter in layer 2
    // 256 Bytes for each Bloom filter in layer 3
    long _bf_size1 = 16 * 1024 * 1024 * 8;
    long _bf_size2 = 64 * 1024 * 8;
    long _bf_size3 = 256 * 8;

    // 256 layer 1 Bloom filters
    // 256 layer 2 Bloom filters associated with one layer 1 BF
    // 256 layer 3 Bloom filters associated with one layer 2 BF
    long _bf_amount1 = 256;
    long _bf_amount2 = 256;
    long _bf_amount3 = 256;

    long _bf_total1 = 256;
    long _bf_total2 = 256 * 256;
    long _bf_total3 = 256 * 256 * 256;

    // Save 256 * 256 * 256 seeds in a Bloom filter in layer 1
    // Save 256 * 256 seeds in a Bloom filter in layer 2
    // Save 256 seeds in a Bloom filter in layer 3
    long _seed_range1 = 256 * 256 * 256;
    long _seed_range2 = 256 * 256;
    long _seed_range3 = 256;

    Layer* _layer[3];
    // Layer a;

    // Mapping configuration
    long _test_num = 10000;
    long _ref_size = 2948627755;

    void genSeedMask();
    void readBase(uint64_t& out, char& base);
    void queryLayer(string&, int, long, long);

   public:
    ShortReadMapper(string ref_path, string read_path, long read_len,
                    long seed_len, long query_shift_amt, long hit_threshold,
                    long ans_margin, long satellate_threshold);

    ~ShortReadMapper();

    void train();
    void map();
};

#endif
