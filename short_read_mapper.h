
#include <cstdlib>
#include <string>

#include "bml_selector.h"
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
    long _bf_size1;
    long _bf_size2;
    long _bf_size3;

    // 256 layer 1 Bloom filters
    // 256 layer 2 Bloom filters associated with one layer 1 BF
    // 256 layer 3 Bloom filters associated with one layer 2 BF
    long _bf_amount1;
    long _bf_amount2;
    long _bf_amount3;

    long _bf_total1;
    long _bf_total2;
    long _bf_total3;

    // Save 256 * 256 * 256 seeds in a Bloom filter in layer 1
    // Save 256 * 256 seeds in a Bloom filter in layer 2
    // Save 256 seeds in a Bloom filter in layer 3
    long _seed_range1;
    long _seed_range2;
    long _seed_range3;

    // Mapping configuration
    long _test_num;
    long _ref_size;

    // Hierarchical Bloom filter layers
    Layer* _layer[3];

    // Full reference sequence
    char* _ref_seq;

    // BML selector
    BMLSelector* _bml_sel;

    // Result
    int _correctly_mapped;
    int _wrongly_mapped;
    int _satellite;
    int _not_mapped;

    void genSeedMask();
    void updateSeed(char&, uint64_t&);
    void updateRefSeq(char&, long);
    int queryLayer(string&, int, long, long);
    string getRefSeqFromLoc(long, int);
    bool isSatellite(int[], int);

   public:
    ShortReadMapper(string&, string&, long, long, long, long, long, long);

    ~ShortReadMapper();

    void trainBF();
    void writeBF();
    void readBF();
    void mapRead();
    void displayResult();
};

#endif
