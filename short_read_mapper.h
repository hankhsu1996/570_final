
#include <cstdlib>
#include <ctime>
#include <string>
#include <unordered_map>

#include "bml_selector.h"
#include "layer.h"

using namespace std;

#ifndef __SHORT_SEQ_MAPPER__
#define __SHORT_SEQ_MAPPER__

class Stopwatch {
   private:
    clock_t _duration;
    clock_t _start_time;

   public:
    void reset() {
        _duration = 0;
        _start_time = 0;
    }
    void start() { _start_time = clock(); }
    void pause() {
        _duration += clock() - _start_time;
        _start_time = 0;
    }
    float getSec() { return (float)_duration / (float)CLOCKS_PER_SEC; }
};

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
    long _satellite_threshold;

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

    // Scoreboard
    int _correctly_mapped;
    int _wrongly_mapped;
    int _satellite;
    int _not_mapped;

    // Total hit count in each layer
    // If hit_cnt > _satellite_threshold, read is satellite
    int _layer1_hit_cnt;
    int _layer2_hit_cnt;
    int _layer3_hit_cnt;

    // Seed count used to ignore satellite when training BF
    unordered_map<uint64_t, int> _seed_cnt;

    // Stopwatch
    Stopwatch* _training_sw;
    Stopwatch* _seeding_sw;
    Stopwatch* _seed_extraction_sw;

    // Private functions
    void genSeedMask();
    void updateSeed(char&, uint64_t&);
    void updateRefSeq(char&, long);
    void initQuery();
    int queryLayer(string&, int, long, long);
    void updateScoreboard(int&, long&, long&);
    string getRefSeqFromLoc(long, int);
    bool isSatellite(int, int[]);

   public:
    ShortReadMapper(string&, string&, long, long, long, long, long, long);
    ~ShortReadMapper();
    void trainBF(bool);
    void writeBF();
    void readBF();
    void mapRead();
    void displayResult();
};

#endif
