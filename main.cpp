#include <iostream>

#include "short_seq_mapper.h"
int main(int argc, char const* argv[]) {
    // File paths
    string ref_path = "../dataset/hg38_short.fa";
    string read_path = "../dataset/single_HSXn_150bp_1f.aln";

    // Configuration
    long read_len = 150;
    long seed_len = 20;
    // When querying, shift the seed by N.
    long query_shift_amt = 1;
    // If more than N seeds are exact matches, then a CML is found.
    long hit_threshold = 37;
    // If the location returned by the Bloom filter - answer < N,
    // then it is a successful mapping.
    long ans_margin = 10;
    // If the Bloom filter has #CMLs > N, then it is considered a satellite DNA.
    long satellite_threshold = 10;

    ShortSeqMapper mapper =
        ShortSeqMapper(ref_path, read_path, read_len, seed_len, query_shift_amt,
                       hit_threshold, ans_margin, satellite_threshold);

    cout << "Start training the Bloom filter" << endl;
    mapper.train();
    // mapper.map();

    return 0;
}