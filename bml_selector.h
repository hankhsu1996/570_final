#include <string>

using namespace std;

#ifndef __BML_SELECTOR__
#define __BML_SELECTOR__

class BMLSelector {
   private:
    // Score setting
    int _match_score = 1;
    int _mismatch_score = -1;
    int _gap_open_score = -3;
    int _gap_extend_score = -2;

    // Current max score
    int _max_score;
    long _map_loc;

   public:
    BMLSelector();
    ~BMLSelector();
    void init();
    void update(const string &, const string &, long);
    int smith_waterman(const string &, const string &, uint32_t &, uint32_t &);
    long getMapLoc();
};

#endif
