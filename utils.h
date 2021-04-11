#include <map>
#include <unordered_map>

using namespace std;

void findAndIncrement(unordered_map<uint64_t, int>& _seed_cnt, uint64_t& seed) {
    unordered_map<uint64_t, int>::iterator it = _seed_cnt.find(seed);
    if (it != _seed_cnt.end()) {
        it->second++;
    }
    else {
        _seed_cnt.insert(make_pair(seed, 1));
    }
}

int getCount(unordered_map<uint64_t, int>& _seed_cnt, uint64_t& seed) {
    unordered_map<uint64_t, int>::const_iterator cnt = _seed_cnt.find(seed);
    if (cnt == _seed_cnt.end()) {
        // Not found, throw error
        cerr << "Seed not found in the _seed_cnt map" << endl;
        exit(1);
    }
    else {
        return cnt->second;
    }
}

void printSeedCnt(unordered_map<uint64_t, int>& _seed_cnt) {
    map<int, int> cnt_map;
    for (const auto& x : _seed_cnt) {
        map<int, int>::iterator it = cnt_map.find(x.second);
        if (it != cnt_map.end()) {
            it->second++;
        }
        else {
            cnt_map.insert(make_pair(x.second, 1));
        }
    }
    for (auto const& x : cnt_map) {
        cout << x.first << ": " << x.second << endl;
    }
}
