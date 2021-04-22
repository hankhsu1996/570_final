#include <cmath>
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

void printHitCnt(int layer, int hit_cnt[], int amount) {
    if (layer == 0) {
        cout << "[index]  ";
        for (int i = 0; i < amount; i++) {
            cout << setw(2) << i % 100 << ' ';
        }
        cout << endl;
    }

    cout << "[layer" << layer << "] ";
    for (int i = 0; i < amount; i++) {
        cout << setw(2) << hit_cnt[i] << ' ';
    }
    cout << endl;
}

int countHitBF(int hit_cnt[], int amount, int threshold) {
    int sum = 0;
    for (int i = 0; i < amount; i++) {
        if (hit_cnt[i] >= threshold) {
            sum += 1;
        }
    }
    return sum;
}

long meanPlusStdev(int arr[], int len, int N) {
    float sum = 0;
    float mean;
    float var = 0;
    for (int i = 0; i < len; i++) {
        sum += arr[i];
    }
    mean = sum / len;
    for (int i = 0; i < len; i++) {
        var += pow(arr[i] - mean, 2);
    }
    return (long)mean + N * (long)sqrt(var / len);
}

int findMax(int arr[], int len) {
    int rv = -1;
    for (int i = 0; i < len; i++) {
        rv = max(rv, arr[i]);
    }
    return rv;
}
