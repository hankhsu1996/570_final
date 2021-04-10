#include "bml_selector.h"

#include <vector>

BMLSelector::BMLSelector() { init(); }

BMLSelector::~BMLSelector() {}

void BMLSelector::init() {
    _max_score = 0;
    _map_loc = 0;
}

void BMLSelector::update(const string &ref_seq, const string &read,
                         long cml_loc) {
    uint32_t temp_end_row = 0;
    uint32_t temp_end_col = 0;
    int new_score = smith_waterman(ref_seq, read, temp_end_row, temp_end_col);

    if (new_score > _max_score) {
        _map_loc = cml_loc + temp_end_col - temp_end_row;
        _max_score = new_score;
    }
}

int BMLSelector::smith_waterman(const string &ref_seq, const string &read,
                                uint32_t &max_end_row, uint32_t &max_end_col) {
    vector<int> align_scorebuffer(1 + ref_seq.size(), 0);
    vector<int> insert_scorebuffer(1 + ref_seq.size(), 0);
    vector<int> delete_scorebuffer(1 + ref_seq.size(), 0);

    int align_diag = 0;
    int insert_diag = 0;
    int delete_diag = 0;

    int temp1 = 0;
    int temp2 = 0;
    int temp3 = 0;

    int temp_align = 0;
    int temp_insert = 0;
    int temp_delete = 0;

    int max_score = 0;

    for (int i = 0; i < read.size(); ++i) {
        if (i == 0) {
            align_diag = 0;
            insert_diag = 0;
            delete_diag = 0;
        }
        else {
            align_diag = -100;
            insert_diag = _gap_open_score + (i - 1) * _gap_extend_score;
            delete_diag = -100;
        }

        align_scorebuffer[0] = -100;
        insert_scorebuffer[0] = _gap_open_score + i * _gap_extend_score;
        delete_scorebuffer[0] = -100;

        for (int j = 0; j < ref_seq.size(); ++j) {
            // alignment
            if (ref_seq[j] == read[i]) {
                temp1 = align_diag + _match_score;
                temp2 = insert_diag + _match_score;
                temp3 = delete_diag + _match_score;
            }
            else {
                temp1 = align_diag + _mismatch_score;
                temp2 = insert_diag + _mismatch_score;
                temp3 = delete_diag + _mismatch_score;
            }

            if (temp1 >= temp2 && temp1 >= temp3 && temp1 > 0) {
                temp_align = temp1;
            }
            else if (temp2 > temp1 && temp2 >= temp3 && temp2 > 0) {
                temp_align = temp2;
            }
            else if (temp3 > temp1 && temp3 > temp2 && temp3 > 0) {
                temp_align = temp3;
            }
            else {
                temp_align = 0;
            }

            // insertion
            temp1 = align_scorebuffer[j + 1] + _gap_open_score;
            temp2 = insert_scorebuffer[j + 1] + _gap_extend_score;
            if (temp1 >= temp2 && temp1 > 0) {
                temp_insert = temp1;
            }
            else if (temp2 > temp1 && temp2 > 0) {
                temp_insert = temp2;
            }
            else {
                temp_insert = 0;
            }

            // deletion
            temp1 = align_scorebuffer[j] + _gap_open_score;
            temp2 = delete_scorebuffer[j] + _gap_extend_score;
            if (temp1 >= temp2 && temp1 > 0) {
                temp_delete = temp1;
            }
            else if (temp2 > temp1 && temp2 > 0) {
                temp_delete = temp2;
            }
            else {
                temp_delete = 0;
            }

            // smith-waterman local alignment
            // final score
            if (temp_align >= temp_insert && temp_align >= temp_delete &&
                temp_align > 0 && temp_align > max_score) {
                max_score = temp_align;
                max_end_row = i;
                max_end_col = j;
            }
            else if (temp_insert > temp_align && temp_insert >= temp_delete &&
                     temp_insert > 0 && temp_insert > max_score) {
                max_score = temp_insert;
                max_end_row = i;
                max_end_col = j;
            }
            else if (temp_delete > temp_align && temp_delete > temp_insert &&
                     temp_delete > 0 && temp_delete > max_score) {
                max_score = temp_delete;
                max_end_row = i;
                max_end_col = j;
            }

            align_diag = align_scorebuffer[j + 1];
            insert_diag = insert_scorebuffer[j + 1];
            delete_diag = delete_scorebuffer[j + 1];

            align_scorebuffer[j + 1] = temp_align;
            insert_scorebuffer[j + 1] = temp_insert;
            delete_scorebuffer[j + 1] = temp_delete;
        }
    }

    return max_score;
}

long BMLSelector::getMapLoc() { return _map_loc; }