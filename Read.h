#ifndef READ_H
#define READ_H

#include "hashtab.h"
#include <string>
#include <vector>
#include <fstream>
#include <bitset>

using namespace std;

const int bitsize = 500;

class modification {
    public:
        modification(short i, short t, short tp) {
            index = i;
            to    = t;
            type  = tp;
        };

        short index;
        short to;
        short type;
};

class modified_read {
    public:
        modified_read(vector<modification> & c, bitset<bitsize> & u, short re)
            :candidate(u) {
                coverage = 0;
                region_edits = re;
                for (int i = 0; i < c.size(); i++)
                    modifications.push_back(modification(c[i]));
            };
        modified_read(bitset<bitsize> & u, short re)
            :candidate(u) {
                coverage = 0;
                region_edits = re;
            };
        ~modified_read() {
        };

        vector<modification> modifications;
        bitset<bitsize> candidate;
        short region_edits;
        double coverage;
};

class Read {
    public:
        Read(const string & h, const unsigned int* s, const string& q, vector<int>& u, const int rl);
        ~Read();

        // string trim(int t);
        string modify(hashtab* all, double homo_coverage, double hete_coverage, double cutoff, double percentile);
        int modify_cc(vector<short> region, vector<int> candidate_subset, hashtab* all, double homo_coverage, double hete_coverage, double cutoff, double percentile);
        vector<short> candidate_region(vector<int> candidate_subset);
        vector<short> candidate_region_chop(vector<int> candidate_subset);
        bool check_trust(modified_read* mr, hashtab* all, vector<short> region, double homo_coverage, double hete_coverage, double cutoff, double percentile);
        bool first_check(hashtab* all, short edit_i, short nt, double homo_coverage, double hete_coverage, double cutoff, double percentile);
        bool first_check(hashtab* all, short edit_i, short nt, modified_read* mr, double homo_coverage, double hete_coverage, double cutoff, double percentile);
        string print_seq();
        string print_qual();
        string print_modified(vector<modification>& mor);
        string print_modified(vector<modification>& mor, int print_nt);

        string header;
        int read_length;
        // int trim_length;
        unsigned int* seq;
        vector<int> quals;
        vector<int> candidate;
        modified_read* consensused_read;

    private:
        bool candidate_intersect(vector<int> candidate_subset, vector<short>& region);
        void candidate_union(vector<int> candidate_subset, vector<short>& region);
        double get_coverage(vector<short> region, hashtab* all);
        // void quality_quicksort(vector<short> & indexes, int left, int right);
    
        const static bool aggressive = true;
        const static short expand_region = 1;
};

#endif

