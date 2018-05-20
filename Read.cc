#include "Read.h"
#include "hashtab.h"
#include <iostream>
#include <math.h>
#include <algorithm>
#include <set>
#include <queue>
#include <assert.h>

// #define TESTING false
// #define _PRINTING true

int hashtab::k;

int k = 0;

Read::Read(const string& h, const unsigned int* s, const string& q, vector<int>& u, const int rl)
    :candidate(u) {
        header = h;
        read_length = rl;
        read_length = rl;
        seq = new unsigned int[read_length];
        for (int i = 0; i < read_length; i++) {
            seq[i] = s[i];
            quals.push_back(q[i]);
        }
        consensused_read = 0;
    }

Read::~Read() {
    delete[] seq;
    if (consensused_read != 0)
        delete consensused_read;
}

string Read::modify(hashtab* all, double homo_coverage, double hete_coverage, double cutoff, double percentile) {
    k = hashtab::k;
    vector< vector<int> > cc_candidate;

    // add first
    cc_candidate.push_back(vector<int>());
    int cc = 0;
    cc_candidate[cc].push_back(candidate[0]);

    for (int i = 1; i < candidate.size(); i ++) {
        // if (candidate[i-1]+hashtab::k-1 < candidate[i]) {
        if (candidate[i] - candidate[i-1] > 1) {
            cc++;
            cc_candidate.push_back(vector<int>());
        }
        cc_candidate[cc].push_back(candidate[i]);
    }

    vector<modification> multi_mors;
    vector<short> chop_region;
    // vector<short> big_region;
    // int chop_modify_code, big_modify_code;
    int chop_modify_code;

    for (cc = 0; cc < cc_candidate.size(); cc++) {
#ifdef _PRINTING
cout << "cc_candidate: " << endl;
for (vector<int>::iterator it = cc_candidate[cc].begin(); it != cc_candidate[cc].end(); it ++) 
    cout << *it << ", ";
cout << endl;
#endif
        chop_region = candidate_region_chop(cc_candidate[cc]);
        // cout << "chop region size: " << chop_region.size() << endl;
        chop_modify_code = modify_cc(chop_region, cc_candidate[cc], all, homo_coverage, hete_coverage, cutoff, percentile);
        // cout << "chop modify code: " << chop_modify_code << endl;
        if (chop_modify_code > 0) {
            continue;
            // return print_modified(multi_mors, cc_candidate[cc].front());
            // try bigger region
            // big_region = candidate_region(cc_candidate[cc]);
            // if (big_region.size() == chop_region.size()) {
            //     // return print_modified(multi_mors, read_length);
            //     continue;
            // } else {
            //     big_modify_code = modify_cc(big_region, cc_candidate[cc], all, homo_coverage, hete_coverage, cutoff);
            //     if (big_modify_code > 0)
            //         // return print_modified(multi_mors, read_length);
            //         continue;
            // }
        }
        // modified!
        for (int c = 0; c < consensused_read->modifications.size(); c++)
            multi_mors.push_back(consensused_read->modifications[c]);
    }

    // modified_read* tmp = consensused_read;
    // consensused_read = new modified_read(multi_mors, tmp->candidate, 0);
    // delete tmp;

    return print_modified(multi_mors);
}


int Read::modify_cc(vector<short> region, vector<int> candidate_subset, hashtab* all, double homo_coverage, double hete_coverage, double cutoff, double percentile) {
    unsigned int max_queue_size = 4000000;

    if (region.size() > 0) 
        sort(region.begin(), region.end());
    else
        return 3;

#ifdef _PRINTING // test
    cout << "candidate kmers: ";
    for (int i = 0; i != candidate_subset.size(); i ++)
        cout << candidate_subset[i] << ", ";
    cout << endl;

    cout << "region: ";
    for (int i = 0; i != region.size(); i ++)
        cout << region[i] << ", ";
    cout << endl;
#endif

    queue<modified_read*> cpq;
    // corrections_compare cpq_comp;

    ////////////////////////////////////////
    // initialize
    ////////////////////////////////////////
    modified_read *cr, *next_cr;
    short edit_i;
    bitset<bitsize> bitcandidate;
    for(int i = 0; i < candidate_subset.size(); i++) {
        if(candidate_subset[i] >= bitsize) {
            cerr << "These reads must be longer than assumed. Increase the variable 'bitsize' in 'Read.h' to the read length." << endl;
            exit(EXIT_FAILURE);
        } else
            bitcandidate.set(candidate_subset[i]);
    }

    double cov = get_coverage(region, all);
#ifdef _PRINTING
    cout << "original coverage: " << cov << endl;
    cout << "region size: " << region.size() << endl;
#endif
    // push one bit correction into queue
    
    // for (short i = 0; i < region.size(); i ++) {
    //     short region_edit;
    //     if (i % 2 == 0)
    //         region_edit = i/2;
    //     else
    //         region_edit = region.size()-1-i/2;
    //     edit_i = region[region_edit];
    //     for (short nt = 0; nt < 4; nt ++) {
    //         if (seq[edit_i] != nt) {
    //             next_cr = new modified_read(bitcandidate, region_edit+1);
    //             next_cr->modifications.push_back(modification(edit_i, nt, 1));

    //             cpq.push(next_cr);
    //         }
    //     }
    // }
    
    for (short region_edit = 0; region_edit < region.size(); region_edit ++) {
        edit_i = region[region_edit];

        for (short nt = 0; nt < 4; nt ++) {
            if (seq[edit_i] != nt) {
                next_cr = new modified_read(bitcandidate, region_edit+1);
                next_cr->modifications.push_back(modification(edit_i, nt, 1));

                cpq.push(next_cr);
            }
        }
    }
    
    // int size = max(1, hashtab::k - region.at(0));
    // for (short region_edit = 0; region_edit < size; region_edit ++) {
    //     edit_i = region[region_edit];
    //     for (short nt = 0; nt < 4; nt ++) {
    //         // cout << "start first check" << endl;
    //         if (first_check(all, edit_i, nt, homo_coverage, hete_coverage, cutoff, percentile)) {
    //             // cout << "first check passed!" << endl;
    //             // cout << "edit_i: " << edit_i << endl;
    //             // cout << "nt: " << nt << endl;
    //             next_cr = new modified_read(bitcandidate, region_edit+1);
    //             next_cr->modifications.push_back(modification(edit_i, nt, 1));

    //             cpq.push(next_cr);
    //         }
    //     }
    // }

    consensused_read = 0;


    signed int candidate_count;

    while(cpq.size() > 0) {    

        /////////////////////////
        // quit if pq is too big
        /////////////////////////
        if(cpq.size() > max_queue_size) {

            if (consensused_read != 0) {
                delete consensused_read;
                consensused_read = 0;
            }
            break;
        }

        cr = cpq.front();
        cpq.pop();

        candidate_count = (signed int)cr->candidate.count();
#ifdef _PRINTING
        cout << "candidate count before check trust: " << candidate_count << endl;
#endif
        /////////////////////////
        // check trust
        /////////////////////////
        // save for later comparison
        if(check_trust(cr, all, region, homo_coverage, hete_coverage, cutoff, percentile)) {

#ifdef _PRINTING
            cout << "check_trust passed, new coverage: " << cr->coverage << endl;
#endif

            // cout << cr->coverage << endl;
            if (cr->coverage > cov && cr->coverage < percentile * hete_coverage) {
                // cout << "yeah! we get this!" << endl;
                consensused_read = cr;
                cov = cr->coverage;
            }
        }
#ifdef _PRINTING
        cout << "candidate count after check trust: " << cr->candidate.count() << endl;
#endif

        // if ((signed int)cr->candidate.count() <= candidate_count) {
        //     for (short region_edit = cr->region_edits; region_edit < min((int)region.size(), cr->region_edits+hashtab::k); region_edit++) {
        //         edit_i = region[region_edit];
        //         for (short nt = 0; nt < 4; nt ++) {
        //             if (first_check(all, edit_i, nt, cr, homo_coverage, hete_coverage, cutoff, percentile)) {
        //                 next_cr = new modified_read(cr->modifications, cr->candidate, region_edit + 1);
        //                 next_cr->modifications.push_back(modification(edit_i, nt, 1));

        //                 cpq.push(next_cr);
        //             }
        //         }
        //     }
        // }

                    

        if ((signed int)cr->candidate.count() <= candidate_count && cr->modifications.size() < 1) {
            // bool cr_added = true;
            for (short region_edit = cr->region_edits; region_edit < region.size(); region_edit++) {
            // for (short i = 0; i < region.size()-(cr->region_edits); i++) {
                // short region_edit;
                // if (i % 2 == 0) 
                //     region_edit = i/2+cr->region_edits;
                // else
                //     region_edit = region.size()-i/2-1;
                edit_i = region[region_edit];

                // cr_added = false;
                for (short nt = 0; nt < 4; nt ++) {
                    if (seq[edit_i] != nt) {
                        // unsigned int* newseq = new unsigned int[read_length];
                        // for (int j = 0; j < read_length; j ++) 
                        //     newseq[j] = seq[j];
                        // for (int c = 0; c < cr->modifications.size(); c ++)
                        //     newseq[cr->modifications[c].index] = cr->modifications[c].to;
                        // int leftmost_kmer = max(0, region_edit - hashtab::k + 1);
                        // int rightmost_kmer = min(read_length-hashtab::k, int(region_edit));
                        // if (all->getcount(&newseq[leftmost_kmer]) <= cutoff && all->getcount(&newseq[rightmost_kmer]) <= cutoff) {
                        //     delete[] newseq;
                        //     continue;
                        // } 
                        // delete[] newseq;

                        next_cr = new modified_read(cr->modifications, cr->candidate, region_edit+1);
                        next_cr->modifications.push_back(modification(edit_i, nt, 1));

                        // cpq.push_back(next_cr);
                        // push_heap(cpq.begin(), cpq.end(), cpq_comp);
                        // cr_added = true;
                        cpq.push(next_cr);
                    }
                }
            }
        }


        if (consensused_read != cr) {
            // cout << "oops, deleted" << endl;
            delete cr;
        }
    }

    // clean up queue
    // for(int i = 0; i < cpq.size(); i++)
    //     delete cpq[i];
    while (!cpq.empty()) {
        delete cpq.front();
        cpq.pop();
    }

    if(consensused_read != 0) {
        // cout << "corrected, modified bits: " << trusted_read->corrections.size() << endl;
        //cerr << header << "\t" << region.size() << "\t" << untrusted_subset.size() << "\t" << nt90 << "\t" << nt99 << "\t" << exp_errors << "\t" << cpq_adds << "\t" << check_count << "\t1\t" << trusted_read->likelihood << endl;
        // change affected kmers to make sure all similiar kmers change to the same direction
        return 0;
    } else {
        return 1;
    }
}


vector<short> Read::candidate_region(vector<int> candidate_subset) {
    vector<short> region;
    if (!candidate_intersect(candidate_subset, region))
        candidate_union(candidate_subset, region);

    short f = region.front();
    short b = region.back();

    if (k-1 >= f) {
        for (short i = f-1; i >= 0; i --)
            region.push_back(i);
    }
    if (read_length-k <= b) {
        for (short i = b+1; i < read_length; i ++)
            region.push_back(i);
    }

    return region;
}


vector<short> Read::candidate_region_chop(vector<int> candidate_subset) {
    // cout << "candidate subset size: " << candidate_subset.size() << endl;
    vector<short> region;
    // cout << "can i get intersect: " << candidate_intersect(candidate_subset, region);
    region.clear();
    if (!candidate_intersect(candidate_subset, region))
        candidate_union(candidate_subset, region);
    // if (!candidate_intersect(candidate_subset, region)) {
    //     return region;
    // }  
    
    // fix front
    int right_leftkmer = candidate_subset.front() - 1;
    if (right_leftkmer < 0) {
        for (int i = region[0]-1; i >= 0; i --) 
            region.push_back(i);
    } else {
        vector<short> front_chop(region);
        region.clear();
        for (int i = 0; i < front_chop.size(); i ++) {
            if (front_chop[i] > right_leftkmer + k - 1)
                region.push_back(front_chop[i]);
        }
    }

    // fix back
    int left_rightkmer = candidate_subset.back() + 1;
    if (left_rightkmer+k-1 >= read_length) {
        for (int i = region.back()+1; i < read_length; i++)
            region.push_back(i);
    } else {
        vector<short> back_chop(region);
        region.clear();
        for (int i = 0; i < back_chop.size(); i ++) {
            if (back_chop[i] < left_rightkmer)
                region.push_back(back_chop[i]);
        }
    }
    // if (left_rightkmer+k-1 < read_length) {
    //     vector<short> back_chop(region);
    //     region.clear();
    //     for (int i = 0; i < back_chop.size(); i++) {
    //         if (back_chop[i] < left_rightkmer)
    //             region.push_back(back_chop[i]);
    //     }
    // } else {
    //     vector<short> chop_region(region);
    //     region.clear();
    //     for (int i = chop_region.back()+1; i < read_length; i++) 
    //         region.push_back(i);
    //     for (int i = 0; i < chop_region.size(); i ++) {
    //         if (chop_region[i] > left_rightkmer-1)
    //             region.push_back(chop_region[i]);
    //     }
    // }

    sort(region.begin(), region.end());

#ifdef _PRINTING
cout << "region: " << endl;
for (vector<short>::iterator it = region.begin(); it != region.end(); it ++)
    cout << *it << ", ";
cout << endl;
#endif

    return region;
}

bool Read::first_check(hashtab* all, short edit_i, short nt, double homo_coverage, double hete_coverage, double cutoff, double percentile) {
    unsigned int* newseq = new unsigned int[read_length];
    for (int j = 0; j < read_length; j ++) 
        newseq[j] = seq[j];
    newseq[edit_i] = nt;

    int kmer_start = max(0, edit_i - hashtab::k + 1);
    double cov = all->getcount(&newseq[kmer_start]);
    bool ret = true;
    if (cov <= cutoff || cov >= percentile * hete_coverage) 
        ret = false;

    delete[] newseq;

    return ret;
}

bool Read::first_check(hashtab* all, short edit_i, short nt, modified_read* mr, double homo_coverage, double hete_coverage, double cutoff, double percentile) {
    unsigned int* newseq = new unsigned int[read_length];
    for (int j = 0; j < read_length; j ++) 
        newseq[j] = seq[j];
    for (int c = 0; c < mr->modifications.size(); c ++)
        newseq[mr->modifications[c].index] = mr->modifications[c].to;
    newseq[edit_i] = nt;

    int kmer_start = max(0, edit_i - hashtab::k + 1);
    double cov = all->getcount(&newseq[kmer_start]);
    bool ret = true;
    if (cov <= cutoff || cov >= percentile * hete_coverage)
        ret = false;

    delete[] newseq;
    
    return ret;
}

bool Read::check_trust(modified_read* mr, hashtab* all, vector<short> region, double homo_coverage, double hete_coverage, double cutoff, double percentile) {
    if (mr->modifications.empty())
        return false;

    unsigned int* newseq = new unsigned int[read_length];
    for (int j = 0; j < read_length; j ++) 
        newseq[j] = seq[j];
    for (int c = 0; c < mr->modifications.size(); c ++)
        newseq[mr->modifications[c].index] = mr->modifications[c].to;

#ifdef _PRINTING    // test
    cout << "modifications: " << endl;
    for (int i = 0; i < mr->modifications.size(); i ++) {
        cout << "idx: " << mr->modifications[i].index << endl;
        cout << "to: " << mr->modifications[i].to << endl;
    }
//     cout << "original sequence: ";
//     for (int i = 0; i < read_length; i ++)
//         cout << seq[i];
//     cout << endl;
//     cout << "new sequence: ";
//     for (int i = 0; i < read_length; i ++)
//         cout << newseq[i];
//     cout << endl;
#endif    // test

    int spos = region.front();
    int epos = region.back();

    int kmer_start = max(0, spos - hashtab::k + 1);
    int kmer_end = min(read_length-hashtab::k, epos);

    double cov = 0;
    // int count = 0;

    // cout << "this is a modification" << endl;
    // cout << "region size: " << region.size() << endl;
    for (int i = kmer_start; i <= kmer_end; i++) {
        // count += 1;
        double ccov = all->getcount(&newseq[i]);
        cov += ccov;

#ifdef _PRINTING
        cout << "mer index: " << i << ", coverage: " << all->getcount(&newseq[i]) << endl;
#endif

        if (ccov <= cutoff || ccov >= percentile * homo_coverage) {
            mr->candidate.set(i);
        } else {
            mr->candidate.reset(i);
        }
    }
    mr->coverage = cov / (kmer_end - kmer_start + 1);
    // cout << "new coverage: " << cov << endl;

    delete[] newseq;
    return (mr->candidate.none());
}


string Read::print_seq() {
    char nts[5] = {'A', 'C', 'G', 'T', 'N'};
    string sseq;
    for (int i = 0; i < read_length; i ++) 
        sseq.push_back(nts[seq[i]]);
    return sseq;
}


string Read::print_qual() {
    string squal;
    for (int i = 0; i < read_length; i ++)
        squal.push_back((char)quals[i]);
    return squal;
}


string Read::print_modified(vector<modification>& mor) {
    return print_modified(mor, read_length);
}


string Read::print_modified(vector<modification>& mor, int print_nt) {
    char nts[5] = {'A', 'C', 'G', 'T', 'N'};
    string sseq;
    vector<int> squals;
    int modify_i;

    for (int i = 0; i < print_nt; i ++) {
        modify_i = -1;
        for (int c = 0; c < mor.size(); c++) {
            if (mor[c].index == i)
                modify_i = c;
        }
        if (modify_i != -1)
            sseq.push_back(nts[mor[modify_i].to]);
        else
            sseq.push_back(nts[seq[i]]);
    }
    return sseq;
}


bool Read::candidate_intersect(vector<int> candidate_subset, vector<short>& region) {
    int start = 0;
    int end = read_length - 1;

    int u;
    for (int i = 0; i < candidate_subset.size(); i ++) {
        u = candidate_subset[i];

        if (start <= u+hashtab::k-1 && u <= end) {
            start = max(start, u);
            end = min(end, u+hashtab::k-1);
        } else {
            return false;
        }
    }

    // intersection is non-empty
    for (short i = start; i <= end; i ++)
        region.push_back(i);

#ifdef _PRINTING
cout << "intersect: " << endl;
for (vector<short>::iterator it = region.begin(); it != region.end(); it ++)
    cout << *it << ", ";
cout << endl;
#endif

    return true;
}

void Read::candidate_union(vector<int> candidate_subset, vector<short>& region) {
    short u;
    set<short> region_set;
    for (int i = 0; i < candidate_subset.size(); i ++) {
        u = candidate_subset[i];

        for (short ui = u; ui < u+hashtab::k; ui++)
            region_set.insert(ui);
    }
    // for (int i = 0; i < region_set.size(); i++) {
    //     region.push_back(region_set[i]);
    // }
    for (set<short>::iterator it = region_set.begin(); it != region_set.end(); it ++)
        region.push_back(*it);
}

double Read::get_coverage(vector<short> region, hashtab* all) {
    int spos = region.front();
    int epos = region.back();

    int kmer_start = max(0, spos-hashtab::k+1);
    int kmer_end   = min(read_length-hashtab::k, epos);

    double cov = 0;

    for (int i = kmer_start; i <= kmer_end; i ++) {
        cov += all->getcount(&seq[i]);
#ifdef _PRINTING
        cout << "mer index: " << i << ", coverage: " << all->getcount(&seq[i]) << endl;
#endif
    }
    return cov / (kmer_end - kmer_start + 1);
}

