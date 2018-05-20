#ifndef hashtab_H
#define hashtab_H

#include <string>
#include <iostream>
#include <vector>
#include <cmath>
#include <numeric>
#include <boost/dynamic_bitset.hpp>

#define PREFIX_SIZE         15

using namespace::std;

struct loc{
    unsigned long long offset;
    unsigned count;
    unsigned cur;
};

class hashtab {
    public:
        hashtab(int _k);
        ~hashtab();

        void add(unsigned kmer[], double count);
        void add(const string &s, double count);
        void add(unsigned long long pre_mer, unsigned long long sur_mer, double count);

        void file_load(istream &mer_in, unsigned long long atgc[]);
        long long unsigned binary_kmer(const string &s);
        long long unsigned binary_rckmer(const string &s);
        string rckmer(const string & s);
        char rc_nt(char ch);
        void iseq_kmer(const string &s, unsigned kmer[]);
        void binary_file_output(char* outf);
        void binary_file_input(char* inf, unsigned long long atgc[]);

        double getcount(unsigned kmer[]);
        double getcount(string &s);
        unsigned int num_kmers();

        static int k;

    private:
        unsigned binary_nt(char ch);
        int count_at(string seq);
        int count_at(unsigned long long seq);

        boost::dynamic_bitset<> bits;
        loc *offsets;
        unsigned long long* surfix_mers;
        double* counts;
        unsigned long long mask;
        unsigned long long number;
};    
    
       
#endif
