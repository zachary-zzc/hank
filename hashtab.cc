#include "hashtab.h"
#include <iostream>
#include <fstream>
#include <cstdlib>

using namespace::std;

hashtab::hashtab(int _k) 
:bits( (unsigned long long int)pow(4.0, PREFIX_SIZE) ){

    if (_k < PREFIX_SIZE) {
        cerr << "k should larger than 17" << endl;
        exit(EXIT_FAILURE);
    }

    k = _k;
    mask = (unsigned long long)pow(4.0, PREFIX_SIZE) - 1;
    offsets = new loc[(unsigned long long)pow(4.0, PREFIX_SIZE)];
    for (unsigned long long i = 0; i <= mask; i++) {
        offsets[i].count = 0; 
        offsets[i].cur = 0;
        offsets[i].offset = 0;
    }
}

    
hashtab::~hashtab() {
    delete[] counts;
    delete[] surfix_mers;
    delete[] offsets;
}

void hashtab::add(unsigned long long pre_mer, unsigned long long sur_mer, double count) {

    bool update = false;
    // for (int i = 0; i < offsets[pre_mer].count; i++) {
    //     if (counts[offsets[pre_mer].offset+i] <= 0) {
    //         surfix_mers[offsets[pre_mer].offset+i] = sur_mer;
    //         counts[offsets[pre_mer].offset+i] = count;
    //         update = true;
    //         break;
    //     }
    // }
    if (offsets[pre_mer].cur <= offsets[pre_mer].count && counts[offsets[pre_mer].offset+offsets[pre_mer].cur] <= 0)  {
        surfix_mers[offsets[pre_mer].offset+offsets[pre_mer].cur] = sur_mer;
        counts[offsets[pre_mer].offset+offsets[pre_mer].cur] = count;
        offsets[pre_mer].cur ++;
        update = true;
    }

    if (!update) {
        cerr << "offset: " << offsets[pre_mer].offset << endl;
        cerr << "count: " << offsets[pre_mer].count << endl;
        cerr << "cur" << offsets[pre_mer].cur << endl;
        for (int i = 0; i < offsets[pre_mer].count; i++) 
            cerr << "i: " << i << ", count: " << counts[offsets[pre_mer].offset+i] << endl;
        cerr << "hash table count error" << endl;
        exit(EXIT_FAILURE);
    }
}


void hashtab::add(const string &s, double count) {
    unsigned long long pre_kmermap = 0, sur_kmermap = 0;
    for (int i = 0; i < PREFIX_SIZE; i++) {
        pre_kmermap <<= 2;
        pre_kmermap |= binary_nt(s[i]);
    }

    for (int i = PREFIX_SIZE; i < k; i++) {
        sur_kmermap <<= 2;
        sur_kmermap |= binary_nt(s[i]);
    }

    add(pre_kmermap, sur_kmermap, count);
}



void hashtab::add(unsigned kmer[], double count){
    // get prefix hash table index
    unsigned long long pre_kmermap = 0, sur_kmermap = 0;
    for (int i = 0; i < PREFIX_SIZE; i++) {
        if (kmer[i] < 4) {
            pre_kmermap <<= 2;
            pre_kmermap |= kmer[i];
        } else {
            // cerr << "Cannot handle N's in hash table" << endl;
            // cerr << "prefix: ";
            // for (int i = 0; i < PREFIX_SIZE; i++)
            //     cerr << kmer[i];
            exit(EXIT_FAILURE);
        }
    }
    
    // save surfix mer and count
    for (int i = PREFIX_SIZE; i < k; i ++) {
        if (kmer[i] < 4) {
            sur_kmermap <<= 2;
            sur_kmermap |= kmer[i];
        } else {
            // cerr << "Cannot handle N's in hash table" << endl;
            // cerr << "surfix: ";
            // for (int i = PREFIX_SIZE; i < k; i++)
            //     cerr << kmer[i];
            exit(EXIT_FAILURE);
        }
    }

    add(pre_kmermap, sur_kmermap, count);
}



double hashtab::getcount(unsigned kmer[]) {
    unsigned long long pre_kmermap = 0, sur_kmermap = 0, idx = -1;
    for (int i = 0; i < PREFIX_SIZE; i++) {
        if (kmer[i] < 4) {
            pre_kmermap <<= 2;
            pre_kmermap |= kmer[i];
        } else {
            // cerr << "Cannot handle N's in reads" << endl;
            // cerr << "prefix: ";
            // for (int i = 0; i < PREFIX_SIZE; i++)
            //     cerr << kmer[i];
            return 0;
        }
    }

    for (int i = PREFIX_SIZE; i < k; i++) {
        if (kmer[i] < 4) {
            sur_kmermap <<= 2;
            sur_kmermap |= kmer[i];
        } else {
            // cerr << "Cannot handle N's in reads" << endl;
            // cerr << "surfix: ";
            // for (int i = PREFIX_SIZE; i < k; i++)
            //     cerr << kmer[i];
            return 0;
        }
    }

    for (int i = 0; i < offsets[pre_kmermap].count; i++) {
        if (surfix_mers[offsets[pre_kmermap].offset+i] == sur_kmermap) {
            idx = offsets[pre_kmermap].offset + i;
            break;
        }
    }
    if (idx == -1) 
        return 0;
    else 
        return counts[idx];
}


double hashtab::getcount(string &s) {
    unsigned* kmer = new unsigned[k];
    for (int i = 0; i < s.length(); i++) 
        kmer[i] = binary_nt(s[i]);
    double count = getcount(kmer);
    delete[] kmer;
    return count;
}
    

void hashtab::file_load(istream &mer_in, unsigned long long atgc[]) {
    string line;
    double count;
    unsigned long long offset = 0;
    
    while (getline(mer_in, line)) {
        if(line[k] != ' ' && line[k] != '\t') {
            cerr << "Kmers are not of expected length: " << k << endl;
            exit(EXIT_FAILURE);
        }
        offsets[binary_kmer(line.substr(0, PREFIX_SIZE))].count += 1;
        offsets[binary_kmer(rckmer(line.substr(0, k)).substr(0, PREFIX_SIZE))].count += 1;
    }
    
    for (unsigned long long i = 0; i <= mask; i++) {
        offsets[i].offset = offset;
        offset += offsets[i].count;
    }

    number = offset;
    surfix_mers = new unsigned long long[number];
    counts = new double[number];
    for (unsigned long long i = 0; i < number; i++) 
        counts[i] = 0;


    // relocate stream to the begining
    mer_in.clear();
    mer_in.seekg(0, mer_in.beg);
    while (getline(mer_in, line)) {
        count = atof(line.substr(k+1).c_str());
       
        add(line.substr(0, k), count);
        add(rckmer(line.substr(0, k)), count);

        if (atgc != NULL) {
            unsigned int at = count_at(line.substr(0, k));
            atgc[0] += at;
            atgc[1] += (k - at);
        }
    }
}


////////////////////////////////////////////////////////////
// binary_file_output
//
// Write hashtab to file in binary format
////////////////////////////////////////////////////////////
void hashtab::binary_file_output(char* outf) {
    unsigned long long mysize = (unsigned long long)bits.size() / 8ULL;
    char* buffer = new char[mysize];
    unsigned int flag = 1;
    for(unsigned long long i = 0; i < mysize; i++) {
        unsigned int temp = 0;
        for(unsigned int j = 0; j < 8; j++) { // read 8 bits from the bitset
            temp <<= 1;
            //unsigned int tmp = i*8 + j;
            //cout << tmp << ",";
            if(bits[i*8 + j])
                temp |= flag;
        }
        buffer[i] = (char)temp;
    }
    ofstream ofs(outf, ios::out | ios::binary);
    ofs.write(buffer, mysize);
    delete[] buffer;
    ofs.close();
}

////////////////////////////////////////////////////////////
// binary_file_input
//
// Read hashtab from file in binary format
////////////////////////////////////////////////////////////
/*
   void hashtab::binary_file_input(char* inf) {
   ifstream ifs(inf, ios::binary);

// get size of file
ifs.seekg(0,ifstream::end);
unsigned long long mysize = ifs.tellg();
ifs.seekg(0);

// allocate memory for file content
char* buffer = new char[mysize];

// read content of ifs
ifs.read (buffer, mysize);

// parse bits
unsigned int flag = 128;
unsigned int temp;
for(unsigned long i = 0; i < mysize; i++) {
temp = (unsigned int)buffer[i];
for(unsigned int j = 0; j < 8; j++) {
if((temp & flag) == flag)
bits.set(i*8 + j);
temp <<= 1;
}
}

delete[] buffer;
}
*/

////////////////////////////////////////////////////////////
// binary_file_input
//
// Read hashtab from file in binary format
////////////////////////////////////////////////////////////
void hashtab::binary_file_input(char* inf, unsigned long long atgc[]) {
    unsigned int flag = 128;
    unsigned int temp;

    ifstream ifs(inf, ios::binary);

    // get size of file
    ifs.seekg(0,ifstream::end);
    unsigned long long mysize = ifs.tellg();
    ifs.seekg(0);

    // allocate memory for file content
    unsigned long long buffersize = 134217728;  // i.e. 4^15 / 8, 16 MB
    if(mysize < buffersize)
        buffersize = mysize;
    char* buffer = new char[buffersize];

    for(unsigned long long b = 0; b < mysize/buffersize; b++) {

        // read content of ifs
        ifs.read (buffer, buffersize);

        // parse bits
        for(unsigned long long i = 0; i < buffersize; i++) {
            temp = (unsigned int)buffer[i];
            for(int j = 0; j < 8; j++) {
                if((temp & flag) == flag) {
                    bits.set((buffersize*b + i)*8 + j);

                    // count gc
                    unsigned int at = count_at((buffersize*b + i)*8 + j);
                    atgc[0] += at;
                    atgc[1] += (k-at);
                }
                temp <<= 1;
            }
        }
    }

    delete[] buffer;
}

////////////////////////////////////////////////////////////
// count_at
//
// Count the A's and T's in the sequence given
////////////////////////////////////////////////////////////
int hashtab::count_at(string seq) {
    int at = 0;
    for(int i = 0; i < seq.size(); i++)
        if(seq[i] == 'A' || seq[i] == 'T')
            at +=  1;
    return at;
}

int hashtab::count_at(unsigned long long seq) {
    int at = 0;
    unsigned long long mask = 3;
    unsigned long long nt;
    for(int i = 0; i < k; i++) {
        nt = seq & mask;
        seq >>= 2;

        if(nt == 0 || nt == 3)
            at++;
    }
    return at;
}


void hashtab::iseq_kmer(const string &s, unsigned kmer[]) {
    for (int i = 0; i < s.length(); i++) 
        kmer[i] = binary_nt(s[i]);
}

//  Convert string  s  to its binary equivalent in  mer .
unsigned long long  hashtab::binary_kmer(const string & s) {
    int  i;
    unsigned long long mer = 0;
    for  (i = 0; i < s.length(); i++) {
        mer <<= 2;
        mer |= binary_nt(s[i]);
    }
    return mer;
}


//  Convert string s to its binary equivalent in mer .
unsigned long long  hashtab::binary_rckmer(const string & s) {
    int  i;
    unsigned long long mer = 0;
    for  (i = s.length()-1; i >= 0; i--) {
        mer <<= 2;
        mer |= 3 - binary_nt(s[i]);
    }
    return mer;
}


//  Return the binary equivalent of  ch .
unsigned hashtab::binary_nt(char ch) {
    switch  (tolower (ch)) {
        case  'a' : return  0;
        case  'c' : return  1;
        case  'g' : return  2;
        case  't' : return  3;
    }
}


string hashtab::rckmer(const string & s) {
    string ss;
    ss.assign(s.rbegin(), s.rend());
    for (int i = 0; i < ss.size(); i++) {
        ss[i] = rc_nt(ss[i]);
    }
    return ss;
}


char hashtab::rc_nt(char ch) {
    switch (tolower(ch)) {
        case 'a' : return 't';
        case 'c' : return 'g';
        case 'g' : return 'c';
        case 't' : return 'a';
    }
}

unsigned int hashtab::num_kmers() {
    return number;
}
