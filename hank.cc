#include "hashtab.h"
#include "Read.h"
#include "edit.h"
#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>
#include <map>
#include <string>
#include <string.h>
#include <cstring>
#include <getopt.h>
#include <omp.h>
#include <cstdlib>
#include <iomanip>
#include <sys/stat.h>
#include "gzstream.h"

// #define TESTING false
// #define _TEST true
// #define _PRINTING true

using namespace std;

////////////////////////////////////////////////////////////
// options
////////////////////////////////////////////////////////////

const static char* myopts = "c:s:e:f:k:m:p:r:zuh";
static struct option long_options[] = {
    {"headers", 0, 0, 1000},
    {"log", 0, 0, 1001},
    {0, 0, 0, 0}
};

struct block {
    unsigned int from;
    unsigned int to;
    double mean_cov;
    short type;
    block(int f, int t)
        :from(f), to(t) {
        };
    ~block(){
    };
};

// -c, homozygous coverage
static double homo_coverage = 0;
// -s, heterozygous coverage
static double hete_coverage = 0;
// -e, trusted/untrusted k-mer cutoff
static double cutoff = 0;
// -f, file of fastq files of reads
// static char* fastqf = NULL;
// -k, kmer size
static int k = 0;
// -m, mer counts
static char* merf = NULL;
// -p, process number
// static int threads = 4;
// -z, zip output files
// static bool zip_output = false;
// -u, output unmodified reads
static bool unmodified_out = false;
// -r, percentile
static double percentile = 1.96;
// --headers, print only normal headers
static bool orig_headers = false;
// --log, output modification log
static bool out_log = false;

// parameters
static bool overwrite_temp = true;
static char* nts = "ACGTN";

struct stats{
    stats() {
        validated    = 0;
        modified     = 0;
        removed      = 0;
        trimmed      = 0;
        trimmed_only = 0;
    }
    unsigned long long validated;
    unsigned long long modified;
    unsigned long long removed;
    unsigned long long trimmed;
    unsigned long long trimmed_only;
};

static void Usage(char* command) {
    /* print to stderr */
    fprintf(stderr, 
            "USAGE: hdyra [options]\n"
            "\n"
            "Consensus sequencing error and heterozygous in fastq file\n"
            "\n"
            "Options: \n"
            " -f <file>     --required\n"
            "   File containing fastq file names, one per line or\n"
            "   two per line for paired end reads.\n"
            " -c <num>      --required\n"
            "   Homozygous coverage\n"
            " -s <num>      --required\n"
            "   Heterozygous coverage\n"
            " -e <num>      --required\n"
            "   Trusted/Untrusted kmer cutoff\n"
            " -k <num>      --required\n"
            "   K-mer size\n"
            " -m <file>     --required\n"
            "   File containing k-mer counts in format `seq\tcount`.\n"
            " -p <num>      <default: 4>\n"
            "   Threads number\n"
            " -r <num>      <default: 1.96>\n"
            "   Repeat scaling threshold, float between 1 and 1.96.\n"
            " -z            <default: false>\n"
            "   Write output files as gzipped.\n"
            " -u            <default: false>\n"
            "   Output unmodified reads.\n"
            " --headers     <default: false>\n"
            "   Output only the original read headers\n"
            " --log         <default: false>\n"
            "   Output a log of all modifications into *.log\n"
            "\n");
    return;
}


////////////////////////////////////////////////////////////
// parse command line
////////////////////////////////////////////////////////////
static void parse_command_line(int argc, char** argv) {
    bool errflg = false;
    int ch;
    optarg = NULL;
    int option_index = 0;
    char* p;

    // parse args
    while (!errflg && ((ch = getopt_long(argc, argv, myopts, long_options, &option_index)) != EOF)) {
        switch(ch) {
            case 'f':
                file_of_fastqf = strdup(optarg);
                break;

            case 'c':
                homo_coverage = double(strtod(optarg, &p));
                if (p == optarg || homo_coverage <= 0) {
                    fprintf(stderr, "Bad homozygous coverage value \"%s\"\n", optarg);
                    errflg = true;
                }
                break;

            case 's':
                hete_coverage = double(strtod(optarg, &p));
                if (p == optarg || hete_coverage <= 0) {
                    fprintf(stderr, "Bad heterozygous coverage value \"%s\"\n", optarg);
                    errflg = true;
                }
                break;

            case 'e':
                cutoff = double(strtod(optarg, &p));
                if (p == optarg || cutoff <= 0) {
                    fprintf(stderr, "Bad trusted/untrusted k-mer cutoff value \"%s\"\n", optarg);
                    errflg = true;
                }
                break;

            case 'k':
                k = int(strtod(optarg, &p));
                if (p == optarg || k <= 2) {
                    fprintf(stderr, "Bad kmer size \"%s\"\n", optarg);
                    errflg = true;
                }
                break;

            case 'r':
                percentile = double(strtod(optarg, &p));
                if (p == optarg || percentile <= 1) {
                    fprintf(stderr, "Bad scaling, should be larger than 1 and smaller than 1.96");
                    errflg = true;
                }
                break;

            case 'm':
                merf = strdup(optarg);
                break;

            case 'p':
                threads = int(strtol(optarg, &p, 10));
                if (p == optarg || threads <= 0) {
                    fprintf(stderr, "Bad number of threads \"%s\"\n", optarg);
                    errflg = true;
                }
                break;

            case 'z':
                zip_output = true;
                break;

            case 1000:
                orig_headers = true;
                break;

            case 1001:
                out_log = true;
                break;

            case 'h':
                Usage(argv[0]);
                exit(EXIT_FAILURE);

            case '?':
                fprintf (stderr, "Unrecognized option -%c\n", optopt);

            default:
                errflg = true;
        }
    }

    if (file_of_fastqf == NULL) {
        cerr << "Must provide a file containing a list of fastq files of reads (-f)" << endl;
        exit(EXIT_FAILURE);
    }

    if (homo_coverage == 0) {
        cerr << "Must provide homozygous coverage (-c)" << endl;
        exit(EXIT_FAILURE);
    }

    if (hete_coverage == 0) {
        cerr << "Must provide heterozygous coverage (-s)" << endl;
        exit(EXIT_FAILURE);
    }

    if (cutoff == 0) {
        cerr << "Must provide trusted/untrusted k-mer cutoff (-e)" << endl;
        exit(EXIT_FAILURE);
    }

    if (k == 0) {
        cerr << "Must provide k-mer size (-k)" << endl;
        exit(EXIT_FAILURE);
    }

    if (merf == NULL) {
        cerr << "Must provide a file of kmer counts (-m)" << endl;
        exit(EXIT_FAILURE);
    }
}


////////////////////////////////////////////////////////////////////////////////
// output_read
//
// Output the given possibly corrected and/or trimmed
// read according to the given options.
////////////////////////////////////////////////////////////////////////////////
static void output_read(ofstream &reads_out, ofstream &corlog_out, int pe_code, string header, string ntseq, string mid, string strqual, string corseq, stats &tstats) {
    
    bool modified = false;

    for (int i = 0; i < corseq.size(); i ++) {
        if (corseq[i] != ntseq[i]) {
            if (corlog_out.good())
                corlog_out << (i+1) << "\t" << corseq[i] << "\t" << ntseq[i] << endl;
            modified = true;
        }
    }
    if (modified)
        tstats.modified ++;

    // update header
    if (!orig_headers) {
        if (modified)
            header += " modified";
        // for trim, add later
        int trimlen = ntseq.size() - corseq.size();
        // if (trimlen == 0 && !modified)
    }
    
    if (!modified) 
        tstats.validated ++;
    
    reads_out << header << endl << corseq << endl << mid << endl << strqual.substr(0, corseq.size()) << endl;
}

// static void classify_kmers(hashtab* all, vector<unsigned int> iseq, double homo_coverage, double hete_coverage, double cutoff, vector<int>& candidate, string header) {
//     candidate.clear();
// 
// #ifdef _PRINTING
//     cout << "get candidate regions" << endl;
// #endif
// 
//     vector<double> icov;
// 
// 
// }

static void classify_kmers(hashtab* all, vector<unsigned int> iseq, double homo_coverage, double hete_coverage, double cutoff, double percentile, vector<int>& candidate, string header) {
    candidate.clear();
    
#ifdef _PRINTING
    cout << "get candidate regions" << endl;
#endif

    // vector<double> icov;
    map<int, double> icov;
    for (int i = 0; i < iseq.size()-k+1; i ++) {
#ifdef _PRINTING
        cout << "kmer index: " << i << ", coverage: " << all->getcount(&iseq[i]) << endl;
#endif
        double cov = all->getcount(&iseq[i]);
        if (cov <= percentile * homo_coverage)
            icov[i] = cov;
        // icov.push_back(all->getcount(&iseq[i]));
    }

    // search locations with significant drop
    vector<int> down_breakpoints;
    vector<int> up_breakpoints;

    for (int idx = 0; idx < icov.size()-1; idx ++) {
    // for (int idx = 0; idx < icov.size()-1; idx ++) {
        if ((icov[idx] - icov[idx+1] >= icov[idx] / 3) && \
            (icov[idx] - icov[idx+1] >= (homo_coverage - hete_coverage) / 3) && \
            (icov[idx] - icov[idx+1] >= cutoff)) {
            if ((icov[idx] <= percentile * homo_coverage) && (icov[idx+1] <= percentile * hete_coverage))
                down_breakpoints.push_back(idx);
            else
                icov[idx+1] = icov[idx];
#ifdef _PRINTING
cout << "kmer cov: " << endl;
for (int i = 0; i < icov.size(); i ++) {
    cout << icov[i] << ", ";
}
cout << endl;
#endif
        }
        if ((icov[idx+1] - icov[idx] >= icov[idx+1] / 3) && \
            (icov[idx+1] - icov[idx] >= (homo_coverage - hete_coverage) / 3) && \
            (icov[idx+1] - icov[idx] >= cutoff)) {
            if ((icov[idx+1] <= percentile * homo_coverage) && (icov[idx] <= percentile * hete_coverage))
                up_breakpoints.push_back(idx);
            else
                icov[idx+1] = icov[idx];
#ifdef _PRINTING
cout << "kmer cov: " << endl;
for (int i = 0; i < icov.size(); i ++) {
    cout << icov[i] << ", ";
}
cout << endl;
#endif
        }
    }

    if (down_breakpoints.size() != 0 && up_breakpoints.size() == 0) {
        up_breakpoints.push_back(icov.size() - 1);
    } else if (up_breakpoints.size() != 0 && down_breakpoints.size() == 0) {
        down_breakpoints.push_back(-1);
    } else if (up_breakpoints.size() != 0 && down_breakpoints.size() != 0) {
        if (down_breakpoints.front() > up_breakpoints.front()) {
            down_breakpoints.insert(down_breakpoints.begin(), -1);
        }
        if (down_breakpoints.back() > up_breakpoints.back()) {
            up_breakpoints.push_back(icov.size() - 1);
        }
    }

#ifdef _PRINTING
cout << "up breakpoints: " << endl;
for (vector<int>::iterator it = up_breakpoints.begin(); it != up_breakpoints.end(); it ++) {
    cout << *it << ", ";
}
cout << endl;
cout << "down breakpoints: " << endl;
for (vector<int>::iterator it = down_breakpoints.begin(); it != down_breakpoints.end(); it ++) {
    cout << *it << ", ";
}
cout << endl;

cout << "icov size: " << icov.size() << endl;
#endif

    int i = 0, j = 0;
    while (i < down_breakpoints.size() && j < up_breakpoints.size()) {
        int down = down_breakpoints[i];
        int up   = up_breakpoints[j];

#ifdef _PRINTING
cout << "down: " << down << endl;
cout << "up: " << up << endl;
#endif
        
        assert(down != up);

        if (up > down) {
            if ((up - down >= k || down == -1 || up == icov.size() - 1) && (up - down < 2 * k)) {
#ifdef _PRINTING
cout << "add to candidate" << endl;
cout << "up: " << up << ", down: " << down << endl;
#endif
                bool repeat = false;
                for (int c = down + 1; c != up + 1; c ++) {
                    if (icov[c] >= hete_coverage * percentile) {
                        repeat = true;
                        break;
                    }
                }
                if (repeat == false) {
                    for (int c = down + 1; c != up + 1; c ++)
                        candidate.push_back(c);
                }
            }
            i += 1;
            j += 1;
        } else {
            j += 1;
        }
        
        // if ((up-down >= k || down == -1 || up == icov.size()-1) && (up - down < 2 * k)) {
//         if (up-down >= k || down == -1 || up == icov.size()-1) {
// #ifdef _PRINTING
// cout << "add to candidate" << endl;
// cout << "up: " << up << ", down: " << down << endl;
// #endif
//             // if (up-down >= k) {
//             for (int c = down+1; c != up+1; c++)
//                 candidate.push_back(c);
//             // }
//             i += 1;
//             j += 1;
//         } else if (down > up)
//             j += 1;
//         else
//             i += 1;
//             j += 1;
    }

#ifdef _PRINTING
cout << "candidates: " << endl;
for (vector<int>::iterator it = candidate.begin(); it != candidate.end(); it ++) 
    cout << *it << ", ";
cout << endl;
#endif
} 


static void consensus_reads(string fqf, int pe_code, hashtab* all, vector<streampos>& starts, vector<unsigned long long>& counts) {
    struct stat info;
    string path_suffix = split(fqf, '/').back();
    string out_dir("."+path_suffix);
    if (stat(out_dir.c_str(), &info) == 0) {
        cerr << "Hidden temporary directory " << out_dir << "already exists and will be used" << endl;
    } else {
        if (mkdir(out_dir.c_str(), S_IRWXU) == -1) {
            cerr << "Failed to create hidden temporary directory " << out_dir << endl;
            exit(EXIT_FAILURE);
        }
    }
    stats * thread_stats = new stats[omp_get_max_threads()];

    unsigned int chunk = 0;
#pragma omp parallel //share(all)
    {
        int tid = omp_get_thread_num();

        ifstream reads_in(fqf.c_str());

        unsigned int tchunk;
        string header, ntseq, mid, strqual, corseq;
        // int trim_length;
        char* nti;
        Read* r;

#pragma omp critical
        tchunk = chunk ++;

        while (tchunk < starts.size()) {
            reads_in.seekg(starts[tchunk]);

            // output
            string toutf(out_dir+"/");
            stringstream tconvert;
            tconvert << tchunk;
            toutf += tconvert.str();

            if (overwrite_temp || stat(toutf.c_str(), &info) == -1) {
                ofstream reads_out(toutf.c_str());

                string tlogf = toutf + ".log";
                ofstream corlog_out;
                if (out_log) 
                    corlog_out.open(tlogf.c_str());

                unsigned long long tcount = 0;
                while (getline(reads_in, header)) {
                    // get sequence
                    getline(reads_in, ntseq);

                    // convert ntseq to stseq
                    vector<unsigned int> iseq;
                    for (int i = 0; i < ntseq.size(); i ++) {
                        nti = strchr(nts, ntseq[i]);
                        iseq.push_back(nti - nts);
                    }
                    
                    // get mid
                    getline(reads_in, mid);
                    // get quality
                    getline(reads_in, strqual);

                    vector<int> candidate;

                    // cout << "read name: " << header << endl;
                    // classify k-mers
#ifdef _PRINTING
                    cout << "read id: " << header << endl;
#endif
                    classify_kmers(all, iseq, homo_coverage, hete_coverage, cutoff, percentile, candidate, header);
                    // very important!!!

                    if (candidate.size() > 0) {
                        // cout << "candidate size: " << candidate.size() << endl;
                        r = new Read(header, &iseq[0], strqual, candidate, iseq.size());
                        // cout << "init new read finished" << endl;
                        corseq = r->modify(all, homo_coverage, hete_coverage, cutoff, percentile);
                        strqual = r->print_qual();
                        output_read(reads_out, corlog_out, pe_code, header, ntseq, mid, strqual, corseq, thread_stats[tid]);
                        delete r;
                    } else {
                        output_read(reads_out, corlog_out, pe_code, header, ntseq, mid, strqual, ntseq, thread_stats[tid]);
                    }

                    if (++tcount == counts[tchunk])
                        break;
                }
                reads_out.close();
            }

#pragma omp critical
            tchunk = chunk++;
        }
        reads_in.close();
    }
    for (int i = 1; i < omp_get_max_threads(); i++) {
        thread_stats[0].validated += thread_stats[i].validated;
        thread_stats[0].modified += thread_stats[i].modified;
    }

    // print stats
    int suffix_index = fqf.rfind(".");
    string outf;
    if (suffix_index == -1) {
        outf = fqf + ".stats.txt";
    } else {
        outf = fqf.substr(0, suffix_index+1) + "stats.txt";
    }
    ofstream stats_out(outf.c_str());
    stats_out << "Validated: " << thread_stats[0].validated << endl;
    stats_out << "Modified: " << thread_stats[0].modified << endl;
    stats_out.close();
}


////////////////////////////////////////////////////////////
// main
////////////////////////////////////////////////////////////
int main(int argc, char** argv) {
    parse_command_line(argc, argv);

    unsigned long long atgc[2] = {0};  // GC counts
    
    hashtab *all = new hashtab(k);
    cout << "homozygous coverage: " << homo_coverage << endl;
    cout << "heterozygous coverage: " << hete_coverage << endl;
    cout << "threads: " << threads << endl;

    string merf_str(merf);
    if (merf_str.substr(merf_str.size()-3) == ".gz") {
        igzstream mer_in(merf);
        all->file_load(mer_in, atgc);
    } else {
        ifstream mer_in(merf);
        all->file_load(mer_in, atgc);
    }

    cout << "get kmer count table finished" << endl;

    vector<string> fastqfs;
    vector<int> pairedend_codes;
    parse_fastq(fastqfs, pairedend_codes);

    string fqf; 
    bool zip;
    // for (vector<string>::iterator it = fastqfs.begin(); it != fastqfs.end(); it ++) {
    for (int f = 0; f != fastqfs.size(); f ++) {
        fqf = fastqfs[f];
        cout << fqf << endl;
        
        // unzipped
        if (fqf.substr(fqf.size()-3) == ".gz") {
            zip = true;
            unzip_fastq(fqf);
        } else
            zip = false;

        // split file
        vector<streampos> starts;
        vector<unsigned long long> counts;
        chunkify_fastq(fqf, starts, counts);

        consensus_reads(fqf, pairedend_codes[f], all, starts, counts);

        // combine
        if (pairedend_codes[f] == 0) {
            combine_output(fqf, string("mor"), unmodified_out);
        }

        if (pairedend_codes[f] == 2) {
            if (!zip)
                combine_output_paired(fastqfs[f-1], fqf, string("mor"), unmodified_out);
            else
                combine_output_paired(fastqfs[f-1].substr(0, fastqfs[f-1].size()-3), fqf, string("mor"), unmodified_out);
        }

        if (zip)
            zip_fastq(fqf);
    }
    delete all;

    return 0;
}

