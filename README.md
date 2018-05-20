# hank
Hank -- Consensus sequencing error and heterozygous in fastq files

Options:
 -f <file>     --required
   File containing fastq file names, one per line or
   two per line for paired end reads.
 -c <num>      --required
   Homozygous coverage
 -s <num>      --required
   Heterozygous coverage
 -e <num>      --required
   Trusted/Untrusted kmer cutoff
 -k <num>      --required
   K-mer size
 -m <file>     --required
   File containing k-mer counts in format `seq  count`.
 -p <num>      <default: 4>
   Threads number
 -r <num>      <default: 1.96>
   Repeat scaling threshold, float between 1 and 1.96.
 -z            <default: false>
   Write output files as gzipped.
 -u            <default: false>
   Output unmodified reads.
 --headers     <default: false>
   Output only the original read headers
 --log         <default: false>
   Output a log of all modifications into *.log

## Kmer count
k-mer frequencies can be fast counted by jellyfish (http://www.cbcb.umd.edu/software/jellyfish/) by the following command lines:

jellyfish count -c ct_size -o output.db -m kmer_size -t threads -s hash_size --both-strands fastq_files
jellyfish dump -c output.db > output.cts

