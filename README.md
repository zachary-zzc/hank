# hank
Hank -- Consensus sequencing error and heterozygous in fastq files

## Kmer count

k-mer frequencies can be fast counted by jellyfish (http://www.cbcb.umd.edu/software/jellyfish/) by the following command lines:

``
jellyfish count -c ct_size -o output.db -m kmer_size -t threads -s hash_size --both-strands fastq_files
``

``
jellyfish dump -c output.db > output.cts
``
