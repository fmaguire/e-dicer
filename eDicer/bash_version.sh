#!/bin/bash

#./jellyfish-2.2.4/bin/jellyfish count -m23 -s 100M -t 10 -C endo_cds.cds -o endo_cds_counts.jf
#
#./jellyfish-2.2.4/bin/jellyfish dump endo_cds_counts.jf > endo_cds_counts.fa 

# get length of query set
name=$(basename $1)

# run kat comp
./KAT/src/kat comp yad1gn_endosymbiont_transcriipts.fasta $1 -o tmp/${name} -m 23 > save || echo "$1 failed" | tee -a edicer_failed_runs 



# parse output for shared kmers and number of kmers in hash2
stats=tmp/${name}.stats

# normalising factors
length_kmers=$(grep -A2 "Total K-mers in:" $stats | grep "Hash 2: " | sed 's/ - Hash 2: //')
length_unique_kmers=$(grep -A2 "Distinct K-mers in:" $stats | grep "Hash 2:" | sed 's/ - Hash 2: //') 

#collisions
total_collisions=$(grep -A2 "Shared K-mers:" $stats | grep "shared found in hash 2:" | sed 's/ - Total shared found in hash 2: //')
unique_collisions=$(grep -A3 "Shared K-mers:" $stats | grep "Distinct shared K-mers:" | sed 's/ - Distinct shared K-mers: //')

all_coll_by_length=$(bc -l <<< "$total_collisions / $length_kmers")
all_coll_by_unique_length=$(bc -l <<< "$total_collisions / $length_unique_kmers")
uniq_coll_by_length=$(bc -l <<< "$unique_collisions / $length_kmers")
uniq_coll_by_unique_length=$(bc -l <<< "$unique_collisions / $length_unique_kmers")

mkdir -p out/yadg1n/$2

echo "\"$name\"\t$length_kmers\t$length_unique_kmers\t$total_collisions\t$unique_collisions\t$all_coll_by_length\t$all_coll_by_unique_length\t$uniq_coll_by_length\t$uniq_coll_by_unique_length" > out/yadg1n/$2/${name}.summary

#> out/$2/${name}_collisions

