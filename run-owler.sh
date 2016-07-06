#! /bin/sh

make -j 4 testing

mkdir -p temp
bin/graphmap-not_release owler2 -t 1 -r test-data/lambda/reads.fastq -d test-data/lambda/reads.fastq -o temp/lambda.mhap
