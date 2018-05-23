#!/bin/bash
# -*- coding: utf-8 -*-

# Code: extract BRAF gene from bam
# Baihan Lin
# Columbia University
# May 2018

grep -o 's3.*bam' case_s3_bam_loc.tumor.ref.txt > s3_tumor.txt
grep -o 's3.*bam' case_s3_bam_loc.normal.ref.txt > s3_normal.txt

# in AWS instance

mkdir tumor
mkdir normal

cd tumor
for f in `cat s3_tumor.txt`
do
	aws s3 cp $f . --dryrun
    samtools view *.bam "Chr7:140719327-140924764" >> tumor_BRAF.txt
    rm *.bam
done
cd ..

cd normal
for f in `cat s3_tumor.txt`
do
	aws s3 cp $f . --dryrun
    samtools view *.bam "Chr7:140719327-140924764" >> normal_BRAF.txt
    rm *.bam
done
cd ..