#!/bin/bash
# -*- coding: utf-8 -*-

# Code: extract BRAF gene from bam
# Baihan Lin
# Columbia University
# May 2018

grep -o 's3.*bam' case_s3_bam_loc.tumor.ref.txt > s3_tumor.txt
grep -o 's3.*bam' case_s3_bam_loc.normal.ref.txt > s3_normal.txt

ssh -i /Users/DoerLBH/Dropbox/git/doerlbh-rabadan.pem ubuntu@ec2-54-152-113-174.compute-1.amazonaws.com

# in AWS instance

mkdir tumor
mkdir normal

cd tumor

for f in `cat s3_tumor.txt`
do
	exists=$(aws s3 ls $f)
	if [ -z "$exists" ]; then
		echo $f >> no_exist.txt
	else
		echo $f >> exist.txt
	fi

	aws s3 cp $f . 
	aws s3 cp $f.bai . 
    samtools view *.bam "chr7:140,719,327-140,924,764" | cut -f 10 >> tumor-BRAF.txt
    samtools view *.bam "chr7:140,719,327-140,924,764" >> tumor-BRAF-full.txt
    rm *.bam*
done

for f in `cat s3_tumor_val.txt`
do
	exists=$(aws s3 ls $f)
	if [ -z "$exists" ]; then
		echo $f >> no_exist.txt
	else
		echo $f >> exist.txt
	fi

	aws s3 cp $f . 
	aws s3 cp $f.bai . 
    samtools view *.bam "chr7:140,719,327-140,924,764" | cut -f 10 >> tumor-BRAF-val.txt
    samtools view *.bam "chr7:140,719,327-140,924,764" >> tumor-BRAF-val-full.txt
    rm *.bam*
done

cd ..

cd normal

for f in `cat s3_normal.txt`
do
	exists=$(aws s3 ls $f)
	if [ -z "$exists" ]; then
		echo $f >> no_exist.txt
	else
		echo $f >> exist.txt
	fi

	aws s3 cp $f .
	aws s3 cp $f.bai . 
    samtools view *.bam "chr7:140,719,327-140,924,764" | cut -f 10 >> normal-BRAF.txt
    samtools view *.bam "chr7:140,719,327-140,924,764" >> normal-BRAF-full.txt
    rm *.bam*
done

for f in `cat s3_normal_val.txt`
do
	exists=$(aws s3 ls $f)
	if [ -z "$exists" ]; then
		echo $f >> no_exist.txt
	else
		echo $f >> exist.txt
	fi

	aws s3 cp $f .
	aws s3 cp $f.bai . 
    samtools view *.bam "chr7:140,719,327-140,924,764" | cut -f 10 >> normal-BRAF-val.txt
    samtools view *.bam "chr7:140,719,327-140,924,764" >> normal-BRAF-val-full.txt
    rm *.bam*
done

cd ..



cd data

sed -e 's/N/0, /g' \
    -e 's/A/1, /g' \
    -e 's/C/2, /g' \
    -e 's/T/3, /g' \
    -e 's/G/4, /g' test-tumor-BRAF.txt > NACTG-test-tumor-BRAF.txt

sed -e 's/N/0, /g' \
    -e 's/A/1, /g' \
    -e 's/C/2, /g' \
    -e 's/T/3, /g' \
    -e 's/G/4, /g' test-normal-BRAF.txt > NACTG-test-normal-BRAF.txt

sed -e 's/N/0, /g' \
    -e 's/A/1, /g' \
    -e 's/C/2, /g' \
    -e 's/T/3, /g' \
    -e 's/G/4, /g' train-normal-BRAF.txt > NACTG-train-normal-BRAF.txt

sed -e 's/N/0, /g' \
    -e 's/A/1, /g' \
    -e 's/C/2, /g' \
    -e 's/T/3, /g' \
    -e 's/G/4, /g' train-tumor-BRAF.txt > NACTG-train-tumor-BRAF.txt


sed "s/.*/1, &/" test-tumor-BRAF.txt > C-test-tumor-BRAF.txt
sed "s/.*/1, &/" train-tumor-BRAF.txt > C-train-tumor-BRAF.txt
sed "s/.*/0, &/" test-normal-BRAF.txt > C-test-normal-BRAF.txt
sed "s/.*/0, &/" train-normal-BRAF.txt > C-train-normal-BRAF.txt

sed "s/.*/, &/" test-tumor-BRAF.txt > t-test-tumor-BRAF.txt
sed "s/.*/, &/" train-tumor-BRAF.txt > t-train-tumor-BRAF.txt
sed "s/.*/, &/" test-normal-BRAF.txt > t-test-normal-BRAF.txt
sed "s/.*/, &/" train-normal-BRAF.txt > t-train-normal-BRAF.txt

paste ./tumor/tumor-BRAF-loc.txt t-train-tumor-BRAF.txt > u-train-tumor-BRAF.txt
paste ./tumor/tumor-BRAF-val-loc.txt t-test-tumor-BRAF.txt > u-test-tumor-BRAF.txt
paste ./normal/normal-BRAF-loc.txt t-train-normal-BRAF.txt > u-train-normal-BRAF.txt
paste ./normal/normal-BRAF-val-loc.txt t-test-normal-BRAF.txt > u-test-normal-BRAF.txt

sed "s/.*/1, &/" u-test-tumor-BRAF.txt | tr -d " \t"  > C-test-tumor-BRAF.txt
sed "s/.*/1, &/" u-train-tumor-BRAF.txt | tr -d " \t" > C-train-tumor-BRAF.txt
sed "s/.*/0, &/" u-test-normal-BRAF.txt | tr -d " \t"  > C-test-normal-BRAF.txt
sed "s/.*/0, &/" u-train-normal-BRAF.txt | tr -d " \t" > C-train-normal-BRAF.txt

cat C-test-tumor-BRAF.txt C-test-normal-BRAF.txt > ref-test-BRAF.csv
cat C-train-tumor-BRAF.txt C-train-normal-BRAF.txt > ref-train-BRAF.csv

rm C-*
rm u-*
rm t-*
rm NACTG-*

GENOME_START = 140719327
GENOME_END = 140924764

mv ref-test-BRAF.csv all-ref-test-BRAF.csv
mv ref-train-BRAF.csv all-ref-train-BRAF.csv

grep "140719[0-2][0-9][0-9]" all-ref-test-BRAF.csv >> not-ref-test-BRAF.csv
grep "1407193[0-1][0-9]" all-ref-test-BRAF.csv >> not-ref-test-BRAF.csv
grep "14071932[0-6]" all-ref-test-BRAF.csv >> not-ref-test-BRAF.csv
grep "1409247[0-9][0-9]" all-ref-test-BRAF.csv >> not-ref-test-BRAF.csv
grep "1409246[7-9][0-9]" all-ref-test-BRAF.csv >> not-ref-test-BRAF.csv
grep "14092466[5-9]" all-ref-test-BRAF.csv >> not-ref-test-BRAF.csv

grep "140719[0-2][0-9][0-9]" all-ref-train-BRAF.csv >> not-ref-train-BRAF.csv
grep "1407193[0-1][0-9]" all-ref-train-BRAF.csv >> not-ref-train-BRAF.csv
grep "14071932[0-6]" all-ref-train-BRAF.csv >> not-ref-train-BRAF.csv
grep "1409247[0-9][0-9]" all-ref-train-BRAF.csv >> not-ref-train-BRAF.csv
grep "1409246[7-9][0-9]" all-ref-train-BRAF.csv >> not-ref-train-BRAF.csv
grep "14092466[5-9]" all-ref-train-BRAF.csv >> not-ref-train-BRAF.csv


# diff -U $(wc -l < all-ref-test-BRAF.csv) all-ref-test-BRAF.csv not-ref-test-BRAF.csv | sed -n 's/^-//p' > ref-test-BRAF.csv
# diff -U $(wc -l < all-ref-train-BRAF.csv) all-ref-train-BRAF.csv not-ref-train-BRAF.csv | sed -n 's/^-//p' > ref-train-BRAF.csv


comm -3 <(sort all-ref-test-BRAF.csv) <(sort not-ref-test-BRAF.csv)  > ref-test-BRAF.csv
comm -3 <(sort all-ref-train-BRAF.csv) <(sort not-ref-train-BRAF.csv)  > ref-train-BRAF.csv


