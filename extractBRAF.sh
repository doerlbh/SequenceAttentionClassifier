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

ssh -i /Users/DoerLBH/Dropbox/git/doerlbh-rabadan.pem ubuntu@ec2-52-205-255-171.compute-1.amazonaws.com

scp -i /Users/DoerLBH/Dropbox/git/doerlbh-rabadan.pem -r ubuntu@ec2-52-205-255-171.compute-1.amazonaws.com:/home/ubuntu/attention/ .

scp -i /Users/DoerLBH/Dropbox/git/doerlbh-rabadan.pem -r models ubuntu@ec2-52-205-255-171.compute-1.amazonaws.com:/home/ubuntu/attention/
scp -i /Users/DoerLBH/Dropbox/git/doerlbh-rabadan.pem -r utils ubuntu@ec2-52-205-255-171.compute-1.amazonaws.com:/home/ubuntu/attention/
scp -i /Users/DoerLBH/Dropbox/git/doerlbh-rabadan.pem *py ubuntu@ec2-52-205-255-171.compute-1.amazonaws.com:/home/ubuntu/attention/
scp -i /Users/DoerLBH/Dropbox/git/doerlbh-rabadan.pem ./data/t*csv ubuntu@ec2-52-205-255-171.compute-1.amazonaws.com:/home/ubuntu/attention/

scp -i /Users/DoerLBH/Dropbox/git/doerlbh-rabadan.pem ubuntu@ec2-52-205-255-171.compute-1.amazonaws.com:/home/ubuntu/attention/*pt* . 

scp *.py bl2681@habanero.rcs.columbia.edu:/rigel/theory/users/bl2681/attention/

ssh -i /Users/DoerLBH/Dropbox/git/doerlbh-rabadan.pem ubuntu@ec2-52-205-255-171.compute-1.amazonaws.com
nohup python adversarial_abblstm.py > adversarial_abblstm.out 2>adversarial_abblstm.err & 
nohup python attn_bi_lstm.py > attn_bi_lstm.out 2>attn_bi_lstm.err & 
nohup python attn_lstm_hierarchical.py > attn_lstm_hierarchical.out 2>attn_lstm_hierarchical.err & 
nohup python ind_rnn_tc.py > ind_rnn_tc.out 2>ind_rnn_tc.err & 
kill 1827
kill 1828
kill 1829
kill 1830

ssh -i /Users/DoerLBH/Dropbox/git/doerlbh-rabadan.pem ubuntu@ec2-52-90-140-248.compute-1.amazonaws.com
nohup python all_attention.py > all_attention.out 2>all_attention.err & 
kill 1545

ssh -i /Users/DoerLBH/Dropbox/git/doerlbh-rabadan.pem ubuntu@ec2-18-206-236-146.compute-1.amazonaws.com
nohup python attn_bi_lstm.py > attn_bi_lstm.out 2>attn_bi_lstm.err & 
nohup python attn_lstm_hierarchical.py > attn_lstm_hierarchical.out 2>attn_lstm_hierarchical.err & 
nohup python adversarial_abblstm.py > adversarial_abblstm.out 2>adversarial_abblstm.err & 
nohup python ind_rnn_tc.py > ind_rnn_tc.out 2>ind_rnn_tc.err & 
nohup python all_attention.py > all_attention.out 2>all_attention.err & 
kill 3831
kill 2843
kill 2844
kill 2845
kill 2846



ssh -i /Users/DoerLBH/Dropbox/git/doerlbh-rabadan.pem ubuntu@ec2-18-206-236-146.compute-1.amazonaws.com

scp -i /Users/DoerLBH/Dropbox/git/doerlbh-rabadan.pem -r samtools-1.8 ubuntu@ec2-18-206-236-146.compute-1.amazonaws.com:/home/ubuntu/

scp -i /Users/DoerLBH/Dropbox/git/doerlbh-rabadan.pem -r htslib-1.8 ubuntu@ec2-18-206-236-146.compute-1.amazonaws.com:/home/ubuntu/

scp -i /Users/DoerLBH/Dropbox/git/doerlbh-rabadan.pem -r  bcftools-1.8/ ubuntu@ec2-18-206-236-146.compute-1.amazonaws.com:/home/ubuntu/

scp -i /Users/DoerLBH/Dropbox/git/doerlbh-rabadan.pem /Users/DoerLBH/Dropbox/git/AttentiveCancerExplorer/data/normal/s3*.txt  ubuntu@ec2-18-206-236-146.compute-1.amazonaws.com:/home/ubuntu/data/normal/

scp -i /Users/DoerLBH/Dropbox/git/doerlbh-rabadan.pem /Users/DoerLBH/Dropbox/git/AttentiveCancerExplorer/data/tumor/s3*.txt  ubuntu@ec2-18-206-236-146.compute-1.amazonaws.com:/home/ubuntu/data/tumor/


aws s3 ls s3://1000genomes/phase3/data/

#!/bin/bash
# -*- coding: utf-8 -*-

# Code: extract BRAF gene from bam
# Baihan Lin
# Columbia University
# May 2018

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

    samtools view *.bam "chr12:25,204,789-25,250,936" >> normal-KRAS-full.txt
    samtools view *.bam "chr9:136,494,433-136,545,786" >> normal-NOTCH1-full.txt
    samtools view *.bam "chr1:114,704,464-114,716,894" >> normal-NRAS-full.txt

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

    samtools view *.bam "chr12:25,204,789-25,250,936" >> normal-KRAS-val-full.txt
    samtools view *.bam "chr9:136,494,433-136,545,786" >> normal-NOTCH1-val-full.txt
    samtools view *.bam "chr1:114,704,464-114,716,894" >> normal-NRAS-val-full.txt

    rm *.bam*
done

cd ..

#!/bin/bash
# -*- coding: utf-8 -*-

# Code: extract BRAF gene from bam
# Baihan Lin
# Columbia University
# May 2018

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

    samtools view *.bam "chr12:25,204,789-25,250,936" >> tumor-KRAS-full.txt
    samtools view *.bam "chr9:136,494,433-136,545,786" >> tumor-NOTCH1-full.txt
    samtools view *.bam "chr1:114,704,464-114,716,894" >> tumor-NRAS-full.txt

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

    samtools view *.bam "chr12:25,204,789-25,250,936" >> tumor-KRAS-val-full.txt
    samtools view *.bam "chr9:136,494,433-136,545,786" >> tumor-NOTCH1-val-full.txt
    samtools view *.bam "chr1:114,704,464-114,716,894" >> tumor-NRAS-val-full.txt

    rm *.bam*
done

cd ..

#!/bin/bash
# -*- coding: utf-8 -*-

# Code: extract BRAF gene from bam
# Baihan Lin
# Columbia University
# May 2018

cd 1000Genomes

for h in `cat HG-train.txt`
do
    aws s3 cp s3://1000genomes/phase3/data/${h}/exome_alignment/ . --recursive --exclude "*" --include "*.mapped*.bam"
    aws s3 cp s3://1000genomes/phase3/data/${h}/exome_alignment/ . --recursive --exclude "*" --include "*.mapped*.bai"
    echo $h >> processed.txt

    samtools view *.bam "12:25,204,789-25,250,936" >> HG-KRAS-full.txt
    samtools view *.bam "9:136,494,433-136,545,786" >> HG-NOTCH1-full.txt
    samtools view *.bam "1:114,704,464-114,716,894" >> HG-NRAS-full.txt

    rm *bam*
done

for h in `cat HG-test.txt`
do
    aws s3 cp s3://1000genomes/phase3/data/${h}/exome_alignment/ . --recursive --exclude "*" --include "*.mapped*.bam"
    aws s3 cp s3://1000genomes/phase3/data/${h}/exome_alignment/ . --recursive --exclude "*" --include "*.mapped*.bai"
    echo $h >> processed.txt

    samtools view *.bam "chr12:25,204,789-25,250,936" >> HG-KRAS-val-full.txt
    samtools view *.bam "chr9:136,494,433-136,545,786" >> HG-NOTCH1-val-full.txt
    samtools view *.bam "chr1:114,704,464-114,716,894" >> HG-NRAS-val-full.txt

    rm *bam*
done

cd ..



#-----------------------------------------------------------------

cd normal 

cut -f 10 normal-KRAS-val-full.txt > normal-read-KRAS-val.txt
cut -f 4 normal-KRAS-val-full.txt > normal-locs-KRAS-val.txt
cut -f 10 normal-KRAS-full.txt > normal-read-KRAS.txt
cut -f 4 normal-KRAS-full.txt > normal-locs-KRAS.txt
cut -f 10 normal-NOTCH1-val-full.txt > normal-read-NOTCH1-val.txt
cut -f 4 normal-NOTCH1-val-full.txt > normal-locs-NOTCH1-val.txt
cut -f 10 normal-NOTCH1-full.txt > normal-read-NOTCH1.txt
cut -f 4 normal-NOTCH1-full.txt > normal-locs-NOTCH1.txt
cut -f 10 normal-NRAS-val-full.txt > normal-read-NRAS-val.txt
cut -f 4 normal-NRAS-val-full.txt > normal-locs-NRAS-val.txt
cut -f 10 normal-NRAS-full.txt > normal-read-NRAS.txt
cut -f 4 normal-NRAS-full.txt > normal-locs-NRAS.txt

sed "s/.*/, &/" normal-read-KRAS-val.txt > t-normal-read-KRAS-val.txt
sed "s/.*/, &/" normal-read-KRAS.txt > t-normal-read-KRAS.txt
sed "s/.*/, &/" normal-read-NOTCH1-val.txt > t-normal-read-NOTCH1-val.txt
sed "s/.*/, &/" normal-read-NOTCH1.txt > t-normal-read-NOTCH1.txt
sed "s/.*/, &/" normal-read-NRAS-val.txt > t-normal-read-NRAS-val.txt
sed "s/.*/, &/" normal-read-NRAS.txt > t-normal-read-NRAS.txt

paste normal-locs-KRAS-val.txt t-normal-read-KRAS-val.txt > u-normal-KRAS-val.txt
paste normal-locs-NOTCH1-val.txt t-normal-read-NOTCH1-val.txt > u-normal-NOTCH1-val.txt
paste normal-locs-NRAS-val.txt t-normal-read-NRAS-val.txt > u-normal-NRAS-val.txt
paste normal-locs-KRAS.txt t-normal-read-KRAS.txt > u-normal-KRAS.txt
paste normal-locs-NOTCH1.txt t-normal-read-NOTCH1.txt > u-normal-NOTCH1.txt
paste normal-locs-NRAS.txt t-normal-read-NRAS.txt > u-normal-NRAS.txt

sed "s/.*/0, &/"  u-normal-KRAS-val.txt > c-normal-KRAS-val.txt
sed "s/.*/0, &/"  u-normal-KRAS.txt > c-normal-KRAS.txt
sed "s/.*/0, &/"  u-normal-NOTCH1-val.txt > c-normal-NOTCH1-val.txt
sed "s/.*/0, &/"  u-normal-NOTCH1.txt > c-normal-NOTCH1.txt
sed "s/.*/0, &/"  u-normal-NRAS-val.txt > c-normal-NRAS-val.txt
sed "s/.*/0, &/"  u-normal-NRAS.txt > c-normal-NRAS.txt

grep "252509[0-9][0-9]" c-normal-KRAS-val.txt >> not-normal-KRAS-val.txt
grep "252508[4-9][0-9]" c-normal-KRAS-val.txt >> not-normal-KRAS-val.txt
grep "2525083[7-9]" c-normal-KRAS-val.txt >> not-normal-KRAS-val.txt
grep "252046[0-9][0-9]" c-normal-KRAS-val.txt >> not-normal-KRAS-val.txt
grep "252047[0-7][0-9]" c-normal-KRAS-val.txt >> not-normal-KRAS-val.txt
grep "2520478[0-8]" c-normal-KRAS-val.txt >> not-normal-KRAS-val.txt

comm -3 <(sort c-normal-KRAS-val.txt) <(sort not-normal-KRAS-val.txt)  > normal-KRAS-val.txt

grep "252509[0-9][0-9]" c-normal-KRAS.txt >> not-normal-KRAS.txt
grep "252508[4-9][0-9]" c-normal-KRAS.txt >> not-normal-KRAS.txt
grep "2525083[7-9]" c-normal-KRAS.txt >> not-normal-KRAS.txt
grep "252046[0-9][0-9]" c-normal-KRAS.txt >> not-normal-KRAS.txt
grep "252047[0-7][0-9]" c-normal-KRAS.txt >> not-normal-KRAS.txt
grep "2520478[0-8]" c-normal-KRAS.txt >> not-normal-KRAS.txt

comm -3 <(sort c-normal-KRAS.txt) <(sort not-normal-KRAS.txt)  > normal-KRAS.txt

grep "1365457[0-9][0-9]" c-normal-NOTCH1-val.txt >> not-normal-NOTCH1-val.txt
grep "13654569[0-9]" c-normal-NOTCH1-val.txt >> not-normal-NOTCH1-val.txt
grep "13654568[7-9]" c-normal-NOTCH1-val.txt >> not-normal-NOTCH1-val.txt
grep "1364943[0-9][0-9]" c-normal-NOTCH1-val.txt >> not-normal-NOTCH1-val.txt
grep "1364944[0-2][0-9]" c-normal-NOTCH1-val.txt >> not-normal-NOTCH1-val.txt
grep "13649443[0-2]" c-normal-NOTCH1-val.txt >> not-normal-NOTCH1-val.txt

comm -3 <(sort c-normal-NOTCH1-val.txt) <(sort not-normal-NOTCH1-val.txt)  > normal-NOTCH1-val.txt

grep "1365457[0-9][0-9]" c-normal-NOTCH1.txt >> not-normal-NOTCH1.txt
grep "13654569[0-9]" c-normal-NOTCH1.txt >> not-normal-NOTCH1.txt
grep "13654568[7-9]" c-normal-NOTCH1.txt >> not-normal-NOTCH1.txt
grep "1364943[0-9][0-9]" c-normal-NOTCH1.txt >> not-normal-NOTCH1.txt
grep "1364944[0-2][0-9]" c-normal-NOTCH1.txt >> not-normal-NOTCH1.txt
grep "13649443[0-2]" c-normal-NOTCH1.txt >> not-normal-NOTCH1.txt

comm -3 <(sort c-normal-NOTCH1.txt) <(sort not-normal-NOTCH1.txt)  > normal-NOTCH1.txt

grep "1147168[0-9][0-9]" c-normal-NRAS-val.txt >> not-normal-NRAS-val.txt
grep "11471679[5-9]" c-normal-NRAS-val.txt >> not-normal-NRAS-val.txt
grep "1147043[0-9][0-9]" c-normal-NRAS-val.txt >> not-normal-NRAS-val.txt
grep "1147044[0-5][0-9]" c-normal-NRAS-val.txt >> not-normal-NRAS-val.txt
grep "11470446[0-3]" c-normal-NRAS-val.txt >> not-normal-NRAS-val.txt

comm -3 <(sort c-normal-NRAS-val.txt) <(sort not-normal-NRAS-val.txt)  > normal-NRAS-val.txt

grep "1147168[0-9][0-9]" c-normal-NRAS.txt >> not-normal-NRAS.txt
grep "11471679[5-9]" c-normal-NRAS.txt >> not-normal-NRAS.txt
grep "1147043[0-9][0-9]" c-normal-NRAS.txt >> not-normal-NRAS.txt
grep "1147044[0-5][0-9]" c-normal-NRAS.txt >> not-normal-NRAS.txt
grep "11470446[0-3]" c-normal-NRAS.txt >> not-normal-NRAS.txt

comm -3 <(sort c-normal-NRAS.txt) <(sort not-normal-NRAS.txt)  > normal-NRAS.txt

rm t-*
rm u-*
rm c-*

shuf -n 1000 normal-KRAS-val.txt > normal-KRAS-test1000.txt
shuf -n 1000 normal-NOTCH1-val.txt > normal-NOTCH1-test1000.txt
shuf -n 1000 normal-NRAS-val.txt > normal-NRAS-test1000.txt

shuf -n 1400000 normal-KRAS.txt > normal-KRAS-train1400000.txt
shuf -n 2000000 normal-NOTCH1.txt > normal-NOTCH1-train2000000.txt
shuf -n 75000 normal-NRAS.txt > normal-NRAS-train75000.txt

# -------------------------- for germline mutations

sed "s/.*/, &/" normal-read-KRAS-val.txt > t-normal-read-KRAS-val.txt
sed "s/.*/, &/" normal-read-KRAS.txt > t-normal-read-KRAS.txt
sed "s/.*/, &/" normal-read-NOTCH1-val.txt > t-normal-read-NOTCH1-val.txt
sed "s/.*/, &/" normal-read-NOTCH1.txt > t-normal-read-NOTCH1.txt
sed "s/.*/, &/" normal-read-NRAS-val.txt > t-normal-read-NRAS-val.txt
sed "s/.*/, &/" normal-read-NRAS.txt > t-normal-read-NRAS.txt

paste normal-locs-KRAS-val.txt t-normal-read-KRAS-val.txt > u-normal-KRAS-val.txt
paste normal-locs-NOTCH1-val.txt t-normal-read-NOTCH1-val.txt > u-normal-NOTCH1-val.txt
paste normal-locs-NRAS-val.txt t-normal-read-NRAS-val.txt > u-normal-NRAS-val.txt
paste normal-locs-KRAS.txt t-normal-read-KRAS.txt > u-normal-KRAS.txt
paste normal-locs-NOTCH1.txt t-normal-read-NOTCH1.txt > u-normal-NOTCH1.txt
paste normal-locs-NRAS.txt t-normal-read-NRAS.txt > u-normal-NRAS.txt

sed "s/.*/1, &/"  u-normal-KRAS-val.txt > c-normal-KRAS-val.txt
sed "s/.*/1, &/"  u-normal-KRAS.txt > c-normal-KRAS.txt
sed "s/.*/1, &/"  u-normal-NOTCH1-val.txt > c-normal-NOTCH1-val.txt
sed "s/.*/1, &/"  u-normal-NOTCH1.txt > c-normal-NOTCH1.txt
sed "s/.*/1, &/"  u-normal-NRAS-val.txt > c-normal-NRAS-val.txt
sed "s/.*/1, &/"  u-normal-NRAS.txt > c-normal-NRAS.txt

rm not*

grep "252509[0-9][0-9]" c-normal-KRAS-val.txt >> not-normal-KRAS-val.txt
grep "252508[4-9][0-9]" c-normal-KRAS-val.txt >> not-normal-KRAS-val.txt
grep "2525083[7-9]" c-normal-KRAS-val.txt >> not-normal-KRAS-val.txt
grep "252046[0-9][0-9]" c-normal-KRAS-val.txt >> not-normal-KRAS-val.txt
grep "252047[0-7][0-9]" c-normal-KRAS-val.txt >> not-normal-KRAS-val.txt
grep "2520478[0-8]" c-normal-KRAS-val.txt >> not-normal-KRAS-val.txt

comm -3 <(sort c-normal-KRAS-val.txt) <(sort not-normal-KRAS-val.txt)  > normal-KRAS-val-g.txt

grep "252509[0-9][0-9]" c-normal-KRAS.txt >> not-normal-KRAS.txt
grep "252508[4-9][0-9]" c-normal-KRAS.txt >> not-normal-KRAS.txt
grep "2525083[7-9]" c-normal-KRAS.txt >> not-normal-KRAS.txt
grep "252046[0-9][0-9]" c-normal-KRAS.txt >> not-normal-KRAS.txt
grep "252047[0-7][0-9]" c-normal-KRAS.txt >> not-normal-KRAS.txt
grep "2520478[0-8]" c-normal-KRAS.txt >> not-normal-KRAS.txt

comm -3 <(sort c-normal-KRAS.txt) <(sort not-normal-KRAS.txt)  > normal-KRAS-g.txt

grep "1365457[0-9][0-9]" c-normal-NOTCH1-val.txt >> not-normal-NOTCH1-val.txt
grep "13654569[0-9]" c-normal-NOTCH1-val.txt >> not-normal-NOTCH1-val.txt
grep "13654568[7-9]" c-normal-NOTCH1-val.txt >> not-normal-NOTCH1-val.txt
grep "1364943[0-9][0-9]" c-normal-NOTCH1-val.txt >> not-normal-NOTCH1-val.txt
grep "1364944[0-2][0-9]" c-normal-NOTCH1-val.txt >> not-normal-NOTCH1-val.txt
grep "13649443[0-2]" c-normal-NOTCH1-val.txt >> not-normal-NOTCH1-val.txt

comm -3 <(sort c-normal-NOTCH1-val.txt) <(sort not-normal-NOTCH1-val.txt)  > normal-NOTCH1-val-g.txt

grep "1365457[0-9][0-9]" c-normal-NOTCH1.txt >> not-normal-NOTCH1.txt
grep "13654569[0-9]" c-normal-NOTCH1.txt >> not-normal-NOTCH1.txt
grep "13654568[7-9]" c-normal-NOTCH1.txt >> not-normal-NOTCH1.txt
grep "1364943[0-9][0-9]" c-normal-NOTCH1.txt >> not-normal-NOTCH1.txt
grep "1364944[0-2][0-9]" c-normal-NOTCH1.txt >> not-normal-NOTCH1.txt
grep "13649443[0-2]" c-normal-NOTCH1.txt >> not-normal-NOTCH1.txt

comm -3 <(sort c-normal-NOTCH1.txt) <(sort not-normal-NOTCH1.txt)  > normal-NOTCH1-g.txt

grep "1147168[0-9][0-9]" c-normal-NRAS-val.txt >> not-normal-NRAS-val.txt
grep "11471679[5-9]" c-normal-NRAS-val.txt >> not-normal-NRAS-val.txt
grep "1147043[0-9][0-9]" c-normal-NRAS-val.txt >> not-normal-NRAS-val.txt
grep "1147044[0-5][0-9]" c-normal-NRAS-val.txt >> not-normal-NRAS-val.txt
grep "11470446[0-3]" c-normal-NRAS-val.txt >> not-normal-NRAS-val.txt

comm -3 <(sort c-normal-NRAS-val.txt) <(sort not-normal-NRAS-val.txt)  > normal-NRAS-val-g.txt

grep "1147168[0-9][0-9]" c-normal-NRAS.txt >> not-normal-NRAS.txt
grep "11471679[5-9]" c-normal-NRAS.txt >> not-normal-NRAS.txt
grep "1147043[0-9][0-9]" c-normal-NRAS.txt >> not-normal-NRAS.txt
grep "1147044[0-5][0-9]" c-normal-NRAS.txt >> not-normal-NRAS.txt
grep "11470446[0-3]" c-normal-NRAS.txt >> not-normal-NRAS.txt

comm -3 <(sort c-normal-NRAS.txt) <(sort not-normal-NRAS.txt)  > normal-NRAS-g.txt

rm t-*
rm u-*
rm c-*

shuf -n 1000 normal-KRAS-val-g.txt > normal-KRAS-test1000-g.txt
shuf -n 1000 normal-NOTCH1-val-g.txt > normal-NOTCH1-test1000-g.txt
shuf -n 1000 normal-NRAS-val-g.txt > normal-NRAS-test1000-g.txt

shuf -n 1400000 normal-KRAS-g.txt > normal-KRAS-train1400000-g.txt
shuf -n 2000000 normal-NOTCH1-g.txt > normal-NOTCH1-train2000000-g.txt
shuf -n 39000 normal-NRAS-g.txt > normal-NRAS-train39000-g.txt

# --------------------------

cd ..
cd tumor

cut -f 10 tumor-KRAS-val-full.txt > tumor-read-KRAS-val.txt
cut -f 4 tumor-KRAS-val-full.txt > tumor-locs-KRAS-val.txt
cut -f 10 tumor-KRAS-full.txt > tumor-read-KRAS.txt
cut -f 4 tumor-KRAS-full.txt > tumor-locs-KRAS.txt
cut -f 10 tumor-NOTCH1-val-full.txt > tumor-read-NOTCH1-val.txt
cut -f 4 tumor-NOTCH1-val-full.txt > tumor-locs-NOTCH1-val.txt
cut -f 10 tumor-NOTCH1-full.txt > tumor-read-NOTCH1.txt
cut -f 4 tumor-NOTCH1-full.txt > tumor-locs-NOTCH1.txt
cut -f 10 tumor-NRAS-val-full.txt > tumor-read-NRAS-val.txt
cut -f 4 tumor-NRAS-val-full.txt > tumor-locs-NRAS-val.txt
cut -f 10 tumor-NRAS-full.txt > tumor-read-NRAS.txt
cut -f 4 tumor-NRAS-full.txt > tumor-locs-NRAS.txt

sed "s/.*/, &/" tumor-read-KRAS-val.txt > t-tumor-read-KRAS-val.txt
sed "s/.*/, &/" tumor-read-KRAS.txt > t-tumor-read-KRAS.txt
sed "s/.*/, &/" tumor-read-NOTCH1-val.txt > t-tumor-read-NOTCH1-val.txt
sed "s/.*/, &/" tumor-read-NOTCH1.txt > t-tumor-read-NOTCH1.txt
sed "s/.*/, &/" tumor-read-NRAS-val.txt > t-tumor-read-NRAS-val.txt
sed "s/.*/, &/" tumor-read-NRAS.txt > t-tumor-read-NRAS.txt

paste tumor-locs-KRAS-val.txt t-tumor-read-KRAS-val.txt > u-tumor-KRAS-val.txt
paste tumor-locs-NOTCH1-val.txt t-tumor-read-NOTCH1-val.txt > u-tumor-NOTCH1-val.txt
paste tumor-locs-NRAS-val.txt t-tumor-read-NRAS-val.txt > u-tumor-NRAS-val.txt
paste tumor-locs-KRAS.txt t-tumor-read-KRAS.txt > u-tumor-KRAS.txt
paste tumor-locs-NOTCH1.txt t-tumor-read-NOTCH1.txt > u-tumor-NOTCH1.txt
paste tumor-locs-NRAS.txt t-tumor-read-NRAS.txt > u-tumor-NRAS.txt

sed "s/.*/1, &/"  u-tumor-KRAS-val.txt > c-tumor-KRAS-val.txt
sed "s/.*/1, &/"  u-tumor-KRAS.txt > c-tumor-KRAS.txt
sed "s/.*/1, &/"  u-tumor-NOTCH1-val.txt > c-tumor-NOTCH1-val.txt
sed "s/.*/1, &/"  u-tumor-NOTCH1.txt > c-tumor-NOTCH1.txt
sed "s/.*/1, &/"  u-tumor-NRAS-val.txt > c-tumor-NRAS-val.txt
sed "s/.*/1, &/"  u-tumor-NRAS.txt > c-tumor-NRAS.txt

grep "252509[0-9][0-9]" c-tumor-KRAS-val.txt >> not-tumor-KRAS-val.txt
grep "252508[4-9][0-9]" c-tumor-KRAS-val.txt >> not-tumor-KRAS-val.txt
grep "2525083[7-9]" c-tumor-KRAS-val.txt >> not-tumor-KRAS-val.txt
grep "252046[0-9][0-9]" c-tumor-KRAS-val.txt >> not-tumor-KRAS-val.txt
grep "252047[0-7][0-9]" c-tumor-KRAS-val.txt >> not-tumor-KRAS-val.txt
grep "2520478[0-8]" c-tumor-KRAS-val.txt >> not-tumor-KRAS-val.txt

comm -3 <(sort c-tumor-KRAS-val.txt) <(sort not-tumor-KRAS-val.txt)  > tumor-KRAS-val.txt

grep "252509[0-9][0-9]" c-tumor-KRAS.txt >> not-tumor-KRAS.txt
grep "252508[4-9][0-9]" c-tumor-KRAS.txt >> not-tumor-KRAS.txt
grep "2525083[7-9]" c-tumor-KRAS.txt >> not-tumor-KRAS.txt
grep "252046[0-9][0-9]" c-tumor-KRAS.txt >> not-tumor-KRAS.txt
grep "252047[0-7][0-9]" c-tumor-KRAS.txt >> not-tumor-KRAS.txt
grep "2520478[0-8]" c-tumor-KRAS.txt >> not-tumor-KRAS.txt

comm -3 <(sort c-tumor-KRAS.txt) <(sort not-tumor-KRAS.txt)  > tumor-KRAS.txt

grep "1365457[0-9][0-9]" c-tumor-NOTCH1-val.txt >> not-tumor-NOTCH1-val.txt
grep "13654569[0-9]" c-tumor-NOTCH1-val.txt >> not-tumor-NOTCH1-val.txt
grep "13654568[7-9]" c-tumor-NOTCH1-val.txt >> not-tumor-NOTCH1-val.txt
grep "1364943[0-9][0-9]" c-tumor-NOTCH1-val.txt >> not-tumor-NOTCH1-val.txt
grep "1364944[0-2][0-9]" c-tumor-NOTCH1-val.txt >> not-tumor-NOTCH1-val.txt
grep "13649443[0-2]" c-tumor-NOTCH1-val.txt >> not-tumor-NOTCH1-val.txt

comm -3 <(sort c-tumor-NOTCH1-val.txt) <(sort not-tumor-NOTCH1-val.txt)  > tumor-NOTCH1-val.txt

grep "1365457[0-9][0-9]" c-tumor-NOTCH1.txt >> not-tumor-NOTCH1.txt
grep "13654569[0-9]" c-tumor-NOTCH1.txt >> not-tumor-NOTCH1.txt
grep "13654568[7-9]" c-tumor-NOTCH1.txt >> not-tumor-NOTCH1.txt
grep "1364943[0-9][0-9]" c-tumor-NOTCH1.txt >> not-tumor-NOTCH1.txt
grep "1364944[0-2][0-9]" c-tumor-NOTCH1.txt >> not-tumor-NOTCH1.txt
grep "13649443[0-2]" c-tumor-NOTCH1.txt >> not-tumor-NOTCH1.txt

comm -3 <(sort c-tumor-NOTCH1.txt) <(sort not-tumor-NOTCH1.txt)  > tumor-NOTCH1.txt

grep "1147168[0-9][0-9]" c-tumor-NRAS-val.txt >> not-tumor-NRAS-val.txt
grep "11471679[5-9]" c-tumor-NRAS-val.txt >> not-tumor-NRAS-val.txt
grep "1147043[0-9][0-9]" c-tumor-NRAS-val.txt >> not-tumor-NRAS-val.txt
grep "1147044[0-5][0-9]" c-tumor-NRAS-val.txt >> not-tumor-NRAS-val.txt
grep "11470446[0-3]" c-tumor-NRAS-val.txt >> not-tumor-NRAS-val.txt

comm -3 <(sort c-tumor-NRAS-val.txt) <(sort not-tumor-NRAS-val.txt)  > tumor-NRAS-val.txt

grep "1147168[0-9][0-9]" c-tumor-NRAS.txt >> not-tumor-NRAS.txt
grep "11471679[5-9]" c-tumor-NRAS.txt >> not-tumor-NRAS.txt
grep "1147043[0-9][0-9]" c-tumor-NRAS.txt >> not-tumor-NRAS.txt
grep "1147044[0-5][0-9]" c-tumor-NRAS.txt >> not-tumor-NRAS.txt
grep "11470446[0-3]" c-tumor-NRAS.txt >> not-tumor-NRAS.txt

comm -3 <(sort c-tumor-NRAS.txt) <(sort not-tumor-NRAS.txt)  > tumor-NRAS.txt

rm t-*
rm u-*
rm c-*

shuf -n 1000 tumor-KRAS-val.txt > tumor-KRAS-test1000.txt
shuf -n 1000 tumor-NOTCH1-val.txt > tumor-NOTCH1-test1000.txt
shuf -n 1000 tumor-NRAS-val.txt > tumor-NRAS-test1000.txt

shuf -n 1400000 tumor-KRAS.txt > tumor-KRAS-train1400000.txt
shuf -n 2000000 tumor-NOTCH1.txt > tumor-NOTCH1-train2000000.txt
shuf -n 75000 tumor-NRAS.txt > tumor-NRAS-train75000.txt
shuf -n 39000 tumor-NRAS.txt > tumor-NRAS-train39000.txt

cd ..

cd 1000Genomes

sort HG-KRAS-full.txt | uniq > s-HG-KRAS-full.txt
sort HG-NOTCH1-full.txt | uniq > s-HG-NOTCH1-full.txt
sort HG-NRAS-full.txt | uniq > s-HG-NRAS-full.txt
sort HG-KRAS-val-full.txt | uniq > s-HG-KRAS-val-full.txt
sort HG-NOTCH1-val-full.txt | uniq > s-HG-NOTCH1-val-full.txt
sort HG-NRAS-val-full.txt | uniq > s-HG-NRAS-val-full.txt

mv s-HG-KRAS-full.txt HG-KRAS-full.txt
mv s-HG-NOTCH1-full.txt HG-NOTCH1-full.txt
mv s-HG-NRAS-full.txt HG-NRAS-full.txt
mv s-HG-KRAS-val-full.txt HG-KRAS-val-full.txt
mv s-HG-NOTCH1-val-full.txt HG-NOTCH1-val-full.txt
mv s-HG-NRAS-val-full.txt HG-NRAS-val-full.txt

cut -f 10 HG-KRAS-val-full.txt > HG-read-KRAS-val.txt
cut -f 4 HG-KRAS-val-full.txt > HG-locs-KRAS-val.txt
cut -f 10 HG-KRAS-full.txt > HG-read-KRAS.txt
cut -f 4 HG-KRAS-full.txt > HG-locs-KRAS.txt
cut -f 10 HG-NOTCH1-val-full.txt > HG-read-NOTCH1-val.txt
cut -f 4 HG-NOTCH1-val-full.txt > HG-locs-NOTCH1-val.txt
cut -f 10 HG-NOTCH1-full.txt > HG-read-NOTCH1.txt
cut -f 4 HG-NOTCH1-full.txt > HG-locs-NOTCH1.txt
cut -f 10 HG-NRAS-val-full.txt > HG-read-NRAS-val.txt
cut -f 4 HG-NRAS-val-full.txt > HG-locs-NRAS-val.txt
cut -f 10 HG-NRAS-full.txt > HG-read-NRAS.txt
cut -f 4 HG-NRAS-full.txt > HG-locs-NRAS.txt

sed "s/.*/, &/" HG-read-KRAS-val.txt > t-HG-read-KRAS-val.txt
sed "s/.*/, &/" HG-read-KRAS.txt > t-HG-read-KRAS.txt
sed "s/.*/, &/" HG-read-NOTCH1-val.txt > t-HG-read-NOTCH1-val.txt
sed "s/.*/, &/" HG-read-NOTCH1.txt > t-HG-read-NOTCH1.txt
sed "s/.*/, &/" HG-read-NRAS-val.txt > t-HG-read-NRAS-val.txt
sed "s/.*/, &/" HG-read-NRAS.txt > t-HG-read-NRAS.txt

paste HG-locs-KRAS-val.txt t-HG-read-KRAS-val.txt > u-HG-KRAS-val.txt
paste HG-locs-NOTCH1-val.txt t-HG-read-NOTCH1-val.txt > u-HG-NOTCH1-val.txt
paste HG-locs-NRAS-val.txt t-HG-read-NRAS-val.txt > u-HG-NRAS-val.txt
paste HG-locs-KRAS.txt t-HG-read-KRAS.txt > u-HG-KRAS.txt
paste HG-locs-NOTCH1.txt t-HG-read-NOTCH1.txt > u-HG-NOTCH1.txt
paste HG-locs-NRAS.txt t-HG-read-NRAS.txt > u-HG-NRAS.txt

sed "s/.*/0, &/"  u-HG-KRAS-val.txt > c-HG-KRAS-val.txt
sed "s/.*/0, &/"  u-HG-KRAS.txt > c-HG-KRAS.txt
sed "s/.*/0, &/"  u-HG-NOTCH1-val.txt > c-HG-NOTCH1-val.txt
sed "s/.*/0, &/"  u-HG-NOTCH1.txt > c-HG-NOTCH1.txt
sed "s/.*/0, &/"  u-HG-NRAS-val.txt > c-HG-NRAS-val.txt
sed "s/.*/0, &/"  u-HG-NRAS.txt > c-HG-NRAS.txt

grep "252509[0-9][0-9]" c-HG-KRAS-val.txt >> not-HG-KRAS-val.txt
grep "252508[4-9][0-9]" c-HG-KRAS-val.txt >> not-HG-KRAS-val.txt
grep "2525083[7-9]" c-HG-KRAS-val.txt >> not-HG-KRAS-val.txt
grep "252046[0-9][0-9]" c-HG-KRAS-val.txt >> not-HG-KRAS-val.txt
grep "252047[0-7][0-9]" c-HG-KRAS-val.txt >> not-HG-KRAS-val.txt
grep "2520478[0-8]" c-HG-KRAS-val.txt >> not-HG-KRAS-val.txt

comm -3 <(sort c-HG-KRAS-val.txt) <(sort not-HG-KRAS-val.txt)  > HG-KRAS-val.txt

grep "252509[0-9][0-9]" c-HG-KRAS.txt >> not-HG-KRAS.txt
grep "252508[4-9][0-9]" c-HG-KRAS.txt >> not-HG-KRAS.txt
grep "2525083[7-9]" c-HG-KRAS.txt >> not-HG-KRAS.txt
grep "252046[0-9][0-9]" c-HG-KRAS.txt >> not-HG-KRAS.txt
grep "252047[0-7][0-9]" c-HG-KRAS.txt >> not-HG-KRAS.txt
grep "2520478[0-8]" c-HG-KRAS.txt >> not-HG-KRAS.txt

comm -3 <(sort c-HG-KRAS.txt) <(sort not-HG-KRAS.txt)  > HG-KRAS.txt

grep "1365457[0-9][0-9]" c-HG-NOTCH1-val.txt >> not-HG-NOTCH1-val.txt
grep "13654569[0-9]" c-HG-NOTCH1-val.txt >> not-HG-NOTCH1-val.txt
grep "13654568[7-9]" c-HG-NOTCH1-val.txt >> not-HG-NOTCH1-val.txt
grep "1364943[0-9][0-9]" c-HG-NOTCH1-val.txt >> not-HG-NOTCH1-val.txt
grep "1364944[0-2][0-9]" c-HG-NOTCH1-val.txt >> not-HG-NOTCH1-val.txt
grep "13649443[0-2]" c-HG-NOTCH1-val.txt >> not-HG-NOTCH1-val.txt

comm -3 <(sort c-HG-NOTCH1-val.txt) <(sort not-HG-NOTCH1-val.txt)  > HG-NOTCH1-val.txt

grep "1365457[0-9][0-9]" c-HG-NOTCH1.txt >> not-HG-NOTCH1.txt
grep "13654569[0-9]" c-HG-NOTCH1.txt >> not-HG-NOTCH1.txt
grep "13654568[7-9]" c-HG-NOTCH1.txt >> not-HG-NOTCH1.txt
grep "1364943[0-9][0-9]" c-HG-NOTCH1.txt >> not-HG-NOTCH1.txt
grep "1364944[0-2][0-9]" c-HG-NOTCH1.txt >> not-HG-NOTCH1.txt
grep "13649443[0-2]" c-HG-NOTCH1.txt >> not-HG-NOTCH1.txt

comm -3 <(sort c-HG-NOTCH1.txt) <(sort not-HG-NOTCH1.txt)  > HG-NOTCH1.txt

grep "1147168[0-9][0-9]" c-HG-NRAS-val.txt >> not-HG-NRAS-val.txt
grep "11471679[5-9]" c-HG-NRAS-val.txt >> not-HG-NRAS-val.txt
grep "1147043[0-9][0-9]" c-HG-NRAS-val.txt >> not-HG-NRAS-val.txt
grep "1147044[0-5][0-9]" c-HG-NRAS-val.txt >> not-HG-NRAS-val.txt
grep "11470446[0-3]" c-HG-NRAS-val.txt >> not-HG-NRAS-val.txt

comm -3 <(sort c-HG-NRAS-val.txt) <(sort not-HG-NRAS-val.txt)  > HG-NRAS-val.txt

grep "1147168[0-9][0-9]" c-HG-NRAS.txt >> not-HG-NRAS.txt
grep "11471679[5-9]" c-HG-NRAS.txt >> not-HG-NRAS.txt
grep "1147043[0-9][0-9]" c-HG-NRAS.txt >> not-HG-NRAS.txt
grep "1147044[0-5][0-9]" c-HG-NRAS.txt >> not-HG-NRAS.txt
grep "11470446[0-3]" c-HG-NRAS.txt >> not-HG-NRAS.txt

comm -3 <(sort c-HG-NRAS.txt) <(sort not-HG-NRAS.txt)  > HG-NRAS.txt

rm t-*
rm u-*
rm c-*

shuf -n 1000 HG-KRAS-val.txt > HG-KRAS-test1000.txt
shuf -n 1000 HG-NOTCH1-val.txt > HG-NOTCH1-test1000.txt
shuf -n 1000 HG-NRAS-val.txt > HG-NRAS-test1000.txt

shuf -n 1400000 HG-KRAS.txt > HG-KRAS-train1400000.txt
shuf -n 2000000 HG-NOTCH1.txt > HG-NOTCH1-train2000000.txt
shuf -n 39000 HG-NRAS.txt > HG-NRAS-train39000.txt

#-------------------------------------------------

cd ..

# ------- somatic mutations -----------

cat ./tumor/tumor-KRAS-train1400000.txt ./normal/normal-KRAS-train1400000.txt  | tr -d " \t"  > ref-somatic-KRAS-train1400000.csv
cat ./tumor/tumor-NOTCH1-train2000000.txt ./normal/normal-NOTCH1-train2000000.txt  | tr -d " \t"  > ref-somatic-NOTCH1-train2000000.csv
cat ./tumor/tumor-NRAS-train75000.txt ./normal/normal-NRAS-train75000.txt  | tr -d " \t"  > ref-somatic-NRAS-train75000.csv

cat ./tumor/tumor-KRAS-test1000.txt ./normal/normal-KRAS-test1000.txt  | tr -d " \t"  > ref-somatic-KRAS-test1000.csv
cat ./tumor/tumor-NOTCH1-test1000.txt ./normal/normal-NOTCH1-test1000.txt  | tr -d " \t"  > ref-somatic-NOTCH1-test1000.csv
cat ./tumor/tumor-NRAS-test1000.txt ./normal/normal-NRAS-test1000.txt  | tr -d " \t"  > ref-somatic-NRAS-test1000.csv

# ------- germline mutations -----------

cat ./1000Genomes/HG-KRAS-train1400000.txt ./normal/normal-KRAS-train1400000-g.txt  | tr -d " \t"  > ref-germline-KRAS-train1400000.csv
cat ./1000Genomes/HG-NOTCH1-train2000000.txt ./normal/normal-NOTCH1-train2000000-g.txt  | tr -d " \t"  > ref-germline-NOTCH1-train2000000.csv
cat ./1000Genomes/HG-NRAS-train39000.txt ./normal/normal-NRAS-train39000-g.txt  | tr -d " \t"  > ref-germline-NRAS-train39000.csv

cat ./1000Genomes/HG-KRAS-test1000.txt ./normal/normal-KRAS-test1000-g.txt  | tr -d " \t"  > ref-germline-KRAS-test1000.csv
cat ./1000Genomes/HG-NOTCH1-test1000.txt ./normal/normal-NOTCH1-test1000-g.txt  | tr -d " \t"  > ref-germline-NOTCH1-test1000.csv
cat ./1000Genomes/HG-NRAS-test1000.txt ./normal/normal-NRAS-test1000-g.txt  | tr -d " \t"  > ref-germline-NRAS-test1000.csv

# ------- full mutations -----------

cat ./tumor/tumor-KRAS-train1400000.txt ./1000Genomes/HG-KRAS-train1400000.txt  | tr -d " \t"  > ref-mutation-KRAS-train1400000.csv
cat ./tumor/tumor-NOTCH1-train2000000.txt ./1000Genomes/HG-NOTCH1-train2000000.txt  | tr -d " \t"  > ref-mutation-NOTCH1-train2000000.csv
cat ./tumor/tumor-NRAS-train39000.txt ./1000Genomes/HG-NRAS-train39000.txt  | tr -d " \t"  > ref-mutation-NRAS-train39000.csv

cat ./tumor/tumor-KRAS-test1000.txt ./1000Genomes/HG-KRAS-test1000.txt  | tr -d " \t"  > ref-mutation-KRAS-test1000.csv
cat ./tumor/tumor-NOTCH1-test1000.txt ./1000Genomes/HG-NOTCH1-test1000.txt  | tr -d " \t"  > ref-mutation-NOTCH1-test1000.csv
cat ./tumor/tumor-NRAS-test1000.txt ./1000Genomes/HG-NRAS-test1000.txt  | tr -d " \t"  > ref-mutation-NRAS-test1000.csv

    # samtools view *.bam "chr12:25,204,789-25,250,936" >> tumor-KRAS-val-full.txt
    # samtools view *.bam "chr9:136,494,433-136,545,786" >> tumor-NOTCH1-val-full.txt
    # samtools view *.bam "chr1:114,704,464-114,716,894" >> tumor-NRAS-val-full.txt

scp -i /Users/DoerLBH/Dropbox/git/doerlbh-rabadan.pem ubuntu@ec2-18-206-236-146.compute-1.amazonaws.com:/home/ubuntu/data/*csv . 

scp -i /Users/DoerLBH/Dropbox/git/doerlbh-rabadan.pem test_SAC.py ubuntu@ec2-18-206-236-146.compute-1.amazonaws.com:/home/ubuntu/data/  

scp -i /Users/DoerLBH/Dropbox/git/doerlbh-rabadan.pem model*.py ubuntu@ec2-54-175-210-181.compute-1.amazonaws.com:/home/ubuntu/attention/  

nohup python model_lr5em2lm1em4.py > model_lr5em2lm1em4.out & 
nohup python model_lr5em2lm1em5.py > model_lr5em2lm1em5.out & 
nohup python model_lr5em2lm1em2.py > model_lr5em2lm1em2.out & 
nohup python model_lr5em2lm1em3.py > model_lr5em2lm1em3.out &
nohup python model_lr1em1lm1em2.py > model_lr1em2lm1em2.out & 
nohup python model_lr1em1lm1em3.py > model_lr1em2lm1em3.out &
nohup python model_lr1em1lm1em4.py > model_lr1em2lm1em4.out &
nohup python model_lr1em1lm1em5.py > model_lr1em2lm1em5.out & 
nohup python model_lr1em1lm1em2.py > model_lr1em1lm1em2.out & 
nohup python model_lr1em1lm1em3.py > model_lr1em1lm1em3.out &
nohup python model_lr1em1lm1em4.py > model_lr1em1lm1em4.out &
nohup python model_lr1em1lm1em5.py > model_lr1em1lm1em5.out & 
15516
15664
15819
15960


ssh -i /Users/DoerLBH/Dropbox/git/doerlbh-rabadan.pem ubuntu@ec2-54-175-210-181.compute-1.amazonaws.com

scp -i /Users/DoerLBH/Dropbox/git/doerlbh-rabadan.pem -r ubuntu@ec2-54-175-210-181.compute-1.amazonaws.com:/home/ubuntu/attention/ . 

cat cv_loss_each_batch_lr1em1lm1em2.txt | cut -f 2 -d '(' | cut -f 1 -d ')' > cv_loss_each_batch_lr1em1lm1em2.txt2
cat cv_loss_each_batch_lr1em1lm1em3.txt | cut -f 2 -d '(' | cut -f 1 -d ')' > cv_loss_each_batch_lr1em1lm1em3.txt2
cat cv_loss_each_batch_lr1em1lm1em4.txt | cut -f 2 -d '(' | cut -f 1 -d ')' > cv_loss_each_batch_lr1em1lm1em4.txt2
cat cv_loss_each_batch_lr1em1lm1em5.txt | cut -f 2 -d '(' | cut -f 1 -d ')' > cv_loss_each_batch_lr1em1lm1em5.txt2
cat cv_loss_each_batch_lr5em2lm1em2.txt | cut -f 2 -d '(' | cut -f 1 -d ')' > cv_loss_each_batch_lr5em2lm1em2.txt2
cat cv_loss_each_batch_lr5em2lm1em3.txt | cut -f 2 -d '(' | cut -f 1 -d ')' > cv_loss_each_batch_lr5em2lm1em3.txt2 
cat cv_loss_each_batch_lr5em2lm1em4.txt | cut -f 2 -d '(' | cut -f 1 -d ')' > cv_loss_each_batch_lr5em2lm1em4.txt2
cat cv_loss_each_batch_lr5em2lm1em5.txt | cut -f 2 -d '(' | cut -f 1 -d ')' > cv_loss_each_batch_lr5em2lm1em5.txt2
cat cv_loss_val_each_batch_lr1em1lm1em2.txt | cut -f 2 -d '(' | cut -f 1 -d ')' > cv_loss_val_each_batch_lr1em1lm1em2.txt2
cat cv_loss_val_each_batch_lr1em1lm1em3.txt | cut -f 2 -d '(' | cut -f 1 -d ')' > cv_loss_val_each_batch_lr1em1lm1em3.txt2
cat cv_loss_val_each_batch_lr1em1lm1em4.txt | cut -f 2 -d '(' | cut -f 1 -d ')' > cv_loss_val_each_batch_lr1em1lm1em4.txt2
cat cv_loss_val_each_batch_lr1em1lm1em5.txt | cut -f 2 -d '(' | cut -f 1 -d ')' > cv_loss_val_each_batch_lr1em1lm1em5.txt2
cat cv_loss_val_each_batch_lr5em2lm1em2.txt | cut -f 2 -d '(' | cut -f 1 -d ')' > cv_loss_val_each_batch_lr5em2lm1em2.txt2
cat cv_loss_val_each_batch_lr5em2lm1em3.txt | cut -f 2 -d '(' | cut -f 1 -d ')' > cv_loss_val_each_batch_lr5em2lm1em3.txt2
cat cv_loss_val_each_batch_lr5em2lm1em4.txt | cut -f 2 -d '(' | cut -f 1 -d ')' > cv_loss_val_each_batch_lr5em2lm1em4.txt2
cat cv_loss_val_each_batch_lr5em2lm1em5.txt | cut -f 2 -d '(' | cut -f 1 -d ')' > cv_loss_val_each_batch_lr5em2lm1em5.txt2

mv cv_loss_each_batch_lr1em1lm1em2.txt2 cv_loss_each_batch_lr1em1lm1em2.txt
mv cv_loss_each_batch_lr1em1lm1em3.txt2 cv_loss_each_batch_lr1em1lm1em3.txt
mv cv_loss_each_batch_lr1em1lm1em4.txt2 cv_loss_each_batch_lr1em1lm1em4.txt
mv cv_loss_each_batch_lr1em1lm1em5.txt2 cv_loss_each_batch_lr1em1lm1em5.txt
mv cv_loss_each_batch_lr5em2lm1em2.txt2 cv_loss_each_batch_lr5em2lm1em2.txt
mv cv_loss_each_batch_lr5em2lm1em3.txt2 cv_loss_each_batch_lr5em2lm1em3.txt
mv cv_loss_each_batch_lr5em2lm1em4.txt2 cv_loss_each_batch_lr5em2lm1em4.txt
mv cv_loss_each_batch_lr5em2lm1em5.txt2 cv_loss_each_batch_lr5em2lm1em5.txt
mv cv_loss_val_each_batch_lr1em1lm1em2.txt2 cv_loss_val_each_batch_lr1em1lm1em2.txt
mv cv_loss_val_each_batch_lr1em1lm1em3.txt2 cv_loss_val_each_batch_lr1em1lm1em3.txt
mv cv_loss_val_each_batch_lr1em1lm1em4.txt2 cv_loss_val_each_batch_lr1em1lm1em4.txt
mv cv_loss_val_each_batch_lr1em1lm1em5.txt2 cv_loss_val_each_batch_lr1em1lm1em5.txt 
mv cv_loss_val_each_batch_lr5em2lm1em2.txt2 cv_loss_val_each_batch_lr5em2lm1em2.txt
mv cv_loss_val_each_batch_lr5em2lm1em3.txt2 cv_loss_val_each_batch_lr5em2lm1em3.txt
mv cv_loss_val_each_batch_lr5em2lm1em4.txt2 cv_loss_val_each_batch_lr5em2lm1em4.txt
mv cv_loss_val_each_batch_lr5em2lm1em5.txt2 cv_loss_val_each_batch_lr5em2lm1em5.txt