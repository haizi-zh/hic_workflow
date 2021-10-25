#!/bin/bash

########
# MIT License
#
# Copyright (c) 2021 Haizi Zheng
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
########

# Setup a workspace to run Juicer over test data: HIC003
#
# More specifically, this script does the following:
#
# 1. Create a directory named juicer_test, and setup the reuqired directory
#    structure
# 2. Download the test dataset, which are FASTQ files
# 3. Download juicertools.jar
# 3. Download other data required by Juicer, such as the reference genome,
#    restriction sites, etc. 
# 4. Run Juicer pipeline, which results in the aligned BAM file, the .hic file,
#    quality control statistics, and domains.

echo "Checking dependencies..."
snakemake --version >/dev/null 2>&1
if [[ $? -ne 0 ]]; then echo "Cannot find snakemake." && exit 1; else echo "snakemake: OK"; fi
bwa >/dev/null 2>&1
if [[ $? -ne 1 ]]; then echo "Cannot find bwa." && exit 1; else echo "bwa: OK"; fi
java -version >/dev/null 2>&1
if [[ $? -ne 0 ]]; then echo "Cannot find java." && exit 1; else echo "java: OK"; fi

set -x

mkdir -p juicer_test
pushd juicer_test

echo "Obtaining scripts ..."
if [[ -e script ]]; then echo "Remove directory script and retry." && exit 1; fi
git_temp=$(mktemp -d)
pushd $git_temp
git clone https://github.com/epifluidlab/juicer4snakemake.git >/dev/null 2>&1
cd juicer4snakemake
git checkout 7ada3795 >/dev/null 2>&1
popd
mv $git_temp/juicer4snakemake/script .
rm -rf $git_temp

echo "Obtaining juicer_tools.jar ..."
curl -L -o "script/juicer_tools.jar" https://s3.amazonaws.com/hicfiles.tc4ga.com/public/juicer/juicer_tools_1.22.01.jar

mkdir -p frag/HIC003
echo "Downloading FASTQ test data ..."
curl -L -o "frag/HIC003/HIC003_S2_L001_R1.fastq.gz" http://juicerawsmirror.s3.amazonaws.com/opt/juicer/work/HIC003/fastq/HIC003_S2_L001_R1_001.fastq.gz
curl -L -o "frag/HIC003/HIC003_S2_L001_R2.fastq.gz" http://juicerawsmirror.s3.amazonaws.com/opt/juicer/work/HIC003/fastq/HIC003_S2_L001_R2_001.fastq.gz

echo "Downloading restriction sites ..."
mkdir -p restriction_sites
curl -L -o "restriction_sites/hg19_DpnII.txt" https://s3.amazonaws.com/juicerawsmirror/opt/juicer/restriction_sites/hg19_DpnII.txt

echo "Downloading reference genome ..."
mkdir -p references
pushd references
curl -L -O https://s3.amazonaws.com/juicerawsmirror/opt/juicer/references/Homo_sapiens_assembly19.fasta
curl -L -O https://s3.amazonaws.com/juicerawsmirror/opt/juicer/references/Homo_sapiens_assembly19.fasta.amb
curl -L -O https://s3.amazonaws.com/juicerawsmirror/opt/juicer/references/Homo_sapiens_assembly19.fasta.ann
curl -L -O https://s3.amazonaws.com/juicerawsmirror/opt/juicer/references/Homo_sapiens_assembly19.fasta.bwt
curl -L -O https://s3.amazonaws.com/juicerawsmirror/opt/juicer/references/Homo_sapiens_assembly19.fasta.pac
curl -L -O https://s3.amazonaws.com/juicerawsmirror/opt/juicer/references/Homo_sapiens_assembly19.fasta.sa
popd

popd
exit 0

snakemake_major_ver=$(snakemake --version | grep -E -o '^[0-9]+')
echo $?
echo $snakemake_major_ver
if [[ ! "7.0" -ge 7 ]]
then
	echo LLAJDF
	exit 1
fi

if [[ $? -ne 0 || ! $snakemake_major_ver -ge 6 ]]
then
	echo Invalid
fi
exit 1

if [[ -e juicer_test ]]
then
	echo "The directory juicer_test already exists. Consider removing it and running the script again."
	exit 1
fi

mkdir -p juicer_test
pushd juicer_test

echo "Obtaining scripts from "

mkdir -p juicer_test/HIC003

echo "Downloading Juicer tools jar..."
curl -L -o "../"

echo "Downloading FASTQ test data..."
curl -L -o "juicer_test/HIC003/HIC003_S2_L001_R1.fastq.gz" http://juicerawsmirror.s3.amazonaws.com/opt/juicer/work/HIC003/fastq/HIC003_S2_L001_R1_001.fastq.gz
curl -L -o "juicer_test/HIC003/HIC003_S2_L001_R2.fastq.gz" http://juicerawsmirror.s3.amazonaws.com/opt/juicer/work/HIC003/fastq/HIC003_S2_L001_R2_001.fastq.gz

