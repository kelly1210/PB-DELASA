#!/bin/bash
prefix="NGS001"
out_dir="../example/"$prefix
fq1="../example/test_r1.fastq.gz"
fq2="../example/test_r2.fastq.gz"
barcode="./barcodes.txt"
config="../lib/config.txt"
sampleinfo="./sample_info.txt"
config_dir="./config"
codeDir="../lib/Libraries/"
libraryinfo="../lib/library_info.txt"
#mkdir $out_dir
${pbdelenv}/bin/python ../bin/monitor/pbdel.py -fq1 $fq1 -fq2 $fq2 \
-inPCR 1 -barcode $barcode -c ${config} -o $out_dir \
-p $prefix -s ${sampleinfo} \
-d ${config_dir} -l ${libraryinfo} -code ${codeDir} -r "./"

