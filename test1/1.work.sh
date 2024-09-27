#!/bin/sh
export PATH="/home/software/anaconda3/bin:$PATH"
prefix="NGS045_20231208_kechuang_S1_L003_nova-III"
out_dir="/home/data1/Analysis/"$prefix

fq1="/home/data1/RawFQ/NGS045_20211120_kechuang_S1_L003_nova-III/NGS045_merge2_R1.fastq.gz"
fq2="/home/data1/RawFQ/NGS045_20211120_kechuang_S1_L003_nova-III/NGS045_merge2_R2.fastq.gz"

barcode="/home/data1/Run_script/${prefix}/barcodes.txt"
#config="/opt/software/PBDEL_pipeline_v2.6.1/lib/config.txt"
config="/home/data1/Run_script/${prefix}/config.txt"
sampleinfo="/home/data1/Run_script/${prefix}/sample_info.txt"
config_dir="/home/data1/Run_script/${prefix}/config"
codeDir="/home/software/PBDEL_pipeline_v3.1/bin/LibDisorder/"
libraryinfo="/home/software/PBDEL_pipeline_v3.1/bin/library/library_info.txt"
#mkdir $out_dir
python /home/software/PBDEL_pipeline_v3.1/bin/monitor/PCR-split.py \
	--fq1 $fq1 --fq2 $fq2 --inPCR 1 --barcode $barcode -c ${config} \
	-o $out_dir -p $prefix -s ${sampleinfo} -d ${config_dir} -l ${libraryinfo} \
	-code ${codeDir} -r "/home/data1/Run_script/${prefix}"


