# #################################################################################
# # Author:      Keke Dong
# # Email:       dongkeke@pharmablock.com
# # Date:        20190401
# # Description: Pipline of DEL analysis
# # Version:     v1.0


import os
import subprocess
import argparse
import time
import re
import sys
sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), '../monitor'))
import process

def run(confg, read1, read2, outfq1, outfq2, njobs, outsummaryfile):
    data_start = time.time()
    print(time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(data_start)) + " prep: starting read prep - trimming")
    
    #fastpdir = confg['fastp_app']
    #numcores = confg['numCores']
    minlength = 15
    adapter1 = confg['adapter1']
    adapter2 = confg['adapter2']
    umiextra = 'False'
    umipos = ''
    umilength = ''
    skiplength = ''

    outdir = os.path.dirname(outfq1)
    readname = re.search(r'(\S+)\.trim', os.path.basename(args.outfq1)).group(1)
    fileprefix = "{}/{}.trim".format(outdir, readname)

    if int(njobs) > 16:
        njobs = "16"

    if read1.endswith(".gz"):
        readnum = subprocess.check_output("zcat {}|head -1|".format(read1) + "perl -lane 'if($_=~/\/(\d)$/){print $1}'", shell=True)
        if readnum.strip() == "1":
            subprocess.check_call("zcat {}".format(read1) + "|perl -lane '$n++;if($n%4==1){$_=~s/\/\d+$//;print}else{print}'"
                + " >{}.temp.R1.fastq".format(fileprefix), shell=True)
            subprocess.check_call("zcat {}".format(read2) + "|perl -lane '$n++;if($n%4==1){$_=~s/\/\d+$//;print}else{print}'"
                + " >{}.temp.R2.fastq".format(fileprefix), shell=True)
            read1 = "{}.temp.R1.fastq".format(fileprefix)
            read2 = "{}.temp.R2.fastq".format(fileprefix)
    else:
        readnum = subprocess.check_output("head -1 {}|".format(read1) + "perl -lane 'if($_=~/\/(\d)$/){print $1}'", shell=True)
        if readnum.strip() == "1":
            subprocess.check_call("less {}".format(read1) + "|perl -lane '$n++;if($n%4==1){$_=~s/\/\d+$//;print}else{print}'"
                + " >{}.temp.R1.fastq".format(fileprefix), shell=True)
            subprocess.check_call("less {}".format(read2) + "|perl -lane '$n++;if($n%4==1){$_=~s/\/\d+$//;print}else{print}'"
                + " >{}.temp.R2.fastq".format(fileprefix), shell=True)
            read1 = "{}.temp.R1.fastq".format(fileprefix)
            read2 = "{}.temp.R2.fastq".format(fileprefix)
    if umiextra == "true":
        cmd = "fastp -i {} -o {} -I {} -O {} --adapter_sequence={} --adapter_sequence_r2={} --thread={} --length_required={} " \
            "--compression=6 --trim_poly_g -q 30 --cut_by_quality3 --correction --umi --umi_loc {} --umi_len {} --umi_skip {} " \
            "-j {}.fastp.json -h {}.fastp.html 2> {}.log".format(read1, outfq1, read2, outfq2, adapter1, adapter2, njobs,
            minlength, umipos, umilength, skiplength, fileprefix, fileprefix, fileprefix)
    else:
        cmd = "fastp -i {} -o {} -I {} -O {} --adapter_sequence={} --adapter_sequence_r2={} --thread={} --length_required={} " \
              "--compression=6 --trim_poly_g -q 30 --cut_tail --correction -j {}.fastp.json -h {}.fastp.html 2> {}.log".format(
              read1, outfq1, read2, outfq2, adapter1, adapter2, njobs, minlength, fileprefix, fileprefix, fileprefix)

    if not os.path.isfile("{}.fastp.json".format(fileprefix)):
        subprocess.check_call(cmd, shell=True)

    if os.path.isfile("{}.temp.R1.fastq".format(fileprefix)):
        os.remove("{}.temp.R1.fastq".format(fileprefix))
        os.remove("{}.temp.R2.fastq".format(fileprefix))
    summary_file = "{}.summary.txt".format(fileprefix)
    process.qc_transfer(fileprefix, summary_file)
    
    print("QC_Summary initiated.\n")
    process.run_qc(summary_file, outsummaryfile, readname, "BatchName")
    
    print(time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(time.time())) + " prep: done")
    print("Total time of trimming and QC: " + str(time.time()-data_start))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='trim the adapter from raw fastq')
    parser.add_argument('--fq1', '-f1', help='fq1 file')
    parser.add_argument('--fq2', '-f2', help='fq2 file')
    parser.add_argument('--config', '-c', help='config file')
    parser.add_argument('--outfq1', '-of1', help='output fq1 file')
    parser.add_argument('--outfq2', '-of2', help='output fq2 file')
    #parser.add_argument('--quality', '-q', help='quality q30')
    parser.add_argument('--cores', '-cores', help='number cores')
    parser.add_argument('--outsummary', '-outsummary', help='output summary file')
    args = parser.parse_args()
    
    cfg = process.resolveConfig(args.config)
    run(cfg, args.fq1, args.fq2, args.outfq1, args.outfq2, args.cores,args.outsummary)
