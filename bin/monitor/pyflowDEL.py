########################################################################
# Description: pyflow pipeline for DEL(library, PCR, cycles)
#########################################################################
#coding=utf-8
import os.path,sys
import ast
from io import StringIO
import subprocess
import pandas as pd
import time
import math
from multiprocessing import Pool
from argparse import ArgumentParser
import process

datawarrior_app = os.environ["datawarrior_app"]
seqtk_software = os.environ["seqtk_software"]
#cutadaptPATH = os.environ['cutadaptPATH']
monitor_dir = sys.path[0] + '/'
pipeline_path = monitor_dir + '../pipeline/'
def update_dependence(self, ori_list, update_list, update_type='add'):
    tmp_includes = set(ori_list)
    if update_type == 'add':
        tmp_includes.update(update_list)
    if update_type == 'remove':
        tmp_includes.remove(update_list)
    return list(tmp_includes)

def workflow():
    data_start = time.time()
    """Main function"""
    args = parse_args()
    #config = args.config
    config = process.resolveConfig(args.config)
    caseFq = args.caseFq
    prefix = args.prefix
    sample_info = args.sampleInfo
    library_info = args.libraryInfo
    output_path = args.output_path

    # sample info
    protein_dict = process.total_sample_info(sample_info)[0]
    concentration_dict = process.total_sample_info(sample_info)[1]
    library_dict = process.total_sample_info(sample_info)[2]
    positive_dict = process.total_sample_info(sample_info)[3] 
    #sampleLists = brc.total_sample_info(sample_info)[8]
    roundNumber_dict = process.total_sample_info(sample_info)[7]
    protein = protein_dict[prefix] # protein name for the sample
    # library info
    library_seq_dict = process.library_info(library_info)[3]
    positive_sample = positive_dict[prefix]
    
    sample_dir = os.path.abspath(output_path)
    batchName = os.path.basename(output_path)
    sample_QC_dir = sample_dir + "/01.QC"
    sample_merge_dir = sample_dir + "/02.Merge"

    run_steps = []
    trimExe = int(config['trimExe'])
    mergeExe = int(config['mergeExe'])
    libraryIDExt = int(config['libraryIDExt'])
    tagfindExt = int(config['tagfindExt'])
    dwarExt = int(config['dwarExt'])
    autoPNGExt = int(config['autoPNGExt'])
    q30 = 30
    
    if trimExe == 1:
        run_steps.extend(['trim'])
    if mergeExe == 1:
        run_steps.extend(['merge'])
    if libraryIDExt == 1:
        run_steps.extend(['libraryID'])
    # 'libraryID': tagfind,dwar
    if tagfindExt == 1:
        run_steps.extend(['tagfind'])
    if dwarExt == 1:
        run_steps.extend(['dwar'])
    if autoPNGExt == 1:
        run_steps.extend(['autoPNG'])
    run_steps = list(set(run_steps))
    
    # case fastq reads 1 and 2
    tfqs = caseFq.split(",")
    tfq1 = os.path.abspath(tfqs[0])
    tfq2 = os.path.abspath(tfqs[1])

    trim_app = pipeline_path + "prep.py"
    qc_summary_app = pipeline_path + "fastq_QCsummary.py"
    cpu = int(config["numCores"])
    if 'trim' in run_steps:
        files = os.listdir(sample_QC_dir)
        file_name = ''
        for file in files:
            if file.endswith("fastq.gz"):
                file_name = file
        if file_name == '':
            case_dir = sample_QC_dir + "/" + prefix
            process.mkdir(case_dir)
            tfq1_prep = "%s/%s.trim.R1.fastq.gz" % (case_dir, prefix)
            tfq2_prep = "%s/%s.trim.R2.fastq.gz" % (case_dir, prefix)

            trim_cmd = " ".join(["python", trim_app, "--fq1", tfq1, "--fq2", tfq2, "--outfq1",
                                 tfq1_prep, "--outfq2", tfq2_prep, "--config", config, "-q", q30, "-cores", str(cpu)])
            tfq1 = tfq1_prep
            tfq2 = tfq2_prep
            subprocess.check_call(trim_cmd, shell=True)

            qc_summary_file = "{}/{}".format(sample_QC_dir, prefix)
            qc_summary_cmd = " ".join(
                ["python", qc_summary_app, "--input", qc_summary_file + '/' + prefix + '.trim.summary.txt', \
                 "--prefix", prefix, "--output", qc_summary_file + '/' + prefix + '.summary.qc.txt', \
                 "--SpecColumn",'SampleID'])
            print(qc_summary_cmd)
            subprocess.check_call(qc_summary_cmd, shell=True)
    # merge
    merge_fa = "%s/%s.combined.fa" % (sample_merge_dir, prefix)
    merge_fastq = "%s/%s.combined.fastq" % (sample_merge_dir, prefix)
    if 'merge' in run_steps:
        merge_cmd = "fastp -i %s -o %s/%s.unmerge.R1.fastq.gz --thread=%s -I %s -O %s/%s.unmerge.R2.fastq.gz -m --merged_out %s" \
                        % (tfq1, sample_merge_dir, prefix, cpu, tfq2, sample_merge_dir, prefix, merge_fastq)
        merge_cmd2 = "%s seq -A %s > %s" % (seqtk_software, merge_fastq, merge_fa)
        if os.path.exists(merge_fastq):
            print('The file {} is existed!'.format(merge_fastq))
        else:
            subprocess.check_call(merge_cmd, shell=True)
        if os.path.exists(merge_fa):
            print('The file {} is existed!'.format(merge_fa))
        else:
            subprocess.check_call(merge_cmd2, shell=True)
        
        unfound_merge_fa = "%s/%s.combined.unfound.lib.fa" % (sample_merge_dir, prefix)
        libraryID_names = library_dict[prefix].split(',')
        if os.path.exists(unfound_merge_fa):
            print('The file {} is existed!'.format(unfound_merge_fa))
        else:
            process.unfoundlib(merge_fa, library_seq_dict, unfound_merge_fa, cpu, libraryID_names)
        constant1 = 'TGACTCCCAAATCGATGTG'
        constant2 = 'CAGACAAGCTTCACCTGC'
        unsplit_out = "%s/%s.combined.unfound.lib.csv" % (sample_merge_dir, prefix)
        lib_out = "%s/%s.combined.unfound.lib.seq.txt" % (sample_merge_dir, prefix)
        if os.path.exists(unsplit_out) and os.path.exists(lib_out):
            print('The file {} is existed!'.format(unsplit_out))
        else:
            process.unfoundlib_out(sample_merge_dir, prefix, unfound_merge_fa, unsplit_out, lib_out, constant1, constant2, 12)

    p=Pool(int(cpu))
    if 'libraryID' in run_steps:
        libraryID_names = library_dict[prefix].split(',') 
        for library_name in libraryID_names:
            #print(library_name)
            #eachlibraryCountSeq(library_info,sample_dir,protein,prefix,library_name,merge_fa,run_steps,positive_sample)
            p.apply_async(eachlibraryCountSeq,args=(library_info,sample_dir,protein,prefix,library_name,merge_fa,run_steps,positive_sample))
    p.close()
    p.join()
# =============================================================================
#         known activity sequence data
# =============================================================================
    positive_lst = []
    for library_name in library_dict[prefix].split(','):
        sample_counts_library_dir = "{}/04.Counts/{}/{}".format(sample_dir,protein,library_name)
        if positive_sample != '0':
            for positive in positive_sample.split(','):
                positive_sample_name = prefix + '.' + positive + '.' + library_name
                splitCountPositive = sample_counts_library_dir + '/' + positive_sample_name + '.split.count.total.txt'
                positive_lst.append(splitCountPositive)
    positive_str = ' '.join(positive_lst)
    

    p=Pool(int(cpu))
    for library_name in library_dict[prefix].split(','):
        p.apply_async(eachlibraryDwar,args=(batchName,config,run_steps,sample_dir,sample_info,library_info,protein,prefix,library_name,positive_str,positive_lst,roundNumber_dict,concentration_dict,merge_fa))
        #eachlibraryDwar(batchName,config,run_steps,sample_dir,sample_info,library_info,protein,prefix,library_name,positive_str,positive_lst,roundNumber_dict,concentration_dict,merge_fa)
    p.close()
    p.join()
  
    print("Total time of PBDEL_pipeline: " + str(time.time() - data_start) + " seconds!")

def mergepositiveCount(single, positive_str, positive_lst, CountFile):
    df_count = pd.read_csv(single, sep='\t', header=0, dtype=object)
    if positive_str != '' and positive_str != ' ' and positive_str != '0':
        for positive in positive_lst:
            positive_df = pd.read_csv(positive, sep='\t', header=0, dtype=object)
            df_count = df_count.append(positive_df, ignore_index=True)
    df_count['Dedup_count'] = df_count['Dedup_count'].astype(int)
    df_count.sort_values(by='Dedup_count', ascending=False)
    df_count.to_csv(CountFile, sep='\t', index=False)

def eachlibraryCountSeq(library_info,sample_dir,protein,prefix,library_name,merge_fa,run_steps,positive_sample):
    # library info
    constant1_dict = process.library_info(library_info)[1]
    cycles_dict = process.library_info(library_info)[2] 
    library_seq_dict = process.library_info(library_info)[3]
    CP_umi_dict = process.library_info(library_info)[4]
    code_smiles_dir = monitor_dir + '../../lib/Libraries/'
    code_smiles = "{}/{}.code.smiles.txt".format(code_smiles_dir, library_name)
    library_seq = library_seq_dict[library_name]
    UMI_N = CP_umi_dict[library_name].count('N')
    constant1 = constant1_dict[library_name]
    cycles = cycles_dict[library_name]
    constant2 = CP_umi_dict[library_name].split('N')[-1]
    sample_split_library_dir = "{}/03.Split/{}/{}".format(sample_dir,protein,library_name)
    sample_counts_library_dir = "{}/04.Counts/{}/{}".format(sample_dir,protein,library_name)
    sample_count_ana_library_dir = "{}/05.Count.analysis/{}/{}".format(sample_dir,protein,library_name)
    sample_dwar_library_dir = "{}/06.Dwar/{}/{}".format(sample_dir,protein,library_name)
    sample_pair_library_dir = "{}/07.Pair/{}/{}".format(sample_dir,protein,library_name)
    process.mkdir(sample_split_library_dir)
    process.mkdir(sample_counts_library_dir)
    process.mkdir(sample_count_ana_library_dir)
    process.mkdir(sample_dwar_library_dir)
    process.mkdir(sample_pair_library_dir)
    
    split_fa = sample_split_library_dir + '/' + prefix + '.' + library_name + '.split.fa'
    # umi extract output
    outUMIseq = sample_split_library_dir + '/' + prefix + '.' + library_name + '.split.seq'
    # count
    outCode = sample_counts_library_dir + '/' + prefix + '.' + library_name + '.split.count'
    outUMIcode = sample_counts_library_dir + '/' + prefix + '.' + library_name + '.split.umi.count'
    
    UMIfile = sample_counts_library_dir + '/' + prefix + '.' + library_name + '.split.umi.count'
    outDedupF = sample_counts_library_dir + '/' + prefix + '.' + library_name + '.split.umi.dedup.count'
    
    splitCountTotal = sample_counts_library_dir + '/' + prefix + '.' + library_name + '.split.count.total.txt'
    splitCountTotalError = sample_counts_library_dir + '/' + prefix + '.' + library_name + '.split.count.error.txt'
    if 'tagfind' in run_steps:
        print('sequence splitting......\n')
        if os.path.exists(split_fa):
            print('The file {} is existed!\n'.format(split_fa))
        else:
            split_app(merge_fa,library_seq,split_fa)
            
        print('umi extracting......\n')
        if os.path.exists(outUMIseq):
            print('The file {} is existed!'.format(outUMIseq))
        else:
            umiExtract_app(split_fa, outUMIseq, constant1, constant2, cycles, UMI_N, library_seq)
        
        print('sequencing counting......\n')
        if os.path.exists(outCode) and os.path.exists(outUMIcode):
            print('The file {} is existed!'.format(outCode))
        else:
            seqCount_app(outUMIseq,outCode,outUMIcode)
        
        print('dedup......\n')
        if os.path.exists(outDedupF):
            print('The file {} is existed!'.format(outDedupF))
        else:
            dedup_app(cycles,UMIfile,outDedupF)
            
        print('tagfinder app starting ......\n')
        if os.path.exists(splitCountTotal):
            print('The file {} is existed!'.format(splitCountTotal))
        else:
            tagfind_app(cycles, code_smiles, outCode, outDedupF, splitCountTotal, splitCountTotalError, prefix, library_name)
        print('tagfinder end')
    if positive_sample != '0':
        for positive in positive_sample.split(','):
            positive_sample_code_smiles = "{}/{}.code.smiles.txt".format(code_smiles_dir, positive)
            positive_sample_name = prefix + '.' + positive + '.' + library_name
            splitCountPositive = sample_counts_library_dir + '/' + positive_sample_name + '.split.count.total.txt'
            splitCountPositiveError = sample_counts_library_dir + '/' + positive_sample_name + '.split.count.error.txt'
            tagfind_app(cycles, positive_sample_code_smiles, outCode, outDedupF, splitCountPositive, splitCountPositiveError, prefix, positive_sample_name)
    return 1

def eachlibraryDwar(batchName,config,run_steps,sample_dir,sample_info,library_info,protein,prefix,library_name,positive_str,positive_lst,roundNumber_dict,concentration_dict,merge_fa):
    
    cycles_dict = process.library_info(library_info)[2] 
    library_size_dict = process.library_info(library_info)[7]  
    # software or args
    dwar_app_BB124 = pipeline_path + 'counts2dwar_BB124.py'
    dwar_app_BB134 = pipeline_path + 'counts2dwar_BB134.py'
    dwar_app_BB234 = pipeline_path + 'counts2dwar_BB234.py'
    code_smiles_dir = monitor_dir + '../../lib/Libraries/'
    code_smiles = "{}/{}.code.smiles.txt".format(code_smiles_dir, library_name)
    
    cycleNM = cyclenm(code_smiles,cycles_dict,library_name)
    sample_counts_library_dir = "{}/04.Counts/{}/{}".format(sample_dir,protein,library_name)
    
    splitCountTotal = sample_counts_library_dir + '/' + prefix + '.' + library_name + '.split.count.total.txt'
    
    sample_split_library_dir = "{}/03.Split/{}/{}".format(sample_dir,protein,library_name)
    sample_count_ana_library_dir = "{}/05.Count.analysis/{}/{}".format(sample_dir,protein,library_name)
    sample_dwar_library_dir = "{}/06.Dwar/{}/{}".format(sample_dir,protein,library_name)

    CountFile = "{}/{}.count.txt".format(sample_count_ana_library_dir, prefix + '.' + library_name)
    # If there are positive samples in a library, merge the data and generate a countfile
    print('Merging count file............\n')
    mergepositiveCount(splitCountTotal, positive_str, positive_lst, CountFile)
    print(CountFile)
    print('cutoff calculation for dedup counts......\n')
    if cycleNM == 2:
        dedup_filter_threshold = process.dedup_filter_2cycle(CountFile)
    else:
        dedup_filter_threshold = process.dedup_filter(CountFile)
    print(dedup_filter_threshold)
    file_Pre_name = [batchName,roundNumber_dict[prefix],protein,
                     concentration_dict[prefix],prefix, library_name, str(library_size_dict[library_name])]
    #print(file_Pre_name)
    if 'dwar' in run_steps:
        FreqCountFile = "{}/{}.count.freq.txt".format(sample_dwar_library_dir, prefix + '.' + library_name)
        DwarFreqCountFile = "{}/{}.count.freq.dwar".format(sample_dwar_library_dir, prefix + '.' + library_name)
        FilterFreqCountFile = "{}/{}.count.freq.filter.txt".format(sample_dwar_library_dir, prefix + '.' + library_name)
        FilterDwarFreqCountFile = "{}/{}.count.freq.filter.dwar".format(sample_dwar_library_dir, prefix + '.' + library_name)
        outSummaryDLP = "{}/{}.count.sum.dlp.txt".format(sample_dwar_library_dir,prefix + '.' + library_name)
        
        sample_dysynthon_filter = '{}/{}.filter.dysynthon.txt'.format(sample_dwar_library_dir,prefix + '.' + library_name)
        sample_allsynthon_filter = '{}/{}.filter.allsynthon.txt'.format(sample_dwar_library_dir,prefix + '.' + library_name)
        outFileFilterPositive = '{}/{}.count.freq.filter.positive.txt'.format(sample_dwar_library_dir,prefix + '.' + library_name)
        outFileFilterPositiveNegative = '{}/{}.count.freq.filter.positive.negative.txt'.format(sample_dwar_library_dir,prefix + '.' + library_name)

        print('compute enrichment (z-score)......\n')
        if os.path.exists(FreqCountFile):
            print('The file {} is existed!'.format(FreqCountFile))
        else:
            enrichment_app(CountFile, FreqCountFile, str(library_size_dict[library_name]))

        #######################################################################
        # If it is an NTC sample, do not modify the count.freq.txt file, otherwise introduce the NTC as a control into the count.freq.txt file to calculate the enrichment factor EF
        ntc_dict = process.total_sample_info(sample_info)[8]
        protein_dict = process.total_sample_info(sample_info)[0]
        proteinName = protein_dict[prefix]
        if 'NTC' in proteinName:
            print('NTC sample!!!')
        else:
            ntc_sample = ntc_dict[prefix]
            ntc_name = '0'
            if ntc_sample != '0':
                ntc_name = protein_dict[ntc_sample]
                ana_dir = sample_dir + '/06.Dwar'
                ntcfile = '{}/{}/{}/{}.{}.count.freq.txt'.format(ana_dir,ntc_name,library_name,ntc_sample,library_name)
                while True:
                    if os.path.exists(ntcfile):
                        print('The paired file {} of case file {} is existed!'.format(ntc_sample,prefix))
                        df_freq = process.merge_nz_score(FreqCountFile, ntcfile)
                        df_freq.to_csv(FreqCountFile,index=None,sep='\t',columns=[
                                'SampleID','LibraryID','Tags','TAG1','TAG1_seq','TAG1_structure','TAG1_FF','TAG1_BB',
                    'TAG2','TAG2_seq','TAG2_structure','TAG2_FF','TAG2_BB',
                    'TAG3','TAG3_seq','TAG3_structure','TAG3_FF','TAG3_BB',
                    'Raw_count','Dedup_count','z-score','Normalized_z-score',
                    'EF','Type'])
                        break
                    else:
                        time.sleep(5)
                        continue
        #######################################################################
        sample_fullpath1 = sample_split_library_dir + '/' + prefix + '.' + library_name + '.split.fa'
        sample_fullpath2 = sample_split_library_dir + '/' + prefix + '.' + library_name + '.temp3.fa'
        
        Combined_counts = process.file_line(merge_fa)/2 
        Split_counts = Combined_counts
        #Coding_counts = Combined_counts
        if not os.path.exists(sample_fullpath1):
            sample_fullpath1 = merge_fa
        Split_counts = process.file_line(sample_fullpath1)/2 
        if not os.path.exists(sample_fullpath2):
            sample_fullpath2 = sample_fullpath1
        try:
            Coding_counts = process.file_line(sample_fullpath2)/2
        except:
            Coding_counts = Combined_counts
        # counts qc list
        counts_qcList = [str(dedup_filter_threshold),str(int(Combined_counts)),str(int(Split_counts)), str(int(Coding_counts))]
        #print(counts_qcList)
        outfileEXCEL = "{}/{}.count.freq.filter.dlp.xlsx".format(sample_dwar_library_dir, prefix + '.' + library_name)
        print(outfileEXCEL)
        print(dedup_filter_threshold)
        print("##############################")
        print("dlp starting ......\n")
        time_start = time.time()
        print("###############################\n")
        dlpList = []
        #try:
        dlpList = process.AggStatistics(FreqCountFile, dedup_filter_threshold, outfileEXCEL, 
                      sample_dysynthon_filter, sample_allsynthon_filter, 
                      FilterFreqCountFile,outFileFilterPositive, outFileFilterPositiveNegative,
                      10)
        #except:
        #    print('Error: The filtered freq file {} is not existed!'.format(FilterFreqCountFile))

        print("Sample {}: dlp summary finished.".format(prefix + '.' + library_name))
        print("dlp total time: {}".format(time.time()-time_start))

        
        agg_count_qc_lst = file_Pre_name + dlpList + counts_qcList
        
        if os.path.exists(outSummaryDLP):
            print('The file {} is existed!'.format(FreqCountFile))
        else:
            CountDLP_app(FreqCountFile,outSummaryDLP,str(agg_count_qc_lst),library_size_dict[library_name])
        print('txt to dwar......\n')
        # txt to dwar file
        txt2dwar_app(FreqCountFile,DwarFreqCountFile,str(dedup_filter_threshold),code_smiles)
        
        # filtered txt to dwar file
        txt2dwar_app(FilterFreqCountFile,FilterDwarFreqCountFile,str(0),code_smiles)
        print('dwar finished!')
        ########################################################################

        if cycleNM == 4:
            process.mkdir(sample_dwar_library_dir + '/BBother/')
            dwar_cmd_124 = " ".join(
                ["python", dwar_app_BB124, "--input", FreqCountFile, \
                 "--output3", sample_dwar_library_dir + '/BBother/' + prefix + '.' + library_name + '.count.freq.filter.BB124.txt', \
                 "--output_dwar4", sample_dwar_library_dir + '/BBother/' + prefix + '.' + library_name + '.count.freq.filter.BB124.dwar'])
            subprocess.check_call(dwar_cmd_124, shell=True)
            dwar_cmd_134 = " ".join(
                ["python", dwar_app_BB134, "--input",FreqCountFile,\
                 "--output3",
                 sample_dwar_library_dir + '/BBother/' + prefix + '.' + library_name + '.count.freq.filter.BB134.txt', \
                 "--output_dwar4",
                 sample_dwar_library_dir + '/BBother/' + prefix + '.' + library_name + '.count.freq.filter.BB134.dwar'])
            subprocess.check_call(dwar_cmd_134, shell=True)
            #print(dwar_cmd_134)
            dwar_cmd_234 = " ".join(
                ["python", dwar_app_BB234, "--input",FreqCountFile, \
                 "--output3",
                 sample_dwar_library_dir + '/BBother/' + prefix + '.' + library_name + '.count.freq.filter.BB234.txt', \
                 "--output_dwar4",
                 sample_dwar_library_dir + '/BBother/' + prefix + '.' + library_name + '.count.freq.filter.BB234.dwar'])
            subprocess.check_call(dwar_cmd_234, shell=True)

def cyclenm(codesmilesfile,cycles_dict,library_name):
    df = pd.read_csv(codesmilesfile,sep='\t',dtype = object)
    df3 = df[df['Cycle']=='3']
    if len(df3)==1 and df3['Structure'].tolist()[0]=='0':
        cycleNM = 2
    else:
        cycles = cycles_dict[library_name]
        cycleNM = cycles.count(',') + 1
    return cycleNM
def split_app(in_combine_fa, libraryID, out_split_fa):
    f_fa = open(out_split_fa, 'w')
    for header, seq in process.read_fasta(in_combine_fa):
        if libraryID in seq:
            f_fa.write(header + '\n')
            f_fa.write(seq + '\n')
        else:
            print('ERROR: libraryID %s mismath seq %s!'%(libraryID, seq))

def umiExtract_app(in_split_fa, umi_extract_seq, constant1, constant2, cycles, UMI_Number, library_seq):
    UMI_Number = int(UMI_Number)
    cycle1, cycle2, cycle3, cycle4, cycle5 = ['','','','','']
    cycleNM = cycles.count(',') + 1
    cycle_str = ''
    if cycleNM == 2:
        cycle1 = cycles.strip().split(',')[0]
        cycle2 = cycles.strip().split(',')[1]
        cycle_str = ['TAG1_seq', 'TAG2_seq', 'UMI']
    elif cycleNM == 3:
        cycle1 = cycles.strip().split(',')[0]
        cycle2 = cycles.strip().split(',')[1]
        cycle3 = cycles.strip().split(',')[2]
        cycle_str = ['TAG1_seq', 'TAG2_seq', 'TAG3_seq', 'UMI']
    elif cycleNM == 4:
        cycle1 = cycles.strip().split(',')[0]
        cycle2 = cycles.strip().split(',')[1]
        cycle3 = cycles.strip().split(',')[2]
        cycle4 = cycles.strip().split(',')[3]
        cycle_str = ['TAG1_seq', 'TAG2_seq', 'TAG3_seq', 'TAG4_seq', 'UMI']
    elif cycleNM == 5:
        cycle1 = cycles.strip().split(',')[0]
        cycle2 = cycles.strip().split(',')[1]
        cycle3 = cycles.strip().split(',')[2]
        cycle4 = cycles.strip().split(',')[3]
        cycle5 = cycles.strip().split(',')[4]
        cycle_str = ['TAG1_seq', 'TAG2_seq', 'TAG3_seq', 'TAG4_seq', 'TAG5_seq', 'UMI']
    else:
        print('Error cycles or beyond the max number of cycles, Please check input cycles...')
    f_prefix = os.path.dirname(umi_extract_seq) + '/' + str(os.path.basename(umi_extract_seq).split('.split')[0])
    
    f_list_seq = open(umi_extract_seq, 'w')
    f_list_seq.write('\t'.join(cycle_str) + '\n')
    
    try:
        cmd = "cutadapt -j 0 -e 0.1 -O 25 --untrimmed-output {}.untrimmed_temp0.fa -g {} -a {}$ -o {}.temp0.fa {} 2> {}.tag.temp0.log" \
            .format(f_prefix, constant1, constant2, f_prefix, in_split_fa, f_prefix)
        print(cmd)
        subprocess.check_call(cmd, shell=True)
    except:
        cmd = "cutadapt -j 1 -e 0.1 -O 25 --untrimmed-output {}.untrimmed_temp0.fa -g {} -a {}$ -o {}.temp0.fa {} 2> {}.tag.temp0.log" \
            .format(f_prefix, constant1, constant2, f_prefix, in_split_fa, f_prefix)
        print(cmd)
        subprocess.check_call(cmd, shell=True)
    
    try:
        cmd = "cutadapt -j 0 -e 0.15 -O 15 --discard-untrimmed -a {}$ -o {}.temp1.fa {}.temp0.fa 2> {}.tag.temp0.log".format( \
           constant2, f_prefix, f_prefix, f_prefix) 
        print(cmd)
        subprocess.check_call(cmd, shell=True)
    except:
        cmd = "cutadapt -j 1 -e 0.15 -O 15 --discard-untrimmed -a {}$ -o {}.temp1.fa {}.temp0.fa 2> {}.tag.temp0.log".format( \
           constant2, f_prefix, f_prefix, f_prefix) 
        print(cmd)
        subprocess.check_call(cmd, shell=True) 
    
    os.remove('{}.temp0.fa'.format(f_prefix))
    fa_temp2 = open('{}.temp2.fa'.format(f_prefix), 'w')
    for header, seq in process.read_fasta('{}.temp1.fa'.format(f_prefix)):
        UMI = seq[-UMI_Number:]
        new_seq = seq[0:len(seq)-UMI_Number]
        new_header = header + ':UMI:' + UMI
        fa_temp2.write(new_header + '\n')
        fa_temp2.write(new_seq + '\n')
    fa_temp2.close()
    os.remove('{}.temp1.fa'.format(f_prefix))
    try:
        cmd = "cutadapt -j 0 -e 0.1 -O 5 --untrimmed-output {}.untrimmed_libraryID.fa -a {} -o {}.temp3.fa {}.temp2.fa 2> {}.tag.temp0.log".\
            format(f_prefix, library_seq, f_prefix, f_prefix, f_prefix)
        print(cmd)
        subprocess.check_call(cmd, shell=True)
    except:
        cmd = "cutadapt -j 1 -e 0.1 -O 5 --untrimmed-output {}.untrimmed_libraryID.fa -a {} -o {}.temp3.fa {}.temp2.fa 2> {}.tag.temp0.log".\
            format(f_prefix, library_seq, f_prefix, f_prefix, f_prefix)
        print(cmd)
        subprocess.check_call(cmd, shell=True)
    os.remove('{}.temp2.fa'.format(f_prefix))
    log = open("{}.tag.temp0.log".format(f_prefix), 'w')

    total_count = 0
    total_count2 = 0
    total_count3 = 0
    for header, seq in process.read_fasta('{}.temp3.fa'.format(f_prefix)):
        total_count += 1
        cycles = cycles.replace(',', '')
        if len(seq) == len(cycles):
            total_count2 += 1
            if cycleNM == 2:
                cycle1_v = seq[0:len(cycle1)]
                cycle2_v = seq[len(cycle1):len(cycle1 + cycle2)]
                if cycle1_v[-2:] == cycle1[-2:] and cycle2_v[-2:] == cycle2[-2:]:
                    total_count3 += 1
                    key_seq = '\t'.join((cycle1_v, cycle2_v, header.split(':')[-1]))  
                    f_list_seq.write(key_seq + '\n')
            elif cycleNM == 3:
                cycle1_v = seq[0:len(cycle1)]
                cycle2_v = seq[len(cycle1):len(cycle1 + cycle2)]
                cycle3_v = seq[len(cycle1 + cycle2):len(seq)]
                if cycle1_v[-2:] == cycle1[-2:] and cycle2_v[-2:] == cycle2[-2:] and cycle3_v[-2:] == cycle3[-2:]:
                    total_count3 += 1
                    key_seq = '\t'.join((cycle1_v, cycle2_v, cycle3_v, header.split(':')[-1]))  
                    f_list_seq.write(key_seq + '\n')
            elif cycleNM == 4:
                cycle1_v = seq[0:len(cycle1)]
                cycle2_v = seq[len(cycle1):len(cycle1 + cycle2)]
                cycle3_v = seq[len(cycle1 + cycle2):len(cycle1 + cycle2 + cycle3)]
                cycle4_v = seq[len(cycle1 + cycle2 + cycle3):len(seq)]
                if cycle1_v[-2:] == cycle1[-2:] and cycle2_v[-2:] == cycle2[-2:] and cycle3_v[-2:] == cycle3[-2:] \
                        and cycle4_v[-2:] == cycle4[-2:]:
                    total_count3 += 1
                    key_seq = '\t'.join((cycle1_v, cycle2_v, cycle3_v, cycle4_v, header.split(':')[-1]))  
                    f_list_seq.write(key_seq + '\n')
            elif cycleNM == 5:
                cycle1_v = seq[0:len(cycle1)]
                cycle2_v = seq[len(cycle1):len(cycle1 + cycle2)]
                cycle3_v = seq[len(cycle1 + cycle2):len(cycle1 + cycle2 + cycle3)]
                cycle4_v = seq[len(cycle1 + cycle2 + cycle3):len(cycle1 + cycle2 + cycle3 + cycle4)]
                cycle5_v = seq[len(cycle1 + cycle2 + cycle3 + cycle4):len(seq)]
                if cycle1_v[-2:] == cycle1[-2:] and cycle2_v[-2:] == cycle2[-2:] and cycle3_v[-2:] == cycle3[-2:] \
                        and cycle4_v[-2:] == cycle4[-2:] and cycle5_v[-2:] == cycle5[-2:]:
                    total_count3 += 1
                    key_seq = '\t'.join((cycle1_v, cycle2_v, cycle3_v, cycle4_v, cycle5_v, header.split(':')[-1]))  
                    f_list_seq.write(key_seq + '\n')
        else:
            log.write("%s\t%s\n" % (seq, len(seq)))

    log.write("Counts after removing librarySeq: %s\n" % (total_count))
    log.write("Counts conforming to specific length: %s\n" % (total_count2))
    log.write("Counts conforming to specific bases: %s\n" % (total_count3))
def seqCount_app(seq_file, count_file, count_umi_file):
    count_f = open(count_file, 'w')
    count_umi_f = open(count_umi_file, 'w')
    seq_f = pd.read_csv(seq_file, sep="\t", header=0)
    #'cycle1\tcycle2\tcycle3\tumi\n'
    list_seq = seq_f.columns.tolist()
    list_seq.remove('UMI')
    list_umi_seq = ['UMI'] + list_seq
    temp0 = seq_f.groupby(list_seq).size()
    temp0_umi = seq_f.groupby(list_umi_seq).size()
    temp0.to_csv(count_f, sep='\t')
    temp0_umi.to_csv(count_umi_f, sep='\t')
def dedup_app(cycles, umifile, dedup_file):
    #cycle1, cycle2, cycle3, cycle4, cycle5 = ['', '', '', '', '']
    cycleNM = cycles.count(',') + 1
    umi_f = ''
    cycle_lst = []
    name_lst = []
    print(cycles)
    print(umifile)
    if cycleNM == 2:
        cycle_lst = ['cycle1', 'cycle2']
        name_lst = ['umis'] + cycle_lst + ['count']
        umi_f = pd.read_csv(umifile, sep="\t", header=0, names=name_lst)
        umi_f['seq'] = umi_f['cycle1'] + umi_f['cycle2']
    elif cycleNM == 3:
        cycle_lst = ['cycle1', 'cycle2', 'cycle3']
        name_lst = ['umis'] + cycle_lst + ['count']
        umi_f = pd.read_csv(umifile, sep="\t", header=0, names=name_lst)
        umi_f['seq'] = umi_f['cycle1'] + umi_f['cycle2'] + umi_f['cycle3']
        #umis = umi_f['umis'].drop_duplicates(keep='first')
    elif cycleNM == 4:
        cycle_lst = ['cycle1', 'cycle2', 'cycle3', 'cycle4']
        name_lst = ['umis'] + cycle_lst + ['count']
        umi_f = pd.read_csv(umifile, sep="\t", header=0, names=name_lst)
        umi_f['seq'] = umi_f['cycle1'] + umi_f['cycle2'] + umi_f['cycle3'] + umi_f['cycle4']
    elif cycleNM == 5:
        cycle_lst = ['cycle1', 'cycle2', 'cycle3', 'cycle4', 'cycle5']
        name_lst = ['umis'] + cycle_lst + ['count']
        umi_f = pd.read_csv(umifile, sep="\t", header=0, names=name_lst)
        umi_f['seq'] = umi_f['cycle1'] + umi_f['cycle2'] + umi_f['cycle3'] + umi_f['cycle4'] + umi_f['cycle5']
    temp0 = pd.DataFrame(umi_f, columns=cycle_lst)
    temp0_umi = temp0.groupby(cycle_lst).size()
    temp0_umi.to_csv(dedup_file, sep='\t')
    
def tagfind_app(cycles, tag_file, raw_count_file, dedup_count_file, out_total_count_file, error_tag_file, readname, Library_name):
    lib = str(Library_name.replace('PBDEL',''))
    cycleNM = cycles.count(',') + 1
    cycle_str = []
    error_tag = open(error_tag_file, 'w')
    codesmile = pd.read_csv(tag_file,sep='\t',header=0,dtype='object')
    codesmile['Cycle_No'] = codesmile['Cycle'] + '.' + codesmile['Number']
    dict_tag = codesmile.set_index(['Cycle_No'])['Tag'].to_dict()
    dict_strcture = codesmile.set_index(['Cycle_No'])['Structure'].to_dict()
    dict_BB = codesmile.set_index(['Cycle_No'])['Reagent'].to_dict()
    dict_tags = codesmile.set_index(['Cycle_No'])["TAG5'-3'"].to_dict()
    
    if 'StructureFF' not in codesmile.columns:
        codesmile['StructureFF'] = '0'
    dict_strctureFF = codesmile.set_index(['Cycle_No'])['StructureFF'].to_dict()
    dict_dedup = {}
    for line in open(dedup_count_file):
        dedup_line = line.strip().split('\t')
        if cycleNM == 2:
            dedup_seq = dedup_line[0] + dedup_line[1]
            dict_dedup[dedup_seq] = dedup_line[2]
        elif cycleNM == 3:
            dedup_seq = dedup_line[0] + dedup_line[1] + dedup_line[2]
            dict_dedup[dedup_seq] = dedup_line[3]
        elif cycleNM == 4:
            dedup_seq = dedup_line[0] + dedup_line[1] + dedup_line[2] + dedup_line[3]
            dict_dedup[dedup_seq] = dedup_line[4]
        elif cycleNM == 5:
            dedup_seq = dedup_line[0] + dedup_line[1] + dedup_line[2] + dedup_line[3] + dedup_line[4]
            dict_dedup[dedup_seq] = dedup_line[5]

    with open(out_total_count_file, 'w') as f:
        if cycleNM == 2:
            cycle_str = ['Tags','TAG1', 'TAG1_seq', 'TAG1_structure', 'TAG1_FF','TAG1_BB','TAG2', 'TAG2_seq', 'TAG2_structure','TAG2_FF','TAG2_BB']
        elif cycleNM == 3:
            cycle_str = ['Tags','TAG1', 'TAG1_seq', 'TAG1_structure', 'TAG1_FF','TAG1_BB','TAG2', 'TAG2_seq', 'TAG2_structure', 'TAG2_FF','TAG2_BB','TAG3', 'TAG3_seq', 'TAG3_structure','TAG3_FF','TAG3_BB']
        elif cycleNM == 4:
            cycle_str = ['Tags','TAG1', 'TAG1_seq', 'TAG1_structure', 'TAG1_FF','TAG1_BB','TAG2', 'TAG2_seq', 'TAG2_structure', 'TAG2_FF','TAG2_BB','TAG3', 'TAG3_seq', \
                         'TAG3_structure', 'TAG3_FF','TAG3_BB','TAG4', 'TAG4_seq', 'TAG4_structure','TAG4_FF','TAG4_BB']
        elif cycleNM == 5:
            cycle_str = ['Tags','TAG1', 'TAG1_seq', 'TAG1_structure', 'TAG1_FF','TAG1_BB','TAG2', 'TAG2_seq', 'TAG2_structure', 'TAG2_FF','TAG2_BB','TAG3', 'TAG3_seq', \
                         'TAG3_structure', 'TAG3_FF','TAG3_BB','TAG4', 'TAG4_seq', 'TAG4_structure', 'TAG4_FF','TAG4_BB','TAG5', 'TAG5_seq', 'TAG5_structure','TAG5_FF','TAG5_BB']
        header_sub = ['SampleID', 'LibraryID'] + cycle_str + ['Raw_count', 'Dedup_count']
        f.write('\t'.join(header_sub) + '\n')
        #print(dict_tag.keys())
        for line in open(raw_count_file):
            if cycleNM == 2:
                (tag1_seq, tag2_seq, raw_count) = line.strip().split('\t')
                seq = tag1_seq + tag2_seq
                if tag1_seq in dict_tag.values() and tag2_seq in dict_tag.values():
                    tag1_key = str(''.join([str(k) for k, v in dict_tag.items() if v == tag1_seq]))
                    tag2_key = str(''.join([str(k) for k, v in dict_tag.items() if v == tag2_seq]))
                    if seq in dict_dedup.keys():
                        dedup_count = dict_dedup[seq]
                    else:
                        dedup_count = 0
                    tags = lib + '_' + dict_tags[tag1_key] + '_' + dict_tags[tag2_key]
                    raw_dwar_str = '\t'.join(
                        (readname, Library_name, tags, tag1_key, tag1_seq, dict_strcture[tag1_key], dict_strctureFF[tag1_key],dict_BB[tag1_key], tag2_key, \
                         tag2_seq, dict_strcture[tag2_key], dict_strctureFF[tag2_key],dict_BB[tag2_key], raw_count, dedup_count))
                    f.write(raw_dwar_str + '\n')
                else:
                    error_tag.write('\t'.join((seq, raw_count)) + '\n')
            elif cycleNM == 3:
                (tag1_seq, tag2_seq, tag3_seq, raw_count) = line.strip().split('\t')
                seq = tag1_seq + tag2_seq + tag3_seq
                if tag1_seq in dict_tag.values() and tag2_seq in dict_tag.values() and tag3_seq in dict_tag.values():
                    tag1_key = str(''.join([str(k) for k, v in dict_tag.items() if v == tag1_seq]))
                    tag2_key = str(''.join([str(k) for k, v in dict_tag.items() if v == tag2_seq]))
                    tag3_key = str(''.join([str(k) for k, v in dict_tag.items() if v == tag3_seq]))
                    #print(tag3_key)
                    if seq in dict_dedup.keys():
                        dedup_count = dict_dedup[seq]
                    else:
                        dedup_count = 0
                    tags = lib + '_' + dict_tags[tag1_key] + '_' + dict_tags[tag2_key] + '_' + dict_tags[tag3_key]
                    raw_dwar_str = '\t'.join((readname, Library_name, tags, tag1_key, tag1_seq, dict_strcture[tag1_key], dict_strctureFF[tag1_key],dict_BB[tag1_key], tag2_key, \
                                          tag2_seq, dict_strcture[tag2_key], dict_strctureFF[tag2_key],dict_BB[tag2_key],tag3_key, tag3_seq, dict_strcture[tag3_key], dict_strctureFF[tag3_key],dict_BB[tag3_key],\
                                          raw_count, dedup_count))
                    f.write(raw_dwar_str + '\n')
                else:
                    error_tag.write('\t'.join((seq, raw_count)) + '\n')
            elif cycleNM == 4:
                (tag1_seq, tag2_seq, tag3_seq, tag4_seq, raw_count) = line.strip().split('\t')
                seq = tag1_seq + tag2_seq + tag3_seq + tag4_seq
                if tag1_seq in dict_tag.values() and tag2_seq in dict_tag.values() and tag3_seq in dict_tag.values() \
                        and tag4_seq in dict_tag.values():
                    tag1_key = str(''.join([str(k) for k, v in dict_tag.items() if v == tag1_seq]))
                    tag2_key = str(''.join([str(k) for k, v in dict_tag.items() if v == tag2_seq]))
                    tag3_key = str(''.join([str(k) for k, v in dict_tag.items() if v == tag3_seq]))
                    tag4_key = str(''.join([str(k) for k, v in dict_tag.items() if v == tag4_seq]))
                    if seq in dict_dedup.keys():
                        dedup_count = dict_dedup[seq]
                    else:
                        dedup_count = 0
                    tags = lib + '_' + dict_tags[tag1_key] + '_' + dict_tags[tag2_key] + '_' + dict_tags[tag3_key] + '_' + dict_tags[tag4_key]
                    raw_dwar_str = '\t'.join((readname, Library_name, tags, tag1_key, tag1_seq, dict_strcture[tag1_key], dict_strctureFF[tag1_key],dict_BB[tag1_key],tag2_key, \
                                          tag2_seq, dict_strcture[tag2_key], dict_strctureFF[tag2_key],dict_BB[tag2_key],tag3_key, tag3_seq, dict_strcture[tag3_key], dict_strctureFF[tag3_key],dict_BB[tag3_key],\
                                        tag4_key, tag4_seq, dict_strcture[tag4_key], dict_strctureFF[tag4_key],dict_BB[tag4_key],raw_count, dedup_count))
                    f.write(raw_dwar_str + '\n')
                else:
                    error_tag.write('\t'.join((seq, raw_count)) + '\n')
            elif cycleNM == 5:
                (tag1_seq, tag2_seq, tag3_seq, tag4_seq, tag5_seq, raw_count) = line.strip().split('\t')
                seq = tag1_seq + tag2_seq + tag3_seq + tag4_seq
                if tag1_seq in dict_tag.values() and tag2_seq in dict_tag.values() and tag3_seq in dict_tag.values() \
                        and tag4_seq in dict_tag.values():
                    tag1_key = str(''.join([str(k) for k, v in dict_tag.items() if v == tag1_seq]))
                    tag2_key = str(''.join([str(k) for k, v in dict_tag.items() if v == tag2_seq]))
                    tag3_key = str(''.join([str(k) for k, v in dict_tag.items() if v == tag3_seq]))
                    tag4_key = str(''.join([str(k) for k, v in dict_tag.items() if v == tag4_seq]))
                    tag5_key = str(''.join([str(k) for k, v in dict_tag.items() if v == tag5_seq]))
                    if seq in dict_dedup.keys():
                        dedup_count = dict_dedup[seq]
                    else:
                        dedup_count = 0
                    tags = lib + '_' + dict_tags[tag1_key] + '_' + dict_tags[tag2_key] + '_' + dict_tags[tag3_key] + '_' + dict_tags[tag4_key] + '_' + dict_tags[tag5_key]
                    raw_dwar_str = '\t'.join((readname, Library_name, tags, tag1_key, tag1_seq, dict_strcture[tag1_key], dict_strctureFF[tag1_key],dict_BB[tag1_key],tag2_key, \
                                          tag2_seq, dict_strcture[tag2_key], dict_strctureFF[tag2_key],dict_BB[tag2_key],tag3_key, tag3_seq, dict_strcture[tag3_key], dict_strctureFF[tag3_key],dict_BB[tag3_key],\
                                        tag4_key, tag4_seq, dict_strcture[tag4_key], dict_strctureFF[tag4_key],dict_BB[tag4_key],tag5_key, tag5_seq, \
                                              dict_strcture[tag5_key], dict_strctureFF[tag5_key],dict_BB[tag5_key],raw_count, dedup_count))
                    f.write(raw_dwar_str + '\n')
                else:
                    error_tag.write('\t'.join((seq, raw_count)) + '\n')

def enrichment_app(count_file, ratio_file, librarysize):
    count_f = pd.read_csv(count_file, sep='\t', dtype=object, header=0)
    count_f['Raw_count'] = count_f['Raw_count'].astype(int)
    count_f['Dedup_count'] = count_f['Dedup_count'].astype(int)
    compounds = count_f.shape[0] # The number of file lines is the number of compound combinations
    count_f['z-score'] = 0
    count_f['Normalized_z-score'] = 0
    df_columns = list(count_f.columns)
    df_empty = pd.DataFrame(columns=df_columns)
    df_empty.loc[0,'Dedup_count']=0
    df_empty = df_empty.fillna(0)
    df_new = count_f.append(df_empty)     
    if len(list(set(count_f['Dedup_count'].tolist()))) > 1:
        if compounds != 0:
            df_new['z-score'] = process.z_score_totallib(df_new['Dedup_count'], 0, librarysize)
            #print(df_new['z-score'])
            df_new['Normalized_z-score'] = df_new['z-score']/(math.sqrt(int(librarysize)))
    print(df_new.columns)
    del count_f
    df_new = df_new.sort_values(by='Dedup_count', ascending=False).round(4)
    df_new.reset_index(drop=True, inplace=True)
    df_new.to_csv(ratio_file, sep='\t', index=False)
def CountDLP_app(CountFreqFile, SummaryCountFile,count_qc_lst,library_size):
    name_lst = ['BatchName','Round','Protein','Protein_CC','SampleID', 'LibraryID', 'LibrarySize',
                'Plane','Line','Dot','DLP','Plane_noline','DLP_noline','Cutoff_dedup_counts','Combined_counts','Split_counts', 
                'Coding_counts', 'Raw_counts','Dedup_counts','Dup_ratio(%)', 
                'Compounds', 'CompoundRatio(%)','Coverage_lib','Raw_counts_Max', 
                'Dedup_counts_Max','Dedup_counts_TOP10_mean']
    df1 = pd.read_csv(CountFreqFile, sep='\t', header=0)
    Raw_counts = df1['Raw_count'].sum()
    Dedup_counts = df1['Dedup_count'].sum()
    Compounds = df1.shape[0] 
    CompoundRatio = '{:.2%}'.format(Compounds/int(library_size))
    Coverage_lib = '{:.2%}'.format(int(Dedup_counts)/int(library_size))
    Raw_counts_Max = df1['Raw_count'].max()
    Dedup_counts_Max = df1['Dedup_count'].max()
    Dedup_counts_TOP10_mean = df1['Dedup_count'].sort_values(ascending=False).head(10).mean()
    Dup_ratio = 0
    if Raw_counts != 0:
        Dup_ratio = '{:.2%}'.format((Raw_counts - Dedup_counts)/Raw_counts)
    file_str = [str(Raw_counts),str(Dedup_counts),str(Dup_ratio),
                str(Compounds), str(CompoundRatio),str(Coverage_lib), 
                str(Raw_counts_Max), str(Dedup_counts_Max), 
                str(Dedup_counts_TOP10_mean)]
    count_qc_lst_list = ast.literal_eval(count_qc_lst)
    with open(SummaryCountFile, 'w') as f:
        f.write('\t'.join(name_lst) + '\n')
        f.write('\t'.join(count_qc_lst_list + file_str) + '\n')
def codesmilesSeq(codesmilesfile,columnNo,count_f):
    taglist_min = []
    taglist_max = []
    codesmile = pd.read_csv(codesmilesfile,sep='\t',header=0,dtype='object')
    codesmile['CycleNum'] = codesmile['Cycle'] + '.' + codesmile['Number']
    #print(set(codesmile['Cycle'].tolist()))
    cycles = list(set(codesmile['Cycle'].tolist()))
    #print(cycles)
    for cycle in cycles:
        
        eachcycle = codesmile[codesmile['Cycle']==cycle]
        taglist_max.append(eachcycle.values[-1].tolist()[-1])
        taglist_min.append(eachcycle.values[0].tolist()[-1])
        
    taglist_min.sort()
    taglist_max.sort()
    
    count_lst = count_f.columns
    if 'TAG1_seq' in count_lst:
        lineminlst_min = ['-','-','-',taglist_min[0],'-','-','-','-',taglist_min[1],'-','-','-','-',taglist_min[2],'-','-','-','-']
        lineminlst_max = ['-','-','-',taglist_max[0],'-','-','-','-',taglist_max[1],'-','-','-','-',taglist_max[2],'-','-','-','-']
        zero_1 = ['-' for i in range(int(columnNo) - 18)]
    else:
        lineminlst_min = ['-',taglist_min[0],'-','-','-',taglist_min[1],'-','-','-',taglist_min[2],'-','-','-']
        lineminlst_max = ['-',taglist_max[0],'-','-','-',taglist_max[1],'-','-','-',taglist_max[2],'-','-','-']
        zero_1 = ['-' for i in range(int(columnNo) - 13)]
    tags_min = '\t'.join(lineminlst_min + zero_1)+'\n'
    tags_max = '\t'.join(lineminlst_max + zero_1)+'\n'
    return tags_min,tags_max  

def txt2dwar_app(ratio_file, dwar_file, dedup_filter_threshold,codesmilesfile):
    count_f = pd.read_csv(ratio_file, sep='\t', dtype=object, header=0)
    count_f['Dedup_count'] = count_f['Dedup_count'].astype(int)
    count_f2 = count_f[count_f['Tags'] != '0']
    columnNo = len(count_f2.columns)
    
    old = StringIO(count_f2.to_csv(sep='\t',index=False))
#    dedup_max = count_f['Dedup_count'].max()
    # dedup_filter_threshold = brc.dedup_filter(ratio_file)

# head
    tag_header = '''<datawarrior-fileinfo>
<version="3.2">
<created="1559631713662">
<rowcount="{0}">
</datawarrior-fileinfo>
'''.format(count_f2.shape[0])

    tag_tail = '''<datawarrior properties>
<axisColumn_2D View_0="Dedup_count">
<axisColumn_2D View_1="Normalized_z-score">
<axisColumn_3D Dedup_counts_0="TAG1">
<axisColumn_3D Dedup_counts_1="TAG2">
<axisColumn_3D Dedup_counts_2="TAG3">
<chartType_2D View="scatter">
<chartType_3D Dedup_counts="scatter">
<chartType_3D Raw_counts="scatter">
<columnWidth_Table_Dedup_count="80">
<columnWidth_Table_z_score="80">
<columnWidth_Table_Library_ID="80">
<columnWidth_Table_Protein_CC_ID="80">
<columnWidth_Table_Raw_count="80">
<columnWidth_Table_TAG1="80">
<columnWidth_Table_TAG1_seq="80">
<columnWidth_Table_TAG1_structure="80">
<columnWidth_Table_TAG1_BB="80">
<columnWidth_Table_TAG2="80">
<columnWidth_Table_TAG2_seq="80">
<columnWidth_Table_TAG2_structure="80">
<columnWidth_Table_TAG2_BB="80">
<columnWidth_Table_TAG3="80">
<columnWidth_Table_TAG3_seq="80">
<columnWidth_Table_TAG3_structure="80">
<columnWidth_Table_TAG3_BB="80">
<detailView="height[Data]=1">
<faceColor3D_3D Dedup_counts="-1250054">
<faceColor3D_3D Raw_counts="-1250054">
<fastRendering_2D View="true">
<fastRendering_3D Dedup_counts="true">
<filter0="#browser#	#disabled#	TAG1	1.001">
<filter1="#double#	TAG1">
<filter21="#double#	Dedup_count{0}{1}{2}{3}">
<filter22="#double#	z_score">
<filter23="#double#	Raw_count">
<filter24="#double#	Normalized_z-score">
<filter25="#string#	Tags">
<filter2="#string#	TAG1_seq">
<filter3="#string#	TAG1_structure">
<filter4="#string#	TAG1_BB">
<filter5="#double#	TAG2">
<filter6="#string#	TAG2_seq">
<filter7="#string#	TAG2_structure">
<filter8="#string#	TAG2_BB">
<filter9="#double#	TAG3">
<filter10="#string#	TAG3_seq">
<filter11="#string#	TAG3_structure">
<filter12="#string#	TAG3_BB">
<filter13="#double#	TAG4">
<filter14="#string#	TAG4_seq">
<filter15="#string#	TAG4_structure">
<filter16="#string#	TAG4_BB">
<filter17="#string#	TAG1_DLP">
<filter18="#string#	TAG2_DLP">
<filter19="#string#	TAG3_DLP">
<filter20="#string#	TAG_DOT">
<filterAnimation0="state=stopped delay=500">
<mainSplitting="0.7227">
<mainView="3D Raw_counts">
<mainViewCount="4">
<mainViewDockInfo0="root">
<mainViewDockInfo1="Table	bottom	0.5">
<mainViewDockInfo2="2D View	center">
<mainViewDockInfo3="3D Raw_counts	right	0.497">
<mainViewName0="Table">
<mainViewName1="2D View">
<mainViewName2="3D Raw_counts">
<mainViewName3="3D Dedup_counts">
<mainViewType0="tableView">
<mainViewType1="2Dview">
<mainViewType2="3Dview">
<mainViewType3="3Dview">
<masterView_3D Raw_counts="3D Dedup_counts">
<rightSplitting="0.62829">
<rotationMatrix_3D Dedup_counts00="0.99021">
<rotationMatrix_3D Dedup_counts01="-0.082155">
<rotationMatrix_3D Dedup_counts02="-0.11283">
<rotationMatrix_3D Dedup_counts10="0.0938">
<rotationMatrix_3D Dedup_counts11="0.99034">
<rotationMatrix_3D Dedup_counts12="0.10212">
<rotationMatrix_3D Dedup_counts20="0.10335">
<rotationMatrix_3D Dedup_counts21="-0.11171">
<rotationMatrix_3D Dedup_counts22="0.98836">
<rowHeight_Table="16">
<showNaNValues_2D View="true">
<showNaNValues_3D Dedup_counts="true">
<showNaNValues_3D Raw_counts="true">
<sizeColumn_3D Dedup_counts="Dedup_count">
<sizeColumn_3D Raw_counts="Raw_count">
</datawarrior properties>

'''.format('\t', dedup_filter_threshold, '\t', count_f2['Dedup_count'].max())
    #生成dwar文件
    tags_min,tags_max = codesmilesSeq(codesmilesfile,columnNo,count_f2)
    #with open(ratio_file, "r") as f:
    old2 = old.read() 
    #print(old)
    with open(dwar_file, 'w') as dwar:
        dwar.seek(0)
        dwar.write(tag_header)
        dwar.write(old2)
        dwar.write(tags_min)
        dwar.write(tags_max)
        dwar.write(tag_tail)
    dwar.close()
    
def parse_args():
    # Parameters to be input.
    parser = ArgumentParser(description='workflow for each sample')
    parser.add_argument("-caseFq", "--caseFq", action="store", dest="caseFq",
                        help="fastq, separate by comma", required=True)
    parser.add_argument("-c", "--config", action="store", dest="config", 
                        help="sample special configure file", required=True)
    parser.add_argument("-o", "--output_path", action="store", dest="output_path", 
                        help="sample output dir", required=True)
    parser.add_argument("-p", "--prefix", action="store", dest="prefix", 
                        help="sample prefix", required=True)
    parser.add_argument("-s", "--sampleInfo", action="store", dest="sampleInfo",
                        help="sampleInfo file", required=True)
    parser.add_argument("-l", "--libraryInfo", action="store", dest="libraryInfo",
                        help="libraryInfo file", required=True)
    args = parser.parse_args()
    tfqs = args.caseFq.split(",")
    if len(tfqs) != 2:
        print("Error: Fastq files are wrong!")
        exit(1)
    return args

if __name__ == '__main__':
    workflow()

    #workflow(args.config, args.caseFq, args.prefix, args.sampleInfo, args.libraryInfo, args.output_path)


            
