#!/opt/software/anaconda3/bin python
#coding=utf-8
import os.path,sys
import subprocess
import pandas as pd
import shutil
import time
from multiprocessing import Pool
from argparse import ArgumentParser
#sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), '../pipeline'))
#from pyflow import WorkflowRunner
import process # 相当于旧版本的brc函数 

#########################################################################
program = os.path.basename(sys.argv[0])
__author__ = "Dong Keke"
__contact__ = "dong_keke@PharmaBlock.com"
__copyright__ = "Copyright 2019, PharmaBlock"
__date__ = "2020/09/01"
__version__ = "2.6"

seqtk_demultiplex = os.environ["seqtk_demultiplex"]

def main():
    start = time.time()
    """Main function"""
    args = parse_args()
    sample_cfg = os.path.abspath(args.config)
    config = process.resolveConfig(sample_cfg)
    inPCR = args.inPCR
    barcode = os.path.abspath(args.barcode)
    read1 = os.path.abspath(args.caseFq1)
    read2 = os.path.abspath(args.caseFq2)
    runscript_path = os.path.abspath(args.runscript_path)
    output_path = os.path.abspath(args.output_path)
    prefix = args.prefix
    sample_info = os.path.abspath(args.samples)
    library_info = os.path.abspath(args.libraryInfo)
    config_dir = os.path.abspath(args.configDir)
    code_smiles_dir = os.path.abspath(args.codeDir)
    
    # generate barcode file for all samples
    generate_barcodes_code(sample_info, barcode)
    
    sample_dir = os.path.abspath(output_path) # Analysis directory
    process.mkdir(sample_dir)
    sample_script_dir = sample_dir + "/00.Script/" # script dir for each sample
    sample_QC_dir = sample_dir + "/01.QC/"
    sample_merge_dir = sample_dir + "/02.Merge/"
    sample_split_dir = sample_dir + "/03.Split/"
    sample_counts_dir = sample_dir + "/04.Counts/"
    sample_count_ana_dir = sample_dir + "/05.Count.analysis/"
    sample_dwar_dir = sample_dir + "/06.Dwar/"
    sample_pair_dir = sample_dir + "/07.Pair/"
    sample_contamination_dir = sample_dir + "/08.Contamination/"
    sample_report_dwar_dir = sample_dir + "/09.report/dlp/"
    #sample_report_PNG_dir = sample_dir + "/09.report/PNG"
    # sample_report_table_dir = sample_dir + "/09.report/table"
    #sample_report_image_dir = sample_dir + "/09.report/image/"
    sample_report_table_dir = sample_dir + "/09.report/table/"
    docking_dir = sample_dir + "/VinaDocking/"

    process.mkdir(sample_dir)
    process.mkdir(sample_script_dir)
    process.mkdir(sample_QC_dir)
    process.mkdir(sample_merge_dir)
    process.mkdir(sample_report_dwar_dir)
    #process.mkdir(sample_report_image_dir)
    process.mkdir(sample_contamination_dir)
    process.mkdir(sample_split_dir)
    process.mkdir(sample_counts_dir)
    process.mkdir(sample_dwar_dir)
    process.mkdir(sample_pair_dir)
    process.mkdir(sample_count_ana_dir)
    process.mkdir(sample_report_table_dir)
    process.mkdir(docking_dir)
    
    runscript_dir = os.path.abspath(runscript_path) # script dir for batch data
    
    monitor_dir = sys.path[0] + '/'
    pipeline_path = monitor_dir + '../pipeline/'
    pyflow_app = monitor_dir + './pyflowDEL.py'
    autoPNG_app = pipeline_path + 'auto_png_pair.py'
    identifyFP_app = pipeline_path + 'IdentifyFP.py'
    QCsummary = pipeline_path + 'QCsummary.py'
    code_smiles_dir = monitor_dir + '../lib/Libraries/'
    ef_plot_app = pipeline_path + 'EF_pair.py'
    offDNA_app = monitor_dir + '../offDNA/offDNAenu.py'
    del_ai_app = monitor_dir + '../offDNA/DEL_predict/DEL_predict.py'
    del_ai_model_app = monitor_dir + '../offDNA/DEL_predict/RF.pkl'
    #del_report_database = monitor_dir + '../DELreport/'
    #web_report_app = monitor_dir + '../DELreport/DELreport_generate.py'
    # 对外报告需要调整，尽量能输出html文件

    #static_del_report =str(config['static_del_report'])
    #del_report_database = str(config['del_report_database'])

#    cpython = config["python"]
    trim_app = pipeline_path + "prep.py"
    #qc_summary_app = config['qc_summary_app']
    #pyflow_app = config["pyflow_app"]
    #fastp_app = config['fastp_app']
    numcores = config['numCores']
    #minlength = config['minLength']
    #autoPNG_app = config['autoPNG_app']
    #seqtk_demultiplex = config['seqtk_demultiplex']
    #report_app = config['report_app']
    #datawarrior_app = config['datawarrior_app']
    #identifyFP_app = config['identifyFP_app']
    #code_smiles_dir = config['code_smiles_dir']
    #files2xlsx_app = config['files2xlsx_app']
    #ef_plot_app = config['ef_plot_app']
    #offDNAcommentsfile = config['offDNAcommentsfile']
    #qc_commentsfile = config['QC_commentsfile']
    #EF_pair_app = config['ef_plot_app']
    #EF_comments = config['EF_comments']

    #offDNA_app = config['offDNA_app']
    #autodockvina_app= config['autodockvina_app']
    #pdb_path= config['pdb_path']
    #protein_lst= config['protein_lst']
    #target_dir = config['target_dir']

    #sar_cutoff = float(config['sar_cutoff'])
    
    sampleNM_MPI = int(config['sampleNM_MPI'])
    multi_numcores = str(int(sampleNM_MPI) * int(numcores))
    #q30 = str(config['q30'])
    single = str(config['single_png'])
    
    #del_ai_app = str(config['del_ai'])
    #del_ai_model_app = str(config['del_ai_model'])
    #web_report_app = str(config['web_report_app'])
    #static_del_report =str(config['static_del_report'])
    #del_report_database = str(config['del_report_database'])
    
    #jobs = numcores*sampleNM_MPI
    app_lst = [ef_plot_app, multi_numcores,autoPNG_app,
               single,identifyFP_app,QCsummary,del_ai_app,del_ai_model_app,offDNA_app]
   
    if inPCR == '1':
        process.check_singlefile(barcode)
        sample_trim_temp_dir = sample_dir + "/01.Trim.temp"
        sample_trim_dir = sample_dir + "/01.Trim/"
        case_prefix = prefix
        process.mkdir(sample_trim_dir)
        process.mkdir(sample_trim_temp_dir)
        
        # 1. trim and qc
        tfq1 = "%s/%s.trim.R1.fastq.gz" % (sample_QC_dir, case_prefix)
        tfq2 = "%s/%s.trim.R2.fastq.gz" % (sample_QC_dir, case_prefix)
        outsummary_f = sample_QC_dir + '/' + case_prefix + '.summary.qc.txt'
        trim_cmd = " ".join(["python", trim_app, "--fq1", read1, "--fq2", read2, "--outfq1", \
                             tfq1, "--outfq2", tfq2, "--config", sample_cfg, "-cores", multi_numcores, \
                                 "--outsummary", outsummary_f])
        print(trim_cmd)
        subprocess.check_call(trim_cmd, shell=True)
        os.system("cp {} {}".format(outsummary_f,sample_report_table_dir))
        
        # 2. split samples
        #seqtk_demultiplex is used to split data based on barcode.
        seqtk_cmd = "{} -1 {} -2 {} -b {} -d {}".format(seqtk_demultiplex, tfq1, tfq2, barcode, sample_trim_temp_dir)
        print(seqtk_cmd)
        subprocess.check_call(seqtk_cmd, shell=True)
        print('split samples finished!')
        # check the barcode file and sample names
        
        # all samples listf
        sample_lst = process.total_sample_info(sample_info)[6]
        df = pd.read_csv(barcode, header=None, sep='\t')
        sample_names_barcode = df.iloc[:, 0].tolist()
        sample_names = process.compareList(sample_lst,sample_names_barcode)
        print(sample_names)
        # 3. sample qc
        sample_run_sh = {}
        p=Pool(int(numcores))
        for sample_name in sample_names:
            #sampleTrim(sample_trim_temp_dir, sample_name,sample_info,sample_trim_dir,numcores)
            p.apply_async(sampleTrim,args=(sample_trim_temp_dir, sample_name,sample_info,sample_trim_dir,numcores))
        p.close()
        p.join()
        
        outfile = sample_report_table_dir + prefix + '.samples.summary.qc.txt'
        process.mergefile(sample_trim_dir,'.summary.qc.txt',outfile)
        
        # 4. generate script for each sample
        for sample_name in sample_names:
            # generate config file for each sample: config_sub
            config_sub = generate_config_code(sample_info, sample_cfg, sample_name, config_dir, code_smiles_dir)
            # generate script for each sample
            outfq1 = sample_trim_dir + '/' + sample_name + '_R1.fastq.gz'
            outfq2 = sample_trim_dir + '/' + sample_name + '_R2.fastq.gz'
            generate_run_code("python", pyflow_app, outfq1, outfq2, config_sub, output_path, sample_name, sample_script_dir, sample_info, library_info)
            sample_run_sh["%s/%s_run.sh" % (sample_script_dir, sample_name)] = "%s/%s_run.out" % (sample_script_dir, sample_name)
        shutil.rmtree(sample_trim_temp_dir)
        generate_script_qc_code(prefix, sample_info, library_info, sample_dir, sample_names, runscript_dir + '/2.qc.sh',app_lst)
        #print(sample_run_sh)
        SampleRunOrder(prefix, sample_run_sh, runscript_dir, sampleNM_MPI)

    elif inPCR == '0':
        generate_run_code("python", pyflow_app, read1, read2, sample_cfg, output_path, prefix, sample_script_dir, sample_info, library_info)
        df = pd.read_csv(barcode, header=None, sep='\t')
        sample_names = df.iloc[:, 0].tolist()
        generate_script_qc_code(prefix, sample_info, library_info, sample_dir, sample_names, runscript_dir + '/2.qc.sh',app_lst)
        SampleRunOrder(prefix, sample_run_sh, runscript_dir, sampleNM_MPI)
    end = time.time()
    print('Task of %s runs %0.2f seconds.' %(args.prefix, end - start))
    
def sampleTrim(sample_trim_temp_dir, sample_name,sample_info,sample_trim_dir,numcores):
    read1 = sample_trim_temp_dir + '/' + sample_name + '_R1.fastq'
    read2 = sample_trim_temp_dir + '/' + sample_name + '_R2.fastq'
    outfq1 = sample_trim_dir + '/' + sample_name + '_R1.fastq.gz'
    outfq2 = sample_trim_dir + '/' + sample_name + '_R2.fastq.gz'

    I5_dict = process.total_sample_info(sample_info)[4][sample_name]
    I3_dict = process.total_sample_info(sample_info)[5][sample_name]

    adapter1 = process.reverse_complement(I3_dict)
    adapter2 = process.reverse_complement(I5_dict)
    sample_prefix = sample_trim_dir + '/' + sample_name
    if int(numcores) > 16:
        numcores = 16
    cmd = "fastp -i {} -o {} -I {} -O {} --adapter_sequence={} --adapter_sequence_r2={} --thread={} --length_required=15 " \
          "--compression=6 -q 30 --trim_poly_g --cut_tail --correction -j {}.fastp.json -h {}.fastp.html 2> {}.log".format( 
                              read1, outfq1, read2, outfq2, adapter1, adapter2, 
                              numcores, sample_prefix,sample_prefix, sample_prefix)
    print(cmd)
    subprocess.check_call(cmd, shell=True)
    
    # get qc file for each sample
    each_sample_summary_file = "{}.summary.txt".format(sample_prefix)
    process.qc_transfer(sample_prefix, each_sample_summary_file)
    sampleQC = "{}.summary.qc.txt".format(sample_prefix)
    process.run_qc(each_sample_summary_file, sampleQC, sample_name, "SampleID")
    
    # script and log file for each sample
    
def generate_barcodes_code(sample_info_f, barcode_file):
    sample_f = pd.read_csv(sample_info_f, sep='\t', dtype=object, header=[0])
    barcode_f = pd.DataFrame(sample_f, columns=['SampleID', 'I5', 'I3'])
    barcode_f.to_csv(barcode_file, sep='\t', header=None, index=None)

def generate_run_code(cpython, pyflow_app, read1, read2, sample_cfg, output_path, samplename, script_dir, sampleInfo, libraryInfo):
    script_file = open("%s/%s_run.sh" % (script_dir, samplename), "w")
    run_header = '''#!/bin/sh'''
    script_file.write(run_header + "\n")
    run_main = '''
{0} {1} --caseFq {2},{3} -c {4} -p {5} -o {6} -s {7} -l {8}
    '''.format(cpython, pyflow_app, read1, read2, sample_cfg, samplename, output_path, sampleInfo, libraryInfo)
    script_file.write(run_main + "\n")
    script_file.close()
    return
#配置时修改以下路径
def generate_script_qc_code(run_id, sample_info, library_info, sample_dir, sample_lst, out_qc_script,app_lst):
    qc_f = open(out_qc_script,'w')
    
    template = """#!/bin/sh
# ../pipeline/EF_pair.py
python {4} \
-i {0} \
-b {1} \
-s {2} \
-l {3} \
-mpi {5}

# ../pipeline/auto_png_pair.py
python {6} \
-d {0} \
-s {2} \
-l {3} \
-mpi {5} \
-single {7}

# ../pipeline/IdentifyFP.py
python {8} \
-i {0}/09.report/table/{1}.lib.summary.qc.txt \
-o {0}/09.report/table/{1}.count.sum.dlp.detail.txt \
-b {1} \
-s {2} \
-l {3}

# ../pipeline/QCsummary.py
python {9} \
-i {0}/09.report/table/ \
-b {1} \
-s {2} \
-l {3} \
-o {0}/09.report/table/{1}.qc.xlsx

# ../offDNA/DEL_predict/DEL_predict.py
python {10} \
-s {2} \
-i1 {0}/09.report/table/{1}.count.sum.dlp.detail.txt \
-i2 {0}/06.Dwar \
-m {11} \
-o {0}/09.report/table/{1}.count.sum.dlp.detail.confirmation.txt

# ../offDNA/offDNAenu.py
python {12} \
-i {0}/09.report/table/{1}.count.sum.dlp.detail.confirmation.txt \
-o1 {0}/09.report/table/ \
-o2 {0}/VinaDocking/ \
-a {0} \
-s {2} \
-n {5}

    """.format(sample_dir, run_id, sample_info,library_info,
    app_lst[0],app_lst[1],app_lst[2],app_lst[3],app_lst[4],app_lst[5],
    app_lst[6],app_lst[7],app_lst[8])
    qc_f.write(template)
    qc_f.close()


def generate_config_code(sample_info_f, common_cfg, samplename, config_dir, code_smiles_dir):
    #print(config_dir)
    process.mkdir(config_dir)
    I5_dict = process.total_sample_info(sample_info_f)[4][samplename]
    I3_dict = process.total_sample_info(sample_info_f)[5][samplename]
    #dedup_method_dict = process.total_sample_info(sample_info_f)[6][samplename]
    #plot_threshold_dict = process.total_sample_info(sample_info_f)[7][samplename]

    config_header = '''[sequence]
CaseSampleID = {0}
I5 = {1}
I3 = {2}
overhang = NA
code_smiles_dir = {3}
'''.format(samplename, I5_dict, I3_dict, code_smiles_dir)
# full path for positive file, but not for code file
    config_sub = "%s/%s.config.txt" % (config_dir, samplename)
    with open(common_cfg, "r") as f:
        old = f.read()
        with open(config_sub, "w") as config:
            config.seek(0)
            config.write(config_header)
            config.write(old)
    return config_sub

def eachSampleRun(sample, shDict):
    time.sleep(5)
    #print "\nRun task %s-%s" %(sample, os.getpid())
    start = time.time()
    sample_cmd = " ".join(["sh", sample, ">", shDict[sample],"2>&1"])
    subprocess.check_call(sample_cmd, shell=True)
    end = time.time()
    print('Task %s runs %0.2f seconds.' %(sample, end - start))
    return 1


def SampleRunOrder(batch, shDict, runscript_dir, sampleNM_MPI=5):
    # sampleNM_MPI: 同时运行的样本个数
    p=Pool(int(sampleNM_MPI))
    for eachSample in shDict.keys():
        p.apply_async(eachSampleRun,args=(eachSample,shDict))
        #hasRun.append(res)
    p.close()
    p.join() # Waiting for all subprocesses done...
    
    start = time.time()
    report_cmd = " ".join(["sh", runscript_dir + '/2.qc.sh', ">", runscript_dir + '/2.qc.out',"2>&1"])
    os.system(report_cmd)
    
    end = time.time()
    print('Task report of %s runs %0.2f seconds.' %(batch, end - start))  
    
def parse_args():
    # Parameters to be input.
    parser = ArgumentParser(description='del analysis pipeline')
    
    parser.add_argument("-fq1", "--caseFq1", action="store", dest="caseFq1",
                        help="The fastq file read1. The format is *.fastq.gz or *.fastq.", required=True, metavar='<string>')
    parser.add_argument("-fq2", "--caseFq2", action="store", dest="caseFq2",
                        help="The fastq file read2. The format is *.fastq.gz or *.fastq.", required=True, metavar='<string>')
    parser.add_argument("-inPCR", action="store", dest="inPCR",
                        help="1 if constant in PCR, 0 or not", required=True, metavar='<int>')
    parser.add_argument("-barcode", action="store", dest="barcode",
                        help="Output barcode file, which containg the barcodes of all samples.", required=True, metavar='<string>')
    parser.add_argument("-c", "--config", action="store", dest="config",
                        help="Sample special configure file", required=True, metavar='<string>')
    parser.add_argument("-d", "--configDir", action="store", dest="configDir",
                        help="Output config directionary for each sample.", required=True, metavar='<string>')
    parser.add_argument("-r", "--runscript_path", action="store", dest="runscript_path",
                        help="Run script directory for each batch data", required=True, metavar='<string>')
    parser.add_argument("-o", "--output_path", action="store", dest="output_path",
                        help="Output directory for each batch data", required=True, metavar='<string>')
    parser.add_argument("-p", "--prefix", action="store", dest="prefix",
                        help="The name of each batch data.", required=True, metavar='<string>')
    parser.add_argument("-s", "--samples", action="store", dest="samples",
                        help="A table for sample information, txt file and 'TAB' as the seperate.", required=True, metavar='<string>')
    parser.add_argument("-l", "--libraryInfo", action="store", dest="libraryInfo",
                        help="A table for library information, txt file and 'TAB' as the seperate.", required=True, metavar='<string>')
    parser.add_argument("-code", "--codeDir", action="store", dest="codeDir",
                        help="TagCode~structure file dir", required=True, metavar='<string>')
    
    return parser.parse_args()

if __name__ == '__main__':
    main()
    
