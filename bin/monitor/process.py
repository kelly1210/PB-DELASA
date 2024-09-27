#!/opt/software/anaconda3/bin python
#coding=utf-8
version = "v2.0"
# add agg sta analysis in the end
from subprocess import Popen, PIPE
import subprocess
import sys,os,re
import glob
import numpy as np
import pandas as pd
from collections import Counter
from collections import defaultdict
import math
from rdkit.Chem import AllChem as ch
from rdkit import DataStructs
from functools import reduce
from rdkit.Chem.rdchem import EditableMol
from rdkit import Chem

def kv_reversal(data_dict):
    # key-value exchange
    kv_list = defaultdict(list)
    for k, v in data_dict.items():
        kv_list[v].append(k)
    kv_list = dict(kv_list)
    return kv_list

def mkdir(path):
    folder = os.path.exists(path)
    if not folder:
        os.makedirs(path)
       
def check_singlefile(singlefile):
    error_flag = 0
    if not os.path.isfile(singlefile):
        print("ERROR: Could not find file %s." % singlefile)
        error_flag = 1
    elif os.path.getsize(singlefile) == 0:
        print("WARNING: File %s size is zero!" % singlefile)
    if error_flag:
        exit(error_flag)
def file_line(filepath):
    count = 0
    for index, line in enumerate(open(filepath, 'r')):
        count += 1
    return count
def resolveConfig(configure_file):
    print(configure_file)
    file = open(configure_file, 'r')
    count_lines = 0
    myDict = {}
    for line in file.readlines():
        line = line.strip()
        count_lines +=1
        if re.match(r'^#.*', line):
            continue
        match_key_value = re.match(r'^\s*([^=\s]+)\s*=\s*(.*)$', line)
        if match_key_value == None:
            continue
        key    = match_key_value.group(1)
        value  = match_key_value.group(2)
        if key == '' or value == '' :
            print("Could not find key or value at line %s.\n" %(count_lines))
            continue
        match_var_value = re.match(r'^\$\((.*)\)(.*)$', value)
        if match_var_value != None:
            variable_key = match_var_value.group(1)
            variable_value = match_var_value.group(2)
            if myDict.has_key(variable_key):
                value = myDict[variable_key] + variable_value
        myDict[key]=value
    return myDict
def read_fasta(file_path):
    """
    Loading FASTA file and return a iterative object
    """
    line = ""
    try:
        fasta_handle = open(file_path, "r")
    except:
        raise IOError("Your input FASTA file is not right!")
    # make sure the file is not empty
    while True:
        line = fasta_handle.readline()
        if line == "":
            return
        if line[0] == ">":
            break
    # when the file is not empty, we try to load FASTA file
    while True:
        if line[0] != ">":
            raise ValueError("Records in Fasta files should start with '>' character")
        title = line[0:].rstrip()
        lines = []
        line = fasta_handle.readline()
        while True:
            if not line:
                break
            if line[0] == ">":
                break
            lines.append(line.rstrip())
            line = fasta_handle.readline()
        yield title, "".join(lines).replace(" ","").replace("\r","")
        if not line:
            return
    fasta_handle.close()

    assert False, "Your input FASTA file have format problem."
def unfoundlib(in_combine_fa, libraryID_all, out_unsplit_fa,cpu, librarynames):
    #libraryID_all dict
    #librarynames (list)
    libs = []
    for lib in librarynames:
        try:
            libs.append(libraryID_all[lib])
        except:
            print('The lib sequence of {} is empty!'.format(lib))
    f_fa = open(out_unsplit_fa, 'w')
    combined_fa_obj = read_fasta(in_combine_fa)
    for header, seq in combined_fa_obj:
        seq_find = [subseq for subseq in libs if subseq in seq]
        if not seq_find:
            f_fa.write(header + '\n')
            f_fa.write(seq + '\n')
def get_count_by_counter(list):
    count = Counter(list)
    count_dict = dict(count)
    return count_dict
def unfoundlib_out(unfound_lib_dir, samplename, unsplit_fa, unsplit_out, lib_out, constant1, constant2, UMI_Number):
    # unfound_lib_dir 为生成的找不到库的seq的序列文件夹，在02.Merge文件夹下
    # samplename 为样本名称
    # unsplit_fa 为后缀为'.combined.unfound.lib.fa'的文件
    UMI_Number = int(UMI_Number)
    #cutadaptDir = config['cutadaptDir']

    cmd = "cutadapt -j 0 -e 0.1 -O 25 -g {} -a {}$ -o {}/{}.temp0.fa {}" \
            .format(constant1, constant2, unfound_lib_dir, samplename, unsplit_fa)
    print(cmd)
    subprocess.check_call(cmd, shell=True)
    cmd = "cutadapt -j 0 -e 0.15 -O 15 \
                    --discard-untrimmed -a {}$ -o {}/{}.temp1.fa {}/{}.temp0.fa ".format( \
       constant2, unfound_lib_dir, samplename, unfound_lib_dir, samplename)
    print(cmd)
    subprocess.check_call(cmd, shell=True)
    lib_file = open(lib_out, 'w')
    lib_list = []
    unsplit_file = open(unsplit_out, 'w')
    cycle_str = ['seq', 'lib', 'UMI']
    unsplit_file.write('\t'.join(cycle_str) + '\n')
    fasta_temp_obj = read_fasta('{}/{}.temp1.fa'.format(unfound_lib_dir, samplename))
    for header, seq in fasta_temp_obj:
        UMI = seq[-UMI_Number:]
        lib = seq[-UMI_Number-10:-UMI_Number]
        if len(seq)-UMI_Number-10 >= 1:
            new_seq = seq[0:len(seq)-UMI_Number-10]
            lib_list.append(lib)
            unsplit_file.write('\t'.join((new_seq, lib, UMI)) + '\n')
    unsplit_file.close()
    lib_count_dict = get_count_by_counter(lib_list)
    for key, value in lib_count_dict.items():
        lib_file.write('\t'.join(([key, str(value)])) + '\n')
    os.remove('{}/{}.temp0.fa'.format(unfound_lib_dir, samplename))
    os.remove('{}/{}.temp1.fa'.format(unfound_lib_dir, samplename))

def total_sample_info(total_info_file):
    # resolve sampleinfo file
    total_info_f = pd.read_csv(total_info_file, sep='\t', dtype=object, header=0)
    total_info_f = total_info_f.dropna(subset=["SampleID"])
    total_info_f.fillna('0')
    #print(total_info_f)
    protein_dict = {} # protein id for each sample
    concentration_dict = {} # protein concentration for each sample
    library_dict = {} # screening librarys for each sample
    positive_dict = {} # names of positive sample for each sample
    I5_dict = {} # I5 for each sample
    I3_dict = {} # I3 for each sample
    #dedup_method_dict = {} # protein id for each sample
    #plot_threshold_dict = {} # protein id for each sample
    roundNumber_dict = {} # protein id for each sample
    ntc_dict = {}
    qpcr_dict = {}
    pair_dict = {}
    pocket_dict = {}
    pocket_diff_dict = {}
    sample_lst = total_info_f['SampleID'].tolist()
    for i in range(total_info_f.shape[0]):
        #print(total_info_f.loc[i, ['Protein']].tolist())
        protein = total_info_f.loc[i, ['Protein']].tolist()[0]
        concentration = total_info_f.loc[i, ['Concentration']].tolist()[0]
        sample = total_info_f.loc[i, ['SampleID']].tolist()[0]
        library = total_info_f.loc[i, ['LibraryID']].tolist()[0]
        positive = total_info_f.loc[i, ['positive_sample']].tolist()[0] 
        I5 = total_info_f.loc[i, ['I5']].tolist()[0]
        I3 = total_info_f.loc[i, ['I3']].tolist()[0]
        #dedup_method = total_info_f.loc[i, ['dedup_method']].tolist()[0]
        #plot_threshold = total_info_f.loc[i, ['plot_threshold']].tolist()[0]
        if 'RoundNM' not in list(total_info_f.keys()):
            print('Error: The sample info file not contain the column RoundNM!')
        if 'NTC' not in list(total_info_f.keys()):
            print('Error: The sample info file not contain the column NTC!')
        if 'QPCR' not in list(total_info_f.keys()):
            print('Error: The sample info file not contain the column QPCR!')
        if 'Pair' not in list(total_info_f.keys()):
            print('Error: The sample info file not contain the column Pair!')
        roundNumber = total_info_f.loc[i, ['RoundNM']].tolist()[0]
        ntc = total_info_f.loc[i, ['NTC']].tolist()[0]
        qpcr = total_info_f.loc[i, ['QPCR']].tolist()[0]
        pair = total_info_f.loc[i, ['Pair']].tolist()[0]
        pocket = total_info_f.loc[i, ['Pocket']].tolist()[0]
        pocket_diff = total_info_f.loc[i, ['Pocket_diff']].tolist()[0]
        protein_dict[sample] = protein 
        concentration_dict[sample] = concentration 
        library_dict[sample] = library 
        positive_dict[sample] = positive
        I5_dict[sample] = I5
        I3_dict[sample] = I3
        #dedup_method_dict[sample] = dedup_method
        #plot_threshold_dict[sample] = plot_threshold
        roundNumber_dict[sample] = roundNumber
        ntc_dict[sample] = ntc
        qpcr_dict[sample] = qpcr
        pair_dict[sample] = pair
        pocket_dict[sample] = pocket
        pocket_diff_dict[sample] = pocket_diff
        
    sampleinfoList = [protein_dict, concentration_dict, library_dict, positive_dict, 
                      I5_dict, I3_dict,sample_lst,roundNumber_dict,ntc_dict,
                      qpcr_dict,pair_dict,pocket_dict,pocket_diff_dict]
    return sampleinfoList
def library_info(library_info_file):
    library_info_f = pd.read_csv(library_info_file, sep='\t', dtype=object, header=[0])
    library_info_f = library_info_f.dropna(subset=["LibraryID"])
    library_info_f.fillna('0')
    PCR1_5_dict = {}
    constant1_dict = {}
    cycles_dict = {}
    library_dict = {}
    CP_umi_dict = {}
    PCR1_3_dict = {}
    positive_dict = {}
    library_size_dict = {}
    library_lst = library_info_f['LibraryID'].tolist()
    for i in range(library_info_f.shape[0]):
        PCR1_5 = library_info_f.loc[i, ['PCR1_5']].tolist()[0]
        constant1 = library_info_f.loc[i, ['constant1']].tolist()[0]
        library = library_info_f.loc[i, ['LibraryID']].tolist()[0]
        cycles = library_info_f.loc[i, ['cycles']].tolist()[0]
        library_seq = library_info_f.loc[i, ['5_3']].tolist()[0]
        CP_umi = library_info_f.loc[i, ['CP_umi']].tolist()[0]
        PCR1_3 = library_info_f.loc[i, ['PCR1_3']].tolist()[0]
        positive_sample = library_info_f.loc[i, ['positive_sample']].tolist()[0]
        librarysize = library_info_f.loc[i, ['LibrarySize']].tolist()[0]
        PCR1_5_dict[library] = PCR1_5
        constant1_dict[library] = constant1
        cycles_dict[library] = cycles
        library_dict[library] = library_seq
        CP_umi_dict[library] = CP_umi
        PCR1_3_dict[library] = PCR1_3
        positive_dict[library] = positive_sample
        library_size_dict[library] = librarysize
    return PCR1_5_dict, constant1_dict, cycles_dict, library_dict, CP_umi_dict, PCR1_3_dict, positive_dict, library_size_dict,library_lst

def reverse_complement(sequence):
    """Reverse Complement Problem
    Find the reverse complement of a DNA string.
    Given: A DNA string Pattern.
    Return: Pattern, the reverse complement of Pattern."""
    complement = ''
    nucleotide_complements = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    for i in range(len(sequence) - 1, -1, -1):
        complement += nucleotide_complements[sequence[i]]
    return complement
def qc_transfer(fileprefix, summary_file):
    qcdict = {}
    head = ""
    qcarray = ["total_reads", "total_bases", "q20_rate", "q30_rate", "gc_content", "corrected_bases",
               "low_quality_reads",
               "too_many_N_reads", "too_short_reads"]
    headarray = []
    # unknown insert size length
    UnknowInsertSizeLen_ratio = 0
    for line in open("{}.fastp.html".format(fileprefix), "r"):
        line = line.strip()
        if UnknowInsertSizeLen_ratio == 0:
            if re.search('Insert size distribution ', line):
                value = re.findall('(\d+\.\d*\%)', line)[0]
                UnknowInsertSizeLen_ratio = value
        else:
            break

    for line in open("{}.fastp.json".format(fileprefix), "r"):
        line = line.strip()
        if re.search("\"before_filtering\"", line):
            head = "raw"
        elif re.search("\"after_filtering\"", line):
            head = "trim"
        elif re.search("filtering_result", line):
            head = "filtered"
        elif re.search("duplication", line):
            head = ""
        if head != "":
            for qc in qcarray:
                if re.search(qc, line):
                    value = re.search(r':\s*(.*?),*$', line).group(1)
                    if float(value) < 1 and float(value) > 0:
                        value = "{:.2f}".format(100.00 * float(value))
                    headarray.append(head + "_" + qc)
                    qcdict[head + "_" + qc] = value
    meanlength = int(int(qcdict["trim_total_bases"]) / int(qcdict["trim_total_reads"]))
    raw_single_reads = int(int(qcdict["raw_total_reads"]) / 2)
    trim_single_reads = int(int(qcdict["trim_total_reads"]) / 2)
    f = open(summary_file, "w")
    f.write('\t'.join([str(raw_single_reads),"raw_single_reads(R1/R2)"])+'\n')
    f.write('\t'.join([str(trim_single_reads),"trim_single_reads(R1/R2)"])+'\n')
    for key in headarray:
        f.write('\t'.join([str(qcdict[key]), str(key)]) + '\n')
    f.write('\t'.join([str(meanlength),"mean_reads_length"])+'\n')
    f.write('\t'.join([str(UnknowInsertSizeLen_ratio),"UnknowInsertSizeLen_ratio"]) +'\n')
    f.close()
def qc_summary(TrimQcFile, casesampleID, ColumnName):
    # ColumnName = 'BatchName'
	with open(TrimQcFile) as f1:
		trim_title = [[]]
		trim_data = []
		for line in f1:
			trim_list = line.split("\t")
			if len(trim_list) != 2:
				continue
			trim_title[0].append(trim_list[1].strip())
			trim_data.append((trim_list[0]))
		trim_total_reads = trim_data[7]
		raw_total_reads = trim_data[2]
		trim_ratio = float(int(trim_total_reads)/int(raw_total_reads) )* 100
		trim_ratio = "%.2f"%trim_ratio
		trim_data.append(trim_ratio)
		trim_title[0].append("trim_ratio(%)")
		trim_title.append(trim_data)
	qc_summary_list = [[ColumnName] + trim_title[0]]
	qc_summary_data = [casesampleID] + trim_title[1]
	qc_summary_list.append(qc_summary_data)
	return (qc_summary_list)

# main run QC.
def run_qc(TrimQcFile, QcSummaryFile, casesampleID, ColumnName):
    if (not os.path.exists(TrimQcFile)):
        print("Error: the file {} is not exists!".format(TrimQcFile))
        sys.exit(1)
    qc_summary_list = qc_summary(TrimQcFile, casesampleID, ColumnName)
    #print(qc_summary_list)
    with open(QcSummaryFile,"w") as f:
        for i in range(len(qc_summary_list)):
            list_str = []
            for data in qc_summary_list[i]:
                data = str(data).strip()
                list_str.append(data)
            line_str = "\t".join(list_str)+"\n"
            f.write(line_str)
    f.close()
def mergefile(filedir,substr,outfile):
    #substr='.summary.qc.txt'
    paths = []
    for root,dirs,files in os.walk(filedir):
        for file in files:
            if substr in file:
                paths.append(os.path.join(root, file))
    #files = os.listdir(filedir)
    #files_new = [filedir+str(i) for i in files if substr in i]
    paths.sort()
    with open(outfile,'w',encoding='gbk') as r:
        n=0
        for i in paths:
            with open(i,'r',encoding='gbk') as f:
                for line in f:
                    if n==0:
                        r.writelines(line)
                        #r.write('\n')
                    else:
                        if 'SampleID' not in line: #跳过含有书名号《的行
                            r.writelines(line)
                            #r.write('\n')
            n+=1
def is_number(s):
    # number is or not
    try:
        float(s)
        return True
    except ValueError:
        pass
    try:
        import unicodedata
        unicodedata.numeric(s)
        return True
    except (TypeError, ValueError):
        pass
    return False
def not_empty(s):
    return s and s.strip()
def compareList(list1,list2):
    new_list1 = list(filter(not_empty, list1))#去除空字符，空格，None
    new_list2 = list(filter(not_empty, list2))
    if sorted(new_list1) == sorted(new_list2):
        return sorted(new_list1)
    else:
        print('Error: The sample in the sample_info file is not same as the barcode file.')

def dedup_filter_2cycle(ratio_file):
    count_f = pd.read_csv(ratio_file, sep='\t', dtype=object, header=0)
    count_f['Dedup_count'] = count_f['Dedup_count'].astype(int)
    dedup_max = count_f['Dedup_count'].max()
    BBmin = 0
    BB1, BB2 = [], []
    dedup_filter_threshold = 1
    # 分析文件中每个循环的BB种类数
    BB1 = set(count_f['TAG1'].tolist())
    BB2 = set(count_f['TAG2'].tolist())
    #BB3 = set(count_f['TAG3'].tolist())
    dedup_max_pre = count_f['Dedup_count'].max()
    # filtering for displaying dot, line, plane
    while (dedup_max >= 1):
        new_dataframe = count_f[count_f['Dedup_count'] >= (dedup_max)]
        BB1 = set(new_dataframe['TAG1'].tolist())
        BB2 = set(new_dataframe['TAG2'].tolist())
        #BB3 = set(new_dataframe['TAG3'].tolist())
        BBmin = min(len(BB1), len(BB2))
        if BBmin < 50:
            dedup_max = dedup_max / 2 + 0.1
        elif BBmin > 100:
            dedup_max = (dedup_max + dedup_max_pre) / 2
            if (dedup_max_pre - dedup_max) <= 1:
                dedup_filter_threshold = dedup_max
                dedup_max = 0.5
            dedup_max_pre = dedup_max
        elif BBmin >= 50:
            dedup_filter_threshold = dedup_max
            dedup_max = 0.5
    return dedup_filter_threshold
def dedup_filter(ratio_file):
    count_f = pd.read_csv(ratio_file, sep='\t', dtype=object, header=0)
    count_f['Dedup_count'] = count_f['Dedup_count'].astype(int)
    dedup_max = count_f['Dedup_count'].max()
    BBmin = 0
    BB1, BB2, BB3 = [], [], []
    #BB_lst = []
    dedup_filter_threshold = 1
    # 分析文件中每个循环的BB种类数
    BB1 = set(count_f['TAG1'].tolist())
    BB2 = set(count_f['TAG2'].tolist())
    BB3 = set(count_f['TAG3'].tolist())
    #BBmin_pre = min(len(BB1), len(BB2), len(BB3))
    dedup_max_pre = count_f['Dedup_count'].max()
    # filtering for displaying dot, line, plane
    while (dedup_max >= 1):
        new_dataframe = count_f[count_f['Dedup_count'] >= (dedup_max)]
        BB1 = set(new_dataframe['TAG1'].tolist())
        BB2 = set(new_dataframe['TAG2'].tolist())
        BB3 = set(new_dataframe['TAG3'].tolist())
        BBmin = min(len(BB1), len(BB2), len(BB3))
        if BBmin < 50:
            dedup_max = dedup_max / 2 + 0.1
        elif BBmin > 100:
            dedup_max = (dedup_max + dedup_max_pre) / 2
            if (dedup_max_pre - dedup_max) <= 1:
                dedup_filter_threshold = dedup_max
                dedup_max = 0.5
            dedup_max_pre = dedup_max
        elif BBmin >= 50:
            dedup_filter_threshold = dedup_max
            dedup_max = 0.5
    return dedup_filter_threshold
def join_r_groups(input_smiles,Starsite='no'):
    # First pass loop over atoms and find the atoms with an AtomMapNum
    str_dict = {}
    try:
        if Starsite=='no':
            pass
        elif Starsite=='yes':
            sites = re.compile(r'\[\*\:\d+\]')
            sites_lst = sites.findall(input_smiles)
            sites_dict = dict(Counter(sites_lst))
            #if max(sites_dict.values())<=2
            keys = []
            for key in sites_dict.keys():
                if sites_dict[key] == 1:
                    keys.append(key)
            singleNm = len(keys)
            atoms_lst = ['[Ge]','[Sn]','[Pb]']
            if singleNm <= 3 and singleNm >=1:
                for i in range(len(keys)):
                    str_i = '.' + atoms_lst[i] + keys[i]
                    input_smiles += (str_i)
                    str_dict[atoms_lst[i]] = keys[i]
        input_mol = Chem.MolFromSmiles(input_smiles)
        #print(Chem.MolToSmiles(input_mol))
        join_dict = defaultdict(list)
        for atom in input_mol.GetAtoms():
            map_num = atom.GetAtomMapNum()
            if map_num > 0:
                join_dict[map_num].append(atom)
        # Second pass, transfer the atom maps to the neighbor atoms
        for idx, atom_list in join_dict.items():
            if len(atom_list) == 2:
                atm_1, atm_2 = atom_list
                nbr_1 = [x.GetOtherAtom(atm_1) for x in atm_1.GetBonds()][0]
                nbr_1.SetAtomMapNum(idx)
                nbr_2 = [x.GetOtherAtom(atm_2) for x in atm_2.GetBonds()][0]
                nbr_2.SetAtomMapNum(idx)
        # Nuke all of the dummy atoms
        new_mol = Chem.DeleteSubstructs(input_mol, Chem.MolFromSmarts('[#0]'))
        #print(Chem.MolToSmiles(new_mol))
        # Third pass - arrange the atoms with AtomMapNum, these will be connected
        bond_join_dict = defaultdict(list)
        for atom in new_mol.GetAtoms():
            map_num = atom.GetAtomMapNum()
            if map_num > 0:
                bond_join_dict[map_num].append(atom.GetIdx())
        #print(bond_join_dict)
        # Make an editable molecule and add bonds between atoms with correspoing AtomMapNum
        em = EditableMol(new_mol)
        for idx, atom_list in bond_join_dict.items():
            if len(atom_list) == 2:
                start_atm, end_atm = atom_list
                em.AddBond(start_atm, end_atm,order=Chem.rdchem.BondType.SINGLE)
        final_mol = em.GetMol()
        # remove the AtomMapNum values
        for atom in final_mol.GetAtoms():
            atom.SetAtomMapNum(0)
        final_mol = Chem.RemoveHs(final_mol)
        final_smiles = Chem.MolToSmiles(final_mol)
        for i in str_dict.keys():
            final_smiles = final_smiles.replace(i,str_dict[i])
        return final_smiles
    except:
        return '0'
def z_score_totallib(x, axis,LibrarySize):
    #axis = 0 表示按列求平均；=1 表示按行求平均
    x = np.array(x,dtype="float32")
    y = np.zeros(int(LibrarySize) - len(x),dtype="float32")
    z=np.append(x,y)
    print([len(z),len(y)])
    xr = np.rollaxis(x, axis=axis)
    xr -= np.mean(z, axis=axis)
    if np.std(z, axis=axis) == 0:
        return np.zeros(len(x))
    else:
        xr /= np.std(z, axis=axis)
        xr[xr<0.0] = 0
        return xr
    
def AggStatistics(file, dedup_count_cutoff, outfileEXCEL, dysythonFile, allsythonFile, outFileFilter,outFileFilterPositive, outFileFilterPositiveNegative, dotNum):
        # 可以指定根据dedup_counts数过滤的cutoff值，file 为未过滤的 ratio file
        # 在流程中可以用过滤后的ratio file
    """
    dotNum 为需要输出的点的个数
    以Tag1=1.314 为例：
首先删除count=1的数据
N: Tag1=1.314 所有组合的个数
SumCount: Tag1=1.314 所有组合的count总和
Max: Tag1=1.314 所有组合的count的最大值
Mean：SumCount除以N
N_Tag2: 当Tag1=1.314, unique Tag2 的个数
N_Tag3: 当Tag1=1.314, unique Tag3 的个数
1.      如果N >= 99% percentile of N (top 1% of N) , Index1=2;
如果95%  <= N < 99% percentile of N, Index1=1;
2.      如果SumCount  >= 99% percentile of SumCount, Index2=2;
如果95% <= SumCount  < 99% percentile of SumCount, Index2=1;
3.      如果Max >= 99% percentile of Max , Index3=2;
如果95% <= Max < 99% percentile of Max, Index3=1;
4.      如果Mean >= 99% percentile of Mean, Index4=2;
如果95% <= Mean < 99% percentile of Mean, Index4=1;
5.      如果N>N_Tag2 or N>N_Tag3, Index5=1; 测序深度很大时，发现此参数意义不大
Sum=Index1, Index2, Index3, Index4, Index5的总和
如果Index5 缺损， 可以根据实际情况直接排除

    """
    # sample_detail = []
# =============================================================================
#     file_prefix = outfileEXCEL.split('/')[-1]
#     sample_detail.append(file_prefix)
# =============================================================================

    writer = pd.ExcelWriter(outfileEXCEL)
    df0 = pd.read_csv(file,sep='\t',header=0,dtype=object)
    df0['Raw_count'] = df0['Raw_count'].astype(int)
    df0['Dedup_count'] = df0['Dedup_count'].astype(int)
    df0 = df0.fillna('0') # 若结构为空，则用“0”来填充
    
    cutoff = 0
    # dedup_count_cutoff = brc.dedup_filter(file)
    
    if int(dedup_count_cutoff) <= 1:
        cutoff = 1.1 # 使dataframe中只保留counts数大于1的数据，这样即使Rawcounts数较大，也会舍去
    else:cutoff = dedup_count_cutoff
    df00 = df0[df0['Dedup_count']>=float(cutoff)] # df00为counts数大于1的dataframe
    #print(df00.columns)
    
    if 'EF' in df00.columns:
        df00['EF'] = df00['EF'].astype(float)
        df00_temp = df00[df00['EF']>1] # df00为counts数大于1的dataframe
        df = pd.DataFrame(df00_temp,columns=['Tags','TAG1','TAG1_structure','TAG1_FF','TAG1_BB',
                                    'TAG2','TAG2_structure','TAG2_FF','TAG2_BB',
                                    'TAG3','TAG3_structure','TAG3_FF','TAG3_BB','Raw_count','Dedup_count',
                                    'z-score','Normalized_z-score','EF','Type'],dtype=object)
    else:
        df = pd.DataFrame(df00,columns=['Tags','TAG1','TAG1_structure','TAG1_FF','TAG1_BB',
                                    'TAG2','TAG2_structure','TAG2_FF','TAG2_BB',
                                    'TAG3','TAG3_structure','TAG3_FF','TAG3_BB','Raw_count','Dedup_count',
                                    'z-score','Normalized_z-score'],dtype=object)
    df['Dedup_count'] = df['Dedup_count'].astype(int)
    df['z-score'] = df['z-score'].astype(float)
    df['Normalized_z-score'] = df['Normalized_z-score'].astype(float)
    #print(df.head())
    if len(df) == 0:
        col_lst = ['Tags', 'TAG1', 'TAG1_structure', 'TAG1_FF', 'TAG1_BB', 'TAG2',
       'TAG2_structure', 'TAG2_FF', 'TAG2_BB', 'TAG3', 'TAG3_structure',
       'TAG3_FF', 'TAG3_BB', 'Raw_count', 'Dedup_count', 'z-score',
       'Normalized_z-score', 'EF', 'Type', 'TAG1_P', 'TAG2_P', 'TAG3_P',
       'TAG1_TAG2_L', 'TAG1_TAG3_L', 'TAG2_TAG3_L', 'TAG_DOT', 'TAG1_P_noline',
       'TAG2_P_noline', 'TAG3_P_noline', 'Y', 'Y0']
        df=pd.DataFrame(columns=col_lst)
        df.to_csv(outFileFilter, sep='\t', index=False) # 过滤后的文件保存为文本文件
        dlpList=['0','0','0','0','0','0']
        return dlpList
    else:
    # =============================================================================
    
        [tag1_sheet,tag1_tags] = tags_excel(df, 'TAG1', 'TAG2', 'TAG3')
        [tag2_sheet,tag2_tags] = tags_excel(df, 'TAG2', 'TAG1', 'TAG3')
        [tag3_sheet,tag3_tags] = tags_excel(df, 'TAG3', 'TAG1', 'TAG2')
        #print(tag1_sheet)
        #print(tag1_sheet.shape[0])
        #print(tag1_sheet.columns)
        # 对mono-synthon的每个tag进行打分（sum赋值）
        if 'EF' in df.columns:
            tag1_add = pd.DataFrame(tag1_sheet[['TAG1','TAG1_structure','TAG1_FF','TAG1_BB','N','SumCount','max','mean','z-score','Normalized_z-score','EF']],dtype=object)
            tag2_add = pd.DataFrame(tag2_sheet[['TAG2','TAG2_structure','TAG2_FF','TAG2_BB','N','SumCount','max','mean','z-score','Normalized_z-score','EF']],dtype=object)
            tag3_add = pd.DataFrame(tag3_sheet[['TAG3','TAG3_structure','TAG3_FF','TAG3_BB','N','SumCount','max','mean','z-score','Normalized_z-score','EF']],dtype=object)
        
            tag1_add = tag1_add.reindex(columns=['TAG1','TAG1_structure','TAG1_FF','TAG1_BB','TAG2','TAG2_structure','TAG2_FF','TAG2_BB','TAG3','TAG3_structure','TAG3_FF','TAG3_BB','N','SumCount','max','mean','z-score','Normalized_z-score','EF'], fill_value='0')
            tag2_add = tag2_add.reindex(columns=['TAG1','TAG1_structure','TAG1_FF','TAG1_BB','TAG2','TAG2_structure','TAG2_FF','TAG2_BB','TAG3','TAG3_structure','TAG3_FF','TAG3_BB','N','SumCount','max','mean','z-score','Normalized_z-score','EF'], fill_value='0')
            tag3_add = tag3_add.reindex(columns=['TAG1','TAG1_structure','TAG1_FF','TAG1_BB','TAG2','TAG2_structure','TAG2_FF','TAG2_BB','TAG3','TAG3_structure','TAG3_FF','TAG3_BB','N','SumCount','max','mean','z-score','Normalized_z-score','EF'], fill_value='0')
        else:
            tag1_add = pd.DataFrame(tag1_sheet[['TAG1','TAG1_structure','TAG1_FF','TAG1_BB','N','SumCount','max','mean','z-score','Normalized_z-score']],dtype=object)
            tag2_add = pd.DataFrame(tag2_sheet[['TAG2','TAG2_structure','TAG2_FF','TAG2_BB','N','SumCount','max','mean','z-score','Normalized_z-score']],dtype=object)
            tag3_add = pd.DataFrame(tag3_sheet[['TAG3','TAG3_structure','TAG3_FF','TAG3_BB','N','SumCount','max','mean','z-score','Normalized_z-score']],dtype=object)
        
            tag1_add = tag1_add.reindex(columns=['TAG1','TAG1_structure','TAG1_FF','TAG1_BB','TAG2','TAG2_structure','TAG2_FF','TAG2_BB','TAG3','TAG3_structure','TAG3_FF','TAG3_BB','N','SumCount','max','mean','z-score','Normalized_z-score'], fill_value='0')
            tag2_add = tag2_add.reindex(columns=['TAG1','TAG1_structure','TAG1_FF','TAG1_BB','TAG2','TAG2_structure','TAG2_FF','TAG2_BB','TAG3','TAG3_structure','TAG3_FF','TAG3_BB','N','SumCount','max','mean','z-score','Normalized_z-score'], fill_value='0')
            tag3_add = tag3_add.reindex(columns=['TAG1','TAG1_structure','TAG1_FF','TAG1_BB','TAG2','TAG2_structure','TAG2_FF','TAG2_BB','TAG3','TAG3_structure','TAG3_FF','TAG3_BB','N','SumCount','max','mean','z-score','Normalized_z-score'], fill_value='0')
    
        #print(tag1_add)
        tags_add = reduce(lambda left,right: pd.DataFrame(pd.merge(left,right,how='outer'),dtype=object), [tag1_add,tag2_add,tag3_add])
        #tags_add.to_csv('test.txt',index=None)
        
        tags_add.fillna(0)
        tags_add[['N','SumCount','max','mean','z-score','Normalized_z-score']] = tags_add[['N','SumCount','max','mean','z-score','Normalized_z-score']].apply(pd.to_numeric)
        if 'EF' in df.columns:
            tags_add[['EF']] = tags_add[['EF']].apply(pd.to_numeric)
        #print(tags_add)
        tag1_sheet_sum = funcNewSumSheet(tag1_sheet,tags_add)
        tag2_sheet_sum = funcNewSumSheet(tag2_sheet,tags_add)
        tag3_sheet_sum = funcNewSumSheet(tag3_sheet,tags_add)
    
        # 计算tag总数
        #planecompounds = tags_add.shape[0]
        
        #tag1_sheet_new = planeLinezscore(tag1_sheet_sum,planecompounds)
        #tag2_sheet_new = planeLinezscore(tag2_sheet_sum,planecompounds)
        #tag3_sheet_new = planeLinezscore(tag3_sheet_sum,planecompounds)
        # print('tag1_sheet_new')
        # print(tag1_sheet_new.keys())
        #print('tag1_sheet_sum: preparing')
        #print(tag1_sheet_sum)
        tag1_sheet_sum.to_excel(writer,'_TAG1', index=None)
        #writer.save()
        tag2_sheet_sum.to_excel(writer,'_TAG2', index=None)
        #writer.save()
        tag3_sheet_sum.to_excel(writer,'_TAG3', index=None)
        #writer.save()
        
        # dy-synthon分析
        #print(df.columns)
        #print(tag1_tags)
        [tag12_sheet,tag12_tags] = dytags_TAG12(tag1_tags, tag2_tags, 'TAG1','TAG2',df) # tag12_tags: top 10 for line;后来改为输出所有线
        [tag13_sheet,tag13_tags] = dytags_TAG12(tag1_tags, tag3_tags, 'TAG1','TAG3',df)
        [tag23_sheet,tag23_tags] = dytags_TAG12(tag2_tags, tag3_tags, 'TAG2','TAG3',df)
        #print(tag12_sheet.keys())
        #print(tag12_tags)
        # 为dy-sython每个二元组合进行打分（sum赋值）
        tag12_add = pd.DataFrame(tag12_sheet,dtype=object)
        tag13_add = pd.DataFrame(tag13_sheet,dtype=object)
        tag23_add = pd.DataFrame(tag23_sheet,dtype=object)
        tag12_add = tag12_add.reindex(columns=['TAG1','TAG1_structure','TAG1_FF','TAG1_BB','TAG2','TAG2_structure','TAG2_FF','TAG2_BB','TAG3','TAG3_structure','TAG3_FF','TAG3_BB','N','SumCount','max','mean','z-score','Normalized_z-score'], fill_value='0')
        tag13_add = tag13_add.reindex(columns=['TAG1','TAG1_structure','TAG1_FF','TAG1_BB','TAG2','TAG2_structure','TAG2_FF','TAG2_BB','TAG3','TAG3_structure','TAG3_FF','TAG3_BB','N','SumCount','max','mean','z-score','Normalized_z-score'], fill_value='0')
        tag23_add = tag23_add.reindex(columns=['TAG1','TAG1_structure','TAG1_FF','TAG1_BB','TAG2','TAG2_structure','TAG2_FF','TAG2_BB','TAG3','TAG3_structure','TAG3_FF','TAG3_BB','N','SumCount','max','mean','z-score','Normalized_z-score'], fill_value='0')
        tag12_add['TAG3'] = '3.0'
        tag13_add['TAG2'] = '2.0'
        tag23_add['TAG1'] = '1.0'
        tagsdy_add = reduce(lambda left,right: pd.DataFrame(pd.merge(left,right,how='outer'),dtype=object), [tag12_add,tag13_add,tag23_add])
        tagsdy_add.fillna(0)
        tagsdy_add[['N','SumCount','max','mean','z-score','Normalized_z-score']] = tagsdy_add[['N','SumCount','max','mean','z-score','Normalized_z-score']].apply(pd.to_numeric)
        #print(tag12_sheet)
        tag12_sheet_sum = funcNewSumSheet(tag12_sheet,tagsdy_add)
        tag13_sheet_sum = funcNewSumSheet(tag13_sheet,tagsdy_add)
        tag23_sheet_sum = funcNewSumSheet(tag23_sheet,tagsdy_add)
        #linecompounds = tagsdy_add.shape[0]
        #tag12_sheet_new = planeLinezscore(tag12_sheet_sum,linecompounds)
        #tag13_sheet_new = planeLinezscore(tag13_sheet_sum,linecompounds)
        #tag23_sheet_new = planeLinezscore(tag23_sheet_sum,linecompounds)
    
        tag12_sheet_sum.to_excel(writer,'_TAG1_TAG2', index=None)
        #writer.save()
        tag13_sheet_sum.to_excel(writer,'_TAG1_TAG3', index=None)
        #writer.save()
        tag23_sheet_sum.to_excel(writer,'_TAG2_TAG3', index=None)
        #writer.save()
        
        #print(tag12_sheet_sum)
    
        # tag1_tags, tag2_tags, tag3_tags 分别为排名前10的tag编号
        
        
        df['TAG1_P'] = str(0)
        df['TAG2_P'] = str(0)
        df['TAG3_P'] = str(0)
        
        df['TAG1_TAG2_L'] = str(0)
        df['TAG1_TAG3_L'] = str(0)
        df['TAG2_TAG3_L'] = str(0)
        
        df['TAG_DOT'] = str(0)
        dedup_mean = float('%.2f'%(df['Dedup_count'].mean()))
        # 判断df是否为空
        if math.isnan(dedup_mean):
            dedup_mean = '0'
        BB1 = tag1_sheet_sum['TAG1'].tolist()
        BB2 = tag2_sheet_sum['TAG2'].tolist()
        BB3 = tag3_sheet_sum['TAG3'].tolist()
    # =============================================================================
    #     BB1_len = len(BB1)
    #     BB2_len = len(BB2)
    #     BB3_len = len(BB3)
    # =============================================================================
        
        # identify big dot
    
        [df000,DotList] = dot_judge(df,dotNum,str(dedup_mean))
        #print(df000)
        # identify plane
        [df01,plane1List,plane1List_noline] = plane_judge(BB1, df000, tag1_sheet_sum, 'TAG1', 'TAG2', 'TAG3')
        [df02,plane2List,plane2List_noline] = plane_judge(BB2, df01, tag2_sheet_sum, 'TAG2', 'TAG1', 'TAG3')
        [df03,plane3List,plane3List_noline] = plane_judge(BB3, df02, tag3_sheet_sum, 'TAG3', 'TAG1', 'TAG2')
        # identify line
        [df04,line12List] = line_judge(df03, tag12_sheet_sum, 'TAG1', 'TAG2', 'TAG1_TAG2_L')
    
        [df05,line13List] = line_judge(df04, tag13_sheet_sum, 'TAG1', 'TAG3', 'TAG1_TAG3_L')
    
        [df06,line23List] = line_judge(df05, tag23_sheet_sum, 'TAG2', 'TAG3', 'TAG2_TAG3_L')
    
        #print('df06')
    
       # Y=1表示点线面阳性信号；Y0=1表示不存在线的面的阳性信号
        df06['Y'] = 0
        df06['Y0'] = 0
        df06['TAG1_P'] = df06['TAG1_P'].astype(int)
        df06['TAG2_P'] = df06['TAG2_P'].astype(int)
        df06['TAG3_P'] = df06['TAG3_P'].astype(int)
        df06['TAG1_TAG2_L'] = df06['TAG1_TAG2_L'].astype(int)
        df06['TAG1_TAG3_L'] = df06['TAG1_TAG3_L'].astype(int)
        df06['TAG2_TAG3_L'] = df06['TAG2_TAG3_L'].astype(int)
        df06['TAG_DOT'] = df06['TAG_DOT'].astype(int)
        
        df06.loc[(df06['TAG1_P'] >= 2) | (df06['TAG2_P'] >=2) | (df06['TAG3_P'] >= 2) | 
                (df06['TAG1_TAG2_L'] >= 2) | (df06['TAG1_TAG3_L'] >= 2) | (df06['TAG2_TAG3_L'] >= 2) |
                (df06['TAG_DOT'] >= 1), 'Y'] = str(1)
    
        df06.loc[(df06['TAG1_P_noline'] == '1') | (df06['TAG2_P_noline'] == '1') | (df06['TAG3_P_noline'] == '1'), 'Y0'] = str(1)
        #print(df06)
        df06.to_csv(outFileFilter, sep='\t', index=False) # 过滤后的文件保存为文本文件
        df06.to_excel(writer,'_count.filter', index=None) # 过滤后的文件也保存在同一个excel表格中
        writer.save()
        #writer.close()
        
        # 输出带有标注的阳性样本
    # =============================================================================
        if 'EF' in df06.columns:
            df07 = pd.DataFrame(df06[['Tags','TAG1','TAG1_structure','TAG1_FF','TAG1_BB',
                                  'TAG2','TAG2_structure','TAG2_FF','TAG2_BB','TAG3',
                                  'TAG3_structure','TAG3_FF','TAG3_BB','Dedup_count',
                                  'z-score','Normalized_z-score','EF','Y']],dtype=object)
        else:
            df07 = pd.DataFrame(df06[['Tags','TAG1','TAG1_structure','TAG1_FF','TAG1_BB',
                                  'TAG2','TAG2_structure','TAG2_FF','TAG2_BB','TAG3',
                                  'TAG3_structure','TAG3_FF','TAG3_BB','Dedup_count',
                                  'z-score','Normalized_z-score','Y']],dtype=object)
        df08 = df07[df07['Y'] == '1']
        df08.to_csv(outFileFilterPositive, sep='\t',index=False)
        df07.to_csv(outFileFilterPositiveNegative, sep='\t',index=False)
    # =============================================================================
        # 若存在大点，则只选dedup counts数排名前20的点TOP20LIST返回的是字符串
        #DotListFilter = ALLLIST(DotList) # 输出所有点
        DotListFilter = TOP20LIST(DotList)
        # 若存在线，则只输出前20条线 调整成输出所有线
        line12ListFilter = ALLLIST(line12List)
        line13ListFilter = ALLLIST(line13List)
        line23ListFilter = ALLLIST(line23List)
        # 若存在面，则只输出前20个面 调整成输出所有面
        plane1ListFilter = ALLLIST(plane1List)
        plane2ListFilter = ALLLIST(plane2List)
        plane3ListFilter = ALLLIST(plane3List)
        
        plane1ListFilterNoline = ALLLIST(plane1List_noline)
        plane2ListFilterNoline = ALLLIST(plane2List_noline)
        plane3ListFilterNoline = ALLLIST(plane3List_noline)
        
        lineList = [line12ListFilter,line13ListFilter,line23ListFilter]
        lineList = [x for x in lineList if x != '0']
        lineList2str = LIST2STR(lineList)
        
        planeList = [plane1ListFilter,plane2ListFilter,plane3ListFilter]
        planeList = [x for x in planeList if x != '0']
        planeList2str = LIST2STR(planeList)
        
        # 面里无线的情况先输出，作为备选项
        planeListNoline = [plane1ListFilterNoline,plane2ListFilterNoline,plane3ListFilterNoline]
        planeListNoline = [x for x in planeListNoline if x != '0']
        planeList2strNoline = LIST2STR(planeListNoline)
        
        DLPSTR = '0'
        Pnoline = '0'
        
        if (planeList2str != '0' and lineList2str != '0' and DotListFilter != '0'):
            DLPSTR = 'DLP'
        if (planeList2str != '0' and lineList2str == '0' and DotListFilter == '0'):
            DLPSTR = 'P'
        if (planeList2str == '0' and lineList2str != '0' and DotListFilter == '0'):
            DLPSTR = 'L'
        if (planeList2str == '0' and lineList2str == '0' and DotListFilter != '0'):
            DLPSTR = 'D'
        if (planeList2str != '0' and lineList2str != '0' and DotListFilter == '0'):
            DLPSTR = 'LP'
        if (planeList2str != '0' and lineList2str == '0' and DotListFilter != '0'):
            DLPSTR = 'DP'
        if (planeList2str == '0' and lineList2str != '0' and DotListFilter != '0'):
            DLPSTR = 'DL'
        if planeList2strNoline != '0':
            Pnoline = '1'
        
        
        dlpList = [planeList2str,lineList2str,DotListFilter,DLPSTR,planeList2strNoline,Pnoline]
            
    # =============================================================================
        # line
        tag12_sheet_sum['TAG3'] = '0'
        tag12_sheet_sum['TAG3_structure'] = '0'
        tag12_sheet_sum['TAG3_FF'] = '0'
        tag12_sheet_sum['TAG3_BB'] = '0'
        
        tag13_sheet_sum['TAG2'] = '0'
        tag13_sheet_sum['TAG2_structure'] = '0'
        tag13_sheet_sum['TAG2_FF'] = '0'
        tag13_sheet_sum['TAG2_BB'] = '0'
        
        tag23_sheet_sum['TAG1'] = '0'
        tag23_sheet_sum['TAG1_structure'] = '0'
        tag23_sheet_sum['TAG1_FF'] = '0'
        tag23_sheet_sum['TAG1_BB'] = '0'
        if 'EF' in tag12_sheet_sum.columns:
            tag12_sheet_sum = tag12_sheet_sum[['TAG1','TAG1_structure','TAG1_FF','TAG1_BB','TAG2','TAG2_structure','TAG2_FF','TAG2_BB','TAG3',
                                           'TAG3_structure','TAG3_FF','TAG3_BB','N','SumCount','max','mean','Sum','z-score','Normalized_z-score','EF']]
            tag13_sheet_sum = tag13_sheet_sum[['TAG1','TAG1_structure','TAG1_FF','TAG1_BB','TAG2','TAG2_structure','TAG2_FF','TAG2_BB','TAG3',
                                           'TAG3_structure','TAG3_FF','TAG3_BB','N','SumCount','max','mean','Sum','z-score','Normalized_z-score','EF']]
            tag23_sheet_sum = tag23_sheet_sum[['TAG1','TAG1_structure','TAG1_FF','TAG1_BB','TAG2','TAG2_structure','TAG2_FF','TAG2_BB','TAG3',
                                           'TAG3_structure','TAG3_FF','TAG3_BB','N','SumCount','max','mean','Sum','z-score','Normalized_z-score','EF']]
        else:
            tag12_sheet_sum = tag12_sheet_sum[['TAG1','TAG1_structure','TAG1_FF','TAG1_BB','TAG2','TAG2_structure','TAG2_FF','TAG2_BB','TAG3',
                                           'TAG3_structure','TAG3_FF','TAG3_BB','N','SumCount','max','mean','Sum','z-score','Normalized_z-score']]
            tag13_sheet_sum = tag13_sheet_sum[['TAG1','TAG1_structure','TAG1_FF','TAG1_BB','TAG2','TAG2_structure','TAG2_FF','TAG2_BB','TAG3',
                                           'TAG3_structure','TAG3_FF','TAG3_BB','N','SumCount','max','mean','Sum','z-score','Normalized_z-score']]
            tag23_sheet_sum = tag23_sheet_sum[['TAG1','TAG1_structure','TAG1_FF','TAG1_BB','TAG2','TAG2_structure','TAG2_FF','TAG2_BB','TAG3',
                                           'TAG3_structure','TAG3_FF','TAG3_BB','N','SumCount','max','mean','Sum','z-score','Normalized_z-score']]
        tag12_sheet_sum = pd.DataFrame(tag12_sheet_sum,dtype=object)
        tag13_sheet_sum = pd.DataFrame(tag13_sheet_sum,dtype=object)
        tag23_sheet_sum = pd.DataFrame(tag23_sheet_sum,dtype=object)
        #print('tag12_sheet_new2')
        #print(tag12_sheet_new.keys())
        
        dysythons_df = reduce(lambda left,right: pd.DataFrame(pd.merge(left,right,how='outer'),dtype=object), [tag12_sheet_sum,tag13_sheet_sum,tag23_sheet_sum])
        #print('dysythons_df')
        #print(dysythons_df.keys())
        dysythons_df['TAG1'].fillna('1.0')
        dysythons_df['TAG2'].fillna('2.0')
        dysythons_df['TAG3'].fillna('3.0')
        dysythons_df.fillna(0)
        dysythons_df_new = dysythons_df.sort_values(by=['Sum','SumCount','N', 'max', 'mean'],ascending=False)
        dysythons_df_new.reset_index(drop=True,inplace=True)
        dysythons_df_new.to_csv(dysythonFile, sep='\t', index=False) # 过滤后的文件保存为文本文件
        
        # plane
        
        try:
            
            tag1_sheet_sum.drop(['N_TAG2','N_TAG3','TAG1_similarity'],axis=1, inplace=True)
            tag2_sheet_sum.drop(['N_TAG1','N_TAG3','TAG2_similarity'],axis=1, inplace=True)
            tag3_sheet_sum.drop(['N_TAG1','N_TAG2','TAG3_similarity'],axis=1, inplace=True)
        except:
            print("['N_TAG2' 'N_TAG3'] not found in axis")
        
        tag1_sheet_sum.reset_index(drop=True, inplace=True)
        tag2_sheet_sum.reset_index(drop=True, inplace=True)
        tag3_sheet_sum.reset_index(drop=True, inplace=True)
        
        # dot
        df06.rename(columns={'Dedup_count':'SumCount'}, inplace=True)
        if 'EF' in df06.columns:
            df_dot = df06[['Tags','TAG1','TAG1_structure','TAG1_FF','TAG1_BB','TAG2','TAG2_structure','TAG2_FF','TAG2_BB','TAG3','TAG3_structure','TAG3_FF','TAG3_BB','SumCount','z-score','Normalized_z-score','EF']]
        else:
            df_dot = df06[['Tags','TAG1','TAG1_structure','TAG1_FF','TAG1_BB','TAG2','TAG2_structure','TAG2_FF','TAG2_BB','TAG3','TAG3_structure','TAG3_FF','TAG3_BB','SumCount','z-score','Normalized_z-score']]
        # merge dlp
        
        
        tag1_sheet_sum = pd.DataFrame(tag1_sheet_sum,dtype=object)
        tag2_sheet_sum = pd.DataFrame(tag2_sheet_sum,dtype=object)
        tag3_sheet_sum = pd.DataFrame(tag3_sheet_sum,dtype=object)
        df_dot = pd.DataFrame(df_dot,dtype=object)
        
        tag12_sheet_sum['TAG3'] = '3.0'
        tag13_sheet_sum['TAG2'] = '2.0'
        tag23_sheet_sum['TAG1'] = '1.0'
        tag1_sheet_sum['TAG2'] = '2.0'
        tag1_sheet_sum['TAG3'] = '3.0'
        tag2_sheet_sum['TAG1'] = '1.0'
        tag2_sheet_sum['TAG3'] = '3.0'
        tag3_sheet_sum['TAG1'] = '1.0'
        tag3_sheet_sum['TAG2'] = '2.0'
        
        df_dot_line_plane = reduce(lambda left,right: pd.DataFrame(pd.merge(left,right,how='outer'),
                                                                   dtype=object),[tag12_sheet_sum,tag13_sheet_sum,tag23_sheet_sum, tag1_sheet_sum,tag2_sheet_sum,tag3_sheet_sum,df_dot])
        
        df_dot_line_plane = df_dot_line_plane.fillna('0')
        df_dot_line_plane.reset_index(drop=True, inplace=True)
        df_dot_line_plane.to_csv(allsythonFile, sep='\t', index=False)
    # =============================================================================
    #     else:
    #         dlpList=['0','0','0','0','0','0']
    #         df_count_filter = pd.DataFrame(columns=['TAG1','TAG1_structure','TAG2','TAG2_structure','TAG3','TAG3_structure','Raw_count','Dedup_count',
    #                                                 'TAG1_P','TAG2_P','TAG3_P','TAG1_TAG2_L',
    #                                                 'TAG1_TAG3_L','TAG2_TAG3_L','TAG_DOT',
    #                                                 'TAG1_P_noline','TAG2_P_noline','TAG3_P_noline','Y','Y0'])
    #         df_count_filter.to_csv(outFileFilter, sep='\t', index=False) # 过滤后的文件保存为文本文件
    # =============================================================================
    
        #print(dlpList)
    
        
        #return sample_detail,dlpList
        return dlpList

def LIST2STR(taglist):
    # 输出列表的前20个元素，并以分号分割，将列表转化为字符串
    tagliststr = '0'
    if taglist:
        tagliststr = ';'.join(taglist)
    else:tagliststr = '0'
    return tagliststr

def TOP20LIST(taglist):
    # 输出列表的前20个元素，并以分号分割，将列表转化为字符串
    if taglist:
        
        if len(taglist) <= 20:
            taglistFilter = ';'.join(taglist)
        else:
            taglistFilter = ';'.join(taglist[:20])
    else:taglistFilter = '0'
    
    return taglistFilter
def ALLLIST(taglist):
    # 输出列表的所有元素，并以分号分割，将列表转化为字符串
    if taglist:
        taglistFilter = ';'.join(taglist)
    else:taglistFilter = '0'
    
    return taglistFilter

def dot_judge(df00,dotNum,dedup_mean):
    # df00 为dedup_counts大于cutoff的dataframe
    # dotNum 为取前多少个点，如20
    tag123_lst = []
    dedup_mean = float(dedup_mean)
    df00.loc[(df00['Dedup_count'] >= 10 * dedup_mean), 'TAG_DOT'] = str(2)  # 判断清晰的点       

    df00.loc[(df00['Dedup_count'] < 10 * dedup_mean) & (df00['Dedup_count'] >= 5 * dedup_mean), 'TAG_DOT'] = str(1)
    dotDF = df00[df00['Dedup_count'] >= 5 * dedup_mean]
    n = 0
    for row in dotDF.iterrows():
        n = n+1
        if n<=int(dotNum):
            str_tag1 = row[1]['TAG1']
            str_tag2 = row[1]['TAG2']
            str_tag3 = row[1]['TAG3']
            str_tag123 = str_tag1 + ',' + str_tag2 + ',' + str_tag3
            tag123_lst.append(str_tag123)
    return df00, tag123_lst

    
def line_judge(df00, tag12_sheet, CombineTAG1, CombineTAG2, TAG12_L):
    # df00 为dedup_counts大于cutoff的dataframe
    # tag12_sheet 为cycle1和cycle2的组合BBlist排序列表
    # CombineTAG1 = 'TAG1'
    # CombineTAG2 = 'TAG2'
    # TAG12_L == 'TAG12_L'
    tag12_lst = []
    #mm_domain = np.max(list(tag12_sheet['Sum'].values)) - np.min(list(tag12_sheet['Sum'].values))
    tag12_sheet_filter = tag12_sheet.loc[(tag12_sheet['N']>=3) & (tag12_sheet['max']>=2)&(tag12_sheet['Sum']>=6)]
    if tag12_sheet_filter.empty:
        pass
    else:
        
        for row in tag12_sheet_filter.iterrows():
            str_tag1 = str(row[1][CombineTAG1])
            str_tag2 = str(row[1][CombineTAG2])
            str_tag12 = str_tag1 + ',' + str_tag2
            # 依据下面的条件可以适当的放松或收紧线的输出的条件
            if (row[1]['Sum'] >= 6.5) and (row[1]['max'] >= 5) and (row[1]['N'] >= 10 or row[1]['SumCount'] >= 16):
                tag12_lst.append(str_tag12)
                df00.loc[(df00[CombineTAG1] == str_tag1) & (df00[CombineTAG2] == str_tag2), TAG12_L] = str(5)
            elif (row[1]['Sum'] >= 8) and (row[1]['max'] >= 3) and row[1]['N'] >= 4:
                tag12_lst.append(str_tag12)
                df00.loc[(df00[CombineTAG1] == str_tag1) & (df00[CombineTAG2] == str_tag2), TAG12_L] = str(4)
            elif (row[1]['Sum'] >= 8) and (row[1]['max'] >= 2) and row[1]['N'] >= 5:
                tag12_lst.append(str_tag12)
                df00.loc[(df00[CombineTAG1] == str_tag1) & (df00[CombineTAG2] == str_tag2), TAG12_L] = str(3)
            elif (row[1]['Sum'] >= 8.5 and row[1]['N'] >= 3):
                tag12_lst.append(str_tag12)
                df00.loc[(df00[CombineTAG1] == str_tag1) & (df00[CombineTAG2] == str_tag2), TAG12_L] = str(2)
            elif (row[1]['Sum'] >= 8):
                df00.loc[(df00[CombineTAG1] == str_tag1) & (df00[CombineTAG2] == str_tag2), TAG12_L] = str(1)
# =============================================================================
#             elif (mm_domain<2.5) and (row[1]['Sum'] >= 8.5):
#                 tag12_lst.append(str_tag12)
#                 df00.loc[(df00[CombineTAG1] == str_tag1) & (df00[CombineTAG2] == str_tag2), TAG12_L] = str(2)
#                 
# =============================================================================
    return df00, tag12_lst

def plane_judge(BB1, df00, tag1_sheet, targetTAG, tag2, tag3):
    # BB1为BB1的tag list
    # df00 为dedup_counts大于cutoff的dataframe
    # tag1_sheet为cycle1的list列表
    # targetTAG为遍历的哪个循环的tag，如‘TAG1’
    # cycleth 为第几个循环的tag，如1
    # tag2 为另外一个循环的tag，如‘TAG2’
    # tag3 为第三个循环的tag,如‘TAG3’
    N_TAG2 = 'N_' + tag2
    N_TAG3 = 'N_' + tag3
    TAG1_P = targetTAG + '_P'
    TAG1_P_noline = targetTAG + '_P_noline'
    df00[TAG1_P_noline]=str(0)
    df00[TAG1_P]=str(0)
    
    TAG1_P_noline_lst = []
    TAG1_P_lst = []
    tag1_N_dict = tag1_sheet.set_index(targetTAG)['N']
    tag1_N_TAG2_dict = tag1_sheet.set_index(targetTAG)[N_TAG2]
    tag1_N_TAG3_dict = tag1_sheet.set_index(targetTAG)[N_TAG3]
    # tag1_SumCount_dict = tag1_sheet.set_index(targetTAG)['SumCount']
    tag1_Max_dict = tag1_sheet.set_index(targetTAG)['max']
    # tag1_Mean_dict = tag1_sheet.set_index(targetTAG)['mean']
    tag1_Sum_dict = tag1_sheet.set_index(targetTAG)['Sum']
    
    #最后打分的最大值和最小值如果太小，说明数据不太能区分
    #mm_domain = np.max(list(tag1_sheet['Sum'].values)) - np.min(list(tag1_sheet['Sum'].values))
    # plane
    for Tag1 in BB1:
        #df00_tag1 = df00[df00[targetTAG] == Tag1]
        tag1_N = tag1_N_dict[Tag1] # 同一个BB1的行数，即
        tag1_N_TAG2 = tag1_N_TAG2_dict[Tag1]
        tag1_N_TAG3 = tag1_N_TAG3_dict[Tag1]
        # tag1_SumCount = tag1_SumCount_dict[Tag1]
        tag1_Max = tag1_Max_dict[Tag1]
        # tag1_Mean = tag1_Mean_dict[Tag1]
        tag1_Sum = tag1_Sum_dict[Tag1]
        df00.loc[df00[targetTAG] == Tag1, TAG1_P_noline] = str(0)
        

        if (tag1_Sum >= 5) and (tag1_N >= 4) and (tag1_Max >=2) and (min(tag1_N_TAG2,tag1_N_TAG3)>=2):
            if tag1_N >= 10 and tag1_Max >=5 and tag1_Sum >=8 and ((tag1_N > tag1_N_TAG2) or (tag1_N > tag1_N_TAG3)):
                df00.loc[df00[targetTAG] == Tag1, TAG1_P] = str(5)
                TAG1_P_lst.append(Tag1)
            
            elif tag1_N >= 5 and tag1_Max >=3 and tag1_Sum >=8 and ((tag1_N > tag1_N_TAG2) or (tag1_N > tag1_N_TAG3)):
                df00.loc[df00[targetTAG] == Tag1, TAG1_P] = str(4)
                TAG1_P_lst.append(Tag1) # 如果要放松条件，则在此条件下可将tag信息加入到列表中
            elif tag1_N >= 5 and tag1_Sum >=8 and ((tag1_N > tag1_N_TAG2) or (tag1_N > tag1_N_TAG3)):
                df00.loc[df00[targetTAG] == Tag1, TAG1_P] = str(3)
                TAG1_P_lst.append(Tag1)
            elif tag1_Sum >=8.5 and ((tag1_N > tag1_N_TAG2) or (tag1_N > tag1_N_TAG3)):
                df00.loc[df00[targetTAG] == Tag1, TAG1_P] = str(2)
                TAG1_P_lst.append(Tag1)
            elif tag1_Sum >=8 and ((tag1_N > tag1_N_TAG2) or (tag1_N > tag1_N_TAG3)):
                df00.loc[df00[targetTAG] == Tag1, TAG1_P] = str(1)
            elif tag1_Sum >=8.5 and ((tag1_N == tag1_N_TAG2) and (tag1_N == tag1_N_TAG3)):
                TAG1_P_noline_lst.append(Tag1)
                df00.loc[df00[targetTAG] == Tag1, TAG1_P_noline] = str(1)

        else:
            df00.loc[df00[targetTAG] == Tag1, TAG1_P] = str(0)
            #df00.loc[df00[targetTAG] == Tag1, TAG1_P_noline] = str(0)
    return df00,TAG1_P_lst,TAG1_P_noline_lst

def dytags_TAG12(tag1, tag2, tag1str,tag2str,df):
    # tag1str = 'TAG1'; tag2str = 'TAG2'
    N = [];SumCount = [];DataMax = [];DataMean = [];TAG1 = [];TAG1_structure = []
    TAG2 = [];TAG2_structure = [];TAG1_BB = []; TAG2_BB = [];TAG1_FF=[];TAG2_FF=[]
    zscore = [];normal_zscore = [];ef=[]
    tag1_structure = tag1str + '_structure'
    tag2_structure = tag2str + '_structure'
    tag1_bb = tag1str + '_BB'
    tag2_bb = tag2str + '_BB'
    tag1_ff = tag1str + '_FF'
    tag2_ff = tag2str + '_FF'
# =============================================================================
#     tag1_seq = tag1str + '_seq'
#     tag2_seq = tag2str + '_seq'
# =============================================================================
    for eachTag1 in tag1:
        for eachTag2 in tag2:
            df2 = df.loc[(df[tag1str] == eachTag1) & (df[tag2str] == eachTag2)]
            if df2.empty:
                continue
            else:               
                TAG1.append(eachTag1)
                TAG2.append(eachTag2)
                if list(set(df2[tag1_structure].tolist())):
                    TAG1_structure.append(list(set(df2[tag1_structure].tolist()))[0])
                    TAG1_BB.append(list(set(df2[tag1_bb].tolist()))[0])
                    TAG1_FF.append(list(set(df2[tag1_ff].tolist()))[0])
                else:
                    TAG1_structure.append('0')
                    TAG1_BB.append('0')
                    TAG1_FF.append('0')
# =============================================================================
#                 if list(set(df2[tag1_seq].tolist())):
#                     TAG1_seq.append(list(set(df2[tag1_seq].tolist()))[0])
#                 else:
#                     TAG1_seq.append('NA')
# =============================================================================
                
                if list(set(df2[tag2_structure].tolist())):
                    TAG2_structure.append(list(set(df2[tag2_structure].tolist()))[0])
                    TAG2_BB.append(list(set(df2[tag2_bb].tolist()))[0])
                    TAG2_FF.append(list(set(df2[tag2_ff].tolist()))[0])
                else:
                    TAG2_structure.append('0')
                    TAG2_BB.append('0')
                    TAG2_FF.append('0')
# =============================================================================
#                 if list(set(df2[tag2_seq].tolist())):
#                     
#                     TAG2_seq.append(list(set(df2[tag2_seq].tolist()))[0])
#                 else:TAG2_seq.append('NA')
# =============================================================================
            
                N.append(int(len(df2)))
                SumCount.append(int(df2['Dedup_count'].sum()))
                DataMax.append(df2['Dedup_count'].max())
                DataMean.append(df2['Dedup_count'].mean())
                zscore.append(df2['z-score'].sum())
                normal_zscore.append(df2['Normalized_z-score'].sum())
                if 'EF' in df2.columns:
                    ef.append(df2['EF'].sum())
    if 'EF' in df.columns:
        sheet1 = pd.DataFrame({tag1str:TAG1,tag1_structure:TAG1_structure,tag1_ff:TAG1_FF,tag1_bb:TAG1_BB,
                           tag2str:TAG2,tag2_structure:TAG2_structure,tag2_ff:TAG2_FF,tag2_bb:TAG2_BB,
                           'N':N,'SumCount':SumCount,'max':DataMax,'mean':DataMean,
                           'z-score':zscore, 'Normalized_z-score':normal_zscore, 'EF':ef},dtype=object)
    else:
        sheet1 = pd.DataFrame({tag1str:TAG1,tag1_structure:TAG1_structure,tag1_ff:TAG1_FF,tag1_bb:TAG1_BB,
                           tag2str:TAG2,tag2_structure:TAG2_structure,tag2_ff:TAG2_FF,tag2_bb:TAG2_BB,
                           'N':N,'SumCount':SumCount,'max':DataMax,'mean':DataMean,
                           'z-score':zscore, 'Normalized_z-score':normal_zscore},dtype=object)
    sheet1['N'] = sheet1['N'].astype('int')
    sheet1['SumCount'] = sheet1['SumCount'].astype('int')
    sheet1['max'] = sheet1['max'].astype('int')
    sheet1['mean'] = np.round(sheet1['mean'].astype('float'),2)
    
    #sheet1_new = funcNewSheet(sheet1)#sheet转换成带有统计信息的sheet
    sheet1_new = sheet1
    tag1tag2 = '_' + tag1str + '_' + tag2str
    #sheet1_new.to_excel(writer,tag1tag2, index=None)
    
    #single_lst = [str('%.1f'%(dict1_N['99%'])),str('%.1f'%(dict1_SumCount['99%'])),
    #              str('%.1f'%(dict1_DataMax['99%'])),str('%.1f'%(dict1_DataMean['99%']))]
    
    # sheet1_new[tag1tag2] = sheet1_new[tag1str] + '|' + sheet1_new[tag1_structure] + ',' + sheet1_new[tag2str] + '|' + sheet1_new[tag2_structure]
    # sheet1_new_tag12 = sheet1_new
    sheet1_new[tag1tag2] = sheet1_new[tag1str].map(str) + ',' + sheet1_new[tag2str].map(str)
    
    
    # print(sheet1_new.keys)
    
    #tags12_top10 = list(sheet1_new.loc[:9,tag1tag2]) # return the top 10 of sum value for the tags list 
    tags12_top10 = list(sheet1_new.loc[:,tag1tag2])
    #return sheet1_new, single_lst, tags12_top10
    try:
        sheet1_new.drop([tag1tag2],axis=1, inplace=True)
    except:
        print('No special columns: %s!'.format(tag1tag2))
    return sheet1_new, tags12_top10
    
def tags_excel(df, tags1,tags2,tags3):
    # tags1='Tag1'
    # tags2='Tag2'
    # tags3='Tag3'
    N = [];SumCount = [];DataMax = [];DataMean = [];N_Tag2 = [];N_Tag3 = [];TAG1 = []
    zscore = [];normal_zscore=[];ef=[]
# =============================================================================
#     TAG1_seq = []
# =============================================================================
    TAG1_structure = [];TAG1_structurebb = [];TAG1_similarity = [];TAG1_FF = []
# =============================================================================
#     seq = tags1 + '_seq'
# =============================================================================
    structure = tags1 + '_structure'
    structurebb = tags1 + '_BB'
    structureff = tags1 + '_FF'
    similarityTAG = tags1 + '_similarity'
    
    for tag1 in list(set(df[tags1].tolist())):
        df1 = df[df[tags1]==tag1]
        TAG1.append(tag1)
        TAG1_similarity.append(1)
        
        if list(set(df1[structure].tolist()))[0]:
            TAG1_structure.append(list(set(df1[structure].tolist()))[0])
            TAG1_structurebb.append(list(set(df1[structurebb].tolist()))[0])
            TAG1_FF.append(list(set(df1[structureff].tolist()))[0])
        else:
            TAG1_structure.append('0')
            TAG1_structurebb.append('0')
            TAG1_FF.append('0')
# =============================================================================
#         TAG1_seq.append(list(set(df1[seq].tolist()))[0])
# =============================================================================
        #print(df1.head())
        N.append(int(len(df1)))
        SumCount.append(int(df1['Dedup_count'].sum()))
        DataMax.append(int(df1['Dedup_count'].max()))
        DataMean.append(df1['Dedup_count'].mean())
        zscore.append(df1['z-score'].sum())
        normal_zscore.append(df1['Normalized_z-score'].sum())
        N_Tag2.append(len(list(set(df1[tags2].tolist()))))
        N_Tag3.append(len(list(set(df1[tags3].tolist()))))
        #print(N_Tag2)
        #print('EF')
        if 'EF' in df1.columns:
            ef.append(df1['EF'].sum())
        #print(ef)
    _tags2 = 'N_' + tags2
    _tags3 = 'N_' + tags3
    if 'EF' in df.columns:
        sheet1 = pd.DataFrame({tags1:TAG1,structure:TAG1_structure,structureff:TAG1_FF,structurebb:TAG1_structurebb,similarityTAG:TAG1_similarity,'N':N,'SumCount':SumCount,'max':DataMax,'mean':DataMean,_tags2:N_Tag2,_tags3:N_Tag3,'z-score':zscore,'Normalized_z-score':normal_zscore,'EF':ef})
    else:
        sheet1 = pd.DataFrame({tags1:TAG1,structure:TAG1_structure,structureff:TAG1_FF,structurebb:TAG1_structurebb,similarityTAG:TAG1_similarity,'N':N,'SumCount':SumCount,'max':DataMax,'mean':DataMean,_tags2:N_Tag2,_tags3:N_Tag3,'z-score':zscore,'Normalized_z-score':normal_zscore})
    
    #sheet1_new = funcNewSheet(sheet1)#sheet转换成带有统计信息的sheet
    sheet1_new = sheet1
    # 求BB相似性
    sheet1_new[similarityTAG] = 0
    sheet1_new = sheet1_new.fillna('0') # 若结构为空，则用“NA”来填充
    preBB = ''
    struct_lst = sheet1_new[structure].tolist()
    if len(struct_lst):
        
        for i in struct_lst:
            if struct_lst[0] == 'NA':
                continue
            elif struct_lst[0] == '0':
                continue
            elif struct_lst[0] == '':
                continue
            elif struct_lst[0] == 'nan':
                continue
            else:
                preBB = struct_lst[0]
                break
        if preBB != '':
            for stru in struct_lst:
                if stru:
                    if stru == 'NA':
                        sheet1_new.loc[sheet1_new[structure] == stru,similarityTAG] = 0
                    elif stru == '0':
                        sheet1_new.loc[sheet1_new[structure] == stru,similarityTAG] = 0
                    elif stru == '':
                        sheet1_new.loc[sheet1_new[structure] == stru,similarityTAG] = 0
                    elif stru == 'nan':
                        sheet1_new.loc[sheet1_new[structure] == stru,similarityTAG] = 0
                    else:
                        sheet1_new.loc[sheet1_new[structure] == stru,similarityTAG] = similarity_func(preBB,stru)

    else:
        print("empty sample for plane")
                
    sheet1_new.reset_index(drop=True, inplace=True)
    #sheet1_new.to_excel(writer,'_'+tags1, index=None)
    
    # tags1_top10 = list(sheet1_new.loc[:9,tags1]) # return the top 10 of sum value for the tags list 
    tags1_top10 = list(sheet1_new.loc[:,tags1]) #尝试输出所有成线的tag，看计算速度如何
    # print(tags1_top10)
    
    #return sheet1_new, single_lst, tags1_top10
    return sheet1_new, tags1_top10
def similarity_func(molA,molB):
    # 通过smiles计算两个化合物的相似性
    molA_ch = ch.MolFromSmiles(molA)
    molA_fp = ch.GetMorganFingerprint(molA_ch, 2)
    molB_ch = ch.MolFromSmiles(molB)
    molB_fp = ch.GetMorganFingerprint(molB_ch, 2)
    sim = DataStructs.TanimotoSimilarity(molA_fp,molB_fp)
    sim_s = '%.4f'%(sim) # 保留4位有效数字
    return sim_s
def funcNewSumSheet(oldsheet,sumsheet):
    # sheet添加特定列，并进行统计学分析
    # sum 范围为0-10
    # 将多个面或线的dataframe合并为一个，记为sumsheet，在此基础上进行打分
    if sumsheet.shape[0] >= 1:
        
        dict1_N = sumsheet['N'].describe(percentiles=[0.25,0.50,0.75,0.95,0.99])
        dict1_SumCount = sumsheet['SumCount'].describe(percentiles=[0.25,0.50,0.75,0.95,0.99])
        dict1_DataMax = sumsheet['max'].describe(percentiles=[0.25,0.50,0.75,0.95,0.99])
        dict1_DataMean = sumsheet['mean'].describe(percentiles=[0.25,0.50,0.75,0.95,0.99])
    
        oldsheet['index1'] = 0
        oldsheet.loc[oldsheet['N'] >= dict1_N['99%'],'index1'] = 2.5
        oldsheet.loc[(oldsheet['N'] >= dict1_N['95%']) & (oldsheet['N'] < dict1_N['99%']), 'index1'] = 2
        oldsheet.loc[(oldsheet['N'] >= dict1_N['75%']) & (oldsheet['N'] < dict1_N['95%']), 'index1'] = 1.5
        oldsheet.loc[(oldsheet['N'] >= dict1_N['50%']) & (oldsheet['N'] < dict1_N['75%']), 'index1'] = 1
        oldsheet.loc[(oldsheet['N'] >= dict1_N['25%']) & (oldsheet['N'] < dict1_N['50%']), 'index1'] = 0.5
    
        oldsheet['index2'] = 0
        oldsheet.loc[oldsheet['SumCount'] >= dict1_SumCount['99%'], 'index2'] = 2.5
        oldsheet.loc[(oldsheet['SumCount'] >= dict1_SumCount['95%']) & (oldsheet['SumCount'] < dict1_SumCount['99%']), 'index2'] = 2
        oldsheet.loc[(oldsheet['SumCount'] >= dict1_SumCount['75%']) & (oldsheet['SumCount'] < dict1_SumCount['95%']), 'index2'] = 1.5
        oldsheet.loc[(oldsheet['SumCount'] >= dict1_SumCount['50%']) & (oldsheet['SumCount'] < dict1_SumCount['75%']), 'index2'] = 1
        oldsheet.loc[(oldsheet['SumCount'] >= dict1_SumCount['25%']) & (oldsheet['SumCount'] < dict1_SumCount['50%']), 'index2'] = 0.5
        
        
        oldsheet['index3'] = 0
        #oldsheet.loc[oldsheet['max'] >= dict1_DataMax['99%'], 'index3'] = 2.5
        oldsheet.loc[oldsheet['max'] >= dict1_DataMax['95%'], 'index3'] = 2
        oldsheet.loc[(oldsheet['max'] >= dict1_DataMax['75%']) & (oldsheet['max'] < dict1_DataMax['95%']), 'index3'] = 1.5
        oldsheet.loc[(oldsheet['max'] >= dict1_DataMax['50%']) & (oldsheet['max'] < dict1_DataMax['75%']), 'index3'] = 1
        oldsheet.loc[(oldsheet['max'] >= dict1_DataMax['25%']) & (oldsheet['max'] < dict1_DataMax['50%']), 'index3'] = 0.5
       
        oldsheet['index4'] = 0
        #oldsheet.loc[oldsheet['mean'] >= dict1_DataMean['99%'], 'index4'] = 2.5
        oldsheet.loc[oldsheet['mean'] >= dict1_DataMean['95%'], 'index4'] = 2
        oldsheet.loc[(oldsheet['mean'] >= dict1_DataMean['75%']) & (oldsheet['mean'] < dict1_DataMean['95%']), 'index4'] = 1.5
        oldsheet.loc[(oldsheet['mean'] >= dict1_DataMean['50%']) & (oldsheet['mean'] < dict1_DataMean['75%']), 'index4'] = 1
        oldsheet.loc[(oldsheet['mean'] >= dict1_DataMean['25%']) & (oldsheet['mean'] < dict1_DataMean['50%']), 'index4'] = 0.5
       
        # sheet1['index5'] = 0
        # sheet1.loc[(sheet1['N'] > sheet1[_tags2]) | (sheet1['N'] > sheet1[_tags3]), 'index5'] = 1
    
        oldsheet['Sum'] = oldsheet['index1'] + oldsheet['index2'] + oldsheet['index3'] + oldsheet['index4']
        # 排序并输出到excel
        newsheet = oldsheet.sort_values(by=['Sum','SumCount','N', 'max', 'mean'],ascending=False)
        newsheet.drop(['index1','index2','index3','index4'],axis=1, inplace=True) # 删除特定列
        newsheet = newsheet.fillna('0') # 若结构为空，则用“NA”来填充
        newsheet.reset_index(drop=True, inplace=True)
    else:
        oldsheet['Sum'] = 0
        newsheet = oldsheet
    
    return newsheet
def merge_nz_score(casefile, ntcfile):
    if os.path.exists(casefile) and os.path.exists(ntcfile):
        case_df = pd.read_csv(casefile, index_col=None, header=0, sep="\t",dtype=object)
        ntc_df = pd.read_csv(ntcfile, index_col=None, header=0, sep="\t", dtype=object)
        ntc_df.rename(columns={"Normalized_z-score": "NTC_Normalized_z-score",'Dedup_count':'NTC_Dedup_count','z-score':'NTC_z-score'}, inplace=True)
        #print(ntc_df.columns)
        #print(case_df.columns)
        ntc_df1 = ntc_df[['Tags','TAG1','TAG1_seq','TAG1_structure','TAG1_FF','TAG1_BB',
                    'TAG2','TAG2_seq','TAG2_structure','TAG2_FF','TAG2_BB',
                    'TAG3','TAG3_seq','TAG3_structure','TAG3_FF','TAG3_BB','NTC_Dedup_count','NTC_z-score',"NTC_Normalized_z-score"]]
        #ntc_df1['NTC_z-score'].astype(float)
        #print(ntc_df1['NTC_z-score'])
        #print([i for i in ntc_df1_lst if float(i) > 0])
        #print(ntc_df1)
        try:
            ntc_df1_lst = list(set(ntc_df1['NTC_z-score'].tolist()))
            ntc_df1_zscore_min = min(float(i) for i in ntc_df1_lst if float(i) > 0)
        except:
            df_case_lst = list(set(case_df['z-score'].tolist()))
            try:
                ntc_df1_zscore_min = min(float(i) for i in df_case_lst if float(i) > 0)
            except:
                ntc_df1_zscore_min = 0
        #print(ntc_df1_zscore_min)
        ntc_case_df = pd.merge(case_df, ntc_df1,on=['Tags','TAG1','TAG1_seq','TAG1_structure','TAG1_FF','TAG1_BB',
                    'TAG2','TAG2_seq','TAG2_structure','TAG2_FF','TAG2_BB',
                    'TAG3','TAG3_seq','TAG3_structure','TAG3_FF','TAG3_BB'],how="left")
        ntc_case_df.fillna(0, inplace=True)
        ntc_case_df["Type"] = ["0"] * len(ntc_case_df)
        #print(ntc_case_df)
        ntc_case_df.loc[(ntc_case_df["TAG1_BB"] != "0") & (ntc_case_df["TAG2_BB"] != "0") & (ntc_case_df["TAG3_BB"] != "0"), "Type"] = "123"
        ntc_case_df.loc[(ntc_case_df["TAG1_BB"] != "0") & (ntc_case_df["TAG2_BB"] != "0") & (ntc_case_df["TAG3_BB"] == "0"), "Type"] = "12"
        ntc_case_df.loc[(ntc_case_df["TAG1_BB"] != "0") & (ntc_case_df["TAG2_BB"] == "0") & (ntc_case_df["TAG3_BB"] != "0"), "Type"] = "13"
        ntc_case_df.loc[(ntc_case_df["TAG1_BB"] == "0") & (ntc_case_df["TAG2_BB"] != "0") & (ntc_case_df["TAG3_BB"] != "0"), "Type"] = "23"
        ntc_case_df.loc[(ntc_case_df["TAG1_BB"] != "0") & (ntc_case_df["TAG2_BB"] == "0") & (ntc_case_df["TAG3_BB"] == "0"), "Type"] = "1"
        ntc_case_df.loc[(ntc_case_df["TAG1_BB"] == "0") & (ntc_case_df["TAG2_BB"] != "0") & (ntc_case_df["TAG3_BB"] == "0"), "Type"] = "2"
        ntc_case_df.loc[(ntc_case_df["TAG1_BB"] == "0") & (ntc_case_df["TAG2_BB"] == "0") & (ntc_case_df["TAG3_BB"] != "0"), "Type"] = "3"
        try:
            ntc_case_df['EF'] = ntc_case_df.apply(lambda row:enrichment_score(row['z-score'], row['NTC_z-score'],ntc_df1_zscore_min), axis=1)
        except:
            ntc_case_df['EF'] = 0
        #print(ntc_case_df)
        return ntc_case_df
    else:
        print("缺少Case或NTC的Normalized_z-score文件！")
def enrichment_score(case,ntc,ntc_min):
    EF = 0
    if float(case) >0 and float(ntc)>0:
        EF = round(float(case)/float(ntc),4)
    elif float(case) >0 and float(ntc) <=0:
        try:
            EF = round(float(case)/float(ntc_min),4)
        except:
            EF = round(float(case),4)
    else:
        EF = 0
    return EF
if __name__ == '__main__':
    ConfigDict = resolveConfig("../config.txt")
    print(ConfigDict["bwaIndex"])
