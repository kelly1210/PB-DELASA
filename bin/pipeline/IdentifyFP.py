# -*- coding: utf-8 -*-
import pandas as pd
import argparse
import os
import sys
sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), '../monitor'))
import process

"""
# identify positive signal based on NTC sample
"""

def Tag2Struc(TagStru_file):
    # process.check_singlefile(TagStru_file)
    df = pd.read_csv(TagStru_file,sep='\t',dtype=object,header=0)
    df.drop_duplicates()
    df = df.dropna(subset=['Cycle'])
    # check file and change file to tagStru dict
    df['TagNumber'] = df['Cycle'] + '.' + df['Number']
    if 'StructureFF' not in list(df.columns):
        df['StructureFF'] = '0'
    tagStru_dict = df.set_index('TagNumber')['Structure'].to_dict()
    tagBB_dict = df.set_index('TagNumber')['Reagent'].to_dict()
    tagFF_dict = df.set_index('TagNumber')['StructureFF'].to_dict()
    return tagStru_dict,tagBB_dict,tagFF_dict

def ntcTAGlst(ntc_lst,libname,summary_df):
    ntc_tag_total = []
    for ntc_sample in ntc_lst:    
        ntc_plane = summary_df.loc[(summary_df['SampleID'] == ntc_sample) & (summary_df['LibraryID'] == libname), 'Plane'].values[0]
        ntc_line = summary_df.loc[(summary_df['SampleID'] == ntc_sample) & (summary_df['LibraryID'] == libname), 'Line'].values[0]
        ntc_tag_total.extend(ntc_plane.strip().split(";"))
        for i in ntc_line.strip().split(";"):
            ntc_tag_total.extend(i.strip().split(","))
        
    ntc_tag_lst_new = list(set(ntc_tag_total))
    if '0' in ntc_tag_lst_new:
        ntc_tag_lst_new.remove('0')
    return ntc_tag_lst_new

def ntcTAGtoBB(ntc_tag_lst_new,tagBB_dict):
    NTCbb_lst = []
    for bbkey in ntc_tag_lst_new:
        NTCbb = tagBB_dict[bbkey]
        NTCbb_lst.extend([NTCbb])
    NTCbb_lst_new = list(set(NTCbb_lst)) 
    if '0' in ntc_tag_lst_new:
        NTCbb_lst_new.remove('0')
    return NTCbb_lst_new

def ntcTAGtoBBline(ntc_tag_lst_new,tagBB_dict):
    NTCbb_lst = []
    NTCbb_lst_new = []
    if ntc_tag_lst_new:
        
        for bbkeyline in ntc_tag_lst_new:
            line = bbkeyline.strip().split(',')
            line_str = ','.join([tagBB_dict[line[0]],tagBB_dict[line[1]]])
            NTCbb_lst.extend([line_str])
        NTCbb_lst_new = list(set(NTCbb_lst)) 
        if '0' in ntc_tag_lst_new:
            NTCbb_lst_new.remove('0')
    return NTCbb_lst_new

def ntcTAGlstLine(ntc_lst,libname,summary_df):
    ntc_tag_total = []
    for ntc_sample in ntc_lst:
        ntc_line = summary_df.loc[(summary_df['SampleID'] == ntc_sample) & (summary_df['LibraryID'] == libname), 'Line'].values[0]
        ntc_tag_total.extend(ntc_line.strip().split(";"))
        
    ntc_tag_lst_new = list(set(ntc_tag_total))
    if '0' in ntc_tag_lst_new:
        ntc_tag_lst_new.remove('0')
    return ntc_tag_lst_new

def ntcTAGtoBBdot(ntc_tag_lst_new,tagBB_dict):
    NTCbb_lst = []
    NTCbb_lst_new = []
    if ntc_tag_lst_new:
        for bbkeydot in ntc_tag_lst_new:
            dot = bbkeydot.strip().split(',')
            dot_str = ','.join([tagBB_dict[dot[0]],tagBB_dict[dot[1]],tagBB_dict[dot[2]]])
            NTCbb_lst.extend([dot_str])
        NTCbb_lst_new = list(set(NTCbb_lst)) 
        if '0' in ntc_tag_lst_new:
            NTCbb_lst_new.remove('0')
    return NTCbb_lst_new

def ntcTAGlstDot(ntc_lst,libname,summary_df):
    # signle for NTC sample
    ntc_tag_total = []
    for ntc_sample in ntc_lst:
        ntc_Dot = summary_df.loc[(summary_df['SampleID'] == ntc_sample) & (summary_df['LibraryID'] == libname), 'Dot'].values[0]
        ntc_tag_total.extend(ntc_Dot.strip().split(";"))
    ntc_tag_lst_new = list(set(ntc_tag_total))
    if '0' in ntc_tag_lst_new:
        ntc_tag_lst_new.remove('0')
    return ntc_tag_lst_new
        
def identifyFP(protein):
    # false positive identified based on sample type
    FP = 'TP'
    if 'NTC' in protein:
        FP = 'FP'
    return FP

def summaryFP(batchName,PreSummaryFile,out_dlp_file,input_SampleInfo,input_LibInfo,code_smiles_dir):
    # identify positive signal
    protein_dict = process.total_sample_info(input_SampleInfo)[0]
    concentration_dict = process.total_sample_info(input_SampleInfo)[1]
    roundNumber_dict = process.total_sample_info(input_SampleInfo)[7]
    sample_lst = process.total_sample_info(input_SampleInfo)[6]
    ntc_dict = process.total_sample_info(input_SampleInfo)[8] 
    lib_dict = process.total_sample_info(input_SampleInfo)[2]
    library_size_dict = process.library_info(input_LibInfo)[7]
    # statistic
    sum_df = pd.read_csv(PreSummaryFile,sep='\t',header=0,dtype=object)
    sum_df.fillna('0')
    
    dlp_f = open(out_dlp_file,'w')
    dlp_name_lst = ['BatchName','Round','Protein','Protein_CC','SampleID','LibraryID','LibrarySize',
                'TAG1','TAG1_structure','TAG1_FF','TAG1_BB','TAG2','TAG2_structure','TAG2_FF','TAG2_BB',
                'TAG3','TAG3_structure','TAG3_FF','TAG3_BB','DLP','FP/TP','Confirmation(1/0)']
    dlp_f.write('\t'.join(dlp_name_lst)+'\n')
    for sample in sample_lst:
        protein = protein_dict[sample]
        ntc_str = ntc_dict[sample] 
        ntc_lst = ntc_str.strip().split(",")
        lib_lst = lib_dict[sample].strip().split(",")
        for lib in lib_lst:
            sample_info_lst = [batchName,roundNumber_dict[sample],protein_dict[sample],concentration_dict[sample],sample,lib,str(library_size_dict[lib])]
            #print(sample_info_lst)
            code_smiles_file = "{}/{}.code.smiles.txt".format(code_smiles_dir, lib)

            tagSTR_dict,tagBB_dict,tagFF_dict = Tag2Struc(code_smiles_file)
            for i in tagSTR_dict.keys():
                if process.is_number(tagSTR_dict[i]):
                    tagSTR_dict[i] = '0'
            for i in tagFF_dict.keys():
                if process.is_number(tagFF_dict[i]):
                    tagFF_dict[i] = '0'
            # Output the string of points, lines, and planes, with planes separated by semicolons.
            plane = '0';line='0';dot='0'
            if sum_df.empty == False:
                if (sample in sum_df['SampleID'].tolist()) and (lib in sum_df['LibraryID'].tolist()):
                    if not sum_df[(sum_df['SampleID'] == sample) & (sum_df['LibraryID'] == lib)].empty:
                        plane = sum_df.loc[(sum_df['SampleID'] == sample) & (sum_df['LibraryID'] == lib), 'Plane'].values[0]
                    
                        line = sum_df.loc[(sum_df['SampleID'] == sample) & (sum_df['LibraryID'] == lib), 'Line'].values[0]
                    
                        dot = sum_df.loc[(sum_df['SampleID'] == sample) & (sum_df['LibraryID'] == lib), 'Dot'].values[0]
                else:
                    print('the sample {0} and the lib {1} are not in the summary file!'.format(sample,lib))
            DLP = '0';FP = '0';conformation = '0' 
            if plane != '0':
                DLP = 'P'
                #print(DLP)
                plane_lst = plane.strip().split(";")
                for plane_single in plane_lst:
                    #print(plane_single)
                    plane_single_bb = tagBB_dict[plane_single]
                    # If there are control samples for the case samples.
                    # print(ntc_str)
                    if ntc_str != '0':
                        # print(ntc_lst)
                        # list of tags for ntc sample
                        ntc_tag_total_lst = ntcTAGlst(ntc_lst,lib,sum_df)
                        # list of bb for ntc sample
                        ntc_bb_total_lst = ntcTAGtoBB(ntc_tag_total_lst,tagBB_dict)
                        # print(ntc_tag_total_lst)
                        # Determine whether the tags corresponding to the planes in the case samples are present in the NTC samples.
                        if plane_single_bb in ntc_bb_total_lst:
                            FP = 'TP'
                            conformation = '1'
                        else:
                            FP = 'TP'
                            conformation = '1'
                    else:
                        FP = 'TP'
                        conformation = '1'
                    if identifyFP(protein) == 'FP':
                        FP = 'FP'
                        conformation = '0'
                    cycle = plane_single.strip().split(".")[0]
                    if cycle == '1':
                        dlp_lst = [plane_single,tagSTR_dict[plane_single],tagFF_dict[plane_single],tagBB_dict[plane_single],'0','0','0','0','0','0','0','0',DLP,FP,conformation]
                        dlp_f.write('\t'.join(sample_info_lst + dlp_lst)+'\n')
                        #print('\t'.join(sample_info_lst + dlp_lst)+'\n')
                    elif cycle == '2':
                        dlp_lst = ['0','0','0','0',plane_single,tagSTR_dict[plane_single],tagFF_dict[plane_single],tagBB_dict[plane_single],'0','0','0','0',DLP,FP,conformation]
                        dlp_f.write('\t'.join(sample_info_lst + dlp_lst)+'\n')
                    elif cycle == '3':
                        
                        dlp_lst = ['0','0','0','0','0','0','0','0',plane_single,tagSTR_dict[plane_single],tagFF_dict[plane_single],tagBB_dict[plane_single],DLP,FP,conformation]
                        dlp_f.write('\t'.join(sample_info_lst + dlp_lst)+'\n')
                    else:
                        continue
            if line != '0':
                # positive signal for line
                DLP = 'L'
                FP = identifyFP(protein) 
                line_lst = line.strip().split(";")
                # print(line_lst)
                
                for line_single in line_lst:
                    line_single_new = ntcTAGtoBBline([line_single],tagBB_dict)
                    if ntc_str != '0':
                        ntc_tag_total_lst = ntcTAGlstLine(ntc_lst,lib,sum_df)
                        ntc_bb_total_lst = ntcTAGtoBBline(ntc_tag_total_lst,tagBB_dict)
                        if line_single_new[0] in ntc_bb_total_lst:
                            FP = 'TP'
                            conformation = '1'
                        else:
                            FP = 'TP'
                            conformation = '1'
                    else:
                        FP = 'TP'
                        conformation = '1'
                        
                    if identifyFP(protein) == 'FP':
                        FP = 'FP'
                        conformation = '0'
                    line_each_tag = line_single.strip().split(",")
                    line_tag1 = line_each_tag[0]
                    line_tag2 = line_each_tag[1]
                    
                    cycle1 = line_tag1.strip().split(".")[0]
                    cycle2 = line_tag2.strip().split(".")[0]
                    
                    if cycle1 == '1' and cycle2 == '2':
                        dlp_lst = [line_tag1,tagSTR_dict[line_tag1],tagFF_dict[line_tag1],tagBB_dict[line_tag1],line_tag2,tagSTR_dict[line_tag2],tagFF_dict[line_tag2],tagBB_dict[line_tag2],'0','0','0','0',DLP,FP,conformation]
                        dlp_f.write('\t'.join(sample_info_lst + dlp_lst)+'\n')
                        print('\t'.join(sample_info_lst + dlp_lst)+'\n')
                    elif cycle1 == '1' and cycle2 == '3':
                        dlp_lst = [line_tag1,tagSTR_dict[line_tag1],tagFF_dict[line_tag1],tagBB_dict[line_tag1],'0','0','0','0',line_tag2,tagSTR_dict[line_tag2],tagFF_dict[line_tag2],tagBB_dict[line_tag2],DLP,FP,conformation]
                        dlp_f.write('\t'.join(sample_info_lst + dlp_lst)+'\n')
                    elif cycle1 == '2' and cycle2 == '3':
                        dlp_lst = ['0','0','0','0',line_tag1,tagSTR_dict[line_tag1],tagFF_dict[line_tag1],tagBB_dict[line_tag1],line_tag2,tagSTR_dict[line_tag2],tagFF_dict[line_tag2],tagBB_dict[line_tag2],DLP,FP,conformation]
                        dlp_f.write('\t'.join(sample_info_lst + dlp_lst)+'\n')
                    else:
                        continue
            if dot != '0':
                # positive signal for dot
                # print(dot)
                DLP = 'D'
                FP = identifyFP(protein)
                dot_lst = dot.strip().split(";")
                
                for dot_single in dot_lst:
                    dot_single_new = ntcTAGtoBBdot([dot_single],tagBB_dict)
                    if ntc_str != '0':
                        ntc_tag_total_lst = ntcTAGlstDot(ntc_lst,lib,sum_df)
                        ntc_bb_total_lst = ntcTAGtoBBdot(ntc_tag_total_lst,tagBB_dict)
                        # Determine whether the tags corresponding to the planes in the case samples are present in the NTC samples, and retain all signals that appear in the NTC samples in this update.
                        if dot_single_new[0] in ntc_bb_total_lst:
                            FP = 'TP'
                            conformation = '1'
                        else:
                            FP = 'TP'
                            conformation = '1'
                    else:
                        FP = 'TP'
                        conformation = '1'
                    if identifyFP(protein) == 'FP':
                        FP = 'FP'
                        conformation = '0'
                    # print(dot_single)
                    dot_each_tag = dot_single.strip().split(",")
                    dot_tag1 = dot_each_tag[0]
                    dot_tag2 = dot_each_tag[1]
                    dot_tag3 = dot_each_tag[2]
                    dlp_lst = [dot_tag1,tagSTR_dict[dot_tag1],tagFF_dict[dot_tag1],tagBB_dict[dot_tag1],dot_tag2,tagSTR_dict[dot_tag2],tagFF_dict[dot_tag2],tagBB_dict[dot_tag2],dot_tag3,tagSTR_dict[dot_tag3],tagFF_dict[dot_tag3],tagBB_dict[dot_tag3],DLP,FP,conformation]
                    dlp_f.write('\t'.join(sample_info_lst + dlp_lst)+'\n')  
    dlp_f.close()

def sameBBgrouped(arr):
    arr_lst = list(set(arr))
    arr_str = '|'.join(arr_lst)
    return arr_str
    
    
def dlp_file_grouped(dlp_f, dlp_out):
    #print(dlp_f)
    df = pd.read_csv(dlp_f, sep='\t', dtype=object)
    df.groupby(['BatchName','Round','Protein','Protein_CC','SampleID','LibraryID','LibrarySize','TAG1_structure',
                'TAG1_BB','TAG2_structure','TAG2_BB','TAG3_structure','TAG3_BB',
                'DLP','FP/TP','Confirmation(1/0)'])[['TAG1','TAG2','TAG3']].agg(sameBBgrouped)
    #print(dlp_out)
    df2 = df.reindex(columns=['BatchName','Round','Protein','Protein_CC','SampleID','LibraryID','LibrarySize',
                'TAG1','TAG1_structure','TAG1_FF','TAG1_BB','TAG2','TAG2_structure','TAG2_FF','TAG2_BB',
                'TAG3','TAG3_structure','TAG3_FF','TAG3_BB','DLP','FP/TP','Confirmation(1/0)'], fill_value='0')
    df2.to_csv(dlp_out, sep='\t',index=None)
    
if __name__ == "__main__":
    #data_start = time.time()
    #print(time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(time.time())))
    parser = argparse.ArgumentParser(description='get deduped count file')
    parser.add_argument('--Input', '-i', help='input dlp summary file')
    parser.add_argument('--OutputSumPath', '-o', help='output summary path dir')
    #parser.add_argument('--OutputTemp', '-o2', help='output temp dlp file')
    parser.add_argument('--Batch', '-b', help='batchName')
    #parser.add_argument('--codeSmiles', '-c', help='codeSmiles file dir')
    parser.add_argument('--SampleInfo', '-s', help='SampleInfo dir')
    parser.add_argument('--LibInfo', '-l', help='LibInfo dir')
    
    args = parser.parse_args()
    
    PreSummaryFile = args.Input # batchName + '.summay.txt'
    batchName = args.Batch
    out_dlp_file = args.OutputSumPath # batchName + '.dlp.summay.txt'
    out_temp_file = out_dlp_file + '.temp'

    code_smiles_dir = sys.path[0] + '/../../lib/Libraries/'
    input_SampleInfo = args.SampleInfo
    input_LibInfo = args.LibInfo
    
    summaryFP(batchName,PreSummaryFile,out_temp_file,input_SampleInfo,input_LibInfo,code_smiles_dir)
    dlp_file_grouped(out_temp_file, out_dlp_file)
    os.remove('{}'.format(out_temp_file))
