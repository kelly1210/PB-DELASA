# -*- coding: utf-8 -*-
"""
QC summary
"""
import pandas as pd
import os,sys
from argparse import ArgumentParser
sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), '../monitor'))
import process

# contaminated sequence
def find_con_lib(con_seq,libinfo):
    if con_seq != '-':
        lib_info_df = pd.read_csv(libinfo,index_col=None,header=0,sep="\t")
        con_df = lib_info_df[lib_info_df["5_3"]==con_seq]
        if len(con_df) != 0:
            con_lib = con_df["LibraryID"].tolist()[0]
        else:
            con_lib = ""
    else:
        con_lib = "-"
    return con_lib
def modified_qc_summary(total_info_file, QC_batch, QC_sample, contaminationFile, QC_lib, sumQC):
    # QC for the batch data
    QC_batch_f = pd.read_csv(QC_batch,sep='\t', header=0)
    # QC for each sample
    QC_sample_f = pd.read_csv(QC_sample,sep='\t', header=0)
    # QC for each library
    QC_lib_f = pd.read_csv(QC_lib,sep='\t', header=0)
    # contamination analysis
    QC_sample_contamination = pd.read_csv(contaminationFile,sep='\t', header =0)
    #print(QC_sample_contamination)
    sumQC_f = open(sumQC,'w')
    column_lst = ['BatchName','RawTotalReads','TrimTotalReads','SampleID','Raw_single_reads(R1/R2)','Trim_single_reads(R1/R2)',
                  'Combined_counts','Split_counts','Coding_counts','Raw_counts',
                  'Dedup_counts','ContaminationSeqMax','ContaminationSeqMaxRatio',
                  'QPCR','QPCR_coverage']
    sumQC_f.write('\t'.join(column_lst) + '\n')
    BatchName = QC_batch_f['BatchName'][0] 
    RawTotalReads = QC_batch_f['raw_single_reads(R1/R2)'][0] 
    TrimTotalReads = QC_batch_f['trim_single_reads(R1/R2)'][0] 

    sample_lst = process.total_sample_info(total_info_file)[8] 
    qpcr_dict = process.total_sample_info(total_info_file)[11] 

    qc_sample_lst = QC_sample_f['SampleID'].tolist()
    qc_lib_lst = QC_lib_f['SampleID'].tolist()
    warnning_info = 'NO'
    sum_Raw_single_reads = 0
    sum_Trim_single_reads = 0
    sum_Combined_counts = 0
    sum_Split_counts = 0
    sum_Coding_counts = 0
    sum_Raw_counts = 0
    sum_Dedup_counts = 0
    sum_seq_counts_ratio = 0
    sum_qpcr = 0
    sum_seq_counts_max = 0
    # print(sample_lst)
    for sample in sample_lst:
        if sample not in qc_sample_lst:
            warnning_info = 'Warnning: Sample information is incomplete: {}.'.format(sample)
        elif sample not in qc_lib_lst:
            print('Warnning: Sample lib information is incomplete: {}.'.format(sample))
        else:
            Raw_single_reads = QC_sample_f[QC_sample_f['SampleID'] == sample]['raw_single_reads(R1/R2)'].tolist()[0]
            sum_Raw_single_reads += int(Raw_single_reads)
            Trim_single_reads = QC_sample_f[QC_sample_f['SampleID'] == sample]['trim_single_reads(R1/R2)'].tolist()[0]
            sum_Trim_single_reads += int(Trim_single_reads)
            #print(sample,QC_lib_f[QC_lib_f['Protein_CC_ID'] == sample]['Combined_counts'].tolist())
            Combined_counts = QC_lib_f[QC_lib_f['SampleID'] == sample]['Combined_counts'].tolist()[0]
            sum_Combined_counts += int(Combined_counts)
            Split_counts = QC_lib_f[QC_lib_f['SampleID'] == sample]['Split_counts'].sum()
            sum_Split_counts += Split_counts
            Coding_counts = QC_lib_f[QC_lib_f['SampleID'] == sample]['Coding_counts'].sum()
            sum_Coding_counts += Coding_counts
            Raw_counts = QC_lib_f[QC_lib_f['SampleID'] == sample]['Raw_counts'].sum()
            sum_Raw_counts += Raw_counts
            Dedup_counts = QC_lib_f[QC_lib_f['SampleID'] == sample]['Dedup_counts'].sum()
            sum_Dedup_counts += Dedup_counts
            # add qpcr
            qpcr = float(qpcr_dict[sample])
            sum_qpcr += qpcr
            if qpcr!=0:
                qpcr_ratio = '{:.2%}'.format(Dedup_counts/qpcr)
            else:
                qpcr_ratio = '0'
            #seq_Raw_reads = QC_sample_contamination['RawReads'][0]
            seq_counts_max = QC_sample_contamination[sample].max()
            seq_counts_ratio = 0
            if Trim_single_reads != 0:
                seq_counts_ratio = '%.4f'%(seq_counts_max / int(Trim_single_reads)) 
            #print(QC_sample_contamination[QC_sample_contamination[sample] == seq_counts_max])
            try: 
                seq_contamin = QC_sample_contamination[QC_sample_contamination[sample] == seq_counts_max]['ContaminSeq'].tolist()[0]
            except:
                seq_contamin = 'NA'
            sum_seq_counts_max += seq_counts_max

            f_sample = [BatchName, str(RawTotalReads), str(TrimTotalReads), sample, str(Raw_single_reads), str(Trim_single_reads),
                        str(Combined_counts),str(Split_counts), str(Coding_counts), str(Raw_counts), str(Dedup_counts),
                        str(seq_contamin),str(seq_counts_ratio),str(qpcr),str(qpcr_ratio)]
            

            sumQC_f.write('\t'.join(f_sample) + '\n')
            
    if sum_qpcr!=0:
        sum_qpcr_cov = '%.2f'%(sum_Dedup_counts/sum_qpcr)
    else:
        sum_qpcr_cov = 0
            
    if sum_Trim_single_reads!=0:
        sum_seq_counts_ratio = float('%.4f'%(sum_seq_counts_max / sum_Trim_single_reads))
    else:
        sum_seq_counts_ratio = 0
        
    f_sum = [str('sum'), str(RawTotalReads), str(TrimTotalReads), BatchName, str(sum_Raw_single_reads),
             str(sum_Trim_single_reads), str(sum_Combined_counts),str(sum_Split_counts), str(sum_Coding_counts),
             str(sum_Raw_counts), str(sum_Dedup_counts),str('-'),str(sum_seq_counts_max),str(sum_qpcr),str(sum_qpcr_cov)]
    f_ratio = [str('ratio(%)'), str(100), str('%.1f'%(TrimTotalReads/RawTotalReads*100)), BatchName, str('%.1f'%(sum_Raw_single_reads/RawTotalReads*100)),
             str('%.1f'%(sum_Trim_single_reads/RawTotalReads*100)), str('%.1f'%(sum_Combined_counts/RawTotalReads*100)),
               str('%.1f'%(sum_Split_counts/RawTotalReads*100)),str('%.1f'%(sum_Coding_counts/RawTotalReads*100)),
               str('%.1f'%(sum_Raw_counts/RawTotalReads*100)), str('%.1f'%(sum_Dedup_counts/RawTotalReads*100)),str('-'),str(sum_seq_counts_ratio),str(sum_qpcr),str(sum_qpcr_cov)]
    sumQC_f.write('\t'.join(f_sum) + '\n')
    sumQC_f.write('\t'.join(f_ratio) + '\n')
    return warnning_info
def DEL_QC_stat(batchname,sum_qc,sam_sum_qc,m_sum_qc,sinfo,libinfo):
    sinfo_df = pd.read_csv(sinfo,index_col=None,header=0,sep="\t")
    sinfo_df1 = sinfo_df[["SampleID","LibraryID"]].copy()
    sinfo_df1["LibraryID"] = sinfo_df1["LibraryID"].apply(lambda x:x.split(',')) 
    sam_sum_df = pd.read_csv(sam_sum_qc,index_col=None,header=0,sep="\t")
    # Q30 for all data 
    sum_qc_df = pd.read_csv(sum_qc,index_col=None,header=0,sep="\t")
    raw_q30 = sum_qc_df.loc[0,"raw_q30_rate"]/100
    trim_q30 = sum_qc_df.loc[0,"trim_q30_rate"]/100
    # Q30 for each sample
    print(sam_sum_df["raw_q30_rate"])
    sam_sum_df["raw_q30_rate"] = sam_sum_df["raw_q30_rate"].astype(float)
    sam_sum_df["trim_q30_rate"] = sam_sum_df["trim_q30_rate"].astype(float)
    sam_sum_df["Raw_Q30"] = sam_sum_df["raw_q30_rate"].apply(lambda x: x/100)
    sam_sum_df["Trim_Q30"] = sam_sum_df["trim_q30_rate"].apply(lambda x: x/100)

    sam_sum_df = pd.merge(sam_sum_df,sinfo_df1,left_on="SampleID",right_on="SampleID")
    sam_sum_df.loc[len(sam_sum_df),"SampleID"] = batchname
    sam_sum_df.loc[len(sam_sum_df)-1,"Raw_Q30"] = raw_q30
    sam_sum_df.loc[len(sam_sum_df)-1,"Trim_Q30"] = trim_q30
    sam_sum_df1 = sam_sum_df[["SampleID","Raw_Q30","Trim_Q30","LibraryID"]]
    # LowQuality_rate
    m_sum_df = pd.read_csv(m_sum_qc,index_col=None,header=0,sep="\t")
    # LowQuality_rate
    RawTotalReads = m_sum_df.loc[0,"RawTotalReads"]
    TrimTotalReads = m_sum_df.loc[0,"TrimTotalReads"]
    LQ_rate = (RawTotalReads-TrimTotalReads)/RawTotalReads

    # sample_index_translated_rate
    sit_rate = m_sum_df.loc[len(m_sum_df)-2,"Raw_single_reads(R1/R2)"]/m_sum_df.loc[len(m_sum_df)-2,"TrimTotalReads"]

    # Lib_translate_rate
    m_sum_df["Lib_translate_rate"] = m_sum_df["Split_counts"]/m_sum_df["Combined_counts"]

    # lib_code_translated_rate
    m_sum_df["Lib_code_translated_rate"] = m_sum_df["Raw_counts"]/m_sum_df["Split_counts"]

    # dup_ratio
    m_sum_df["Dup_ratio"] = (m_sum_df["Raw_counts"]-m_sum_df["Dedup_counts"])/m_sum_df["Raw_counts"]

    #sequcing depth
    m_sum_df["Sequcing_Depth"] = m_sum_df["Raw_counts"]/m_sum_df["QPCR"]

    # NGS_ratio
    m_sum_df["NGS_ratio"] = m_sum_df["Dedup_counts"]/m_sum_df["QPCR"]

    # summary
    qc_df = m_sum_df[["SampleID","Raw_single_reads(R1/R2)","Lib_translate_rate","Lib_code_translated_rate",\
                      "Dup_ratio","NGS_ratio","Sequcing_Depth","ContaminationSeqMax","ContaminationSeqMaxRatio","QPCR"]].copy()
    qc_df["QPCR"] = qc_df["QPCR"].apply(lambda x: "{:.3E}".format(x))

    con_ratio = qc_df.loc[len(qc_df)-1,"ContaminationSeqMaxRatio"]
    qc_df1 = qc_df.loc[0:len(qc_df)-2,:]
    qc_df1 = pd.merge(qc_df1,sam_sum_df1,on="SampleID")
    qc_df1_len = len(qc_df1)-1
    qc_df1["LowQuality_rate"] = ["-"]*qc_df1_len + [LQ_rate]
    qc_df1["Sample_index_translated_rate"] = ["-"]*qc_df1_len+[sit_rate]
    qc_df1.loc[len(qc_df1)-1,"ContaminationSeqMaxRatio"] = con_ratio
    qc_df1["Contamination_Librarys"] = qc_df1.loc[:len(qc_df1)-2,"ContaminationSeqMax"].apply(lambda x: find_con_lib(x,libinfo))
    qc_df1["Violations"] = [0]*len(qc_df1)
    for i in range(len(qc_df1)-1):
        n = 0
        if qc_df1.loc[i,"Raw_Q30"] < 0.75:
            n+=1
        if qc_df1.loc[i,"Trim_Q30"] < 0.85:
            n+=1
        if qc_df1.loc[i,"Raw_single_reads(R1/R2)"] < 1000:
            n+=1
        if qc_df1.loc[i,"Lib_translate_rate"] < 0.9:
            n+=1
        if qc_df1.loc[i,"Lib_code_translated_rate"] < 0.8:
            n+=1
        if qc_df1.loc[i,"Dup_ratio"] < 0.5 :
            n+=1
        if qc_df1.loc[i,"Dup_ratio"] > 0.9 :
            n+=1
        if qc_df1.loc[i,"NGS_ratio"] > 1:
            n+=1
        if qc_df1.loc[i,"NGS_ratio"] < 0.1:
            n+=1
        if qc_df1.loc[i,"Sequcing_Depth"] < 1.5:
            n+=1
        qc_df1.loc[i,"Violations"] = n
    n1 = 0
    if qc_df1.loc[len(qc_df1)-1,"Raw_Q30"] < 0.75:
        n1+=1
    if qc_df1.loc[len(qc_df1)-1,"Trim_Q30"] < 0.85:
        n1+=1
    if qc_df1.loc[len(qc_df1)-1,"LowQuality_rate"] > 0.8:
        n1+=1
    if qc_df1.loc[len(qc_df1)-1,"Sample_index_translated_rate"] < 0.9:
        n1+=1
    if qc_df1.loc[len(qc_df1)-1,"Lib_translate_rate"] < 0.9:
        n1+=1
    if qc_df1.loc[len(qc_df1)-1,"Lib_code_translated_rate"] < 0.8:
        n1+=1
    if qc_df1.loc[len(qc_df1)-1,"Dup_ratio"] < 0.5:
        n1+=1
    if qc_df1.loc[len(qc_df1)-1,"Dup_ratio"] > 0.9:
        n1+=1
    if qc_df1.loc[len(qc_df1)-1,"Sequcing_Depth"] < 1.5:
        n1+=1
    if qc_df1.loc[len(qc_df1)-1,"NGS_ratio"] > 1:
        n1+=1
    if qc_df1.loc[len(qc_df1)-1,"NGS_ratio"] < 0.1:
        n1+=1
    if qc_df1.loc[len(qc_df1)-1,"ContaminationSeqMaxRatio"] > 0.05:
        n1+=1
    qc_df1.loc[len(qc_df1)-1,"Violations"] = n1

    qc_df1.loc[len(qc_df1)-1,"ContaminationSeqMax"] = "-"

    qc_df1["Raw_Q30"] = qc_df1["Raw_Q30"].apply(lambda x: format(x,'.2%'))
    qc_df1["Trim_Q30"] = qc_df1["Trim_Q30"].apply(lambda x: format(x,'.2%'))
    qc_df1.loc[len(qc_df1)-1,"LowQuality_rate"] = format(qc_df1.loc[len(qc_df1)-1,"LowQuality_rate"],".2%")
    qc_df1.loc[len(qc_df1)-1,"Sample_index_translated_rate"] = format(qc_df1.loc[len(qc_df1)-1,"Sample_index_translated_rate"],".2%")
    qc_df1["Lib_translate_rate"] = qc_df1["Lib_translate_rate"].apply(lambda x: format(x,'.2%'))
    qc_df1["Lib_code_translated_rate"] = qc_df1["Lib_code_translated_rate"].apply(lambda x: format(x,'.2%'))
    qc_df1["Dup_ratio"] = qc_df1["Dup_ratio"].apply(lambda x: format(x,'.2%'))
    qc_df1["NGS_ratio"] = qc_df1["NGS_ratio"].apply(lambda x: format(x,'.2%'))
    qc_df1["Sequcing_Depth"] = qc_df1["Sequcing_Depth"].apply(lambda x:format(x,'.2'))
    qc_df1["ContaminationSeqMaxRatio"] = qc_df1["ContaminationSeqMaxRatio"].apply(lambda x: format(x,'.2%'))
    qc_df2 = qc_df1[["SampleID","Raw_Q30","Trim_Q30","LowQuality_rate","Sample_index_translated_rate",\
        "Raw_single_reads(R1/R2)","Lib_translate_rate","Lib_code_translated_rate","Dup_ratio","NGS_ratio","Sequcing_Depth","ContaminationSeqMax",\
        "ContaminationSeqMaxRatio","Contamination_Librarys","Violations","QPCR"]]
    return qc_df2
def raw_reads_qc1(input_QC1_sum):
    QC_batch_f = pd.read_csv(input_QC1_sum, sep='\t', header=0)
    RawTotalReads = QC_batch_f['raw_single_reads(R1/R2)'][0]  
    return RawTotalReads

def ContaminationTable(sample_lst,input_QC1_sum,ana_dir, contamin_file):
    # contamin_file = table_dir + '/' + runid + '.contamination.txt'
    contamin_dir = ana_dir + '/02.Merge/'
    contamin_out = ana_dir + '/08.Contamination/'
    process.mkdir(contamin_out)
    seq_lst=[]
    for sample in sample_lst:
        contamin_f_name = contamin_dir + sample + '.combined.unfound.lib.seq.txt'
        contamin_out_f_name = contamin_out + sample + '.combined.unfound.lib.seq.sort.txt'
        contamin_f = pd.read_csv(contamin_f_name,sep='\t',header=None,names=['ContaminSeq','Reads'])
        contamin_f = contamin_f[contamin_f['ContaminSeq'] != '']
        contamin_f = contamin_f.sort_values(by='Reads',ascending=False)
        contamin_f = contamin_f.reset_index(drop=True)
        # print(contamin_out_f_name)
        # print(contamin_f.loc[:3,'ContaminSeq'])
        contamin_f.to_csv(contamin_out_f_name,sep='\t',index=None)
        # seq_dict = contamin_f.loc[:5].groupby('ContaminSeq').Reads.apply(list).to_dict()
        seq_lst.extend(list(contamin_f.loc[:3,'ContaminSeq'].tolist()))
    # print(len(list(set(seq_lst))))
    seq_lst = list(set(seq_lst))
    if '' in seq_lst:
        seq_lst.remove('')
    raw_single_reads = raw_reads_qc1(input_QC1_sum)

    name_lst = []
    name_lst.extend(['ContaminSeq'])
    name_lst.extend(sample_lst)
    name_lst.extend(['RawReads','SumContamin','SumRatioContamin'])

    df_out = pd.DataFrame(columns=name_lst)
    df_out.loc[:, 'SumContamin'] = 0

    for sample in sample_lst:
        contamin_f_name = contamin_dir + sample + '.combined.unfound.lib.seq.txt'
        contamin_f = pd.read_csv(contamin_f_name, sep='\t', header=None, names=['ContaminSeq', 'Reads'])
        lines = 0
        #print(seq_lst)
        for seq in seq_lst:
            df_out.loc[lines, 'ContaminSeq'] = seq
            # print(contamin_f.loc[contamin_f['ContaminSeq'] == seq,'Reads'].tolist())
            contamin_seq_lst = contamin_f.loc[contamin_f['ContaminSeq'] == seq,'Reads'].tolist()
            if contamin_seq_lst:
                df_out.loc[lines, sample] = contamin_f.loc[contamin_f['ContaminSeq'] == seq,'Reads'].tolist()[0]
            else:
                df_out.loc[lines, sample] = 0
            lines += 1
    df_out['RawReads'] = 0
    df_out['SumContamin'] = df_out.drop('ContaminSeq',axis=1).apply(lambda x: x.sum(), axis=1)
    df_out.loc[:,'RawReads'] = int(raw_single_reads)
    df_out['SumRatioContamin'] = df_out['SumContamin']/df_out['RawReads']
    # df_out = df_out[df_out['SumRatioContamin'] >= 0.01]
    df_out['SumRatioContamin'] = df_out['SumRatioContamin'].map(lambda x:format(x,".2%"))
    df_out = df_out.sort_values(by='SumRatioContamin',ascending=False)
    df_out = df_out.reset_index(drop=True)
    #print(df_out)

    df_out.loc[:10].to_csv(contamin_file,sep='\t',index=0)
    # return contamin_file
def severalCSV2xlsx(ana_dir,batchname,sinfo,libinfo,out_file):
    # out_file = {batchname}.qc.xlsx
    qc_comment = sys.path[0] + '/../../lib/QC_comments.txt'
    qc_com_df = pd.read_csv(qc_comment,index_col=None,header=0,sep="\t")

    sheet_file=['summary.qc','samples.summary.qc','modified.summary.qc','contamination','lib.summary.qc','dlp.detail','QC','Comments']
    writer = pd.ExcelWriter(out_file)
    files = os.listdir(ana_dir)
    file1 = batchname + '.summary.qc.txt'
    file2 = batchname + '.samples.summary.qc.txt'
    file3 = batchname + '.modified.summary.qc.txt'
    file4 = batchname + '.contamination.txt'
    file5 = batchname + '.lib.summary.qc.txt'
    file6 = batchname + '.count.sum.dlp.detail.txt'
    
    sample_lst = process.total_sample_info(sinfo)[6]
    # contamination file
    ContaminationTable(sample_lst,ana_dir + '/' + file1,ana_dir + '/../../', ana_dir + '/' + file4)
    # modified.summary.qc.txt
    modified_qc_summary(sinfo, ana_dir + '/' + file1, ana_dir + '/' + file2, ana_dir + '/' + file4, ana_dir + '/' + file5, ana_dir + '/' + file3)
    
    if file1 in files:
        df = pd.read_csv(ana_dir + '/' + file1,sep = "\t",header=0)
        df.to_excel(writer, sheet_name='summary.qc')
    if file2 in files:
        df = pd.read_csv(ana_dir + '/' + file2,sep = "\t",header=0)
        df.to_excel(writer, sheet_name='samples.summary.qc')
    if file3 in files:
        df = pd.read_csv(ana_dir + '/' + file3,sep = "\t",header=0)
        df.to_excel(writer, sheet_name='modified.summary.qc')
    if file4 in files:
        df = pd.read_csv(ana_dir + '/' + file4,sep = "\t",header=0)
        df.to_excel(writer, sheet_name='contamination')
    if file5 in files:
        df = pd.read_csv(ana_dir + '/' + file5,sep = "\t",header=0)
        df.to_excel(writer, sheet_name='lib.summary.qc')
    if file6 in files:
        df = pd.read_csv(ana_dir + '/' + file6,sep = "\t",header=0)
        df.to_excel(writer, sheet_name='dlp.detail')
        df1 = DEL_QC_stat(batchname,ana_dir+'/'+file1,ana_dir+'/'+file2,ana_dir+'/'+file3,sinfo,libinfo)
        df1.to_excel(writer, sheet_name='QC')
        qc_com_df.to_excel(writer,sheet_name="Comments")
    else:
        for var in sheet_file:
            df = pd.DataFrame()
            df.to_excel(writer, sheet_name=var)
    writer.close()
        
if __name__ == "__main__":
    #data_start = time.time()

    parser = ArgumentParser(description="Generate InnerReport")
    parser.add_argument('--input', '-i', help='report dir')
    #parser.add_argument('--input2', '-i2', help='QC comments context')
    parser.add_argument('--batchid', '-b', help='runid')
    parser.add_argument('--sample_info','-s',help='sample info')
    parser.add_argument('--library_info','-l',help='library info')
    parser.add_argument('--output', '-o', help='SampleInfo dir')
    args = parser.parse_args()

    severalCSV2xlsx(args.input,args.batchid,args.sample_info,args.library_info,args.output)
    
