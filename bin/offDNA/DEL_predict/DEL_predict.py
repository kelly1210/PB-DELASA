#! python3
# -*- encoding: utf-8 -*-
'''
@File    :  DEL_predict_v2.0.py
@Time    :  2023/05/31 17:44:57
@Autor   :  Meng_xiangfei
@Version :  12.0
Contact  :  meng_xiangfei@pharmablock.com
DESC     :  End-to-end prediction of DEL signals is primarily achieved by utilizing RF algorithms to build models, and the predicted signals are annotated based on new Confirmation labeling logic.
            
'''

from joblib import load
import pandas as pd
from functools import reduce
import os,warnings,re,argparse,sys
sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), '../../monitor'))
import process
warnings.filterwarnings("ignore")

def data_deal(count_file,tag1,tag2):
    count_df = pd.read_csv(count_file, index_col=None, header=0, sep="\t")
    count_df1 = count_df[[tag1,tag2,"Dedup_count"]]
    total_n = len(count_df)
    total_dedup_count = count_df["Dedup_count"].sum()
    count_df2 = count_df1.groupby([tag1,tag2])["Dedup_count"].sum()
    count_df2_max = count_df1.groupby([tag1,tag2])["Dedup_count"].max()
    count_df2_max1 = count_df2_max.reset_index()
    count_df2_max1.rename(columns={'Dedup_count':'Max'},inplace=True)
    count_df2_1 = count_df2.reset_index()
    count_df3 = count_df1.value_counts([tag1,tag2])
    count_df3_1 = count_df3.reset_index()
    count_df3_1.rename(columns = {0:"N"},inplace=True)

    rcount_df = reduce(lambda left,right: pd.merge(left,right,on=[tag1,tag2],how='inner'),[count_df2_1,\
                        count_df2_max1,count_df3_1])
    
    rcount_df = rcount_df[rcount_df["Dedup_count"] != 0]
    rcount_df['Line'] = [0]*len(rcount_df)
    rcount_df["Total_N"] = [total_n]*len(rcount_df)
    rcount_df["Total_Dedup_count"] = [total_dedup_count]*len(rcount_df)
    rcount_df["Mean_Dedup"] = rcount_df["Dedup_count"] / rcount_df["N"]
    rcount_df["Dedup/Total_Dedup"] = rcount_df["Dedup_count"] / rcount_df["Total_Dedup_count"]
    rcount_df["N/Total_N"] = rcount_df["N"] / rcount_df["Total_N"]
    rcount_df["Max/Mean_Dedup"] = rcount_df["Max"]/rcount_df["Mean_Dedup"]
    rcount_df["Mean_fold"] = rcount_df["Mean_Dedup"]*rcount_df["Total_N"]/rcount_df["Total_Dedup_count"]
    rcount_df["Dedup/Total_Dedup"].apply(lambda x: format(x,'.7f'))
    rcount_df["N/Total_N"].apply(lambda x: format(x, '.7f'))
    rcount_df["Max/Mean_Dedup"].apply(lambda x: format(x,'.7f'))
    rcount_df.rename(columns={tag1:"X1",tag2:"X2"},inplace=True)
    return rcount_df

def batch_data_deal(sinfo,batchpath,outpath):
    filter_counts = []
    total_rdf = pd.DataFrame()
    protein_dict = process.total_sample_info(sinfo)[0]
    for root,dirs,files in os.walk(batchpath):
        for file in files:
            if "count.freq.filter.txt" in file:
                filter_counts.append(os.path.join(root,file))
    for file in filter_counts:
        for i in [["TAG1","TAG2"],["TAG1","TAG3"],["TAG2","TAG3"]]:
            try:
                tmp_rdf = data_deal(file,i[0],i[1])
                sample_id_lib_id = file.split("/")[-1].replace(".count.freq.filter.txt","")
                sample_id = sample_id_lib_id.strip().split('.')[0]
                protein = protein_dict[sample_id]
                lib_id = sample_id_lib_id.strip().split('.')[1]
                tmp_rdf["SampleID"] = [sample_id]*len(tmp_rdf)
                tmp_rdf["LibraryID"] = [lib_id]*len(tmp_rdf)
                tmp_rdf["Protein"] = [protein]*len(tmp_rdf)
            except:
                continue
            total_rdf = pd.concat([total_rdf,tmp_rdf])
        
    total_rdf1 = total_rdf[total_rdf["N"]!=1]
    total_rdf1["Dedup/Total_Dedup"].apply(lambda x: format(x,'.7f'))
    total_rdf1["N/Total_N"].apply(lambda x: format(x, '.7f'))
    return total_rdf1

def predict_signal(sinfo_file,dlp_file,batch_dwar,model,out_file):
    batchname = re.sub(".count.sum.dlp.detail.txt","",dlp_file.split("/")[-1])
    signal_df =pd.read_csv(dlp_file,index_col=None,header=0,sep="\t")
    signal_df["TAG1"] = signal_df["TAG1"].apply(lambda x: "%.3f"%x)
    signal_df["TAG2"] = signal_df["TAG2"].apply(lambda x: "%.3f"%x)
    signal_df["TAG3"] = signal_df["TAG3"].apply(lambda x: "%.3f"%x)
    
    ntc_signal_df = signal_df[signal_df["Protein"].str.contains("NTC")]
    ntc_signal_df.reset_index(inplace=True)

    # TP: line,plane
    tp_case_signal_df = signal_df[(~signal_df["Protein"].str.contains("NTC")) & (signal_df["FP/TP"]=="TP") & (signal_df["DLP"]!="D")]
    # TP: Dot
    tp_case_signal_D_df = signal_df[(~signal_df["Protein"].str.contains("NTC")) & (signal_df["FP/TP"]=="TP") & (signal_df["DLP"]=="D")]
    # FP
    fp_case_signal_df = signal_df[(~signal_df["Protein"].str.contains("NTC")) & (signal_df["FP/TP"]=="FP")]

    tp_case_signal_df["Confirmation(1/0)"] = 0
    tp_case_signal_df.reset_index(inplace=True)

    format_signal_df = batch_data_deal(sinfo_file,batch_dwar,batchname)
    format_signal_df["X1"] = format_signal_df["X1"].apply(lambda x: "%.3f"%x)
    format_signal_df["X2"] = format_signal_df["X2"].apply(lambda x: "%.3f"%x)
    

    DEL_X = format_signal_df[["Dedup_count","N","Max","Mean_Dedup","Dedup/Total_Dedup","N/Total_N","Max/Mean_Dedup","Mean_fold"]]
    clf = load(model)
    line_predict = clf.predict(DEL_X)
    format_signal_df["Predict"] = line_predict
    ntc_format_signal_df = format_signal_df[(format_signal_df["Protein"].str.contains("NTC"))&(format_signal_df["Predict"]==1)]
    ntc_format_signal_df.reset_index(inplace=True)

    case_format_signal_df = format_signal_df[(~format_signal_df["Protein"].str.contains("NTC"))&(format_signal_df["Predict"]==1)]
    case_format_signal_df.reset_index(inplace=True)

    # NTC signal(dot,line)
    for i in range(len(ntc_format_signal_df)):
        sampleid = ntc_format_signal_df.loc[i,"SampleID"]
        libid = ntc_format_signal_df.loc[i,"LibraryID"]
        t1 = ntc_format_signal_df.loc[i,"X1"]
        t2 = ntc_format_signal_df.loc[i,"X2"]
        ntc_signal_df.loc[ntc_signal_df[(ntc_signal_df["SampleID"]==sampleid)&(ntc_signal_df["LibraryID"]==libid)&\
            (ntc_signal_df["TAG1"]==t1)&(ntc_signal_df["TAG2"]==t2)].index,"Confirmation(1/0)"]=1
        ntc_signal_df.loc[ntc_signal_df[(ntc_signal_df["SampleID"]==sampleid)&(ntc_signal_df["LibraryID"]==libid)&\
            (ntc_signal_df["TAG1"]==t1)&(ntc_signal_df["TAG3"]==t2)].index,["Confirmation(1/0)"]]=1
        ntc_signal_df.loc[ntc_signal_df[(ntc_signal_df["SampleID"]==sampleid)&(ntc_signal_df["LibraryID"]==libid)&\
            (ntc_signal_df["TAG2"]==t1)&(ntc_signal_df["TAG3"]==t2)].index,["Confirmation(1/0)"]]=1

    # NTC signal(plane)
    ntc_format_signal_df1 = ntc_format_signal_df[ntc_format_signal_df["X1"].str.contains("1.",regex = False)]
    ntc_format_signal_df2_1 = pd.DataFrame(ntc_format_signal_df[ntc_format_signal_df["X1"].str.contains("2.",regex = False)])
    ntc_format_signal_df2_2 = pd.DataFrame(ntc_format_signal_df[ntc_format_signal_df["X2"].str.contains("2.",regex = False)])
    ntc_format_signal_df2_2.rename(columns={"X1":"X2","X2":"X1"},inplace=True)

    ntc_format_signal_df3 = ntc_format_signal_df[ntc_format_signal_df["X2"].str.contains("3.",regex = False)]
    ntc_format_signal_df3.rename(columns={"X1":"X2","X2":"X1"},inplace=True)

    ntc_format_signal_df4 = pd.concat([ntc_format_signal_df1,ntc_format_signal_df2_1,ntc_format_signal_df2_2,ntc_format_signal_df3])
    ntc_format_signal_df4.reset_index(inplace=True)
    
    ntc_planes_df = pd.DataFrame(ntc_format_signal_df4.groupby(["SampleID","LibraryID","X1"])["Predict"].sum())
    ntc_planes_df = ntc_planes_df[ntc_planes_df["Predict"] >= 2]
    ntc_planes_df.reset_index(inplace=True)
    
    for i in range(len(ntc_planes_df)):
        sampleid = ntc_planes_df.loc[i,"SampleID"]
        libid = ntc_planes_df.loc[i,"LibraryID"]
        t1 = ntc_planes_df.loc[i,"X1"]
        ntc_signal_df.loc[ntc_signal_df[(ntc_signal_df["SampleID"]==sampleid)&(ntc_signal_df["LibraryID"]==libid)&\
            (ntc_signal_df["DLP"]=="P")&(ntc_signal_df["TAG1"]==t1)].index,"Confirmation(1/0)"]=1
        ntc_signal_df.loc[ntc_signal_df[(ntc_signal_df["SampleID"]==sampleid)&(ntc_signal_df["LibraryID"]==libid)&\
            (ntc_signal_df["DLP"]=="P")&(ntc_signal_df["TAG2"]==t1)].index,"Confirmation(1/0)"]=1
        ntc_signal_df.loc[ntc_signal_df[(ntc_signal_df["SampleID"]==sampleid)&(ntc_signal_df["LibraryID"]==libid)&\
            (ntc_signal_df["DLP"]=="P")&(ntc_signal_df["TAG3"]==t1)].index,"Confirmation(1/0)"]=1
    ntc_signal_df[ntc_signal_df["DLP"]=="D"]["Confirmation(1/0)"]=1
    ntc_signal_df.reset_index(inplace=True)

    # case signal(dot,line)
    for i in range(len(case_format_signal_df)):
        sampleid = case_format_signal_df.loc[i,"SampleID"]
        libid = case_format_signal_df.loc[i,"LibraryID"]
        t1 = case_format_signal_df.loc[i,"X1"]
        t2 = case_format_signal_df.loc[i,"X2"]
        tp_case_signal_df.loc[tp_case_signal_df[(tp_case_signal_df["SampleID"]==sampleid)&(tp_case_signal_df["LibraryID"]==libid)&\
            (tp_case_signal_df["TAG1"]==t1)&(tp_case_signal_df["TAG2"]==t2)].index,"Confirmation(1/0)"]=1
        tp_case_signal_df.loc[tp_case_signal_df[(tp_case_signal_df["SampleID"]==sampleid)&(tp_case_signal_df["LibraryID"]==libid)&\
            (tp_case_signal_df["TAG1"]==t1)&(tp_case_signal_df["TAG3"]==t2)].index,["Confirmation(1/0)"]]=1
        tp_case_signal_df.loc[tp_case_signal_df[(tp_case_signal_df["SampleID"]==sampleid)&(tp_case_signal_df["LibraryID"]==libid)&\
            (tp_case_signal_df["TAG2"]==t1)&(tp_case_signal_df["TAG3"]==t2)].index,["Confirmation(1/0)"]]=1

    # case signal(plane)
    case_format_signal_df1 = case_format_signal_df[case_format_signal_df["X1"].str.contains("1.",regex = False)]
    case_format_signal_df2_1 = pd.DataFrame(case_format_signal_df[case_format_signal_df["X1"].str.contains("2.",regex = False)])
    case_format_signal_df2_2 = pd.DataFrame(case_format_signal_df[case_format_signal_df["X2"].str.contains("2.",regex = False)])
    case_format_signal_df2_2.rename(columns={"X1":"X2","X2":"X1"},inplace=True)

    case_format_signal_df3 = case_format_signal_df[case_format_signal_df["X2"].str.contains("3.",regex = False)]
    case_format_signal_df3.rename(columns={"X1":"X2","X2":"X1"},inplace=True)

    case_format_signal_df4 = pd.concat([case_format_signal_df1,case_format_signal_df2_1,case_format_signal_df2_2,case_format_signal_df3])
    case_format_signal_df4.reset_index(inplace=True)
    
    case_planes_df = pd.DataFrame(case_format_signal_df4.groupby(["SampleID","LibraryID","X1"])["Predict"].sum())
    case_planes_df = case_planes_df[case_planes_df["Predict"] >= 2]
    case_planes_df.reset_index(inplace=True)
    
    for i in range(len(case_planes_df)):
        sampleid = case_planes_df.loc[i,"SampleID"]
        libid = case_planes_df.loc[i,"LibraryID"]
        t1 = case_planes_df.loc[i,"X1"]
        tp_case_signal_df.loc[tp_case_signal_df[(tp_case_signal_df["SampleID"]==sampleid)&(tp_case_signal_df["LibraryID"]==libid)&\
            (tp_case_signal_df["DLP"]=="P")&(tp_case_signal_df["TAG1"]==t1)].index,"Confirmation(1/0)"]=1
        tp_case_signal_df.loc[tp_case_signal_df[(tp_case_signal_df["SampleID"]==sampleid)&(tp_case_signal_df["LibraryID"]==libid)&\
            (tp_case_signal_df["DLP"]=="P")&(tp_case_signal_df["TAG2"]==t1)].index,"Confirmation(1/0)"]=1
        tp_case_signal_df.loc[tp_case_signal_df[(tp_case_signal_df["SampleID"]==sampleid)&(tp_case_signal_df["LibraryID"]==libid)&\
            (tp_case_signal_df["DLP"]=="P")&(tp_case_signal_df["TAG3"]==t1)].index,"Confirmation(1/0)"]=1
    tp_case_signal_df.reset_index(inplace=True)

    # maked
    ntc_signal_confir_df = ntc_signal_df[ntc_signal_df["Confirmation(1/0)"]==1]
    ntc_signal_confir_df.reset_index(inplace=True,drop=True)
    signal_predict_df = pd.concat([ntc_signal_df,tp_case_signal_df,tp_case_signal_D_df,fp_case_signal_df])
    signal_predict_df.drop(["index","level_0"],axis=1,inplace=True)
    signal_predict_df.to_csv(out_file,index=None,header=True,sep="\t")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-s","--saminfo",help="sample info file")
    parser.add_argument("-i1","--dlp",help="dlp file")
    parser.add_argument("-i2","--dwar",help="dwar path")
    parser.add_argument("-m","--model",help="model path")
    parser.add_argument("-o","--out",help="out file")
    args = parser.parse_args()
    predict_signal(args.saminfo,args.dlp,args.dwar,args.model,args.out)