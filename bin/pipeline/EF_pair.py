#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
@Desc: enrichment analysis based paired samples
"""

import os, warnings, argparse, sys, time
import pandas as pd
import shutil
from io import StringIO
sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), '../monitor'))
import process
warnings.filterwarnings("ignore")

def enrichment_score(case,ntc):
    EF = 0
    if float(case) >0 and float(ntc)>0:
        EF = float(case)/float(ntc)
    else:
        EF = 0
    return EF
        
def merge_pair_score(casefile, ntcfile,ntc_name,ana_dir, proteinName_pair,sample, lib, pair_name_lst):
    #pair_dict = {}
    print(casefile,ntcfile,ntc_name,ana_dir,proteinName_pair,sample,lib,pair_name_lst)
    if len(pair_name_lst) != 0:
        # when no paired sample, pair_name_lst=['0']
        if len(pair_name_lst) ==1 and pair_name_lst[0] == '0':
            # merge case and ntc sample
            if os.path.exists(casefile):
                case_df = pd.read_csv(casefile, index_col=None, header=0, sep="\t",dtype=object)
                case_df.drop(['Raw_count','TAG1_P','TAG2_P','TAG3_P','TAG1_TAG2_L','TAG1_TAG3_L','TAG2_TAG3_L','TAG_DOT','TAG1_P_noline','TAG2_P_noline','TAG3_P_noline','Normalized_z-score','Y','Y0'], axis=1, inplace=True)
                case_name2 = '(' + sample + ')'
                if 'EF' not in case_df.columns:
                    case_df['EF'] = case_df['z-score']
                case_df.rename(columns={'Dedup_count':'Dedup_count'+case_name2,'z-score':'z-score'+case_name2,'EF':'EF'+case_name2}, inplace=True)
                if os.path.exists(ntcfile):
                    ntc_df = pd.read_csv(ntcfile, index_col=None, header=0, sep="\t", dtype=object)
                    ntc_df.drop(['Raw_count','TAG1_P','TAG2_P','TAG3_P','TAG1_TAG2_L','TAG1_TAG3_L','TAG2_TAG3_L','TAG_DOT','TAG1_P_noline','TAG2_P_noline','TAG3_P_noline','Normalized_z-score','Y','Y0'], axis=1, inplace=True)
                    ntc_name2 = '(' + ntc_name + ')'
                    ntc_df.rename(columns={'Dedup_count':'Dedup_count'+ntc_name2,'z-score':'z-score'+ntc_name2}, inplace=True)
                    ntc_df1 = ntc_df[['Tags','TAG1','TAG1_structure','TAG1_FF','TAG1_BB',
                                'TAG2','TAG2_structure','TAG2_FF','TAG2_BB',
                                'TAG3','TAG3_structure','TAG3_FF','TAG3_BB','Dedup_count'+ntc_name2,'z-score'+ntc_name2]]
                    ntc_case_df = pd.merge(case_df, ntc_df1,on=['Tags','TAG1','TAG1_structure','TAG1_FF','TAG1_BB',
                                'TAG2','TAG2_structure','TAG2_FF','TAG2_BB',
                                'TAG3','TAG3_structure','TAG3_FF','TAG3_BB'],how="outer")
                    ntc_case_df.fillna(0,inplace=True)
                    ntc_case_df["Type"] = ["0"] * len(ntc_case_df)
                    ntc_case_df.loc[(ntc_case_df["TAG1_BB"] != "0") & (ntc_case_df["TAG2_BB"] != "0") & (ntc_case_df["TAG3_BB"] != "0"), "Type"] = "123"
                    ntc_case_df.loc[(ntc_case_df["TAG1_BB"] != "0") & (ntc_case_df["TAG2_BB"] != "0") & (ntc_case_df["TAG3_BB"] == "0"), "Type"] = "12"
                    ntc_case_df.loc[(ntc_case_df["TAG1_BB"] != "0") & (ntc_case_df["TAG2_BB"] == "0") & (ntc_case_df["TAG3_BB"] != "0"), "Type"] = "13"
                    ntc_case_df.loc[(ntc_case_df["TAG1_BB"] == "0") & (ntc_case_df["TAG2_BB"] != "0") & (ntc_case_df["TAG3_BB"] != "0"), "Type"] = "23"
                    ntc_case_df.loc[(ntc_case_df["TAG1_BB"] != "0") & (ntc_case_df["TAG2_BB"] == "0") & (ntc_case_df["TAG3_BB"] == "0"), "Type"] = "1"
                    ntc_case_df.loc[(ntc_case_df["TAG1_BB"] == "0") & (ntc_case_df["TAG2_BB"] != "0") & (ntc_case_df["TAG3_BB"] == "0"), "Type"] = "2"
                    ntc_case_df.loc[(ntc_case_df["TAG1_BB"] == "0") & (ntc_case_df["TAG2_BB"] == "0") & (ntc_case_df["TAG3_BB"] != "0"), "Type"] = "3"
                    ntc_case_df['EF'+case_name2] = ntc_case_df.apply(lambda row:enrichment_score(row['z-score'+case_name2], row['z-score'+ntc_name2]), axis=1)
                    print('insert_libid')
                    ntc_case_df.insert(loc=0,column='LibraryID',value=lib)
                    ntc_case_df.insert(loc=0,column='SampleID',value=sample)
                    print(ntc_case_df)
                    
                    new_columns=['SampleID','LibraryID','Tags','TAG1','TAG1_structure','TAG1_FF','TAG1_BB',
                                'TAG2','TAG2_structure','TAG2_FF','TAG2_BB',
                                'TAG3','TAG3_structure','TAG3_FF','TAG3_BB',
                                'Dedup_count'+ntc_name2,'z-score'+ntc_name2,
                                'Dedup_count'+case_name2,'z-score'+case_name2,'EF'+case_name2,
                                'Type']
                    
                    ntc_case_df_new = ntc_case_df[new_columns]
                    return ntc_case_df_new
                else:
                    case_df.fillna(0,inplace=True)
                    case_df["Type"] = ["0"] * len(case_df)
                    case_df.loc[(case_df["TAG1_BB"] != "0") & (case_df["TAG2_BB"] != "0") & (case_df["TAG3_BB"] != "0"), "Type"] = "123"
                    case_df.loc[(case_df["TAG1_BB"] != "0") & (case_df["TAG2_BB"] != "0") & (case_df["TAG3_BB"] == "0"), "Type"] = "12"
                    case_df.loc[(case_df["TAG1_BB"] != "0") & (case_df["TAG2_BB"] == "0") & (case_df["TAG3_BB"] != "0"), "Type"] = "13"
                    case_df.loc[(case_df["TAG1_BB"] == "0") & (case_df["TAG2_BB"] != "0") & (case_df["TAG3_BB"] != "0"), "Type"] = "23"
                    case_df.loc[(case_df["TAG1_BB"] != "0") & (case_df["TAG2_BB"] == "0") & (case_df["TAG3_BB"] == "0"), "Type"] = "1"
                    case_df.loc[(case_df["TAG1_BB"] == "0") & (case_df["TAG2_BB"] != "0") & (case_df["TAG3_BB"] == "0"), "Type"] = "2"
                    case_df.loc[(case_df["TAG1_BB"] == "0") & (case_df["TAG2_BB"] == "0") & (case_df["TAG3_BB"] != "0"), "Type"] = "3"
                    case_df.insert(loc=0,column='LibraryID',value=lib)
                    case_df.insert(loc=0,column='SampleID',value=sample)
                    
                    new_columns=['SampleID','LibraryID','Tags','TAG1','TAG1_structure','TAG1_FF','TAG1_BB',
                                'TAG2','TAG2_structure','TAG2_FF','TAG2_BB',
                                'TAG3','TAG3_structure','TAG3_FF','TAG3_BB',
                                'Dedup_count'+case_name2,'z-score'+case_name2,'EF'+case_name2,
                                'Type']
                    
                    case_df_new = case_df[new_columns]
                    return case_df_new
        else:
            if os.path.exists(casefile) and os.path.exists(ntcfile):
                case_df = pd.read_csv(casefile, index_col=None, header=0, sep="\t",dtype=object)
                case_df.drop(['Raw_count','TAG1_P','TAG2_P','TAG3_P','TAG1_TAG2_L','TAG1_TAG3_L','TAG2_TAG3_L','TAG_DOT','TAG1_P_noline','TAG2_P_noline','TAG3_P_noline','Normalized_z-score','Y','Y0'], axis=1, inplace=True)
                case_name2 = '(' + sample + ')'
                case_df.rename(columns={'Dedup_count':'Dedup_count'+case_name2,'z-score':'z-score'+case_name2,'EF':'EF'+case_name2}, inplace=True)
                
                ntc_df = pd.read_csv(ntcfile, index_col=None, header=0, sep="\t", dtype=object)
                #print(ntc_df.columns)
                ntc_df.drop(['Raw_count','TAG1_P','TAG2_P','TAG3_P','TAG1_TAG2_L','TAG1_TAG3_L','TAG2_TAG3_L','TAG_DOT','TAG1_P_noline','TAG2_P_noline','TAG3_P_noline','Normalized_z-score','Y','Y0'], axis=1, inplace=True)
                ntc_name2 = '(' + ntc_name + ')'
                ntc_df.rename(columns={'Dedup_count':'Dedup_count'+ntc_name2,'z-score':'z-score'+ntc_name2}, inplace=True)
                ntc_df1 = ntc_df[['Tags','TAG1','TAG1_structure','TAG1_FF','TAG1_BB',
                            'TAG2','TAG2_structure','TAG2_FF','TAG2_BB',
                            'TAG3','TAG3_structure','TAG3_FF','TAG3_BB','Dedup_count'+ntc_name2,'z-score'+ntc_name2]]
                try:
                    ntc_case_df = pd.merge(case_df, ntc_df1,on=['Tags','TAG1','TAG1_structure','TAG1_FF','TAG1_BB',
                            'TAG2','TAG2_structure','TAG2_FF','TAG2_BB',
                            'TAG3','TAG3_structure','TAG3_FF','TAG3_BB'],how="outer")
                except:
                    ntc_case_df = pd.DataFrame(columns=list(case_df.columns)+['Dedup_count'+ntc_name2,'z-score'+ntc_name2])

                print('pair_name_lst: {}'.format(pair_name_lst))
                for pair_name in pair_name_lst:
                    pair_path = '{}/{}/{}/{}.{}.count.freq.filter.txt'.format(ana_dir,proteinName_pair,lib,pair_name,lib)
                    if os.path.exists(pair_path):
                        pair_df = pd.read_csv(pair_path, index_col=None, header=0, sep="\t",dtype=object)
                        pair_df.drop(['Raw_count','Normalized_z-score','TAG1_P','TAG2_P','TAG3_P','TAG1_TAG2_L','TAG1_TAG3_L','TAG2_TAG3_L','TAG_DOT','TAG1_P_noline','TAG2_P_noline','TAG3_P_noline','Y','Y0'], axis=1, inplace=True)
                        pair_name2 = '(' + pair_name + ')'
                        pair_df.rename(columns={'Dedup_count':'Dedup_count'+pair_name2,'z-score':'z-score'+pair_name2,"EF": "EF" + pair_name2}, inplace=True)
                        try:
                            ntc_case_df = pd.merge(ntc_case_df,pair_df,on=['Tags','TAG1','TAG1_structure','TAG1_FF','TAG1_BB',
                            'TAG2','TAG2_structure','TAG2_FF','TAG2_BB',
                            'TAG3','TAG3_structure','TAG3_FF','TAG3_BB','Type'],how="outer")
                        except:
                            ntc_case_df=pd.DataFrame(columns=list(ntc_case_df.columns)+['Dedup_count'+pair_name2,'z-score'+pair_name2,"EF" + pair_name2])
                        
                        ntc_case_df.fillna(0, inplace=True)
                        ntc_case_df["Type"] = ["0"] * len(ntc_case_df)
                        ntc_case_df.loc[(ntc_case_df["TAG1_BB"] != "0") & (ntc_case_df["TAG2_BB"] != "0") & (ntc_case_df["TAG3_BB"] != "0"), "Type"] = "123"
                        ntc_case_df.loc[(ntc_case_df["TAG1_BB"] != "0") & (ntc_case_df["TAG2_BB"] != "0") & (ntc_case_df["TAG3_BB"] == "0"), "Type"] = "12"
                        ntc_case_df.loc[(ntc_case_df["TAG1_BB"] != "0") & (ntc_case_df["TAG2_BB"] == "0") & (ntc_case_df["TAG3_BB"] != "0"), "Type"] = "13"
                        ntc_case_df.loc[(ntc_case_df["TAG1_BB"] == "0") & (ntc_case_df["TAG2_BB"] != "0") & (ntc_case_df["TAG3_BB"] != "0"), "Type"] = "23"
                        ntc_case_df.loc[(ntc_case_df["TAG1_BB"] != "0") & (ntc_case_df["TAG2_BB"] == "0") & (ntc_case_df["TAG3_BB"] == "0"), "Type"] = "1"
                        ntc_case_df.loc[(ntc_case_df["TAG1_BB"] == "0") & (ntc_case_df["TAG2_BB"] != "0") & (ntc_case_df["TAG3_BB"] == "0"), "Type"] = "2"
                        ntc_case_df.loc[(ntc_case_df["TAG1_BB"] == "0") & (ntc_case_df["TAG2_BB"] == "0") & (ntc_case_df["TAG3_BB"] != "0"), "Type"] = "3"
                        ntc_case_df['EF'] = ntc_case_df.apply(lambda row:enrichment_score(row['z-score'], row['NTC_z-score']), axis=1)
                        print('insert column libid')
                        ntc_case_df.insert(loc=0,column='LibraryID',value=lib)
                        ntc_case_df.insert(loc=0,column='SampleID',value=sample)
                        print(ntc_case_df)
                        new_columns=['SampleID','LibraryID','Tags','TAG1','TAG1_structure','TAG1_FF','TAG1_BB',
                                    'TAG2','TAG2_structure','TAG2_FF','TAG2_BB',
                                    'TAG3','TAG3_structure','TAG3_FF','TAG3_BB',
                                    'Dedup_count'+ntc_name2,'z-score'+ntc_name2,
                                    'Dedup_count'+case_name2,'z-score'+case_name2,'EF'+case_name2,
                                    'Dedup_count'+pair_name2,'z-score'+pair_name2,"EF" + pair_name2,
                                    'Type']
                        #print(ntc_case_df.columns)
                        ntc_case_df = ntc_case_df[new_columns]
                        print(ntc_case_df)
                return ntc_case_df
            else:
                print("缺少Case或NTC的Normalized_z-score文件, {0} {1}！".format(casefile,ntcfile))

def txt2dwar_pair(ratio_file, dwar_file, dedup_filter_threshold,codesmilesfile,pair_name_lst,ntc_lst,sample):
    # dwar
    count_f = pd.read_csv(ratio_file, sep='\t', dtype=object, header=0)
    count_f2 = count_f[count_f['Tags'] != '0']
    columnNo = len(count_f2.columns)
    old = StringIO(count_f2.to_csv(sep='\t',index=False))

# head
    tag_header = '''<datawarrior-fileinfo>
<version="3.2">
<created="1559631713662">
<rowcount="{0}">
</datawarrior-fileinfo>
'''.format(count_f2.shape[0])
    if ntc_lst !='1' and ntc_lst!='0':
        if pair_name_lst[0] != '1' and pair_name_lst[0] != '0':
            tag_tail = '''<datawarrior properties>
<axisColumn_2D EF_0="Tags">
<axisColumn_2D EF_1="EF({1})">
<axisColumn_{1}_0="TAG1">
<axisColumn_{1}_1="TAG2">
<axisColumn_{1}_2="TAG3">
<chartType_2D EF="scatter">
<chartType_{0}="scatter">
<chartType_{1}="scatter">
<chartType_{2}="scatter">
<colorColumn_2D EF="EF({1})">
<colorListMode_2D EF="HSBLong">
<color_2D EF_0="-65536">
<color_2D EF_1="-16776961">
<columnWidth_Table_Dedup_count({0})="210">
<columnWidth_Table_Dedup_count({1})="210">
<columnWidth_Table_Dedup_count({2})="210">
<columnWidth_Table_EF({1})="210">
<columnWidth_Table_EF({2})="210">
<columnWidth_Table_z-score({0})="210">
<columnWidth_Table_z-score({1})="210">
<columnWidth_Table_z-score({2})="210">
<columnWidth_Table_Library_ID="80">
<columnWidth_Table_Protein_CC_ID="80">
<columnWidth_Table_TAG1="80">
<columnWidth_Table_TAG1_BB="80">
<columnWidth_Table_TAG1_FF="80">
<columnWidth_Table_TAG2="80">
<columnWidth_Table_TAG2_BB="80">
<columnWidth_Table_TAG2_BB="80">
<columnWidth_Table_TAG3="80">
<columnWidth_Table_TAG3_BB="80">
<columnWidth_Table_TAG3_FF="80">
<columnWidth_Table_Tags="80">
<columnWidth_Table_Type="80">
<detailView="height[Data]=1">
<faceColor3D_{0}="-1250054">
<faceColor3D_{1}="-1250054">
<faceColor3D_{2}="-1250054">
<fastRendering_2D EF="true">
<fastRendering_{0}="true">
<fastRendering_{1}="true">
<fastRendering_{2}="true">
<filter0="#category#	SampleID	-">
<filter1="#string#	Tags">
<filter2="#double#	TAG1">
<filter3="#string#	TAG1_BB">
<filter4="#double#	TAG2">
<filter5="#string#	TAG2_BB">
<filter6="#double#	TAG3">
<filter7="#string#	TAG3_BB">
<filter8="#double#	TAG4">
<filter9="#string#	TAG4_BB">
<filter10="#double#	Dedup_count({1})">
<filter11="#double#	z_score({1})">
<filter12="#double#	EF({1})">
<filter13="#double#	Dedup_count({0})">
<filter14="#double#	z_score({0})">
<filter15="#double#	EF({0})">
<filter16="#double#	Dedup_count({2})">
<filter17="#double#	z_score({2})">
<filter18="#double#	EF({2})">
<mainSplitting="0.72255">
<mainView="Table">
<mainViewCount="5">
<mainViewDockInfo0="root">
<mainViewDockInfo1="Table	center">
<mainViewDockInfo2="2D EF	bottom	0.503">
<mainViewDockInfo3="{0}	right	0.32">
<mainViewDockInfo4="{1}	right	0.50">
<mainViewName0="Table">
<mainViewName1="2D EF">
<mainViewName2="{0}">
<mainViewName3="{1}">
<mainViewName4="{2}">
<mainViewType0="tableView">
<mainViewType1="2Dview">
<mainViewType2="3Dview">
<mainViewType3="3Dview">
<mainViewType4="3Dview">
<masterView_{0}="{1}">
<masterView_{2}="{1}">
<rightSplitting="0.62984">
<rotationMatrix_{1}00="0.99733">
<rotationMatrix_{1}01="-0.066506">
<rotationMatrix_{1}02="-0.030133">
<rotationMatrix_{1}10="0.07021">
<rotationMatrix_{1}11="0.9868">
<rotationMatrix_{1}12="0.14591">
<rotationMatrix_{1}20="0.02003">
<rotationMatrix_{1}21="-0.14764">
<rotationMatrix_{1}22="0.98885">
<rowHeight_Table="16">
<showNaNValues_2D EF="true">
<showNaNValues_{0}="true">
<showNaNValues_{1}="true">
<showNaNValues_{2}="true">
<sizeColumn_{0}="z-score({0})">
<sizeColumn_{1}="z-score({1})">
<sizeColumn_{2}="z-score({2})">
</datawarrior properties>
    
    '''.format(ntc_lst,sample,pair_name_lst[0])
            #print(pair_name_lst[0])
            #print(tag_tail)
        else:
            tag_tail = '''<datawarrior properties>
<axisColumn_2D EF_0="Tags">
<axisColumn_2D EF_1="EF({1})">
<axisColumn_{1}_0="TAG1">
<axisColumn_{1}_1="TAG2">
<axisColumn_{1}_2="TAG3">
<chartType_2D EF="scatter">
<chartType_{0}="scatter">
<chartType_{1}="scatter">
<colorColumn_2D EF="EF({1})">
<colorListMode_2D EF="HSBLong">
<color_2D EF_0="-65536">
<color_2D EF_1="-16776961">
<columnWidth_Table_Dedup_count({0})="210">
<columnWidth_Table_Dedup_count({1})="210">
<columnWidth_Table_EF({1})="210">
<columnWidth_Table_z-score({0})="210">
<columnWidth_Table_z-score({1})="210
<columnWidth_Table_Library_ID="80">
<columnWidth_Table_Protein_CC_ID="80">
<columnWidth_Table_TAG1="80">
<columnWidth_Table_TAG1_BB="80">
<columnWidth_Table_TAG1_FF="80">
<columnWidth_Table_TAG2="80">
<columnWidth_Table_TAG2_BB="80">
<columnWidth_Table_TAG2_BB="80">
<columnWidth_Table_TAG3="80">
<columnWidth_Table_TAG3_BB="80">
<columnWidth_Table_TAG3_FF="80">
<columnWidth_Table_Tags="80">
<columnWidth_Table_Type="80">
<detailView="height[Data]=1">
<faceColor3D_{0}="-1250054">
<faceColor3D_{1}="-1250054">
<fastRendering_2D EF="true">
<fastRendering_{0}="true">
<fastRendering_{1}="true">
<filter0="#category#	SampleID	-">
<filter1="#string#	Tags">
<filter2="#double#	TAG1">
<filter3="#string#	TAG1_BB">
<filter4="#double#	TAG2">
<filter5="#string#	TAG2_BB">
<filter6="#double#	TAG3">
<filter7="#string#	TAG3_BB">
<filter8="#double#	TAG4">
<filter9="#string#	TAG4_BB">
<filter10="#double#	Dedup_count({1})">
<filter11="#double#	z_score({1})">
<filter12="#double#	EF({1})">
<filter13="#double#	Dedup_count({0})">
<filter14="#double#	z_score({0})">
<filter15="#double#	EF({0})">
<mainSplitting="0.72255">
<mainView="Table">
<mainViewCount="5">
<mainViewDockInfo0="root">
<mainViewDockInfo1="Table	center">
<mainViewDockInfo2="2D EF	bottom	0.503">
<mainViewDockInfo3="{0}	right	0.50">
<mainViewName0="Table">
<mainViewName1="2D EF">
<mainViewName2="{0}">
<mainViewName3="{1}">
<mainViewType0="tableView">
<mainViewType1="2Dview">
<mainViewType2="3Dview">
<mainViewType3="3Dview">
<mainViewType4="3Dview">
<masterView_{0}="{1}">
<rightSplitting="0.62984">
<rotationMatrix_{1}00="0.99733">
<rotationMatrix_{1}01="-0.066506">
<rotationMatrix_{1}02="-0.030133">
<rotationMatrix_{1}10="0.07021">
<rotationMatrix_{1}11="0.9868">
<rotationMatrix_{1}12="0.14591">
<rotationMatrix_{1}20="0.02003">
<rotationMatrix_{1}21="-0.14764">
<rotationMatrix_{1}22="0.98885">
<rowHeight_Table="16">
<showNaNValues_2D EF="true">
<showNaNValues_{0}="true">
<showNaNValues_{1}="true">
<sizeColumn_{0}="z-score({0})">
<sizeColumn_{1}="z-score({1})">
</datawarrior properties>
    
    '''.format(ntc_lst,sample)
    else:
        tag_tail = '''<datawarrior properties>
<axisColumn_2D EF_0="Tags">
<axisColumn_2D EF_1="EF({0})">
<axisColumn_{0}_0="TAG1">
<axisColumn_{0}_1="TAG2">
<axisColumn_{0}_2="TAG3">
<chartType_2D EF="scatter">
<chartType_{0}="scatter">
<colorColumn_2D EF="EF({0})">
<colorListMode_2D EF="HSBLong">
<color_2D EF_0="-65536">
<color_2D EF_1="-16776961">
<columnWidth_Table_Dedup_count({0})="210">
<columnWidth_Table_EF({0})="210">
<columnWidth_Table_z-score({0})="210">
<columnWidth_Table_Library_ID="80">
<columnWidth_Table_Protein_CC_ID="80">
<columnWidth_Table_TAG1="80">
<columnWidth_Table_TAG1_BB="80">
<columnWidth_Table_TAG1_FF="80">
<columnWidth_Table_TAG2="80">
<columnWidth_Table_TAG2_BB="80">
<columnWidth_Table_TAG2_BB="80">
<columnWidth_Table_TAG3="80">
<columnWidth_Table_TAG3_BB="80">
<columnWidth_Table_TAG3_FF="80">
<columnWidth_Table_Tags="80">
<columnWidth_Table_Type="80">
<detailView="height[Data]=1">
<faceColor3D_{0}="-1250054">
<fastRendering_2D EF="true">
<fastRendering_{0}="true">
<filter0="#category#	SampleID	-">
<filter1="#string#	Tags">
<filter2="#double#	TAG1">
<filter3="#string#	TAG1_BB">
<filter4="#double#	TAG2">
<filter5="#string#	TAG2_BB">
<filter6="#double#	TAG3">
<filter7="#string#	TAG3_BB">
<filter8="#double#	TAG4">
<filter9="#string#	TAG4_BB">
<filter10="#double#	Dedup_count({0})">
<filter11="#double#	z_score({0})">
<filter12="#double#	EF({0})">
<mainSplitting="0.72255">
<mainView="Table">
<mainViewCount="5">
<mainViewDockInfo0="root">
<mainViewDockInfo1="Table	center">
<mainViewDockInfo2="2D EF	bottom	0.503">
<mainViewDockInfo3="{0}	right	0.50">
<mainViewName0="Table">
<mainViewName1="2D EF">
<mainViewName2="{0}">
<mainViewType0="tableView">
<mainViewType1="2Dview">
<mainViewType2="3Dview">
<mainViewType3="3Dview">
<mainViewType4="3Dview">
<rightSplitting="0.62984">
<rotationMatrix_{0}00="0.99733">
<rotationMatrix_{0}01="-0.066506">
<rotationMatrix_{0}02="-0.030133">
<rotationMatrix_{0}10="0.07021">
<rotationMatrix_{0}11="0.9868">
<rotationMatrix_{0}12="0.14591">
<rotationMatrix_{0}20="0.02003">
<rotationMatrix_{0}21="-0.14764">
<rotationMatrix_{0}22="0.98885">
<rowHeight_Table="16">
<showNaNValues_2D EF="true">
<showNaNValues_{0}="true">
<sizeColumn_{0}="z-score({0})">
</datawarrior properties>
    
    '''.format(sample)
    # dwar
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
def codesmilesSeq(codesmilesfile,columnNo,count_f):
    taglist_min = []
    taglist_max = []
    count_f = count_f.fillna(0)
    codesmile = pd.read_csv(codesmilesfile,sep='\t',header=0,dtype='object')
    codesmile['CycleNum'] = codesmile['Cycle'] + '.' + codesmile['Number']
    cycles = list(set(codesmile['Cycle'].tolist()))
    for cycle in cycles:
        
        eachcycle = codesmile[codesmile['Cycle']==cycle]
        taglist_max.append(eachcycle.values[-1].tolist()[-1])
        taglist_min.append(eachcycle.values[0].tolist()[-1])
        
    taglist_min.sort()
    taglist_max.sort()
    
    count_lst = count_f.columns
    
    dedup_col = [i for i in count_lst if 'Dedup_count' in i]
    zscore_col = [i for i in count_lst if 'z-score' in i]
    EF_col = [i for i in count_lst if 'EF' in i]
    min_1 = []
    max_1 = []
    
    if len(dedup_col) == 3 and len(zscore_col)==3 and len(EF_col) == 2:
        count_f[zscore_col[0]] = count_f[zscore_col[0]].astype(float)
        count_f[zscore_col[1]] = count_f[zscore_col[1]].astype(float)
        count_f[zscore_col[2]] = count_f[zscore_col[2]].astype(float)
        count_f[EF_col[0]] = count_f[EF_col[0]].astype(float)
        count_f[EF_col[1]] = count_f[EF_col[1]].astype(float)
        count_f[dedup_col[0]] = count_f[dedup_col[0]].astype(int)
        count_f[dedup_col[1]] = count_f[dedup_col[1]].astype(int)
        count_f[dedup_col[2]] = count_f[dedup_col[2]].astype(int)
        
        dedup_lst_min = min([count_f[dedup_col[0]].min(),count_f[dedup_col[1]].min(),count_f[dedup_col[2]].min()])
        dedup_lst_max = max([count_f[dedup_col[0]].max(),count_f[dedup_col[1]].max(),count_f[dedup_col[2]].max()])
        
        z_score_lst_min = min([count_f[zscore_col[0]].min(),count_f[zscore_col[1]].min(),count_f[zscore_col[2]].min()])
        z_score_lst_max = max([count_f[zscore_col[0]].max(),count_f[zscore_col[1]].max(),count_f[zscore_col[2]].max()])
        #print(z_score_lst_max)
        
        ef_min = min([count_f[EF_col[0]].min(),count_f[EF_col[1]].min()])
        ef_max = min([count_f[EF_col[0]].max(),count_f[EF_col[1]].max()])
        min_1 = [dedup_lst_min,z_score_lst_min,dedup_lst_min,z_score_lst_min,ef_min,dedup_lst_min,z_score_lst_min,ef_min,0]
        max_1 = [dedup_lst_max,z_score_lst_max,dedup_lst_max,z_score_lst_max,ef_max,dedup_lst_max,z_score_lst_max,ef_max,0]
    elif len(dedup_col) == 2 and len(zscore_col)==2 and len(EF_col) == 1:
        count_f[zscore_col[0]] = count_f[zscore_col[0]].astype(float)
        count_f[zscore_col[1]] = count_f[zscore_col[1]].astype(float)
        count_f[dedup_col[0]] = count_f[dedup_col[0]].astype(int)
        count_f[dedup_col[1]] = count_f[dedup_col[1]].astype(int)
        dedup_lst_min = min([count_f[dedup_col[0]].min(),count_f[dedup_col[1]].min()])
        dedup_lst_max = max([count_f[dedup_col[0]].max(),count_f[dedup_col[1]].max()])
        
        z_score_lst_min = min([count_f[zscore_col[0]].min(),count_f[zscore_col[1]].min()])
        z_score_lst_max = max([count_f[zscore_col[0]].max(),count_f[zscore_col[1]].max()])
        
        ef_min = count_f[EF_col[0]].min()
        ef_max = count_f[EF_col[0]].max()
        min_1 = [dedup_lst_min,z_score_lst_min,dedup_lst_min,z_score_lst_min,ef_min,0]
        max_1 = [dedup_lst_max,z_score_lst_max,dedup_lst_max,z_score_lst_max,ef_max,0]
    
    else:
        print('ERROR: This is not the pair dwar!!!')

    min_1 = [str(i) for i in min_1]
    max_1 = [str(i) for i in max_1]
    
    if 'TAG1_seq' in count_lst:
        lineminlst_min = ['-','-','-',taglist_min[0],'-','-','-','-',taglist_min[1],'-','-','-','-',taglist_min[2],'-','-','-','-']
        lineminlst_max = ['-','-','-',taglist_max[0],'-','-','-','-',taglist_max[1],'-','-','-','-',taglist_max[2],'-','-','-','-']
    else:
        lineminlst_min = ['-','-','-',taglist_min[0],'-','-','-',taglist_min[1],'-','-','-',taglist_min[2],'-','-','-']
        lineminlst_max = ['-','-','-',taglist_max[0],'-','-','-',taglist_max[1],'-','-','-',taglist_max[2],'-','-','-']
    lineminlst_min = [str(i) for i in lineminlst_min]
    lineminlst_max = [str(i) for i in lineminlst_max]
    
    tags_min = '\t'.join(lineminlst_min + min_1)+'\n'
    tags_max = '\t'.join(lineminlst_max + max_1)+'\n'
    
    return tags_min,tags_max  
def main(casefile, ntcfile, ntc_name, ana_dir, proteinName_pair,sample, lib, pair_name_lst,temp_file, dwar_file,codesmilesfile):
	rdf = merge_pair_score(casefile, ntcfile,ntc_name,ana_dir, proteinName_pair,sample, lib, pair_name_lst)
	rdf.to_csv(temp_file, header=True, index=None, sep="\t")
	txt2dwar_pair(temp_file, dwar_file,1,codesmilesfile,pair_name_lst,ntc_name,sample)

if __name__ == "__main__":
    data_start = time.time()
    print(time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(time.time())))
    parser = argparse.ArgumentParser(description='get deduped count file')
    parser.add_argument('--InputPath', '-i', help='input count file')
    parser.add_argument('--Batch', '-b', help='batchName')
    parser.add_argument('--SampleInfo', '-s', help='SampleInfo dir')
    parser.add_argument('--LibInfo', '-l', help='LibInfo dir')
    parser.add_argument('--mpi', '-mpi', help='sample mpi numcores')
    #parser.add_argument('--CodeSmiles', '-c', help='codesmiles dir')
    args = parser.parse_args()
    
    codesmiles_dir = sys.path[0] + '/../../lib/Libraries/'
    ana_dir = args.InputPath + '/06.Dwar'
    pair_dir = args.InputPath + '/07.Pair'
    report_dir = args.InputPath + '/09.report/dlp'
    
    input_SampleInfo = args.SampleInfo
    input_LibInfo = args.LibInfo

    protein_dict = process.total_sample_info(input_SampleInfo)[0]  
    library_dict = process.total_sample_info(input_SampleInfo)[2]  
    sample_lst_all = process.total_sample_info(input_SampleInfo)[8]  
    ntc_dict = process.total_sample_info(input_SampleInfo)[10]  
    pair_dict = process.total_sample_info(input_SampleInfo)[12]  
    sample_lst_0 = []
    pair_list = []
    print(pair_dict)
    for key, value in pair_dict.items():
        if value != '0':
            sample_lst_0.extend([key]) 
            pair_list.extend([value])
    sample_lst = list(set(sample_lst_0))
    pair_list = list(set(pair_list))
    print(sample_lst)
    library_size_dict = process.library_info(input_LibInfo)[7]
    sample_lst_0_ntc=[]
    ntc_0_lst = []
    for key, value in ntc_dict.items():
        if value != '0':
            sample_lst_0_ntc.extend([key])
        else:
            ntc_0_lst.extend([key])
    sample_list_ntc = [i for i in sample_lst_0_ntc if i not in pair_list]
    if len(sample_list_ntc)>0:
        sample_lst.extend(sample_list_ntc)
    case_no_ntc=[]
    for i in ntc_0_lst:
        if 'NTC' not in protein_dict[i]:
            case_no_ntc.extend([i])
    case_no_ntc=[i for i in case_no_ntc if i not in pair_list+sample_list_ntc]
    if len(case_no_ntc)>0:
        sample_lst.extend(case_no_ntc)
    sample_lst = list(set(sample_lst))
    
    if len(sample_lst) != 0:
        # 07.pair
        for sample in sample_lst:
            proteinName = protein_dict[sample]
            if 'NTC' in proteinName:
                print('NTC sample!!!')
            else:
                ntc_sample = str(ntc_dict[sample])
                ntc_name = '0'
                print(ntc_name)
                if ntc_sample != '0':
                    ntc_name = protein_dict[ntc_sample]
                lib_lst = library_dict[sample].strip().split(",")
                # paired sample,such as sample adding known inhibitors
                pair_name = pair_dict[sample]
                # if no ntc sample, pair_name_lst=['0']
                pair_name_lst = list(set(pair_name.strip().split(',')))
                if pair_name_lst[0] in protein_dict.keys():
                    if pair_name_lst[0] != '0':
                        proteinName_pair=''
                        try:
                            proteinName_pair = protein_dict[pair_name_lst[0]]
                        except:
                            print('No paired sample!')
                        if proteinName == proteinName_pair:
                            print('Paired: {0} is paired with {1}!'.format(proteinName,proteinName_pair))
                        else:
                            print('Warnning: {0} is not paired with {1}!'.format(proteinName,proteinName_pair))

                for lib in lib_lst:
                    print(sample+'.'+lib)
                    codesmilesfile = codesmiles_dir + '/' + lib + '.code.smiles.txt'
                    ntcfile_filter = '{}/{}/{}/{}.{}.count.freq.filter.txt'.format(ana_dir,ntc_name,lib,ntc_sample,lib)
                    casefile_filter = '{}/{}/{}/{}.{}.count.freq.filter.txt'.format(ana_dir,proteinName,lib,sample,lib)
                    temp_file = '{}/{}/{}/{}.{}.pair.txt'.format(pair_dir,proteinName,lib,sample,lib)
                    dwar_file = '{}/{}/{}/{}.{}.pair.ef.dwar'.format(pair_dir,proteinName,lib,sample,lib)
                    ef_file_dir = '{}/{}/{}'.format(report_dir,proteinName,lib)

                    if os.path.exists(casefile_filter) and os.path.exists(ntcfile_filter):
                        main(casefile_filter, ntcfile_filter, ntc_sample,ana_dir, proteinName,sample, lib, pair_name_lst, temp_file,dwar_file,codesmilesfile)
                        if not os.path.exists(dwar_file):
                            print('Pair dwar {} is not exist!!!'.format(dwar_file))
                        else:
                            process.mkdir(ef_file_dir)
                            shutil.copy(dwar_file, ef_file_dir)
                    else:
                        print('The sample {} or the ntc smaple {} is not exists.'.format(proteinName+'.'+lib,ntc_name+'.'+lib))
        

