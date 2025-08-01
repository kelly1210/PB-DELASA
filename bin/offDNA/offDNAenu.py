# -*- coding: utf-8 -*-
import logging
import pandas as pd
import numpy as np
import time
import argparse
import os
import sys
import subprocess
import multiprocessing as mp
from multiprocessing import Pool
from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem
from rdkit.Chem import QED,Descriptors, rdMolDescriptors
sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), '../monitor'))
import process
import clusfps_deal
import mce18

import warnings
warnings.filterwarnings("ignore")

VinaProg = os.environ['VinaProg']

def func(i):
    logging.info(f'  第{i}个任务')
    return i  
def call_back(res):  
    logging.info(f'\t\t\tfunc函数的返回值--->{res}')

def offDNAmain(dlpfile,offDNAtempfile,sample_dir,dotNumber):
    # dotNumber
    print('Reading dlp detailed file......')
    df_dlp = pd.read_csv(dlpfile,sep='\t',dtype='object')
    df_dlp = df_dlp.applymap((lambda x: "".join(x.split()) if type(x) is str else x)) 
    try:
        df_dlp.replace("0.000","0",inplace=True)
    except:
        print('The Format of dlp detailed file is correct!')
    
    df_temp = df_dlp[df_dlp['FP/TP']=='TP']
    df_dlp2 = df_temp[df_temp['Confirmation(1/0)']=='1']
    print('Positive Hits number: {}.'.format(len(df_dlp2)))
    # positive sample
    fp_sample = list(set(df_dlp2['SampleID'].tolist()))
    
    offDNAsum_f_lst = ['BatchName','Round','Protein','Protein_CC','SampleID','LibraryID','LibrarySize','Tags',
                'TAG1','TAG1_structure','TAG1_FF','TAG1_BB','TAG2','TAG2_structure','TAG2_FF','TAG2_BB',
                'TAG3','TAG3_structure','TAG3_FF','TAG3_BB','Compound','DLP','Plane','Line','SARScore','Dedup_count','z-score','Normalized_z-score','EF']
    dot_df = pd.DataFrame(columns=offDNAsum_f_lst)

    print('Positive sample number: {}.'.format(len(fp_sample)))
    for eachsample in fp_sample:
        df_dlp3 = df_dlp2[df_dlp2['SampleID'] == eachsample]
        BatchName = df_dlp3['BatchName'].tolist()[0]
        Round = df_dlp3['Round'].tolist()[0]
        Protein = df_dlp3['Protein'].tolist()[0]
        Protein_CC = df_dlp3['Protein_CC'].tolist()[0]
        fp_lib = list(set(df_dlp3['LibraryID'].tolist()))
        print('Positive lib number for the sample {0}: {1}'.format(eachsample,len(fp_lib)))
        for eachlib in fp_lib:
            print(eachsample,eachlib)
            df_dlp4 = df_dlp3[df_dlp3['LibraryID'] == eachlib]
            LibrarySize = df_dlp4['LibrarySize'].tolist()[0]
            
            FilterFreqCountFile = "{}/06.Dwar/{}/{}/{}.count.freq.filter.txt".format(sample_dir,Protein,eachlib, eachsample + '.' + eachlib)
            df_count = pd.read_csv(FilterFreqCountFile,sep='\t',dtype='object')
            df_count['Dedup_count'] = df_count['Dedup_count'].astype(int)
            if 'EF' in df_count.columns:
                df_count['EF'] = df_count['EF'].astype(float)
            else:
                df_count['EF'] = 1
            df_dlp_dot = df_dlp4[df_dlp4['DLP']=='D']
            df_dlp_line = df_dlp4[df_dlp4['DLP']=='L']
            df_dlp_plane = df_dlp4[df_dlp4['DLP']=='P']
            print([len(df_dlp_dot),len(df_dlp_line),len(df_dlp_plane)])
            linetag = []
            planetag = []
            
            planetag = [i.strip() for i in list(set(df_dlp_plane['TAG1'].tolist() + df_dlp_plane['TAG2'].tolist() + df_dlp_plane['TAG3'].tolist()))]
            if len(planetag) > 0 and '0' in planetag:
                planetag.remove('0')
            planetag = list(set(planetag))
            if len(df_dlp_line) != 0:
        
                for index,row in df_dlp_line.iterrows():
                    if row['TAG1'].strip() == '0' and row['TAG2'].strip()!='0' and row['TAG3'].strip()!='0':
                        linetag.append(','.join([row['TAG2'].strip(),row['TAG3'].strip()]))
                    if row['TAG2'].strip() == '0' and row['TAG1'].strip()!='0' and row['TAG3'].strip()!='0':
                        linetag.append(','.join([row['TAG1'].strip(),row['TAG3'].strip()]))
                    if row['TAG3'].strip() == '0' and row['TAG1'].strip()!='0' and row['TAG2'].strip()!='0':
                        linetag.append(','.join([row['TAG1'].strip(),row['TAG2'].strip()]))
                # For lines, list the compounds where lines intersect and label them with LL, and finally list the compounds without intersection
                linetag = list(set(linetag))      
                if len(linetag) > 0 and '0' in linetag:
                    linetag.remove('0')
                if len(linetag) > 0 and np.nan in linetag:
                    linetag.remove(np.nan)
                dot_df = linedlp(dot_df,df_count,df_dlp_line,dotNumber,planetag,linetag)
            if len(df_dlp_plane) != 0:
                dot_df = planedlp(dot_df,df_count,df_dlp_plane,dotNumber)
                
            if len(df_dlp_dot) != 0:
                dot_df = dotdlp(df_count,dot_df,df_dlp_dot,planetag,linetag)
            dot_df['BatchName'] = BatchName
            dot_df['SampleID'] = dot_df['SampleID'].fillna(eachsample)
            dot_df['LibraryID'] = dot_df['LibraryID'].fillna(eachlib)
            dot_df['Round'] = dot_df['Round'].fillna(Round)
            dot_df['Protein'] = dot_df['Protein'].fillna(Protein)
            dot_df['Protein_CC'] = dot_df['Protein_CC'].fillna(Protein_CC)
            dot_df['LibrarySize'] = dot_df['LibrarySize'].fillna(LibrarySize)

    dot_df = dot_df.reindex(columns=offDNAsum_f_lst, fill_value='0')
    dot_df.fillna('0')
    dot_df.to_csv(offDNAtempfile, sep='\t', index=False)
    print('dataframe lines: ',len(dot_df))
    print(time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(time.time())))

def lp(dot_new,dot_tag1_bb,dot_tag2_bb,dot_tag3_bb,planetag,linetag):
    plane = []
    if len(planetag)>0:
        if dot_tag1_bb!='0':
            if dot_tag1_bb in planetag:
                plane.append(dot_tag1_bb)
        if dot_tag2_bb!='0':
            if dot_tag2_bb in planetag:
                plane.append(dot_tag2_bb)
        if dot_tag3_bb!='0':
            if dot_tag3_bb in planetag:
                plane.append(dot_tag3_bb)
    dot_new['Plane'] = ';'.join(plane)
    
    line = []
    dot_tag12,dot_tag13,dot_tag23=['0','0','0']
    if dot_tag1_bb != '0' and dot_tag2_bb != '0':
        
        dot_tag12 = ','.join([dot_tag1_bb,dot_tag2_bb])
    if dot_tag1_bb != '0' and dot_tag3_bb != '0':
        dot_tag13 = ','.join([dot_tag1_bb,dot_tag3_bb])
    if dot_tag2_bb != '0' and dot_tag3_bb != '0':
        dot_tag23 = ','.join([dot_tag2_bb,dot_tag3_bb])
    #print(dot_tag12,dot_tag13,dot_tag23)
    
    if len(linetag)>0:
        if dot_tag12!='0':
            if dot_tag12 in linetag:
                line.append(dot_tag12)
        if dot_tag13!='0':
            if dot_tag13 in linetag:
                line.append(dot_tag13)
        if dot_tag23!='0':
            if dot_tag23 in linetag:
                line.append(dot_tag23)
    dot_new['Line'] = ';'.join(line)
    #print(linetag,line)
    #print(planetag,plane)
    
    return dot_new,plane,line
    
def dedupdot(dot_df,newsheet_pre,newsheet):
    #newsheet = pd.DataFrame(colums=newsheet_pre.columns)
    for index, row in newsheet_pre.iterrows():
        tag1 = row['TAG1'].strip()
        tag2 = row['TAG2'].strip()
        tag3 = row['TAG3'].strip()
        df_len = len(dot_df[(dot_df['TAG1']==tag1) & (dot_df['TAG2']==tag2) &(dot_df['TAG3']==tag3)])
        if df_len==0:
            newsheet = newsheet.append(row)
        else:
            continue
    return newsheet

def dotdlp(df_count,dot_df,df_dlp_dot,planetag,linetag):
    for index, row in df_dlp_dot.iterrows():
        dot_tag1_bb = str(row['TAG1'].strip())
        dot_tag2_bb = str(row['TAG2'].strip())
        dot_tag3_bb = str(row['TAG3'].strip())
        
        dot_new = df_count[(df_count['TAG1']==dot_tag1_bb) & 
                      (df_count['TAG2']==dot_tag2_bb) &
                      (df_count['TAG3']==dot_tag3_bb)]
        if len(dot_new) > 0:
            dot_new['DLP'] = 'D'
            newsheet = pd.DataFrame(columns=dot_new.columns)
            newsheet = dedupdot(dot_df,dot_new,newsheet)
            dot_df=dot_df.append(newsheet,ignore_index=True)
    return dot_df

def linedlp(dot_df,df_count,df_dlp_line,dotNumber,planetag,linetag):
    print(len(df_dlp_line))
    for index, row in df_dlp_line.iterrows():
        line_tag1_bb = str(row['TAG1'].strip())
        line_tag2_bb = str(row['TAG2'].strip())
        line_tag3_bb = str(row['TAG3'].strip())
        print(dot_df)
        print([line_tag1_bb,line_tag2_bb,line_tag3_bb])
        if line_tag1_bb == '0' and line_tag2_bb != '0' and line_tag3_bb != '0':
            line_new = df_count[(df_count['TAG2']==line_tag2_bb) &
                  (df_count['TAG3']==line_tag3_bb)]
            if len(line_new) > 0:
                line_new['SARScore'] = 0
                line_new['DLP'] = 'DL'
                line_new,plane,line = lp(line_new,line_tag1_bb,line_tag2_bb,line_tag3_bb,planetag,linetag)                
                
                if len(plane)>=1:
                    line_new['DLP'] = 'DLP'
                if len(line)>=2:
                    line_new['DLP'] = 'DLL'
                
                bb1_dict = {}
                for i in line_new['TAG1'].tolist():
                    bb1_dict[i] = line_new[line_new['TAG1']==i]['TAG1_structure'].tolist()[0]
                
                bb1_sim_dict = bbsmatrix(bb1_dict)
                for i in line_new['TAG1'].tolist():
                    line_new.loc[line_new['TAG1']==i,'SARScore'] = bb1_sim_dict[i]
                newsheet_pre = line_new.sort_values(by=['Dedup_count'],ascending=False)
                dot_df=dot_df.append(newsheet_pre[:int(dotNumber)],ignore_index=True)
        if line_tag2_bb == '0' and line_tag1_bb != '0' and line_tag3_bb != '0':
            line_new = df_count[(df_count['TAG1']==line_tag1_bb) &
                  (df_count['TAG3']==line_tag3_bb)]
            if len(line_new) > 0:
                line_new['DLP'] = 'DL'
                line_new,plane,line = lp(line_new,line_tag1_bb,line_tag2_bb,line_tag3_bb,planetag,linetag)
                
                if len(plane)>=1:
                    line_new['DLP'] = 'DLP'
                if len(line)>=2:
                    line_new['DLP'] = 'DLL'
                bb2_dict = {}
                for i in line_new['TAG2'].tolist():
                    bb2_dict[i] = line_new[line_new['TAG2']==i]['TAG2_structure'].tolist()[0]
                bb2_sim_dict = bbsmatrix(bb2_dict)
                for i in line_new['TAG2'].tolist():
                    line_new.loc[line_new['TAG2']==i,'SARScore'] = bb2_sim_dict[i]
                newsheet_pre = line_new.sort_values(by=['Dedup_count'],ascending=False)
                dot_df=dot_df.append(newsheet_pre[:int(dotNumber)],ignore_index=True)
        if line_tag3_bb == '0' and line_tag1_bb != '0' and line_tag2_bb != '0':
            line_new = df_count[(df_count['TAG2']==line_tag2_bb) &
                  (df_count['TAG1']==line_tag1_bb)]
            if len(line_new) > 0:
                line_new['DLP'] = 'DL'
                line_new,plane,line = lp(line_new,line_tag1_bb,line_tag2_bb,line_tag3_bb,planetag,linetag)
                if len(plane)>=1:
                    line_new['DLP'] = 'DLP'
                if len(line)>=2:
                    line_new['DLP'] = 'DLL'
                bb3_dict = {}
                for i in line_new['TAG3'].tolist():
                    bb3_dict[i] = line_new[line_new['TAG3']==i]['TAG3_structure'].tolist()[0]
                bb3_sim_dict = bbsmatrix(bb3_dict)
                for i in line_new['TAG3'].tolist():
                    line_new.loc[line_new['TAG3']==i,'SARScore'] = bb3_sim_dict[i]
                newsheet_pre = line_new.sort_values(by=['Dedup_count'],ascending=False)
                dot_df=dot_df.append(newsheet_pre[:int(dotNumber)],ignore_index=True)
    return dot_df

def planedlp(dot_df,df_count,df_dlp_plane,dotNumber):
    dot_df_new = pd.DataFrame(dot_df)
    for index, row in df_dlp_plane.iterrows():
        plane_tag1_bb = str(row['TAG1'])
        plane_tag2_bb = str(row['TAG2'])
        plane_tag3_bb = str(row['TAG3'])
        if plane_tag1_bb == '0' and plane_tag2_bb == '0' and plane_tag3_bb != '0':
            plane_new = df_count[(df_count['TAG3']==plane_tag3_bb)]
            if len(plane_new) >0:
                try:
                    plane_new.loc[:,'DLP'] = 'DP'
                    plane_new.loc[:,'Plane'] = plane_tag3_bb
                except:
                    print('No TAG3 {}!!!'.format(plane_tag3_bb))
                newsheet_pre = plane_new.sort_values(by=['Dedup_count'],ascending=False)
                dot_df_new=dot_df_new.append(newsheet_pre[:int(dotNumber)],ignore_index=True)
        if plane_tag1_bb == '0' and plane_tag3_bb == '0' and plane_tag2_bb != '0':
            plane_new = df_count[(df_count['TAG2']==plane_tag2_bb)]
            if len(plane_new) >0:
                try:
                    plane_new.loc[:,'DLP'] = 'DP'
                    plane_new.loc[:,'Plane'] = plane_tag2_bb
                except:
                    print('No TAG2 {}!!!'.format(plane_tag2_bb))
                newsheet_pre = plane_new.sort_values(by=['Dedup_count'],ascending=False)
                dot_df_new=dot_df_new.append(newsheet_pre[:int(dotNumber)],ignore_index=True)
        if plane_tag3_bb == '0' and plane_tag2_bb == '0' and plane_tag1_bb != '0':
            plane_new = df_count[(df_count['TAG1']==plane_tag1_bb)]
            if len(plane_new) >0:
                try:
                    plane_new.loc[:,'DLP'] = 'DP'
                    plane_new.loc[:,'Plane'] = plane_tag1_bb
                except:
                    print('No TAG1 {}!!!'.format(plane_tag1_bb))
                    
                newsheet_pre = plane_new.sort_values(by=['Dedup_count'],ascending=False)
                dot_df_new=dot_df_new.append(newsheet_pre[:int(dotNumber)],ignore_index=True)
    return dot_df_new

def Pfizer_MPO_score(mw,logp,hbd,arom):
    if mw < 400:
        mw_score = 1
    elif mw <= 600:
        mw_score = 3-1.0/200*mw
    else:
        mw_score = 0
    if logp < -1:
        logp_score = -1
    elif logp <= 1:
        logp_score = (logp + 1)/2.0
    elif logp < 3:
        logp_score = 1
    elif logp < 5:
        logp_score = (logp - 1)/2.0
    else:
        logp_score = 0  
    if arom <= 2:
        arc_score = 1
    elif arom == 3:
        arc_score = 0.5
    elif arom >= 4:
        arc_score = 0
    if hbd <=2:
        hbd_score = 1
    elif hbd == 3:
        hbd_score = 0.5
    else:
        hbd_score = 0
    mpo_score = mw_score * (logp_score + arc_score + hbd_score)
    return mpo_score

def bbsmatrix(bbdict):
    for index in bbdict.keys():
        bbmol = Chem.MolFromSmiles(bbdict[index])
        try:
            bbfp = AllChem.GetMorganFingerprint(bbmol,2)
            bbdict[index] = bbfp
        except:
            bbdict[index] = '0'
    bbsize = len(bbdict)
    table_sim = {}
    for iindex in bbdict.keys():
        ifp = bbdict[iindex]
        similarity = 0
        for jindex in bbdict.keys():
            jfp = bbdict[jindex]
            try:
                similarity += DataStructs.DiceSimilarity(ifp,jfp)
            except:
                continue
        similarity = float('%.4f'%(similarity / bbsize))
        table_sim[iindex] = similarity
    sorted(table_sim.items(),key=lambda i:i[1],reverse=True)
    return table_sim
        
def FP_morgan(mol1,mol2):
    fp1 = AllChem.GetMorganFingerprint(mol1,4)
    fp2 = AllChem.GetMorganFingerprint(mol2,4)
    sim = DataStructs.TanimotoSimilarity(fp1,fp2)
    novelty = round(1-sim,2)
    return novelty
def target_compound(targetfile):
    df1 = pd.read_csv(targetfile,sep='\t',header=0,dtype=object)
    df2 = df1[df1['Smiles']!='0']
    df2['Smiles'] = df2['Smiles'].apply(lambda x: x.strip().split('.')[0])
    bb_dict = df2.set_index('Entry Number')['Smiles']
    target_dict = df2.set_index('Entry Number')['Target']
    return bb_dict,target_dict

class enumerationMPI():
    def __init__(self):
        self.manager = mp.Manager
        self.typedict = self.manager().dict()
        self.compdict = self.manager().dict()
        self.mwdict = self.manager().dict()
        self.logpdict = self.manager().dict()
        self.hbddict = self.manager().dict()
        self.hbadict = self.manager().dict()
        self.psadict = self.manager().dict()
        self.qeddict = self.manager().dict()
        self.fsp3dict = self.manager().dict()
        self.mce18dict = self.manager().dict()
        self.vinadict = self.manager().dict()
        self.mpodict = self.manager().dict()
    def enumeration(self,ffdict,vina_lst,mpinum):
        mpinum = int(mpinum)
        p=mp.Pool(mpinum)
        for ff in ffdict.keys():
            p.apply_async(self.enumerationMap,args=(ff,ffdict,vina_lst))
        p.close()
        p.join()  
    def enumerationMap(self,ff,ffdict,vina_lst):
        temp_smiles = ffdict[ff]
        ffs = temp_smiles.strip().split('.')
        join_smiles = '0'
        self.compdict[ff] = '0'
        self.typedict[ff] = '0'
        if len(ffs) ==3 and ffs[0]!='0' and ffs[1]!='0' and ffs[2]!='0':
            self.typedict[ff] = '123'
            try:
                join_smiles = process.join_r_groups(temp_smiles,Starsite='no')
                self.compdict[ff] = join_smiles
            except:
                print('Error: Cannot enumeration!!!')
        if len(ffs) ==3 and ffs[0]!='0' and ffs[1]!='0' and ffs[2]=='0':
            self.typedict[ff] = '12'
            try:
                join_smiles = process.join_r_groups('.'.join(ffs[0:2]),Starsite='no')
                self.compdict[ff] = join_smiles
            except:
                print('Error: Cannot enumeration!!!')
        if len(ffs) ==3 and ffs[0]!='0' and ffs[1]=='0':
            self.typedict[ff] = '1'
            if ffs[2] != '0':
                self.typedict[ff] = '13'
            try:
                join_smiles = process.join_r_groups('.'.join(ffs[0:1]),Starsite='no')
                self.compdict[ff] = join_smiles
            except:
                print('Error: Cannot enumeration!!!')
        if len(ffs) ==3 and ffs[0]=='0' and ffs[1]!='0':
            self.typedict[ff] = '2'
            if ffs[2] != '0':
                self.typedict[ff] = '23'
        if len(ffs) ==3 and ffs[0]=='0' and ffs[1]=='0' and ffs[2]!='0':
            self.typedict[ff] = '3'
        pdbpath = vina_lst[0];proteinfile = vina_lst[1];proteinName = vina_lst[2];outputdocking = vina_lst[3]
        mw,logp,hbd,hba,psa,qed,fsp3,mce18score,vinadockingscore,mposcore = self.properties(ff,join_smiles,pdbpath,proteinfile,proteinName,outputdocking)
        self.mwdict[ff] = mw
        self.logpdict[ff] = logp
        self.hbddict[ff] = hbd
        self.hbadict[ff] = hba
        self.psadict[ff] = psa
        self.qeddict[ff] = qed
        self.fsp3dict[ff] = fsp3
        self.mce18dict[ff] = mce18score
        self.vinadict[ff] = vinadockingscore
        self.mpodict[ff] = mposcore
    def properties(self,ff,join_smiles,pdbpath,proteinfile,proteinName,outputdocking):
        try:
            compoundmol = Chem.MolFromSmiles(join_smiles)
        except:
            print(join_smiles)
        mw = 0;logp = 0;hbd = 0; arom = 0
        psa=0;qed=0;fsp3=0;mce18score=0;vinadockingscore=0;mposcore=0
        try:
            mw =Descriptors.MolWt(compoundmol)
            mw = round(mw,2)
        except:
            print('The enumerated compound {} error in cacluating mw!'.format(ff))
        try:
            logp =Descriptors.MolLogP(compoundmol)
            logp = round(logp,2)
            #logpdict[ff] = logp
        except:
            print('The enumerated compound {} error in cacluating logp!'.format(ff))
        try:
            hbd = rdMolDescriptors.CalcNumLipinskiHBD(compoundmol)
            hbd = round(hbd,2)
            #hbddict[ff] = hbd
        except:
            print('The enumerated compound {} error in cacluating hbd!'.format(ff))
        try:
            hba = rdMolDescriptors.CalcNumLipinskiHBA(compoundmol)
            hba = round(hba,2)
            #hbadict[ff] = hba
        except:
            print('The enumerated compound {} error in cacluating hba!'.format(ff))
        try:
            psa =Descriptors.TPSA(compoundmol)
            psa = round(psa,2)
            #psadict[ff] = psa
        except:
            print('The enumerated compound {} error in cacluating psa!'.format(ff))
        try:
            qed = QED.qed(compoundmol)
            qed = str(round(qed,2))
        except:
            print('The enumerated compound {} error in cacluating qed!'.format(ff))
        try:
            fsp3 = Descriptors.FractionCSP3(compoundmol)
            fsp3 = round(fsp3,2)
        except:
            print('The enumerated compound {} error in cacluating fsp3!'.format(ff))
        try:
            mce18score = mce18.MCE18(compoundmol).CalculateMCE18()
            mce18score = str(round(mce18score,2))
        except:
            print('The enumerated compound {} error in cacluating mce18!'.format(ff))
        pocketNum=1
        try:
            print('Starting vina docking: {}!!!'.format(ff))
            vinadockingscore = vinadocking(pdbpath,proteinfile,proteinName,compoundmol,pocketNum,ff,outputdocking)
            vinadockingscore = (round(vinadockingscore,2))
        except:
            print('The enumerated compound {} error in cacluating docking score!'.format(ff))
        if arom == 'error':
            mposcore = 'Error'
        else:
            try:
                mposcore = round(Pfizer_MPO_score(mw,logp,hbd,arom),2)
            except:
                mposcore = 'Error'
        return mw,logp,hbd,hba,psa,qed,fsp3,mce18score,vinadockingscore,mposcore
class noveltySimilarity():
    def __init__(self):
        self.manager = mp.Manager
        self.novelty_dict = self.manager().dict()
    def noveltySim(self,keys,bb_dict,inputmol):
        mol2 = Chem.MolFromSmiles(bb_dict[keys])
        self.novelty_dict[keys] = FP_morgan(inputmol,mol2)
    def flow(self,cpu,bb_dict,inputmol):
        p=Pool(int(cpu))
        for keys in bb_dict.keys():
             p.apply_async(self.noveltySim,args = (keys,bb_dict,inputmol))
        p.close()
        p.join()
        return self.novelty_dict
        
def zscorePocket(zscore1,zscore2):
    if float(zscore1) >= float(zscore2):
        return 'Inner'
    else:
        return 'Unknown'
def mechanism_identify(args,args_diff):
    mechanism_lst = []
    if len(args) == 1:
        return 'Inner'
    if len(args) == 2:
        if float(args[0]) >= float(args[1]):
            mechanism_lst.extend('A')
    if len(args) == 3:
        if (float(args[0]) >= float(args[1])) or (float(args[0]) >= float(args[2])):
            mechanism_lst.extend('A')
    if len(args_diff) == 1:
        if (float(args[0]) >= float(args_diff[0])):
            mechanism_lst.extend('B')
    if len(args_diff) == 2:
        if (float(args[0]) >= float(args_diff[0])) or (float(args[0]) >= float(args_diff[1])):
            mechanism_lst.extend('B')
    if len(mechanism_lst) == 0:
        return 'Other'
    else:
        return ','.join(mechanism_lst)
    
def dict_delzero(dict1):
    for k in list(dict1.keys()):
        if dict1[k] == '0':
            del dict1[k]
    return dict1
def dlp_fun(x):
    '''
    DLL: the intersection point on two lines
    DLP：the compound on the line and plane
    DL: the compound on the line
    DP: the compound on the plane
    D: the compound with bigger counts
    '''
    if x=='DLL':
        return 5
    if x=='DLP':
        return 4
    if x=='DL':
        return 3
    if x=='DP':
        return 2
    if x=='D':
        return 1

def eachcluster(sample_lst,dlp_str,pocket_sample,pocket_diff_sample,df2,samples_alllst,case_sample,ntcsample,sample_dir,ntc_protein,protein,cpu,outputdocking):
    lib_path = sys.path[0] + '/../../lib/'
    proteinfile = lib_path + 'proteinlist.txt'
    targetDir = lib_path + 'inhibitors/'
    pdbpath = lib_path + 'pdb/'
    if len(sample_lst)!=0:
        for i in range(0,len(sample_lst)):
            dlp_str.extend(['DLP('+sample_lst[i]+')','Plane('+sample_lst[i]+')','Line('+sample_lst[i]+')',
      'Dedup_count('+sample_lst[i]+')','z-score('+sample_lst[i]+')','EF('+sample_lst[i]+')','Repeat('+sample_lst[i]+')'])
    dlp_str.extend(['ClusterID'])
    if pocket_sample != '0':
        pocket_lst = pocket_sample.strip().split(',')
        for m in pocket_lst:
            dlp_str.extend(['Pocket(' + m +')'])
    pocket_diff_lst = []
    if pocket_diff_sample != '0':
        pocket_diff_lst = pocket_diff_sample.strip().split(',')
        for m in pocket_diff_lst:
            dlp_str.extend(['Pocket(' + m +')'])
    dlp_str.extend(['Mechanism'])
    x,y = df2.shape  
    df3 = df2.drop_duplicates(subset=['uniqtags'], keep='first').reset_index(drop=True)
    prop_str = ['Type','MW','logp','HBD','HBA','TPSA','QED','MCE18','Fsp3','VinaScore','MPO_score']
    novelty_str = ['Novelty','Novelty_smiles','Entry Number']
    #print(df3.head())
    new_column = ['BatchName','Round','Protein','Protein_CC','SampleID',
  'LibraryID','LibrarySize','Tags','TAG1','TAG1_structure',
  'TAG1_FF','TAG1_BB','TAG2','TAG2_structure','TAG2_FF',
  'TAG2_BB','TAG3','TAG3_structure','TAG3_FF','TAG3_BB','Compound','DLP','Plane','Line',
  'Dedup_count','z-score','EF','uniqtags','uniqstructure','uniqid','BBs']+dlp_str+prop_str+novelty_str
    df4 = pd.concat([df3,pd.DataFrame(np.zeros((x,len(dlp_str+prop_str+novelty_str))),columns=dlp_str+prop_str+novelty_str)],axis=1)
    df4['DLP'] = df4['DLP'].apply(lambda x: dlp_fun(x))
    df40 = df4[new_column]
    df40['Compound'] = '0'
    df40.dropna(subset=['uniqtags'],axis=0,inplace=True)
    ntc_protein_dir = '{}/06.Dwar/{}/'.format(sample_dir,ntc_protein)
    ntc_dedup_count_dict = {}
    ntc_zscore_dict = {}
    if os.path.exists(ntc_protein_dir):
        libs = os.listdir(ntc_protein_dir)
        for eachlib in libs:
            FilterFreqCountFile = "{}/06.Dwar/{}/{}/{}.count.freq.filter.txt".format(sample_dir,ntc_protein,eachlib, ntcsample + '.' + eachlib)
            df_ntc = pd.read_csv(FilterFreqCountFile,sep='\t',dtype='object')
            df_ntc['LibraryID'] = eachlib
            df_ntc['uniqstructure'] = df_ntc['LibraryID'] + '_' + df_ntc['TAG1_structure']+ '_'+ df_ntc['TAG2_structure']+ '_' + df_ntc['TAG3_structure']
            df_ntc = df_ntc.drop_duplicates(subset=['uniqstructure'], keep='first')
            ntc_dedup_count_dict = df_ntc[['uniqstructure', 'Dedup_count']].set_index("uniqstructure").to_dict(orient='dict')['Dedup_count']
            ntc_zscore_dict = df_ntc[['uniqstructure', 'z-score']].set_index("uniqstructure").to_dict(orient='dict')['z-score']
        
    for j in samples_alllst:
        if j in set(df40['SampleID'].tolist()):
            df41 = df40[df40['SampleID']==j]
            bbs_dict = df41['uniqstructure'].value_counts().to_dict()
            df41 = df41.drop_duplicates(subset=['uniqstructure'], keep='first')
            dlp_dict = df41[['uniqstructure','DLP']].set_index("uniqstructure").to_dict(orient='dict')['DLP']
            plane_dict = df41[['uniqstructure', 'Plane']].set_index("uniqstructure").to_dict(orient='dict')['Plane']
            line_dict = df41[['uniqstructure', 'Line']].set_index("uniqstructure").to_dict(orient='dict')['Line']
            dedup_dict = df41[['uniqstructure', 'Dedup_count']].set_index("uniqstructure").to_dict(orient='dict')['Dedup_count']
            zscore_dict = df41[['uniqstructure', 'z-score']].set_index("uniqstructure").to_dict(orient='dict')['z-score']
            ef_dict = df41[['uniqstructure', 'EF']].set_index("uniqstructure").to_dict(orient='dict')['EF']
            
            for i in bbs_dict.keys():
                df40.loc[df40['uniqstructure']==i,'Repeat(' + j +')'] = bbs_dict[i]
                df40.loc[df40['uniqstructure']==i,'DLP(' + j +')'] = dlp_dict[i]
                df40.loc[df40['uniqstructure']==i,'Plane(' + j +')'] = plane_dict[i]
                df40.loc[df40['uniqstructure']==i,'Line(' + j +')'] = line_dict[i]
                df40.loc[df40['uniqstructure']==i,'Dedup_count(' + j +')'] = dedup_dict[i]
                df40.loc[df40['uniqstructure']==i,'z-score(' + j +')'] = zscore_dict[i]
                df40.loc[df40['uniqstructure']==i,'EF(' + j +')'] = ef_dict[i]
                if len(ntc_dedup_count_dict)!=0:
                    if i in ntc_dedup_count_dict.keys():
                        df40.loc[df40['uniqstructure']==i,'Dedup_count(' + ntcsample +')'] = ntc_dedup_count_dict[i]
                if len(ntc_zscore_dict)!=0:
                    if i in ntc_zscore_dict.keys():
                        df40.loc[df40['uniqstructure']==i,'z-score(' + ntcsample +')'] = ntc_zscore_dict[i]
    df40 = df40.drop_duplicates(subset=['uniqstructure'], keep='first').reset_index(drop=True)
    if pocket_sample != '0':
        pocket_lst = pocket_sample.strip().split(',')
        for m in pocket_lst:
            df40['Pocket(' + m +')'] = 'unKnown'
            df40['Pocket(' + m +')'] = df40[['z-score(' + case_sample +')','z-score(' + m +')']].apply(lambda x: zscorePocket(x['z-score(' + case_sample +')'],x['z-score(' + m +')']),axis=1)
        z_score_col = ['z-score({})'.format(case_sample)] + ['z-score({})'.format(i) for i in pocket_lst]
        pocket_diff = ['z-score({})'.format(i) for i in pocket_diff_lst]
        if (len(z_score_col) == 2) and (len(pocket_diff)==0):
            df40['Mechanism'] = df40[z_score_col].apply(lambda x: mechanism_identify([x[z_score_col[0]],x[z_score_col[1]]],[]),axis=1)
        if len(z_score_col) == 3 and (len(pocket_diff)==0):
            df40['Mechanism'] = df40[z_score_col].apply(lambda x: mechanism_identify([x[z_score_col[0]],x[z_score_col[1]],x[z_score_col[2]]],[]),axis=1)
        if len(z_score_col) == 3 and (len(pocket_diff)==1):
            df40['Mechanism'] = df40[z_score_col+pocket_diff].apply(lambda x: mechanism_identify([x[z_score_col[0]],x[z_score_col[1]],x[z_score_col[2]]],[x[pocket_diff[0]]]),axis=1)
        if len(z_score_col) == 2 and (len(pocket_diff)==1):
            df40['Mechanism'] = df40[z_score_col+pocket_diff].apply(lambda x: mechanism_identify([x[z_score_col[0]],x[z_score_col[1]]],[x[pocket_diff[0]]]),axis=1)
        if len(z_score_col) == 2 and (len(pocket_diff)==2):
            df40['Mechanism'] = df40[z_score_col+pocket_diff].apply(lambda x: mechanism_identify([x[z_score_col[0]],x[z_score_col[1]]],[x[pocket_diff[0]],x[pocket_diff[1]]]),axis=1)
        if len(z_score_col) == 3 and (len(pocket_diff)==2):
            df40['Mechanism'] = df40[z_score_col+pocket_diff].apply(lambda x: mechanism_identify([x[z_score_col[0]],x[z_score_col[1]],x[z_score_col[2]]],[x[pocket_diff[0]],x[pocket_diff[1]]]),axis=1)
    df40['uniqid_ffs'] = df40['TAG1_FF'].apply(lambda x: str(x).strip().split('.')[0]) + '.' + df40['TAG2_FF'].apply(lambda x: str(x).strip().split('.')[0]) + '.' + df40['TAG3_FF'].apply(lambda x: str(x).strip().split('.')[0])
    ffs_dict = df40[['BBs', 'uniqid_ffs']].set_index("BBs").to_dict(orient='dict')['uniqid_ffs']
    type_dict = {}
    comp_dict = {};mwdict = {};logpdict = {};hbddict = {};hbadict = {};psadict={}
    qeddict = {};fsp3dict = {};mce18dict = {};vinadict = {};mpodict = {}
    vina_lst = [pdbpath,proteinfile,protein,outputdocking + '/' + protein + '/']
    # enumeration and physichemical properties calculation
    enu = enumerationMPI()
    enu.enumeration(ffs_dict,vina_lst, cpu)
    type_dict = enu.typedict
    comp_dict = enu.compdict
    mwdict = enu.mwdict
    logpdict = enu.logpdict; hbddict = enu.hbddict; hbadict = enu.hbadict
    psadict = enu.psadict; qeddict = enu.qeddict; fsp3dict = enu.fsp3dict; mce18dict = enu.mce18dict
    vinadict = enu.vinadict; mpodict = enu.mpodict
    comp_dict2 = dict_delzero(comp_dict)
    if os.path.exists(targetDir + '/' + protein + '.txt'):
        bb_dict,target_dict = target_compound(targetDir + '/' + protein + '.txt')
        for m in comp_dict2.keys():
            inputsmiles = comp_dict2[m]
            inputmol = Chem.MolFromSmiles(inputsmiles)
            novelty_dict = {}
            novel = noveltySimilarity()
            novel.flow(cpu,bb_dict,inputmol)
            novelty_dict = novel.novelty_dict
            # novelty
            min_key = min(novelty_dict, key=lambda x: novelty_dict[x])
            df40.loc[df40['BBs']==m,'Entry Number'] = min_key
            df40.loc[df40['BBs']==m,'Novelty'] = novelty_dict[min_key]
            df40.loc[df40['BBs']==m,'Novelty_smiles'] = bb_dict[min_key]
    for m in comp_dict.keys():
        inputsmiles = comp_dict[m]
        df40.loc[df40['BBs']==m,'Type'] = type_dict[m]
        df40.loc[df40['BBs']==m,'Compound'] = inputsmiles
        df40.loc[df40['BBs']==m,'MW'] = round(float(mwdict[m]),2)
        df40.loc[df40['BBs']==m,'logp'] = round(float(logpdict[m]),2)
        df40.loc[df40['BBs']==m,'HBD'] = round(float(hbddict[m]),2)
        df40.loc[df40['BBs']==m,'HBA'] = round(float(hbadict[m]),2)
        df40.loc[df40['BBs']==m,'TPSA'] = round(float(psadict[m]),2)
        df40.loc[df40['BBs']==m,'QED'] = round(float(qeddict[m]),2)
        df40.loc[df40['BBs']==m,'MCE18'] = round(float(mce18dict[m]),2)
        df40.loc[df40['BBs']==m,'Fsp3'] = round(float(fsp3dict[m]),2)
        df40.loc[df40['BBs']==m,'VinaScore'] = round(float(vinadict[m]),2)
        df40.loc[df40['BBs']==m,'MPO_score'] = round(float(mpodict[m]),2)
    df41 = df40.drop(columns = ['ClusterID']).reset_index(drop=True)
    df40_0 = df41[df41['Compound']=='0'].reset_index()
    df40_1 = df41[df41['Compound']!='0'].reset_index()
    temp0 = clusfps_deal.clusterDF(df40_1,'Compound','moecfp','Butina',cutoff=0.7,nclusters=1,radius=2,clusterIDname='ClusterID')

    df40 = df40_0.append(temp0,ignore_index=True)

    df40 = df40[new_column]
    df40 = df40[df40['Type'].isin(['123','12','1'])]
    df40['TAG1'] = df40['TAG1'].astype(str)
    df40['TAG2'] = df40['TAG2'].astype(str)
    df40['TAG3'] = df40['TAG3'].astype(str)
    df40 = df40.astype({'ClusterID':int,'Dedup_count('+case_sample+')':int,'MPO_score':float,'Repeat('+case_sample+')':int})
    mechanisms = list(set(df40['Mechanism'].tolist()))
    if len(mechanisms) >0:
        for eachmechanism in mechanisms:
            temp_df = df40[df40['Mechanism']==eachmechanism]
            for eachcutcolumn in ['ClusterID']:
                cluster_lst = list(set(temp_df[eachcutcolumn].tolist()))
                cluster_lst.sort()
                cluster_dict = dict(zip(cluster_lst,list(range(1,len(cluster_lst)+1))))
                for cluster in cluster_lst:
                    df40.loc[(df40['Mechanism']==eachmechanism)&(df40[eachcutcolumn]==cluster),eachcutcolumn] = cluster_dict[cluster]
    return df40
    
def clustering(offDNAfile,outfilepath,commentfile,sample_dir,cpu,sampleinfo,outputdocking):
    # duplicated based on BBs and marked as Repeat
    # sample_dir: Analysis,06.dwar
    df1 = pd.read_csv(offDNAfile,sep='\t',header=0,dtype=object)
    # statistic analysis and compound clustring
    protein_dict = process.total_sample_info(sampleinfo)[0]
    ntc_dict = process.total_sample_info(sampleinfo)[8] 
    pair_dict = process.total_sample_info(sampleinfo)[10]
    pocket_dict = process.total_sample_info(sampleinfo)[11]
    pocket_diff_dict = process.total_sample_info(sampleinfo)[12]
    protein_sample_dict = process.kv_reversal(protein_dict)
    sample_lst_0=[]
    for key, value in pair_dict.items():
        if value != '0':
            sample_lst_0.extend([key])  
    sample_lst_0 = list(set(sample_lst_0))
    for protein in protein_sample_dict.keys():
        if 'NTC' not in protein:
            outfile = outfilepath + '/' + protein + '_offDNA.txt'
            outfile_excel = outfilepath + '/' + protein + '_offDNA.xlsx'
            outfile_error = outfilepath + '/' + protein + '_offDNA_error.txt'
            outfile_excel_error = outfilepath + '/' + protein + '_offDNA_error.xlsx'
            samples_alllst = protein_sample_dict[protein]
            case_sample_lst = list(set(sample_lst_0)&set(samples_alllst))
            if protein in set(df1['Protein'].tolist()):
                df2 = df1[df1['Protein']==protein]
                df2['uniqid'] = df2['SampleID'] + '_' + df2['LibraryID'] + '_' + df2['TAG1_structure']+ '_'+ df2['TAG2_structure']+ '_' + df2['TAG3_structure']
                df2['BBs'] = df2['LibraryID'] + '_' + df2['TAG1_BB']+ '_'+ df2['TAG2_BB']+ '_' + df2['TAG3_BB']
                df2['uniqtags'] = df2['SampleID'] + '_' + df2['Tags']
                df2['uniqstructure'] = df2['LibraryID'] + '_' + df2['TAG1_structure']+ '_'+ df2['TAG2_structure']+ '_' + df2['TAG3_structure']
                print(set(df2['SampleID'].tolist()))
                # Duplicating based on structures in the library
                # clustering
                if len(case_sample_lst) == 0:
                    print("Warning: Just one condition of NTC!!!")
                    sample_lst = samples_alllst
                    dlp_str = []
                    pocket_sample = '0'
                    pocket_diff_sample='0'
                    case_sample = sample_lst[0]
                    ntcsample = ntc_dict[case_sample]
                    if ntcsample != '0':
                        ntc_protein = protein_dict[ntcsample]
                        dlp_str.extend(['Dedup_count('+ntcsample+')','z-score('+ntcsample+')'])
                    else:
                        ntc_protein = '0'
                    df40 = eachcluster(sample_lst,dlp_str,pocket_sample,pocket_diff_sample,df2,samples_alllst,case_sample,ntcsample,sample_dir,ntc_protein,protein,cpu,outputdocking)
                    df40 = df40.drop(columns=['SampleID','uniqid','uniqtags','uniqstructure','BBs','DLP','Plane','Line','Dedup_count','z-score','EF'])
                    df40 = df40.sort_values(by=['EF('+case_sample+')', 'Dedup_count('+case_sample+')', 'MPO_score', 'ClusterID','LibraryID'],ascending=[False,False,False,True,True])
                    df40_0 = df40[df40['Compound']=='0']
                    df40_1 = df40[df40['Compound']!='0']
                    df40_1.to_csv(outfile,sep='\t',index=False)
                    severalCSV2xlsx(outfile,commentfile,outfile_excel)
                    df40_0.to_csv(outfile_error,sep='\t',index=False)
                    severalCSV2xlsx(outfile_error,commentfile,outfile_excel_error)
                elif len(case_sample_lst) == 1:
                    case_sample = case_sample_lst[0]
                    ntcsample = ntc_dict[case_sample]
                    print('ntcsample: {}'.format(ntcsample))
                    if ntcsample != '0':
                        ntc_protein = protein_dict[ntcsample]
                    else:
                        ntc_protein = '0'
                    print('ntc_protein: {}'.format(ntc_protein))
                    pocket_sample = pocket_dict[case_sample]
                    pocket_diff_sample = pocket_diff_dict[case_sample]
                    sample_lst = [i for i in samples_alllst]
                    sample_lst.remove(case_sample)
                    dlp_str = []
                    if ntcsample != '0':
                        dlp_str.extend(['Dedup_count('+ntcsample+')','z-score('+ntcsample+')'])
                    dlp_str.extend(['DLP('+case_sample+')','Plane('+case_sample+')','Line('+case_sample+')',
                      'Dedup_count('+case_sample+')','z-score('+case_sample+')','EF('+case_sample+')','Repeat('+case_sample+')'])
                    print(dlp_str)
                    df40 = eachcluster(sample_lst,dlp_str,pocket_sample,pocket_diff_sample,df2,samples_alllst,case_sample,ntcsample,sample_dir,ntc_protein,protein,cpu,outputdocking)
                    df40 = df40.drop(columns=['SampleID','uniqid','uniqtags','uniqstructure','BBs','DLP','Plane','Line','Dedup_count','z-score','EF'])
                    df40 = df40.sort_values(by=['EF('+case_sample+')', 'Dedup_count('+case_sample+')', 'MPO_score', 'ClusterID','LibraryID'],ascending=[False,False,False,True,True])
                    
                    df40_0 = df40[df40['Compound']=='0']
                    df40_1 = df40[df40['Compound']!='0']
                    df40_1.to_csv(outfile,sep='\t',index=False)
                    severalCSV2xlsx(outfile,commentfile,outfile_excel)
                    df40_0.to_csv(outfile_error,sep='\t',index=False)
                    severalCSV2xlsx(outfile_error,commentfile,outfile_excel_error)
                    
def vinadocking(pdbpath,PROpdbfile,proteinName,drgmol,pocketNum,BBs,outputpath):
    # pdbid is lowercase
    vina_score = 0
    pro_df = pd.read_csv(PROpdbfile,sep='\t',dtype='object',header=0)
    pdbid = pro_df[pro_df['Protein']==proteinName]['PDBID'].tolist()[0]
    if pdbid:
        pdbfile = pdbpath + '/' + pdbid + '.pdb'
        outputpathbb = outputpath + '/' + BBs   
        process.mkdir(outputpath)
        process.mkdir(outputpathbb)
        drgfile = outputpathbb + '/' + BBs + '.pdb'
        if not os.path.exists(drgfile):
            m3d = Chem.AddHs(drgmol)
            AllChem.EmbedMolecule(m3d, randomSeed=3)
            AllChem.MMFFOptimizeMolecule(m3d)
            Chem.MolToPDBFile(m3d,drgfile)
            #try:
            vina_cmd = " ".join(['perl', VinaProg+ '/AutoBlindDock.pl', pdbfile, drgfile, str(pocketNum), outputpathbb])
            print(vina_cmd)
            try:
                subprocess.check_call(vina_cmd, shell=True)
            except:
                print('vina scoring error!!!')
        vinafiles=os.listdir(outputpathbb)
        vinaresult_dir = []
        for file in vinafiles:
            if (".txt" not in file) & (".gz" not in file) & (".pdb" not in file) & (".mol2" not in file):
                vinaresult_dir.append(file)
        eachdrg = outputpathbb + '/' + vinaresult_dir[0]
        configtxt = eachdrg + '/' + 'config.txt'
        try:
            config_df=pd.read_csv(configtxt,sep=" ",header=0)
            vina_score = config_df['score'].tolist()[0]
        except:
            vina_score=0
    else:
        vina_score = 0
    return vina_score
    
def severalCSV2xlsx(input_file,commentsfile,out_file):
    writer = pd.ExcelWriter(out_file)
    df = pd.read_csv(input_file,sep = "\t",header=0,encoding="gbk")
    df.to_excel(writer, sheet_name='Hits')
    df = pd.read_csv(commentsfile,sep = "\t",header=0,encoding="gbk")
    df.to_excel(writer, sheet_name='Comments')
    writer.close()
    
def duplication(input_file,out_file):
    df=pd.read_csv(input_file,header=0,sep='\t',dtype=object)
    df.fillna('0')
    df2=pd.DataFrame(columns=df.columns,dtype=object)
    for index, row in df.iterrows():
        samplename = row['SampleID']
        libname = row['LibraryID']
        tag1 = str(row['TAG1'])
        tag2 = str(row['TAG2'])
        tag3 = str(row['TAG3'])
        df3 = df2[(df2['SampleID']==samplename)&(df2['LibraryID']==libname)&(df2['TAG1']==tag1)&(df2['TAG2']==tag2)&(df2['TAG3']==tag3)]
        if len(df3)==0:
            #print(row)
            df4=pd.DataFrame(columns=df.columns)
            df4.loc[len(df4)] = row.tolist()
            df2 = df2.append(df4)
        if len(df3)==1:
            dlptype = row['DLP']
            planetype = row['Plane']
            linetype = row['Line']
            sarscore = float(row['SARScore'])
            dlpscore = 0
            df3_dlpscore = 0
            df3_dlp = df3['DLP'].tolist()[0]
            df3_plane = df3['Plane'].tolist()[0]
            df3_line = df3['Line'].tolist()[0]
            df3_sarscore = float(df3['SARScore'].tolist()[0])
            if row['DLP'] == 'DLL':
                dlpscore = 5
            if row['DLP'] == 'DLP':
                dlpscore = 4
            if row['DLP'] == 'DL':
                dlpscore = 3
            if row['DLP'] == 'DP':
                dlpscore = 2
            if row['DLP'] == 'D':
                dlpscore = 1
            
            if df3_dlp == 'DLL':
                df3_dlpscore = 5
            if df3_dlp == 'DLP':
                df3_dlpscore = 4
            if df3_dlp == 'DL':
                df3_dlpscore = 3
            if df3_dlp == 'DP':
                df3_dlpscore = 2
            if df3_dlp == 'D':
                df3_dlpscore = 1
            
            if dlpscore > df3_dlpscore:
                df2[(df2['SampleID']==samplename)&(df2['LibraryID']==libname)&(df2['TAG1']==tag1)&(df2['TAG2']==tag2)&(df2['TAG3']==tag3)].loc[:,'DLP'] = dlptype
            
            plane_lst = set([planetype,df3_plane])
            if len(plane_lst) >0:
                if '0' in plane_lst:
                    plane_lst.remove('0')
            line_lst = set([linetype,df3_line])
            if len(line_lst) >0:
                if '0' in line_lst:
                    line_lst.remove('0')
                if np.nan in line_lst:
                    line_lst.remove(np.nan)
            
            if sarscore > df3_sarscore:
                df2[(df2['SampleID']==samplename)&(df2['LibraryID']==libname)&(df2['TAG1']==tag1)&(df2['TAG2']==tag2)&(df2['TAG3']==tag3)].loc[:,'SARScore'] = str(sarscore)
            df2['Plane'][(df2['SampleID']==samplename)&(df2['LibraryID']==libname)&(df2['TAG1']==tag1)&(df2['TAG2']==tag2)&(df2['TAG3']==tag3)] = ';'.join(str(i) for i in plane_lst)
            df2['Line'][(df2['SampleID']==samplename)&(df2['LibraryID']==libname)&(df2['TAG1']==tag1)&(df2['TAG2']==tag2)&(df2['TAG3']==tag3)] = ';'.join(str(i) for i in line_lst)
            if len(line_lst)>=2:
                df2['DLP'][(df2['SampleID']==samplename)&(df2['LibraryID']==libname)&(df2['TAG1']==tag1)&(df2['TAG2']==tag2)&(df2['TAG3']==tag3)] = 'DLL'
    df2.to_csv(out_file, sep='\t', index=False)        
     
def parse_args():
    """Parses arguments from cmd"""
    parser = argparse.ArgumentParser(description='get optimized compounds')
    parser.add_argument('-i',"--Input",
                        help='input fp dlp summary file',
                        type=str, required=True, metavar='<string>')
    parser.add_argument('-o1',"--outfilepath",
                        help='output off-DNA summary dir',
                        type=str, required=True, metavar='<string>')
    parser.add_argument('-o2',"--vinapath",
                        help='output vina score dir',
                        type=str, required=True, metavar='<string>')
    parser.add_argument('-a',"--Analysis",
                        help='analysis path, and the directory of "06.dwar" will be used.',
                        type=str, required=True, metavar='<string>')
    parser.add_argument('-s',"--SampleInfo",
                        help='SampleInfo file',
                        type=str, required=True, metavar='<string>')
    parser.add_argument('-n',"--njobs",
                        help='the number of jobs for novelty calculation',
                        required=True)
    return parser.parse_args()
    
def main():
    """Main function"""
    args = parse_args()
    outfilepath = args.outfilepath
    print('Start generating offDNA compound...'+'\n')
    OutputtempDNA = outfilepath + '/temp0.txt'
    offDNAmain(args.Input,OutputtempDNA,args.Analysis,'20')
    print('Start drop duplication...'+'\n')
    Outputtemp2 = outfilepath + '/temp1.txt'
    duplication(OutputtempDNA, Outputtemp2)
    print('Start calculating phisicalchemical properties...'+'\n')
    
    commentfile = sys.path[0] + '/../../lib/comments.txt'
    outputdocking = args.vinapath
    process.mkdir(outputdocking)   
    clustering(Outputtemp2,outfilepath,commentfile,args.Analysis,args.njobs,args.SampleInfo,outputdocking)
    os.remove(OutputtempDNA)
    os.remove(Outputtemp2)
    print('offDNA end!')

if __name__ == "__main__":
    main()
    
    
    
