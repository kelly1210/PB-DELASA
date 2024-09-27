#! python3
# -*- encoding: utf-8 -*-
'''
@File    :  clusfps_deal.py
@Time    :  2023/04/28 11:11:50
@Autor   :  Meng_xiangfei
@Version :  1.0
Contact  :  meng_xiangfei@pharmablock.com
DESC     :  根据相似性阈值或者指定簇的数量对smiles进行聚类分析,input、compound列名、output为必需参数；
            默认按照相似性分簇，即algorithm为b，对应Butina算法，使用morgan指纹，radius为2，cutoff值(1-similarity)为0.3，可以根据需要调整上述参数；
            若按照指定簇数分簇，则设置algorithm为m，对应Murtagh算法（默认使用Wards层次聚类算法），另外需要设置簇数nclusts；
'''
from rdkit import Chem
from rdkit.Chem import AllChem, MACCSkeys
#from rdkit.SimDivFilters import rdSimDivPickers
from rdkit.Chem.Fingerprints import FingerprintMols
from rdkit import DataStructs
from rdkit.ML.Cluster import Butina,Murtagh
from rdkit.ML.Cluster import ClusterUtils
import pandas as pd
import argparse
#from PyBioMed.PyMolecule.fingerprint import CalculateAtomPairsFingerprint,CalculateDaylightFingerprint
#from PyBioMed.PyGetMol.Getmol import ReadMolFromSmile
#from PyBioMed.PyMolecule.moe import GetMOE
#from PyBioMed.PyMolecule.fingerprint import CalculateEstateFingerprint, CalculateAtomPairsFingerprint, \
#    CalculateDaylightFingerprint, CalculateECFP2Fingerprint, CalculateECFP4Fingerprint, CalculateECFP6Fingerprint, \
#    CalculateFCFP2Fingerprint, CalculateFCFP4Fingerprint, CalculateFCFP6Fingerprint, CalculateFP2Fingerprint, \
#    CalculateMACCSFingerprint, CalculatePharm2D2pointFingerprint, CalculatePharm2D3pointFingerprint, \
#    CalculateTopologicalTorsionFingerprint, CalculateGhoseCrippenFingerprint, CalculatePubChemFingerprint, \
#    CalculateFP3Fingerprint, CalculateFP4Fingerprint
#/home/data2/dong_keke/del_ai/bin/feature_module_pybiomed_rdkit_muti_prosess.py
# 分子指纹
fpOpdict = {'tp':'Topological Fingerprints','mc':'MACCS Keys','mo':'Morgan Fingerprints'}
# 聚类算法,Butina需要根据cutoff对Compound进行聚类分析,Murtagh需要根据指定的簇数量对Compound进行聚类分析
algOpdict = {'b':'Butina','m':'Murtagh'}

class ChemParse(object):
    def __init__( self ,input):
        self.source = ''
        self.fps = ''
        self.input = input
        
    def sdf_reader(self):
        self.source = Chem.SDMolSupplier(self.input)
    
    # 新增读取csv格式（tab作为分隔符）smiles文件
    def smi_reader(self,column):
        smi_df = pd.read_csv(self.input,index_col=None,header=0,sep="\t",encoding='gb18030')
        smis = smi_df[column].tolist()
        #print(smis)
        self.source = [Chem.MolFromSmiles(smi) for smi in smis]
        return smi_df
    def df_reader(self,column):
        #smi_df = pd.read_csv(self.input,index_col=None,header=0,sep="\t")
        smis = self.input[column].tolist()
        #print(smis)
        self.source = [Chem.MolFromSmiles(smi) for smi in smis]
        #return self.input
    def smiList_reader(self):
        self.source = [Chem.MolFromSmiles(smi) for smi in self.input]

    def get_fps(self,ftype,radius):
        if self.source == '' or self.source == None:
            pass
        else:
            if ftype == 'mo':
                fps=[AllChem.GetMorganFingerprint(m,radius,useFeatures=True) for m in self.source]
            elif ftype == 'moecfp':
                fps=[AllChem.GetMorganFingerprint(m,radius,useFeatures=False) for m in self.source]
            elif ftype == 'tp':
                fps = [FingerprintMols.FingerprintMol(m) for m in self.source]
            elif ftype == 'mc':
                fps = [MACCSkeys.GenMACCSKeys(m) for m in self.source]
            else:
                print('Error fingertype!')
            # 注释掉的这几种分子指纹计算相似性时不能使用rdkit中自带的相似性方法计算
                
        self.fps = fps

    def clusterOutput(self, output, cdict):
        sdfout = Chem.SDWriter(output)
        
        for index, m in enumerate(self.source):
            classid = cdict[index][0]
            isCentroid = cdict[index][1]
            m.SetProp("class",str(classid))
            m.SetProp("isCentroid",isCentroid)
            sdfout.write(m)
        sdfout.close()
        

class Fingerprint_Cluster(object):
    def __init__( self, fps):
        self.fplist = fps
        self.dist = []
        self.cdict = {} 
        self.clustdict = {}

    #generate the distance matrix
    def distance_matrix(self):
        self.dist = []
        nfps = len(self.fplist)
        for i in range(1,nfps):
            sims = DataStructs.BulkTanimotoSimilarity(self.fplist[i],self.fplist[:i])
            self.dist.extend([1-x for x in sims])
        
    #generate cluster dict as {1:[1,2,3],2:[4,5,6]...}
    def cluster_dict(self,algorithm, cutoff=0.7, ncluster=1):
        if algorithm == 'Butina':
            self.ClusterFps_Butina(self.dist,len(self.fplist),cutoff)
        elif algorithm == 'Murtagh':
            self.ClusterFps_Murtagh(self.dist,len(self.fplist),ncluster)
        
    def ClusterFps_Butina(self, dists, nfps,cutoff):
        self.cdict = {}
        cs = Butina.ClusterData(dists,nfps,cutoff,isDistData=True)
        for index, eachcs in enumerate(cs): 
            self.clustdict[index+1] = eachcs
            for eachid in eachcs:
                self.cdict[eachid] = [index+1]
                if eachid == eachcs[0]:
                    self.cdict[eachid].append("true")
                else:
                    self.cdict[eachid].append("flase")

    def ClusterFps_Murtagh(self, dists, nfps, ncluster):
        self.cdict = {}
        cs= Murtagh.ClusterData(dists,len(self.fplist),Murtagh.WARDS,isDistData=1)
        
        """
        # 不同层次聚类算法
        if method == 'Wards':
            cs = Murtagh.ClusterData(dists,len(self.fplist),Murtagh.WARDS,isDistData=1)
        elif method == 'SLINK':
            cs = Murtagh.ClusterData(dists,len(self.fplist),Murtagh.SLINK,isDistData=1)
        elif method == 'CLINK':
            cs = Murtagh.ClusterData(dists,len(self.fplist),Murtagh.CLINK,isDistData=1)
        elif method == 'UPGMA':
            cs = Murtagh.ClusterData(dists,len(self.fplist),Murtagh.UPGMA,isDistData=1)
        """
        splitClusts=ClusterUtils.SplitIntoNClusters(cs[0],ncluster)
        #centroids = [ClusterUtils.FindClusterCentroidFromDists(x,dists) for x in splitClusts]
        for index, cluster in enumerate(splitClusts):
            children = cluster.GetPoints() 
            pts = [x.GetData() for x in children]  
            self.clustdict[index+1] = pts 
            for pt in pts:
                self.cdict[pt] = [index + 1] 
                if pt == pts[0]:
                    self.cdict[pt].append("true")
                else:
                    self.cdict[pt].append("flase")
def clusterDF(smi_df,column,fp_type,clu_method,cutoff=0.3,nclusters=1,radius=2,clusterIDname='ClusterID'):
    # 直接读取dataframe
    # column = 'Compound'
    smiparse = ChemParse(smi_df)
    smiparse.df_reader(column)
    smiparse.get_fps(fp_type,radius)
    fpCluster = Fingerprint_Cluster(smiparse.fps)
    fpCluster.distance_matrix()
    fpCluster.cluster_dict(clu_method,cutoff,nclusters)
    for c in fpCluster.clustdict:
        for id in fpCluster.clustdict[c]:
            smi_df.loc[id,clusterIDname] = str(c)
    return smi_df
    #smi_df.to_csv(output,index=None,header=True,sep="\t")

def main(smi_file,column,fp_type,clu_method,output,cutoff=0.3,nclusters=1,radius=2):
    smiparse = ChemParse(smi_file)
    smi_df = smiparse.smi_reader(column)
    smiparse.get_fps(fp_type,radius)
    fpCluster = Fingerprint_Cluster(smiparse.fps)
    fpCluster.distance_matrix()
    fpCluster.cluster_dict(clu_method,cutoff,nclusters)
    for c in fpCluster.clustdict:
        for id in fpCluster.clustdict[c]:
            smi_df.loc[id,"ClusterID"] = str(c)
    smi_df.to_csv(output,index=None,header=True,sep="\t")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--smi',dest='input',help='intput smiles file',required=True)
    parser.add_argument('--column',dest='column',help='column name of smiles',default='Smiles',required=True)
    parser.add_argument('--fp',dest='fp',help='fingerprint type: tp,mc,mo (Topological Fingerprints, MACCS Keys, Morgan Fingerprints), default is mo', default='mo')
    parser.add_argument('--algorithm',dest='algorithm',help='cluster algorithm :b,m (Butina, Murtagh), default is b', default='b')
    parser.add_argument('--out',dest='output',help='output sdf file',required=True)
    parser.add_argument('--cutoff',dest='cutoff',help='distThresh(0-1),elements within this range of each other are considered to be neighbors, needed for Butina cluster algorithm, default is 0.3', type=float, default=0.3)
    parser.add_argument('--nclusts',dest='nclusts',help='number of clusters, needed for Murtagh cluster algorithm, default is 1',type=int, default=1)
    parser.add_argument('--radius',dest='radius',help=' the radius of the Morgan fingerprint, default is 2',type=int, default=2)
    args = parser.parse_args()
    args.algorithm = algOpdict[args.algorithm]
    main(args.input,args.column,args.fp,args.algorithm,args.output,args.cutoff,args.nclusts,args.radius)