import os
import argparse
import sys
import pandas as pd
import subprocess
import time
import shutil
from multiprocessing import Pool
sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), '../monitor'))
#sys.path.append('/opt/software/tagFinder-inhouse/bin/monitor')
import process

datawarrior_app = os.environ["datawarrior_app"]

def PNGtemplete(sample_id, dwar_path, dwar_file, png_file, dwam_file):
    """
    :param sample_id: 'DEL008_PBDEL0001'
    :param dwar_path: '{}/06.Dwar/P013/DEL008/PBDEL0001.count.freq.filter.dwar'
    :param dwar_file: 'PBDEL0001.count.freq.filter.dwar'
    :param png_file: '{}/06.Dwar/P013/DEL008/test.png'
    :param png_file: '{}/06.Dwar/P013/DEL008/test.dwam'
    :return:
    """

    template = """<macro name="{0}">
<task name="openFile">
fileName={1}
</task>
<task name="selectWindow">
viewName={2}
</task>
<task name="selectView">
viewName=3D Raw_counts
</task>
<task name="createViewImage">
fileName={3}
viewName=3D Raw_counts
allViews=true
width=1456
format=png
keepAspect=true
dpi=150
transparentBG=false
height=692
target=file
</task>
<task name="selectView">
viewName=Table
</task>
<task name="closeWindow">
viewName={2}
</task>
</macro>""".format(sample_id, dwar_path, dwar_file, png_file)
    with open(dwam_file, 'w') as dwam:
        dwam.write(template)


def autoPNG(input_path, sample_info, library_info, njobs,single):  

    """Automatically generate filtered PNG images, store them in the report folder, 
    and copy a copy to the dwar folder in the report for easy comparative analysis.
    """
    protein_dict = process.total_sample_info(sample_info)[0]
    library_dict = process.total_sample_info(sample_info)[2]
    ntc_dict = process.total_sample_info(sample_info)[8]
    
    sample_lst = process.total_sample_info(sample_info)[6]
    dlp_dir = input_path + "/09.report/dlp/"
    process.mkdir(dlp_dir)
    
    process.mkdir(input_path + '/09.report/PNG_pair/')
    process.mkdir(input_path + "/09.report/PNG_pair/dwam")
    
    # 生成dlpSummaryFile
    paths = input_path.strip().split('/')
    paths_new = [i for i in paths if i!='']
    dlpSummaryFile = input_path + '/09.report/table/'+paths_new[-1] + '.lib.summary.qc.txt'
    process.mergefile(input_path + '/06.Dwar/','.count.sum.dlp.txt',dlpSummaryFile)
    for sample in sample_lst:
        print(sample)
        # The names of all libraries corresponding to samples with protein concentration.
        libraryID_names = library_dict[sample].split(',')
        #p=Pool(int(njobs))
        for library_name in libraryID_names:
            #p.apply_async(singleLibPNG,args=(sample,library_name,input_path,protein_dict,single))
            singleLibPNG(sample,library_name,input_path,protein_dict,single)
        #p.close()
        #p.join()
        
        time.sleep(40)
        
        #p=Pool(int(njobs))
        for library_name in libraryID_names:
    #        p.apply_async(singleLibCP,args=(sample,library_name,input_path,protein_dict,ntc_dict,dlpSummaryFile,single))
            singleLibCP(sample,library_name,input_path,protein_dict,ntc_dict,dlpSummaryFile,single)
        #p.close()
        #p.join()
        
        #shutil.rmtree(input_path + "/09.report/PNG/dwam) # remove dwam dir

def singleLibPNG(sample,library_name,input_path,protein_dict,single='NO'):
    sample_id = sample + '.' + library_name
    protein = protein_dict[sample]
    
    dwar_file_single = sample_id + '.count.freq.filter.dwar'
    dwar_path_single = input_path + "/06.Dwar/" + protein + '/' + library_name + '/' + dwar_file_single
    dwar_file = sample_id + '.pair.ef.dwar'
    dwar_path = input_path + "/07.Pair/" + protein + '/' + library_name + '/' + dwar_file
    dwam_file = input_path + "/09.report/PNG_pair/dwam/" + sample_id + '.pair.dwam'
    png_file = input_path + "/09.report/PNG_pair/" + sample_id + '.pair.png'
    # paired png file
    if os.path.exists(dwar_path):
        if os.path.exists(dwam_file) and os.path.exists(png_file): 
            print('The file {} is existed!'.format(png_file))
        else:
            PNGtemplete(sample_id, dwar_path, dwar_file, png_file, dwam_file)
            cmd = "{0} {1}".format(datawarrior_app, dwam_file)
            print(cmd)
            subprocess.check_call(cmd, shell=True)
    else:
        print('The file {} is not exists!'.format(dwar_file))
    # single png file,The PNGs for the default paired samples have been generated, but the PNGs for individual samples will not be generated.
    if single == 'YES':
        process.mkdir(input_path + '/09.report/PNG/')
        process.mkdir(input_path + "/09.report/PNG/dwam")
        dwam_file_single = input_path + "/09.report/PNG/dwam/" + sample_id + '.dwam'
        png_file_single = input_path + "/09.report/PNG/" + sample_id + '.png'
        if os.path.exists(dwar_path_single):
            if os.path.exists(dwam_file_single) and os.path.exists(png_file_single):
                print('The file {} is existed!'.format(png_file_single))
            else:
                PNGtemplete(sample_id, dwar_path_single, dwar_file_single, png_file_single, dwam_file_single)
                cmd = "{0} {1}".format(datawarrior_app, dwam_file_single)
                subprocess.check_call(cmd, shell=True)
        else:
            print('The file {} is not exists!'.format(dwar_file))
    return 1
        
def singleLibCP(sample,library_name,input_path,protein_dict,ntc_dict,dlpSummaryFile,single='NO'):
    sample_id = sample + '.' + library_name
    protein = protein_dict[sample]
    dwar_file_single = sample_id + '.count.freq.filter.dwar'
    dwar_path_single = input_path + "/06.Dwar/" + protein + '/' + library_name + '/' + dwar_file_single
    #dwam_file_single = input_path + "/09.report/PNG/dwam/" + sample_id + '.dwam'
    png_file_single = input_path + "/09.report/PNG/" + sample_id + '.png'
    
    dwar_file = sample_id + '.pair.ef.dwar'
    dwar_path = input_path + "/07.Pair/" + protein + '/' + library_name + '/' + dwar_file
    filterDLP = input_path + "/06.Dwar/" + protein + '/' + library_name + '/' + sample_id + '.count.freq.filter.dlp.xlsx'
    #dwam_file = input_path + "/09.report/PNG_pair/dwam/" + sample_id + '.pair.dwam'
    png_file = input_path + "/09.report/PNG_pair/" + sample_id + '.pair.png'
    if os.path.exists(dwar_path_single):
        if 'NTC' not in protein:
            df = pd.read_csv(dlpSummaryFile,sep='\t',header=0, dtype=object)
            df1 = df[(df['Protein'] == protein) & (df['LibraryID'] == library_name) & (df['DLP']!='0')]
            sampleList = list(set(df1['SampleID'].tolist()))
            if len(sampleList)!=0:
                report_dwar_path = input_path + "/09.report/dlp/" + protein + '/' + library_name
                process.mkdir(report_dwar_path)
                if os.path.exists(png_file):
                    shutil.copy(png_file, report_dwar_path)
                if os.path.exists(dwar_path):
                    shutil.copy(dwar_path, report_dwar_path)
                if os.path.exists(filterDLP):
                    shutil.copy(filterDLP, report_dwar_path)
                if single=='YES':
                    if os.path.exists(png_file_single):
                        shutil.copy(png_file_single, report_dwar_path)
                    if os.path.exists(dwar_path_single):
                        shutil.copy(dwar_path_single, report_dwar_path)

                ntc_lst = ntc_dict[sample].strip().split(',')
                if ntc_lst[0] !='0':
                    for ntcsample in ntc_lst:
                        
                        ntc_dwar_path = input_path + "/06.Dwar/" + protein_dict[ntcsample] + '/' + library_name + '/' + ntcsample + '.' + library_name + '.count.freq.filter.dwar'
                        ntc_filterDLP = input_path + "/06.Dwar/" + protein_dict[ntcsample] + '/' + library_name + '/' + ntcsample + '.' + library_name + '.count.freq.filter.dlp.xlsx'
                        
                        ntc_png_file = input_path + "/09.report/PNG/" + ntcsample + '.' + library_name + '.png'
                        if os.path.exists(ntc_filterDLP):
                            if not os.path.exists('{}/{}'.format(report_dwar_path,ntcsample + '.' + library_name + '.count.freq.filter.dlp.xlsx')):
                                shutil.copy(ntc_filterDLP, report_dwar_path)
                        if single=='YES':
                            if os.path.exists(ntc_png_file):
                                if not os.path.exists('{}/{}'.format(report_dwar_path,ntcsample + '.' + library_name + '.png')):
                                    shutil.copy(ntc_png_file, report_dwar_path)
                            if os.path.exists(ntc_dwar_path):
                                if not os.path.exists('{}/{}'.format(report_dwar_path,ntcsample + '.' + library_name + '.count.freq.filter.dwar')):
                                    shutil.copy(ntc_dwar_path, report_dwar_path)
                        
                else:
                    print('There is no ntc sample for {}!!!'.format(sample_id))
            else:
                print('The sample {} is no dlp!!!'.format(sample_id))
                    
    else:
        print('The file {} is not exists!'.format(dwar_file))
    return 1


if __name__ == "__main__":
    data_start = time.time()
    #print(time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(time.time())))
    parser = argparse.ArgumentParser(description='Auto PNGfile')
    parser.add_argument('--dir', '-d', dest="InputPath", help="Input path")
    parser.add_argument('--s', '-s', dest="SampleInfo", help="Sample information file")
    parser.add_argument('--l', '-l', dest="LibraryInfo", help="Library information file")
    #parser.add_argument("--datawarrior", "-w", dest="datawarrior", help="datawarrior path")
    parser.add_argument('--mpi', '-mpi', dest="mpi", help='sample mpi numcores')
    parser.add_argument('--single', '-single', dest="single", help='output single png file or not;YES or NO, default NO')
    args = parser.parse_args()
    autoPNG(args.InputPath, args.SampleInfo, args.LibraryInfo, args.mpi, args.single)