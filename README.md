## PB-DELASA
version 1.0  
data 2024-0103  
PB-DELASA is a pipeline for DEL analysis that designed on N12 UMI.  
Installation requirements for PB-DELASA:  
### 1. enviroment  
conda create -n pbdel python==3.12  
conda activate pbdel  
### 2. Installation  
conda install -c bioconda fastp  
conda install -c bioconda pandaseq  
conda install conda-forge/label/cf202003::openbabel 3.0 
conda install rdkit -c rdkit  
conda install -c bioconda cutadapt  
conda install -c anaconda scikit-learn  
pip install openpyxl  
conda install anaconda::joblib  
wget https://github.com/jameslz/fastx-utils/raw/master/seqtk_demultiplex  
chmod +x seqtk_demultiplex  
### 3. seqtk  
git clone https://github.com/lh3/seqtk.git  
cd seqtk  
make  
### 4. vina  
wget https://github.com/ccsb-scripps/AutoDock-Vina/releases/download/v1.2.5/vina_1.2.5_linux_x86_64  
### 5. mgltools  
(https://ccsb.scripps.edu/download/532/)  
tar zxvf mgltools_x86_64Linux2_1.5.7.tar.gz  
cd mgltools_x86_64Linux2_1.5.7  
bash install.sh OR ./install.sh -d ../mgltools-bin  
source ./initMGLtools.sh  
export PATH=/data/software/mgltools-bin/bin:$PATH  
### 6. CB-dock(https://cadd.labshare.cn/cb-dock/php/manual.php)  
tar -xvzf CB-Dock.tar.gz  
cd CB-Dock  
./setup.sh [path to mgltools]/bin/python [path to vina]/vina  
export VinaProg= [path_to_cb-dock]/CB-Dock/prog/  
### 7. Datawarrior(https://openmolecules.org/datawarrior/download.html)  
tar zxvf datawarrior610.tar.gz  
cd datawarrior_linux  
install based on readme.txt and install.sh file(root)  
### 8. path  
export seqtk_demultiplex=[path_to_seqtk_demultiplex]/seqtk_demultiplex  
export datawarrior_app=[path_to_datawarrior]/datawarrior  
export seqtk_software=[path_to_seqtk]/seqtk  
export VinaProg=[path_to_vina]/CB-Dock/prog/  
export cutadaptPATH=[path_to_cutadapt]/  
### 9. usage  
pbdel.py

                [-h] -fq1 <string> -fq2 <string> -inPCR <int> -barcode <string>

               -c <string> -o <string> -p <string> -s <string> -d <string> 
               
               -l <string> -code <string> -r <string>

optional arguments:

  -h, --help            
  
                        show this help message and exit
  
  -fq1 <string>, --caseFq1 <string>
  
                        The fastq file read1. The format is *.fastq.gz or
                        
                        *.fastq.
                        
  -fq2 <string>, --caseFq2 <string>
  
                        The fastq file read2. The format is *.fastq.gz or
                        
                        *.fastq.
                        
  -inPCR <int>         
  
                        1 if constant in PCR, 0 or not
  
  -barcode <string>     Output barcode file, which containg the barcodes of
  
                        all samples.
                        
  -c <string>, --config <string>
  
                        Sample special configure file
                        
  -d <string>, --configDir <string>
  
                        Output config directionary for each sample.
                        
  -r <string>, --runscript_path <string>
  
                        Run script directory for each batch data
                        
  -o <string>, --output_path <string>
  
                        Output directory for each batch data
                        
  -p <string>, --prefix <string>
  
                        The name of each batch data.
                        
  -s <string>, --samples <string>
  
                        A table for sample information, txt file and 'TAB' as
                        
                        the seperate.
                        
  -l <string>, --libraryInfo <string>
  
                        A table for library information, txt file and 'TAB' as
                        
                        the seperate.
                        
  -code <string>, --codeDir <string>
  
                        TagCode~structure file dir
#### Example 1:
inputfile: fastq file, code file, config file, sample_info.txt, library_info.txt  
outfile: outdir(dwar file, png file, html)  
process: decoding and off-DNA optimizing  
nohup ./example1.sh &  
The result is in ./example/NGS001:  
├── 00.Script  # script for each sample  
├── 01.QC  # qc for fastq file  
├── 01.Trim  # generate trimed fastq file   
├── 02.Merge  # generate merged fasta file for paired-end fastq file  
├── 03.Split  # split sequence based on sample and library  
├── 04.Counts  # sequences for each libary of each sample  
├── 05.Count.analysis  # combinated sequence counts  
├── 06.Dwar  # generate dwar file for visualization  
├── 07.Pair  # analysis for paired samples, such as case and ntc samples  
├── 08.Contamination  # contatmination sequence  
├── 09.report  # tables and png files  
└── VinaDocking  # docking result  

#### Example 2:

process: off-DNA recommend  

nohup ./example2.sh &  

#### New data  
Raw sequencing data storage: ./example/*_r1.fastq.gz, ./example/*_r2.fastq.gz  

To add new samples or modify existing libraries, update the following files:  
1. ./lib/library_info.txt**: Contains metadata for each library, including:  
   - LibraryID  
   - LibrarySize  
   - PCR1_5 (indicates constant sequence presence on PCR1_5; write 0 if absent)  
   - constant1 sequence  
   - Cycle composition format  
   - Sequences in 5→3 orientation  
   - CP_umi (UMI-containing closing primer sequences with structure notation)  
   - PCR1_3 (constant sequence presence on PCR1_3)  
   - positive_sample (indicates if library contains positive controls)  

2. ./lib/Libraries/*.code.smiles.txt  
   - Template files (e.g., libtest.code.smiles.txt) defining BB-sequence mappings  
   - Format according to provided examples  

3. ./test1/sample_info.txt  
   - Sample metadata fields:  
     *ExperimentID, Protein name, Protein Concentration, Screening Round number, SampleID, barcode in I5, barcode in I3, Library ID, positive_sample, dedup_method, plot_threshold, NTC, QPCR, Pair, Pocket, Pocket_diff, Client*  
   - Default parameters may be used  

4. ./test1/barcodes.txt**  
   - Automatically generated barcode manifest file based on "sample_info.txt" 


