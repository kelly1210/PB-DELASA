import time
import argparse
import os 
import sys
sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), '../monitor'))

#QC Summary of nopaired sample.
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
    print("##############################")
    print("QC_Summary initiated.\n")
    data_start = time.time()
    if (not os.path.exists(TrimQcFile)):
        print("Error: the file {} is not exists!".format(TrimQcFile))
        sys.exit(1)
    qc_summary_list = qc_summary(TrimQcFile, casesampleID, ColumnName)
    with open(QcSummaryFile,"w") as f:
        for i in range(len(qc_summary_list)):
            list_str = []
            for data in qc_summary_list[i]:
                data = str(data).strip()
                list_str.append(data)
                line_str = "\t".join(list_str)+"\n"
                f.write(line_str)
    f.close()
    print(time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(time.time())) + " QC Summary: done. output QC Summary: " + QcSummaryFile)
    print ("Total time of nopair qc summary: " + str(time.time()-data_start))

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description = 'summary qc results')
    parser.add_argument('--input','-i',dest="TrimQCfile",help = "TrimQCfile")
    parser.add_argument('--prefix', '-p', dest="Prefix", help="CasesampleID")
    parser.add_argument('--output','-o',dest="QcSummaryFile",help = "qc summary output")
    parser.add_argument('--SpecColumn','-s',dest="SpecColumn",help = "special column name")
    
    args = parser.parse_args()
    run_qc(args.TrimQCfile, args.QcSummaryFile, args.Prefix, args.SpecColumn)

