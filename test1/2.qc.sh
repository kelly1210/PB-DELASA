#!/bin/sh
# ../pipeline/EF_pair.py
#python /home/data2/dong_keke/software/PB-DELASA_v/bin/monitor/../pipeline/EF_pair.py -i /home/data2/dong_keke/software/PB-DELASA_v/example/NGS001 -b NGS001 -s /home/data2/dong_keke/software/PB-DELASA_v/test1/sample_info.txt -l /home/data2/dong_keke/software/PB-DELASA_v/lib/library_info.txt -mpi 80

# ../pipeline/auto_png_pair.py
#python /home/data2/dong_keke/software/PB-DELASA_v/bin/monitor/../pipeline/auto_png_pair.py -d /home/data2/dong_keke/software/PB-DELASA_v/example/NGS001 -s /home/data2/dong_keke/software/PB-DELASA_v/test1/sample_info.txt -l /home/data2/dong_keke/software/PB-DELASA_v/lib/library_info.txt -mpi 80 -single NO

# ../pipeline/IdentifyFP.py
python /home/data2/dong_keke/software/PB-DELASA_v/bin/monitor/../pipeline/IdentifyFP.py -i /home/data2/dong_keke/software/PB-DELASA_v/example/NGS001/09.report/table/NGS001.lib.summary.qc.txt -o /home/data2/dong_keke/software/PB-DELASA_v/example/NGS001/09.report/table/NGS001.count.sum.dlp.detail.txt -b NGS001 -s /home/data2/dong_keke/software/PB-DELASA_v/test1/sample_info.txt -l /home/data2/dong_keke/software/PB-DELASA_v/lib/library_info.txt

# ../pipeline/QCsummary.py
python /home/data2/dong_keke/software/PB-DELASA_v/bin/monitor/../pipeline/QCsummary.py -i /home/data2/dong_keke/software/PB-DELASA_v/example/NGS001/09.report/table/ -b NGS001 -s /home/data2/dong_keke/software/PB-DELASA_v/test1/sample_info.txt -l /home/data2/dong_keke/software/PB-DELASA_v/lib/library_info.txt -o /home/data2/dong_keke/software/PB-DELASA_v/example/NGS001/09.report/table/NGS001.qc.xlsx

# ../offDNA/DEL_predict/DEL_predict.py
python /home/data2/dong_keke/software/PB-DELASA_v/bin/monitor/../offDNA/DEL_predict/DEL_predict.py -s /home/data2/dong_keke/software/PB-DELASA_v/test1/sample_info.txt -i1 /home/data2/dong_keke/software/PB-DELASA_v/example/NGS001/09.report/table/NGS001.count.sum.dlp.detail.txt -i2 /home/data2/dong_keke/software/PB-DELASA_v/example/NGS001/06.Dwar -m /home/data2/dong_keke/software/PB-DELASA_v/bin/monitor/../offDNA/DEL_predict/RF.pkl -o /home/data2/dong_keke/software/PB-DELASA_v/example/NGS001/09.report/table/NGS001.count.sum.dlp.detail.confirmation.txt

# ../offDNA/offDNAenu.py
python /home/data2/dong_keke/software/PB-DELASA_v/bin/monitor/../offDNA/offDNAenu.py -i /home/data2/dong_keke/software/PB-DELASA_v/example/NGS001/09.report/table/NGS001.count.sum.dlp.detail.confirmation.txt -o1 /home/data2/dong_keke/software/PB-DELASA_v/example/NGS001/09.report/table/ -o2 /home/data2/dong_keke/software/PB-DELASA_v/example/NGS001/VinaDocking/ -a /home/data2/dong_keke/software/PB-DELASA_v/example/NGS001 -s /home/data2/dong_keke/software/PB-DELASA_v/test1/sample_info.txt -n 80

    
