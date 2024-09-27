#!/bin/sh
prefix="NGS002"
out_dir="../"$prefix

python ../bin/offDNA/offDNAenu_v10.py \
-i1 ${out_dir}/09.report/table/test.txt \
-i2 ../lib/comments_20210725.txt -o1 ${outdir}/09.report/table/ \
-o2 ${outdir}/VinaDocking/ \
-a ${outdir} \
-s ./sample_info.txt \
-VinaProg /home/data2/dong_keke/software/CB-Dock/prog/AutoBlindDock.pl \
-pdbPath ./ \
-ProLst ./proteinlist.txt \
-targetDir ../bin/offDNA/TargetInhibitors -cpu 10

