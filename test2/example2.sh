#!/bin/sh
prefix="NGS001"
out_dir="../"$prefix

python ../bin/offDNA/offDNAenu_v10.py \
-i1 ${out_dir}/09.report/table/test.txt \
-i2 ../lib/comments.txt -o1 ${outdir}/09.report/table/ \
-o2 ${outdir}/VinaDocking/ \
-a ${outdir} \
-s ./sample_info.txt \
-VinaProg ${CBDOCK}/AutoBlindDock.pl \
-pdbPath ./ \
-ProLst ./proteinlist.txt \
-targetDir ../bin/offDNA/TargetInhibitors -cpu 10

