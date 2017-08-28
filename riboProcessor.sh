#!/usr/bin/bash
export NM=$1
source $2
cd ${NM}
module add igmm/apps/R/3.3.3 igmm/apps/BEDTools/2.26.0 igmm/apps/samtools/1.2 igmm/apps/texlive/2015
export DESC="Not described"
if [ -f desc ]
 then
  export DESC=`head -n1 desc`
fi

# Nuclear
grep -v chrM ${NM}.${GENOME}.riboCounts.bed | grep '+' | perl -walne '$s=$F[1]-16;$e=$F[1]+15;$c=0;if($F[4]){$c=$F[4];};next if($s<0);print join"\t",($F[0],$s,$e,$c,$c,"+");' > ${NM}.nuc.plusStrandFrag.bed
grep -v chrM ${NM}.${GENOME}.riboCounts.bed | grep '-' | perl -walne '$s=$F[2]-15;$e=$F[2]+16;$c=0;if($F[4]){$c=$F[4];};next if($s<0);print join"\t",($F[0],$s,$e,$c,$c,"-");' > ${NM}.nuc.minusStrandFrag.bed
# chrM
grep chrM ${NM}.${GENOME}.riboCounts.bed | grep '+' | perl -walne '$s=$F[1]-15;$e=$F[1]+16;$c=0;if($F[4]){$c=$F[4];};next if($s<0);print join"\t",($F[0],$s,$e,$c,$c,"+");' > ${NM}.chrM.plusStrandFrag.bed
grep chrM ${NM}.${GENOME}.riboCounts.bed | grep '-' | perl -walne '$s=$F[2]-16;$e=$F[2]+15;$c=0;if($F[4]){$c=$F[4];};next if($s<0);print join"\t",($F[0],$s,$e,$c,$c,"-");' > ${NM}.chrM.minusStrandFrag.bed


bedtools getfasta -fi ${GENOMEFASTA} -bed ${NM}.nuc.plusStrandFrag.bed -name -tab -s -fo ${NM}_nuc.plusStrandFrag_faTab 2> /dev/null
bedtools getfasta -fi ${GENOMEFASTA} -bed ${NM}.nuc.minusStrandFrag.bed -name -tab -s -fo ${NM}_nuc.minusStrandFrag_faTab 2> /dev/null

bedtools getfasta -fi ${GENOMEFASTA} -bed ${NM}.chrM.plusStrandFrag.bed -name -tab -s -fo ${NM}_chrM_plusStrandFrag_faTab 2> /dev/null
bedtools getfasta -fi ${GENOMEFASTA} -bed ${NM}.chrM.minusStrandFrag.bed -name -tab -s -fo ${NM}_chrM_minusStrandFrag_faTab 2> /dev/null

Rscript --vanilla ${SRC}/faTabLogo.R ${NM}_chrM_minusStrandFrag_faTab chrM_rev ${GENOSHORT}ChrMrev 1 0
Rscript --vanilla ${SRC}/faTabLogo.R ${NM}_chrM_plusStrandFrag_faTab chrM_for ${GENOSHORT}ChrMfor 1 0
rm ${NM}_chrM_plusStrandFrag_faTab ${NM}_chrM_minusStrandFrag_faTab
rm ${NM}_nuc.plusStrandFrag_faTab ${NM}_nuc.minusStrandFrag_faTab
rm ${NM}.nuc.plusStrandFrag.bed ${NM}.nuc.minusStrandFrag.bed ${NM}.chrM.plusStrandFrag.bed ${NM}.chrM.minusStrandFrag.bed

###
# Read end biases
grep -v '=' ${NM}.postTrim.fastq | grep -v N | grep -v E | perl -walne 'if($F[0]){if($F[0]=~/^[ACGT]+$/){if(length($F[0])>30){$o=substr($F[0],0,10);print"read\t$o"}}}' > ${NM}_startN10.faTab
grep -v '=' ${NM}.postTrim.fastq | grep -v N | grep -v E | perl -walne 'if($F[0]){if($F[0]=~/^[ACGT]+$/){if(length($F[0])>30){$o=substr($F[0],-10,10);print"read\t$o"}}}' > ${NM}_endN10.faTab

Rscript --vanilla ${SRC}/faTabLogo.R ${NM}_startN10.faTab readStarts equal 0 1
Rscript --vanilla ${SRC}/faTabLogo.R ${NM}_endN10.faTab readEnds equal 0 1
rm ${NM}_startN10.faTab
rm ${NM}_endN10.faTab
###

###
# Read length distributions
grep -v '=' ${NM}.postTrim.fastq | grep -v N | grep -v E | perl -walne 'if($F[0]){if($F[0]=~/^[ACGT]+$/){print length($F[0]);}}' | perl -w ${SRC}/freqdist.pl 0 > ${NM}.trimmedSeqDist
samtools idxstats ${NM}.${GENOME}.q30.sorted.bam > ${NM}.idxstats
samtools view ${NM}.${GENOME}.q30.sorted.bam | awk '{print$10}' | perl -walne 'if($F[0]){if($F[0]=~/^[ACGT]+$/){print length($F[0]);}}' | perl -w ${SRC}/freqdist.pl 0 > ${NM}.alnSeqDist
samtools view ${NM}.${GENOME}.q30.sorted.bam chrM | awk '{print$10}' | perl -walne 'if($F[0]){if($F[0]=~/^[ACGT]+$/){print length($F[0]);}}' | perl -w ${SRC}/freqdist.pl 0 > ${NM}.alnSeqDistchrM
Rscript --vanilla ${SRC}/summaryNumbers.R
###

###
# Produce triRate analyses
Rscript --vanilla ${SRC}/triRateVecCalculation.R
###

###
# Mitochondrial regional ANALYSIS
grep chrM ${NM}.toLoad > ${NM}.chrM.toLoad
Rscript --vanilla ${SRC}/mcircosPlotAuto.R
# For some reason the higher density sampling breaks on Eddie
#convert -density 1200 -flatten -background white mCircos.svg mCircos.jpg
convert mCircos.svg mCircos.jpg
###
# Compile report
cp ${SRC}/combined.tex ./
pdflatex combined.tex
cp combined.pdf ${NM}_summary.pdf
