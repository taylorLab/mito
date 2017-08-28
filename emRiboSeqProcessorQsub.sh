#!/usr/bin/bash 
# Martin Taylor, June 2015, 
# updated 15th April 2016
# updated 21st July 2016 - use new mm9 genome with chrM and rDNA masks and including _random regions
# updated 3rd Nov 2016 to fix corner turning minus-strand transform, add optional start points and refine qsub options.
# updated 24 Feb 2017 to generalise for alternate genomes and aternate coordinate transforms, e.g. hyden-seq
module add igmm/apps/BEDTools/2.26.0 igmm/apps/bowtie/2.2.6 igmm/apps/samtools/1.2 igmm/apps/python/2.7.10
export NM=$1
export FA=$2
export PA=$3
export STARTSTEP=$4
source ${PA}

if [ ! -d "${SCRATCH}/${NM}" ]
   then
       mkdir -p "${SCRATCH}/${NM}"
fi

cd "${SCRATCH}/${NM}"


if [ ! ${STARTSTEP} ]
then
    export STARTSTEP=1
fi

if (( ${STARTSTEP} == 1 ))
then
    # Stage fastq in..
    if [ ! -f ${FA} ]
    then
	if [ -L ${NM}.fastq.gz ]
    	then
        	rm ${NM}.fastq.gz
	fi
	qsub -q staging -l h_rt=01:00:00 -N ${NM}_stage -o ${LOG}/${NM}_stage.o -e ${LOG}/${NM}_stage.e -cwd -b y -V -r yes -notify '~/mstsrc/qsync.sh ${MASTER}/${FA} ./ ; ln -s ${FA} ${NM}.fastq.gz'
    fi
    # Clip residual 3' P1 (blue) adaptors (with the size selection in the protocol this has a tiny impact but 
    # is potentially useful in some settings).
    qsub -l h_rt=04:00:00 -l h_vmem=8g -N ${NM}_ca -hold_jid ${NM}_stage -o ${LOG}/${NM}_ca.o -e ${LOG}/${NM}_ca.e -cwd -b y -V 'cutadapt --format=fastq --overlap=12 -a ATCACCGACTGCCCATAGAGAGG ${NM}.fastq.gz --minimum-length=30 --output=${NM}.trimmed.fastq --untrimmed-output=${NM}.untrimmed.fastq; sleep 10; cat ${NM}.trimmed.fastq ${NM}.untrimmed.fastq > ${NM}.postTrim.fastq'
 
    # Align
qsub -l h_rt=12:00:00 -pe sharedmem 6 -l h_vmem=16g -hold_jid ${NM}_ca -N ${NM}_bt -o ${LOG}/${NM}_bt.o -e ${LOG}/${NM}_bt.e -cwd -b y -V 'bowtie2 -p 8 -x ${GENOMEPATH}/${GENOME} -U ${NM}.postTrim.fastq -S ${NMBASE}.sam > ${NMBASE}.sam.log'

qsub -l h_rt=12:00:00 -pe sharedmem 6 -l h_vmem=16g -hold_jid ${NM}_ca -N ${NM}_btCirc -o ${LOG}/${NM}_btCirc.e -e ${LOG}/${NM}_btCirc.e -cwd -b y -V 'bowtie2 -p 8 -x ${GENOMEPATH}/${GENOMECIRC} -U ${NM}.postTrim.fastq -S ${NMBASE}.chrMcirc.sam > ${NMBASE}.sam.chrMcirc.log'

fi

if (( ${STARTSTEP} < 3 ))
then
    # Get chrM circularisation junction mapping reads passing q30.
    qsub -l h_rt=05:00:00 -l h_vmem=10g -hold_jid ${NM}_btCirc -N ${NM}_circPipe -o ${LOG}/${NM}_circPipe.o -e ${LOG}/${NM}_circPipe.e -cwd -b y -V 'samtools view -b -q 30 -S ${NMBASE}.chrMcirc.sam > ${NMBASE}.chrMcirc.q30.bam; samtools sort ${NMBASE}.chrMcirc.q30.bam ${NMBASE}.chrMcirc.q30.sorted ; samtools index ${NMBASE}.chrMcirc.q30.sorted.bam ; samtools view -b -q 30 ${NMBASE}.chrMcirc.q30.sorted.bam chrMendJoin > ${NMBASE}.chrMcirc.q30.sorted.chrMjunctionOnly.bam ; bamToBed -i ${NMBASE}.chrMcirc.q30.sorted.chrMjunctionOnly.bam > ${NMBASE}.chrMcirc.q30.bed ; cat ${NMBASE}.chrMcirc.q30.bed | perl ${SRC}/emRiboToRiboBed.pl | perl ${SRC}/emRiboToRiboBedCirc.pl > ${NMBASE}.chrMcirc.ribo.bed'

# Filter on min map quality 30
qsub -l h_rt=05:00:00 -l h_vmem=10g -hold_jid ${NM}_bt -N ${NM}_linPipe -o ${LOG}/${NM}_linPipe.o -e ${LOG}/${NM}_linPipe.e -cwd -b y -V 'samtools view -b -q 30 -S ${NMBASE}.sam > ${NMBASE}.q30.bam ; samtools sort ${NMBASE}.q30.bam ${NMBASE}.q30.sorted ; samtools index ${NMBASE}.q30.sorted.bam ; samtools idxstats ${NMBASE}.q30.sorted.bam > ${NMBASE}.q30.sorted.bam.idxstats ; bamToBed -i ${NMBASE}.q30.sorted.bam | perl ${SRC}/emRiboToRiboBed.pl > ${NMBASE}.ribo.bed '
# combine the standard and circularisation mapping bed files. Order output for sweep algorithms.
fi


if (( ${STARTSTEP} < 4 ))
then
    # Combine both sets of BED file, count signal at 5' end, off-set and flip strands and do coordinate transforms
    # for circularisation, oh and calculate the trinucleotide context for each signal site.
    qsub -l h_rt=04:00:00 -l h_vmem=8g -hold_jid ${NM}_linPipe,${NM}_circPipe -N ${NM}.merge -o ${LOG}/${NM}_merge.o -e ${LOG}/${NM}_merge.e -cwd -b y -V 'perl -w ${SRC}/stripReadFromBed.pl ${NMBASE}.chrMcirc.ribo.bed ${NMBASE}.ribo.bed | sort -k1,1 -k2,2n -k6 | perl -w ${SRC}/cleanMm9ChrMBedCoordinatesForCircularisation.pl 0 > ${NMBASE}.riboMerge.bed ; bedtools groupby -i ${NMBASE}.riboMerge.bed -g 1,3,6 -c 2 -o count > ${NMBASE}.riboCounts.bedish ; cat ${NMBASE}.riboCounts.bedish | perl $SRC/emRiboToRiboCountsBed.pl > ${NMBASE}.riboCounts.bed ; cat ${NMBASE}.riboCounts.bed | perl -w ${SRC}/cleanMm9ChrMBedCoordinatesForCircularisation.pl 1 > ${NMBASE}.riboCounts.triNuc.bed ;bedtools getfasta -fi ${GENOMEFASTA} -bed ${NMBASE}.riboCounts.triNuc.bed -tab -s -name -fo ${NMBASE}.riboCounts.triNuc;  perl ${SRC}/triNucToLoad.pl ${NMBASE}.riboCounts.triNuc > ${NM}.toLoad'
fi

if (( ${STARTSTEP} < 5 ))
then
    # Clean-up and de-stage (if we have the main output file - if not, keep for post-mortem).
    qsub -q staging -l h_rt=01:00:00 -hold_jid ${NM}.merge -N ${NM}_destage -o ${LOG}/${NM}_destage.o -e ${LOG}/${NM}_destage.e -cwd -b y -V -r yes -notify '~/mstsrc/qsync.sh ${NM}.toLoad ${NMBASE}.q30.sorted.bam ${NMBASE}.q30.sorted.bam.bai ${NMBASE}.q30.sorted.bam.idxstats ${MASTER}/${NM}/'
    #qsub -N ${NM}_destageRM -hold_jid ${NM}_destage -o ${LOG}/${NM}_destageRM.o -e ${LOG}/${NM}_destageRM.e -cwd -b y -V 'if [ -n ${NM}.toLoad ]; then rm ../${NM}/${NM}* ../${NM}/*.fastq.gz ; fi'
fi
