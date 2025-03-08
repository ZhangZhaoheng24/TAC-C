#BSUB –J call
#BSUB -n 5
#BSUB -q standardA
#BSUB -o call.out
#BSUB -e call.err

# call peak
macs2 callpeak -t ${sp}_rep1.uni.dedup.sorted2.bam \
    -g $genome_size -q 0.01 -f BAM --nomodel --extsize 150 --shift -75 --call-summits --nolambda \
    --outdir ./${sp}_TAC_rep1_peak -n ${sp}_TAC_rep1 2> ${sp}_TAC-C_rep1.macs2.log
macs2 callpeak -t ${sp}_rep2.uni.dedup.sorted2.bam \
    -g $genome_size -q 0.01 -f BAM --nomodel --extsize 150 --shift -75 --call-summits --nolambda \
    --outdir ./${sp}_TAC_rep2_peak -n ${sp}_TAC_rep2 2> ${sp}_TAC-C_rep2.macs2.log
cat ${sp}_TAC_rep1_peak/${sp}_TAC_rep1_peaks.narrowPeak |awk 'BEGIN{FS=OFS="\t"}{print $1,$2,$3}' |sort -k1,1 -k2,2n |uniq > ${sp}_TAC_rep1_peak/${sp}_1.bed
cat ${sp}_TAC_rep2_peak/${sp}_TAC_rep2_peaks.narrowPeak |awk 'BEGIN{FS=OFS="\t"}{print $1,$2,$3}' |sort -k1,1 -k2,2n |uniq > ${sp}_TAC_rep2_peak/${sp}_2.bed
bedtools intersect -a ${sp}_TAC_rep1_peak/${sp}_1.bed -b ${sp}_TAC_rep2_peak/${sp}_2.bed -wa -wb |awk 'BEGIN{FS=OFS="\t"}{print $1,$2,$3"\n"$4,$5,$6}' |sort -k1,1 -k2,2n |uniq|bedtools merge -i - -d 1 > ${sp}_TAC_peak_merge.bed
#call loop

python digest_genome.py -r dpnii -o ${sp}_DpnII.bed ${sp}.genome.fa

hichipper -o ./${sp}_TAC_rep1_loop \
    -ii ${sp}_rep1.all.validPairs2 \
    -rf ${sp}_DpnII.bed -p ${sp}_TAC_rep1_peak/${sp}_TAC_rep1_peaks.narrowPeak \
    --skip-diffloop --basic-qc call
hichipper -o ./${sp}_TAC_rep2_loop \
    -ii validPairs/${sp}_rep2.all.validPairs2 \
    -rf 00.Genome/${sp}_DpnII.bed -p ${sp}_TAC_rep2_peak/${sp}_TAC_rep2_peaks.narrowPeak \
    --skip-diffloop --basic-qc call 
cat ./${sp}_TAC_rep1_loop/one.filt.intra.loop_counts.bedpe |awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$4,$5,$6,$8}' |sort -k1,1 -k2,2n |uniq > ./${sp}_TAC_rep1_loop/${sp}_TAC_rep1.bedgraph
cat ./${sp}_TAC_rep2_loop/one.filt.intra.loop_counts.bedpe |awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$4,$5,$6,$8}' |sort -k1,1 -k2,2n |uniq > ./${sp}_TAC_rep2_loop/${sp}_TAC_rep2.bedgraph
hicMergeLoops \
    -i ./${sp}_TAC_rep1_loop/${sp}_TAC_rep1.bedgraph ./${sp}_TAC_rep2_loop/${sp}_TAC_rep2.bedgraph \
    -o ${sp}_10kb_merged_loop.bedpe -r 10000