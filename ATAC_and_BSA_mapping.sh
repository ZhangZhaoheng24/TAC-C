#ATAC BSA
set -x

mkdir $name
cd $name

fastp -i $rawdata_dir/${name}_1.fq.gz \
    -I $rawdata_dir/${name}_2.fq.gz \
    -o ./${name}.clean.R1.fq.gz -O ./${name}.clean.R2.fq.gz \
    --detect_adapter_for_pe -w 8 --compression 9 \
    -h ./${name}.html -j ./${name}.json
#R1
bwa mem -t 8 -a /public-supool/home/zhangzh/proj/TAC/data/00.Genome/bwa/${genome} \
    ./${name}.clean.R1.fq.gz |\
     samtools view -@ 8 -bS -F 3328 - > ${name}.R1.bam
samtools sort -@ 8 -m 2500M -n -T /public-supool/home/zhangzh/tmp/${name}.R1 \
    -o ${name}.R1.sorted.bam \
    ${name}.R1.bam
#R2
bwa mem -t 8 -a /public-supool/home/zhangzh/proj/TAC/data/00.Genome/bwa/${genome} \
    ./${name}.clean.R2.fq.gz |\
     samtools view -@ 8 -bS -F 3328 - > ${name}.R2.bam
samtools sort -@ 8 -m 2500M -n -T /public-supool/home/zhangzh/tmp/${name}.R2 \
    -o ${name}.R2.sorted.bam \
    ${name}.R2.bam
#merge
python /public-supool/home/xmliu/software/hicpro/HiC-Pro_3.1.0/scripts/mergeSAM.py \
    -q 5 -t -v -f ${name}.R1.sorted.bam -r ${name}.R2.sorted.bam -o ${name}.merge.bam
java -Djava.io.tmpdir=/public-supool/home/zhangzh/tmp -XX:ParallelGCThreads=8 -Xmx40g -jar \
    /public-supool/home/xmliu/software/picard/picard.jar \
    MarkDuplicates -I ${name}.merge.bam -O ${name}.merge.picard.bam \
    -M ${name}.merge.rmdup.metrics.txt -REMOVE_DUPLICATES true

samtools sort -@ 8 ${name}.merge.picard.bam -o ${name}.merge.picard.sort.bam
samtools index -@ 8 -c ${name}.merge.picard.sort.bam
