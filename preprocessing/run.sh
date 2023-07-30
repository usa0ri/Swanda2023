# !/bin/bash

exp_metadata=/home/data/exp_metadata.csv

while IFS=, read -a col
do
    smpl=${col[2]}
    readfile=${col[3]}
    org=${col[10]}

    if [ "${smpl:0:3}" != "exp" ]; then
        continue
    fi

    echo $smpl
    echo $readfile
    echo $org

    cutadapt -j 15 \
    	-g "^GGG" -a "A{10}" -n 2\
        -m 15 \
    	-o /home/result/trimmed_${readfile}.fastq.gz \
    	/home/result/trimmed5_umi_${readfile}.fastq.gz > log_cutadapt_trimmed_${readfile}.txt

    if [ "$org" = "Mus musculus" ]; then
        ref_rRNA_path=/home/ref/Mus_musculus_yuanhui/rRNA/rRNA_merged_RNA_central_ENS
    elif [ "$org" = "Homo sapiens" ]; then
        ref_rRNA_path=/home/ref/Homo_sapiens_yuanhui/ref_rRNA/rRNA_NCBI_ENS_merged
    fi
    # echo $ref_rRNA_path

    bowtie --quiet \
		-q \
		-v 0 \
		--norc \
		-p 12 \
		-S \
		--sam-nohead \
		--un /home/result/rRNA_left_${readfile}.fastq \
		-q $ref_rRNA_path \
		<(zcat /home/result/trimmed_${readfile}.fastq.gz) \
		| awk 'BEGIN{FS="\t"}{if($2==0){print}}' \
		>> /home/result/rRNA_align_${readfile}.fastq
    gzip /home/result/rRNA_left_${readfile}.fastq
    gzip /home/result/rRNA_align_${readfile}.fastq

    if [ "$org" = "Mus musculus" ]; then
        ref_STAR_path=/home/ref/Mus_musculus_109_saori/STAR_ref
    elif [ "$org" = "Homo sapiens" ]; then
        ref_STAR_path=/home/ref/Homo_sapiens_109_saori/STAR_ref
    fi

    STAR --runThreadN 12 \
        --genomeDir $ref_STAR_path \
        --readFilesIn /home/result/rRNA_left_${readfile}.fastq.gz \
        --outSAMtype BAM SortedByCoordinate  \
        --readFilesCommand zcat \
        --runDirPerm All_RWX \
        --outFileNamePrefix /home/result/STAR_align_$readfile\_ \
        --outSAMattributes All \
        --outFilterScoreMinOverLread 0 \
        --outFilterMatchNminOverLread 0 \
        --outBAMsortingBinsN 200

    file=/home/result/dedup_STAR_align_${readfile}.bam
    samtools view -H $file > /home/result/header
    samtools view $file \
        | grep -P "^\S+\s0\s" \
        | grep -P "NH:i:1\b" \
        | grep -E -w 'NM:i:0|NM:i:1|NM:i:2' \
        | cat /home/result/header -| samtools view -bS ->/home/result/uniq_STAR_align_${readfile}.bam
    samtools index /home/result/uniq_STAR_align_${readfile}.bam

    # rm -rf /home/result/rRNA_left_${readfile}.fastq.gz
    # rm -rf /home/result/rRNA_align_${readfile}.fastq.gz
    # rm -rf /home/result/trimmed_${readfile}.fastq.gz

done < $exp_metadata