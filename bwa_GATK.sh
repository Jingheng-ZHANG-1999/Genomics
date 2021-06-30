#bwa index ${fa_hg38}
#gatk CreateSequenceDictionary -R ${fa_hg38}  -O /zs32/home/kyuan/reference.and.annotations/UCSC.Human.hg38/hg38.dict
#gatk BwaMemIndexImageCreator -I ${fa_hg38}  -O /zs32/home/kyuan/reference.and.annotations/UCSC.Human.hg38/hg38.img
##check order   samtools view -H ERR1039527_RG.bam | grep SQ | cut -f 2 | awk '{ sub(/^SN:/, ""); print;}'

bs="ERR1039527"
fq1="/zs32/home/kyuan/project/ZhangJingheng/fastq/${bs}_1.fastq.gz"
fq2="/zs32/home/kyuan/project/ZhangJingheng/fastq/${bs}_2.fastq.gz"
trim_galore --gzip -o trimmed/ --paired ${fq1} ${fq2}

fa_hg38="/zs32/home/kyuan/reference.and.annotations/UCSC.Human.hg38/hg38.fa"
bwa mem -t 8 -M  ${fa_hg38} trimmed/${bs}_1_val_1.fq.gz trimmed/${bs}_2_val_2.fq.gz | samtools view -Sb - > bam/${bs}.bam
&& echo "** bwa mapping done **"
gatk AddOrReplaceReadGroups -I ERR1039527.bam -O ERR1039527_RG.bam -ID 1 -LB lib1 -PL illumina -PU unit1 -SM gbm --SO coordinate

gatk ReorderSam \
          -I /zs32/home/kyuan/project/ZhangJingheng/sam/SRR6914922_markedup.bam \
          -O /zs32/home/kyuan/project/ZhangJingheng/sam/SRR6914922_reordered.bam \
          -R /zs32/home/kyuan/reference.and.annotations/UCSC.Human.hg38/hg38.fa
gatk SortSam \



###samtools sort -@ 12 -O bam  -o ${bs}_sorted.bam --reference ${fa_hg38} ${bs}.bam

gatk MarkDuplicates -I SRR6914922_sorted.bam -O SRR6914922_markedup.bam -M markdup_metrics.txt

rm SRR6914922.bam
rm
###

###get contig.sam
gatk FindBreakpointEvidenceSpark \
     -I SRR6914922_reordered.bam \
     --aligner-index-image /zs32/home/kyuan/reference.and.annotations/UCSC.Human.hg38/hg38.img \
     --kmers-to-ignore /zs32/home/kyuan/reference.and.annotations/UCSC.Human.hg38/hg38_kmers_to_ignore.txt \
     -O SRR6914922_assemblies.sam

###Structural Variation

 gatk StructuralVariationDiscoveryPipelineSpark \
     -I SRR6914922_reordered.bam \
     -R /zs32/home/kyuan/reference.and.annotations/UCSC.Human.hg38/hg38.2bit \
    --contig-sam-file SRR6914922_contig.sam \
     --aligner-index-image /zs32/home/kyuan/reference.and.annotations/UCSC.Human.hg38/hg38.img \
     --kmers-to-ignore /zs32/home/kyuan/reference.and.annotations/UCSC.Human.hg38/hg38_kmers_to_ignore.txt \
     -O SRR6914922_structural_variants.vcf

 gatk SvDiscoverFromLocalAssemblyContigAlignmentsSpark \
     -I SRR6914922_markedup.bam \
     -R /zs32/home/kyuan/reference.and.annotations/UCSC.Human.hg38/hg38.2bit \
     -O /zs32/home/kyuan/project/ZhangJingheng/vcf

###Variations
gatk HaplotypeCaller -R ${fa_hg38} --emit-ref-confidence GVCF -I SRR6914922_markedup.bam -O SRR6914922_WXS.g.vcf 