# Code of denovo, QC and blast

## 1. Trinity denovo assemble genome

~~~R
#merge straight
vim roly_poly_strait.txt

strait strait1 strait-1_FRAS192303306-1r_1.clean.fq.gz strait-1_FRAS192303306-1r_2.clean.fq.gz
strait strait2 strait-2_FRAS192303307-1r_1.clean.fq.gz strait-2_FRAS192303307-1r_2.clean.fq.gz
strait strait4 strait-3_FRAS192303308-1r_1.clean.fq.gz strait-3_FRAS192303308-1r_2.clean.fq.gz

Trinity --seqType fq --max_memory 50G  \
        --left strait-1_FRAS192303306-1r_1.clean.fq.gz,strait-2_FRAS192303307-1r_1.clean.fq.gz,strait-3_FRAS192303308-1r_1.clean.fq.gz \
        --right strait-1_FRAS192303306-1r_2.clean.fq.gz,strait-2_FRAS192303307-1r_2.clean.fq.gz,strait-3_FRAS192303308-1r_2.clean.fq.gz \
        --output trinity_strait_all_merge \
        --samples_file roly_poly_strait.txt \
        --CPU 10

#merge curly
vim roly_poly_curl.txt

curl curl_1 curl-1_FRAS192303309-1b_1.clean.fq.gz curl-1_FRAS192303309-1b_2.clean.fq.gz
curl curl_3 curl-2_FRAS192303310-1r_1.clean.fq.gz curl-2_FRAS192303310-1r_2.clean.fq.gz
curl curl_4 curl-3_FRAS192303311-1r_1.clean.fq.gz curl-3_FRAS192303311-1r_2.clean.fq.gz

Trinity --seqType fq --max_memory 50G  \
        --left curl-1_FRAS192303309-1b_1.clean.fq.gz,curl-2_FRAS192303310-1r_1.clean.fq.gz,curl-3_FRAS192303311-1r_1.clean.fq.gz \
        --right curl-1_FRAS192303309-1b_2.clean.fq.gz,curl-2_FRAS192303310-1r_2.clean.fq.gz,curl-3_FRAS192303311-1r_2.clean.fq.gz \
        --output trinity_curl_all_merge \
        --samples_file roly_poly_curl.txt \
        --CPU 20

~~~

## 2. Summary and re-alignment

~~~R
/home/zhaolei/miniconda3/envs/trinity/bin/TrinityStats.pl /mnt/data/userdata/zhaolei/project/denovo_roly/trinity_strait_all_merge/Trinity.fasta > Trinitystats.log

/home/zhaolei/miniconda3/envs/trinity/bin/TrinityStats.pl /mnt/data/userdata/zhaolei/project/denovo_roly/trinity_curl_all_merge/Trinity.fasta > Trinitystats.log

bowtie2-build Trinity.fasta Trinity_fasta_index

bowtie2 -p 20 -x /mnt/data/userdata/zhaolei/project/denovo_roly/trinity_strait_all_merge/Trinity_fasta_index \
-1 strait-1_FRAS192303306-1r_1.clean.fq.gz -2 strait-1_FRAS192303306-1r_2.clean.fq.gz | samtools sort -O bam -@ 20 -o - > strait1.bam

bowtie2 -p 20 -x /mnt/data/userdata/zhaolei/project/denovo_roly/trinity_strait_all_merge/Trinity_fasta_index \
-1 strait-2_FRAS192303307-1r_1.clean.fq.gz -2 strait-2_FRAS192303307-1r_2.clean.fq.gz | samtools sort -O bam -@ 20 -o - > strait2.bam

bowtie2 -p 20 -x /mnt/data/userdata/zhaolei/project/denovo_roly/trinity_strait_all_merge/Trinity_fasta_index \
-1 strait-3_FRAS192303308-1r_1.clean.fq.gz -2 strait-3_FRAS192303308-1r_2.clean.fq.gz | samtools sort -O bam -@ 20 -o - > strait3.bam


bowtie2 -p 20 -x /mnt/data/userdata/zhaolei/project/denovo_roly/trinity_curl_all_merge/Trinity_fasta_index \
-1 curl-1_FRAS192303309-1b_1.clean.fq.gz -2 curl-1_FRAS192303309-1b_2.clean.fq.gz | samtools sort -O bam -@ 20 -o - > curl1.bam

bowtie2 -p 20 -x /mnt/data/userdata/zhaolei/project/denovo_roly/trinity_curl_all_merge/Trinity_fasta_index \
-1 curl-2_FRAS192303311-1r_1.clean.fq.gz -2 curl-2_FRAS192303311-1r_2.clean.fq.gz | samtools sort -O bam -@ 20 -o - > curl2.bam

bowtie2 -p 20 -x /mnt/data/userdata/zhaolei/project/denovo_roly/trinity_curl_all_merge/Trinity_fasta_index \
-1 curl-3_FRAS192303310-1r_1.clean.fq.gz -2 curl-3_FRAS192303310-1r_2.clean.fq.gz | samtools sort -O bam -@ 20 -o - > curl3.bam
~~~

## 3. Homology quantification between species

~~~R
nucmer --prefix Drosophila_strait \
/mnt/data/userdata/zhaolei/program/refer/Drosophila/Drosophila_melanogaster/UCSC/dm6/Sequence/WholeGenomeFasta/genome.fa \
/mnt/data/userdata/zhaolei/project/denovo_roly/trinity_strait_all_merge/Trinity.fasta

nucmer --prefix Drosophila_curl \
/mnt/data/userdata/zhaolei/program/refer/Drosophila/Drosophila_melanogaster/UCSC/dm6/Sequence/WholeGenomeFasta/genome.fa \
/mnt/data/userdata/zhaolei/project/denovo_roly/trinity_curl_all_merge/Trinity.fasta

dnadiff -d Drosophila_strait.delta
dnadiff -d Drosophila_curl.delta

show-coords -r Drosophila_strait.delta > Drosophila_strait_anno.coord
show-coords -r Drosophila_curl.delta > Drosophila_curl_anno.coord


nucmer --prefix mouse_mm10_strait \
/mnt/data/public_data/reference/Mus_musculus_UCSC/UCSC/mm10/Sequence/WholeGenomeFasta/genome.fa \
/mnt/data/userdata/zhaolei/project/denovo_roly/trinity_strait_all_merge/Trinity.fasta

nucmer --prefix mouse_mm10_curl \
/mnt/data/public_data/reference/Mus_musculus_UCSC/UCSC/mm10/Sequence/WholeGenomeFasta/genome.fa \
/mnt/data/userdata/zhaolei/project/denovo_roly/trinity_curl_all_merge/Trinity.fasta

dnadiff -d mouse_mm10_strait.delta
dnadiff -d mouse_mm10_curl.delta

show-coords -r mouse_mm10_strait.delta > mouse_mm10_strait_anno.coord
show-coords -r mouse_mm10_curl.delta > mouse_mm10_curl_anno.coord


nucmer --prefix human_hg19_strait \
/mnt/data/sequencedata/hg19/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa \
/mnt/data/userdata/zhaolei/project/denovo_roly/trinity_strait_all_merge/Trinity.fasta

nucmer --prefix human_hg19_curl \
/mnt/data/sequencedata/hg19/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa \
/mnt/data/userdata/zhaolei/project/denovo_roly/trinity_curl_all_merge/Trinity.fasta

dnadiff -d human_hg19_strait.delta
dnadiff -d human_hg19_curl.delta

show-coords -r human_hg19_strait.delta > human_hg19_strait_anno.coord
show-coords -r human_hg19_curl.delta > human_hg19_curl_anno.coord
~~~

