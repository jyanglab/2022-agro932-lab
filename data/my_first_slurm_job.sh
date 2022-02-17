#!/bin/bash -l
#SBATCH -D /common/jyanglab/jyang21/courses/2022-agro932-lab
#SBATCH -o /common/jyanglab/jyang21/courses/2022-agro932-lab/slurm-log/steve-stdout-%j.txt
#SBATCH -e /common/jyanglab/jyang21/courses/2022-agro932-lab/slurm-log/steve-stderr-%j.txt
#SBATCH -J theta
#SBATCH -t 1:00:00
#SBATCH --mail-user=your_email_address@gmail.com
#SBATCH --mail-type=END #email if ends
#SBATCH --mail-type=FAIL #email if fails

set -e
set -u
# insert your script here

module load bwa samtools
mkdir largedata/lab5/
cp data/Zea_mays.B73_RefGen_v4.dna.chromosome.Mt.fa largedata/lab5

# simulate 20 individuals
cd largedata/lab5
for i in {1..20}
do
   wgsim Zea_mays.B73_RefGen_v4.dna.chromosome.Mt.fa -e 0 -d 500 -N 50000 -1 100 -2 100 -r 0.1  -R 0 -X 0 l$i.read1.fq l$i.read2.fq
done

# alignment
module load bwa samtools bcftools
# index the reference genome
bwa index Zea_mays.B73_RefGen_v4.dna.chromosome.Mt.fa

# using bwa mem to align the reads to the reference genome 
for i in {1..20}; do bwa mem Zea_mays.B73_RefGen_v4.dna.chromosome.Mt.fa l$i.read1.fq l$i.read2.fq | samtools view -bSh - > l$i.bam; done
# sort
for i in *.bam; do samtools sort $i -o sorted_$i; done
# index them
for i in sorted*.bam; do samtools index $i; done


### index the genome assembly
samtools faidx Zea_mays.B73_RefGen_v4.dna.chromosome.Mt.fa
### Run `mpileup` to generate VCF format
ls sorted_l*bam > bamlist.txt
samtools mpileup -g -f Zea_mays.B73_RefGen_v4.dna.chromosome.Mt.fa -b bamlist.txt > myraw.bcf
bcftools call myraw.bcf -cv -Ob -o snps.bcf

### Extract allele frequency at each position
bcftools query -f '%CHROM %POS %AF1\n' snps.bcf > frq.txt
bcftools query -f '%CHROM %POS %REF %ALT [\t%GT]\n' snps.bcf > geno.txt

