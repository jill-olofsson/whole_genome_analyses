## Step 0: prepare for alignment 
	# define location reference
REF=/home/$USER/allo_refgenome
	# define location to data
data=/home/$USER/cleaned_reads

mkdir map
## Step 1: Map in Bowtie2
cat samples.txt | while read line ; do bowtie2 -x ${REF}/ASEM_C4_v1.0 --no-unal -1  ${data}/"$line".1.fq.gz -2 ${data}/"$line".2.fq.gz -p 4 -S map/"$line"_bowtie_ASEM_C4_v1.sam ; done

## Step 2: convert sam to bam
cat samples.txt | while read line ; do samtools view -S -b -h map/"$line"_bowtie_ASEM_C4_v1.sam > map/"$line"_bowtie_ASEM_C4_v1.bam ; done

## Step 3: only keeping uniquely aligned reads in proper pairs
mkdir map/1_proper_pairs
cat samples.txt | while read line ; do samtools view -h -bq -f2 map/"$line"_bowtie_ASEM_C4_v1.bam > map/1_proper_pairs/"$line"_bowtie_ASEM_C4_v1_unique_f2.bam ; done

## Step 4: Sort and index bam-file
cat samples.txt | while read line ;  do samtools sort map/1_proper_pairs/"$line"_bowtie_ASEM_C4_v1_unique_f2.bam > map/1_proper_pairs/"$line"_bowtie_ASEM_C4_v1_unique_f2_sort.bam ; done
	# index
cat samples.txt | while read line ; do samtools index map/1_proper_pairs/"$line"_bowtie_ASEM_C4_v1_unique_f2_sort.bam ; done

## Step 5: call SNPs with mpileup and bcftools
mkdir map/2_SNPs
cat samples.txt | while read line ; do samtools mpileup -q 20 -B -uvIf ${REF}/ASEM_C4_v1.0.fasta map/1_proper_pairs/"$line"_bowtie_ASEM_C4_v1_unique_f2_sort.bam | bcftools call -cv -O v - | bcftools filter -i 'DP>=2 & DP<=50 & QUAL>=20' - > map/2_SNPs/"$line"_bowtie_ASEM_C4_v1_unique_f2_sort_covfilt_allsite.vcf ; done

