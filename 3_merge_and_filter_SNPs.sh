#Step 0: prep
	# define location reference
REF=/home/$USER/allo_refgenome
	# define location to data
data=/home/$USER/cleaned_reads

##Step 1: merge all the vcf files
	# compress vcf files and index them
cat samples.txt | while read line ; do bgzip map/2_SNPs/"$line"_bowtie_ASEM_C4_v1_unique_f2_sort_covfilt_allsite.vcf ; done
cat samples.txt | while read line ; do bcftools index map/2_SNPs/"$line"_bowtie_ASEM_C4_v1_unique_f2_sort_covfilt_allsite.vcf.gz ; done
	# make a list of all the compressed vcf files
ls map/2_SNPs/*_bowtie_ASEM_C4_v1_unique_f2_sort_covfilt_allsite.vcf.gz > vcf_files.txt
	# merge the files and zip the resulting vcf file to save storage space
bcftools merge -m id -l vcf_files.txt > merged_ASEM_C4_v1_unique_f2_covfilt_allsites.vcf
bgzip merged_plate01-12_ASEM_C4_v1_unique_f2_covfilt_allsites.vcf
bcftools index merged_plate01-12_ASEM_C4_v1_unique_f2_covfilt_allsites.vcf

# Step 2: mask TEs
bedtools intersect -a merged_plate01-12_ASEM_C4_v1_unique_f2_covfilt_allsites.vcf.gz -b ${REF}/ASEM_C4_v1.0.fasta.out.TE.gff -v -header -wa > merged_plate01-12_ASEM_C4_v1_unique_f2_covfilt_allsites_TEmask.vcf
bgzip merged_plate01-12_ASEM_C4_v1_unique_f2_covfilt_allsites_TEmask.vcf
bcftools index merged_plate01-12_ASEM_C4_v1_unique_f2_covfilt_allsites_TEmask.vcf.gz

#Step 3: only keep chr's
vcftools --gzvcf merged_plate01-12_ASEM_C4_v1_unique_f2_covfilt_allsites_TEmask.vcf.gz --chr Chromosome_1 --chr Chromosome_2 --chr Chromosome_3 --chr Chromosome_4 --chr Chromosome_5 --chr Chromosome_6 --chr Chromosome_7 --chr Chromosome_8 --chr Chromosome_9 --recode --out merged_plate01-12_ASEM_C4_v1_unique_f2_covfilt_allsites_TEmask_chr
bgzip merged_plate01-12_ASEM_C4_v1_unique_f2_covfilt_allsites_TEmask_chr.recode.vcf
bcftools index merged_plate01-12_ASEM_C4_v1_unique_f2_covfilt_allsites_TEmask_chr.recode.vcf.gz

# Step 4: calculate missingness per individual
vcftools --gzvcf merged_plate01-12_ASEM_C4_v1_unique_f2_covfilt_allsites_TEmask_chr.recode.vcf.gz --missing-indv --out merged_plate01-12_ASEM_C4_v1_unique_f2_covfilt_allsites_TEmask_chr
