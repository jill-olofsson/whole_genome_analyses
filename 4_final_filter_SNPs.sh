#Step 0: prep
	# define location reference
REF=/home/$USER/allo_refgenome
	# define location to data
data=/home/$USER/cleaned_reads
	# definle location of file with individuals to filter away
remove=/home/$USER/path_to_file

# Step 1: filter away samples with > 99% missing data
vcftools --gzvcf merged_plate01-12_ASEM_C4_v1_unique_f2_covfilt_allsites_TEmask_chr.recode.vcf.gz --remove $remove/rm_99.txt --recode --out merged_plate01-12_ASEM_C4_v1_unique_f2_covfilt_allsites_TEmask_chr_max99
bgzip merged_plate01-12_ASEM_C4_v1_unique_f2_covfilt_allsites_TEmask_chr_max99.recode.vcf
bcftools index merged_plate01-12_ASEM_C4_v1_unique_f2_covfilt_allsites_TEmask_chr_max99.recode.vcf.gz

# Step 2: only keep SNPs
vcftools --gzvcf merged_plate01-12_ASEM_C4_v1_unique_f2_covfilt_allsites_TEmask_chr_max99.recode.vcf.gz --mac 1 --recode --out merged_plate01-12_ASEM_C4_v1_unique_f2_covfilt_allsites_TEmask_chr_max99_SNP
bgzip merged_plate01-12_ASEM_C4_v1_unique_f2_covfilt_allsites_TEmask_chr_max99_SNP.recode.vcf
bcftools index merged_plate01-12_ASEM_C4_v1_unique_f2_covfilt_allsites_TEmask_chr_max99_SNP.recode.vcf.gz

# Step 3: filter away sites with more than 90 % missing data and those with minor allele freq <0.05 (mac 56)
vcftools --gzvcf merged_plate01-12_ASEM_C4_v1_unique_f2_covfilt_allsites_TEmask_chr_max99_SNP.recode.vcf.gz --mac 56 --max-missing 0.1 --recode --out merged_plate01-12_ASEM_C4_v1_unique_f2_covfilt_allsites_TEmask_chr_max99_SNP_mac56_missing90
bgzip merged_plate01-12_ASEM_C4_v1_unique_f2_covfilt_allsites_TEmask_chr_max99_SNP_mac56_missing90.recode.vcf
bcftools index merged_plate01-12_ASEM_C4_v1_unique_f2_covfilt_allsites_TEmask_chr_max99_SNP_mac56_missing90.recode.vcf.gz
	# resulting file provided in files/vcf

# Step 4: thin the resulting file to only keep SNPs at least 1000 bp apart
vcftools --gzvcf merged_plate01-12_ASEM_C4_v1_unique_f2_covfilt_allsites_TEmask_chr_max99_SNP_mac56_missing90.recode.vcf.gz --thin 1000 --recode --out merged_plate01-12_ASEM_C4_v1_unique_f2_covfilt_allsites_TEmask_chr_max99_SNP_mac56_missing90_thin1k
	# resulting file provided in files/vcf

