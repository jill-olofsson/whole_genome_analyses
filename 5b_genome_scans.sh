# This is the genome scan analyses

	# for the genome scan only individuals assigned to nuclear clade II and IIIa were used

# Step 0: define the locations of all programs and files needed
files=/home/$USER/files/genome_scans

# Fst scans
	# Step 1: Filter the vcf file only keeping individuals assigned to nuclear clade II and IIIa
		# setting the mac filter to 21 gives an approximate maf of 0.05
vcftools --gzvcf merged_plate01-12_ASEM_C4_v1_unique_f2_covfilt_allsites_TEmask_chr_max99_SNP.recode.vcf.gz --max-missing 0.1 --mac 21 --weir-fst-pop ${files}/ntII_samples.txt --weir-fst-pop ${files}/ntIIIa_samples.txt --fst-window-size 500000 --fst-window-step 10000 --out ntII_ntIIIa_500kb10kb_90missing_mac21
			# resulting vcf file provided in files/geneome_scan

	# Step 2: Calculate null distributions 
		# write an empty file to store all the distributions in
touch ntIIntIIIa_fst_null_dist.txt
		# samle var number of variants 1000 times
		# calc fst on a site by site bases
			#remove the sampled file
		# turn negative values to 0
			# remove the fst file
		# calculate the mean of the values in the file
			# remove the fixed fst file
		# add the distribution of mean to the file to store in 
		# remove the file with the mean fst distribution
		# get the header
grep "^#" ntII_ntIIIa_500kb10kb_90missing_mac21.recode.vcf > ntIIntIIIa_header

for var in $(seq 1 198);
  do
  echo $var > ntIIntIIIa_mean_fst_dist_1.txt
  for rep in $(seq 1 1000) ;
    do grep "^Chromosome" ntII_ntIIIa_500kb10kb_90missing_mac21.recode.vcf | shuf -n $var | cat ntIIntIIIa_header - > ntIIntIIIa_1.vcf
       vcftools --vcf ntIIntIIIa_1.vcf --weir-fst-pop ${files}/ntII_samples.txt --weir-fst-pop ${files}/ntIIIa_samples.txt --out ntIIntIIIa_1
       rm ntIIntIIIa_1.vcf
       sed 's/-.*/0/' ntIIntIIIa_1.weir.fst > ntIIntIIIa_fixed_1.weir.fst
       rm ntIIntIIIa_1.weir.fst
       awk '{ total += $3 } END { print sprintf("%.9f",total/NR) }' ntIIntIIIa_fixed_1.weir.fst >> ntIIntIIIa_mean_fst_dist_1.txt
       sed 's/\t/\n/g' ntIIntIIIa_mean_fst_dist_1.txt > temp && mv temp ntIIntIIIa_mean_fst_dist_1.txt
       rm ntIIntIIIa_fixed_1.weir.fst;
  done
  paste ntIIntIIIa_fst_null_dist.txt ntIIntIIIa_mean_fst_dist_1.txt > temp && mv temp ntIIntIIIa_fst_null_dist.txt
  rm ntIIntIIIa_mean_fst_dist_1.txt;
done

	# Step 3: Find Fst outliers and plot
		# this is done in python using script: 500kb10kb_allchr_90missing_mac21_ntII-ntIIIa_20191203.ipynb
		# all necessary files are provided in the the folder genome_scans/Fst

	# Step 4: make the bean plots for comparisions of background to the DE genes
		# calculate the needed Fst using vcftools
vcftools --gzvcf merged_plate01-12_ASEM_C4_v1_unique_f2_covfilt_allsites_TEmask_chr_max99_SNP.recode.vcf.gz --max-missing 0.1 --mac 21 --weir-fst-pop ${files}/ntII_samples.txt --weir-fst-pop ${files}/ntIIIa_samples.txt --fst-window-size 100000 --fst-window-step 10000 --out ntII_ntIIIa_100kb10kb_90missing_mac21
		# take out every 10 line for each chr to make over-lapping windows
cat ${files}/chr.txt | while read line ; do grep "$line" ntII_ntIIIa_100kb10kb_90missing_mac21.windowed.weir.fst | sed -n '0~10p' >> 100kb_windows2_fst.txt ; done 
		# Plot the results in R
R
library(beanplot)
			# import data
tab<-read.table("100kb_windows2_fst.txt", header=F, dec=".", sep="\t")
colnames(tab)<-c("CHROM","BIN_START","BIN_END", "N_VARIANTS","WEIGHTED_FST","MEAN_FST")
			# turn all negative values to 0
tab$MEAN_FST <- ifelse(tab$MEAN_FST < 0, 0, tab$MEAN_FST)
			# define the values for the windows including the DE genes
genes<-data.frame(genes=c("pck","ppdk","asp-at"),fst=c(0.563888,0.281476,0.738131))
			# plot
beanplot(tab$MEAN_FST, col="grey80", beanlines="median", beanlinewd=0.5, horizontal = TRUE, xlab="Mean FST")
for(i in 1:length(genes[,1])){
  abline(v=genes[i,2], col="red")
  text(x=genes[i,2], y=0.75, as.factor(genes[i,1]))
}
abline(v=quantile(tab$MEAN_FST, probs = c(0.05,0.95)), col="blue")
abline(v=quantile(tab$MEAN_FST, probs = c(0.01,0.98)), col="green")

###########

# dxy scans
	# Step 1: calculate the allele frequency at each site in each of the two nuclear clades using vcftools
vcftools --gzvcf merged_plate01-12_ASEM_C4_v1_unique_f2_covfilt_allsites_TEmask_chr_max99.recode.vcf.gz --max-missing 0.1 --keep ${files}/ntII_samples.txt --freq --out ntII_covfilt_allsites_TEmask_chr_max99_90missing
		#
awk '{OFS="\t"}NF==5' ntII_covfilt_allsites_TEmask_chr_max99_90missing.frq > ntII_covfilt_allsites_TEmask_chr_max99_90missing_5col.frq
awk '{OFS="\t"}NF==6' ntII_covfilt_allsites_TEmask_chr_max99_90missing.frq > ntII_covfilt_allsites_TEmask_chr_max99_90missing_6col.frq
awk '{OFS="\t"}NF==7' ntII_covfilt_allsites_TEmask_chr_max99_90missing.frq > ntII_covfilt_allsites_TEmask_chr_max99_90missing_7col.frq
awk '{OFS="\t"}NF==8' ntII_covfilt_allsites_TEmask_chr_max99_90missing.frq > ntII_covfilt_allsites_TEmask_chr_max99_90missing_8col.frq
		#
vcftools --gzvcf merged_plate01-12_ASEM_C4_v1_unique_f2_covfilt_allsites_TEmask_chr_max99.recode.vcf.gz --max-missing 0.1 --keep ${files}/ntIIIa_samples.txt --freq --out ntIIIa_covfilt_allsites_TEmask_chr_max99_90missing
		#
awk '{OFS="\t"}NF==5' ntIIIa_covfilt_allsites_TEmask_chr_max99_90missing.frq > ntIIIa_covfilt_allsites_TEmask_chr_max99_90missing_5col.frq
awk '{OFS="\t"}NF==6' ntIIIa_covfilt_allsites_TEmask_chr_max99_90missing.frq > ntIIIa_covfilt_allsites_TEmask_chr_max99_90missing_6col.frq
awk '{OFS="\t"}NF==7' ntIIIa_covfilt_allsites_TEmask_chr_max99_90missing.frq > ntIIIa_covfilt_allsites_TEmask_chr_max99_90missing_7col.frq
awk '{OFS="\t"}NF==8' ntIIIa_covfilt_allsites_TEmask_chr_max99_90missing.frq > ntIIIa_covfilt_allsites_TEmask_chr_max99_90missing_8col.frq

	# Step 2: Filter the frequency files and merge them 
		# this is done in python using script: dxy_filter_and_merge.ipynb

	# Step 3: calculate the dxy
		# this is done in python using script: dxy_calc_ntIIvsntIIIa.ipynb

	# Step 4: Calculate null distributions
touch ntII-ntIIIa_dxy_null.txt
for var in $(seq 1 100);
  do
  echo $var > ntII-ntIIIa.txt
  grep "^Chromosome" ntII-ntIIIa_filtered_dxy.frq | awk '{print $NF}' | shuf > ntII-ntIIIa_null.txt
  cat ntII-ntIIIa.txt ntII-ntIIIa_null.txt > temp && mv temp ntII-ntIIIa.txt
  rm ntII-ntIIIa_null.txt
  paste ntII-ntIIIa_dxy_null.txt ntII-ntIIIa.txt > temp && mv temp ntII-ntIIIa_dxy_null.txt
  rm ntII-ntIIIa.txt;
done
cut -f2- ntII-ntIIIa_dxy_null.txt > temp && mv temp ntII-ntIIIa_dxy_null.txt
awk -F'\t' '{print NF; exit}' ntII-ntIIIa_dxy_null.txt		#
awk 'BEGIN {OFS="\t"}; {print $1, $2}' ntII-ntIIIa_filtered_dxy.frq | paste - ntII-ntIIIa_dxy_null.txt > ntII-ntIIIa_dxy_null_comb.txt
		# in R
R
tab<-read.table("ntII-ntIIIa_dxy_null_comb.txt", header =T, sep="\t", dec=".")

chr1<-tab[tab[,1]=="Chromosome_1",]
chr1_slide<-seq(1, 76910000, by=10000)
chr1_dxy_slide<-data.frame()
for(i in 1:length(chr1_slide)){
  x<-chr1[chr1[,2]>chr1_slide[i] & chr1[,2]<(chr1_slide[i]+500000),]
  for(n in 1:50){
    if(length(x[,1])==0){
      chr1_dxy_slide[i,n]<-"nan"
    } else {
      chr1_dxy_slide[i,n]<-mean(x[,n+2])
    }
  }
}
chr2<-tab[tab[,1]=="Chromosome_2",]
chr2_slide<-seq(1, 88290000, by=10000)
chr2_dxy_slide<-data.frame()
for(i in 1:length(chr2_slide)){
  x<-chr2[chr2[,2]>chr2_slide[i] & chr2[,2]<(chr2_slide[i]+500000),]
  for(n in 1:50){
    if(length(x[,1])==0){
      chr2_dxy_slide[i,n]<-"nan"
    } else {
      chr2_dxy_slide[i,n]<-mean(x[,n+2])
    }
  }
}
chr3<-tab[tab[,1]=="Chromosome_3",]
chr3_slide<-seq(1, 89510000, by=10000)
chr3_dxy_slide<-data.frame()
for(i in 1:length(chr3_slide)){
  x<-chr3[chr3[,2]>chr3_slide[i] & chr3[,2]<(chr3_slide[i]+500000),]
  for(n in 1:50){
    if(length(x[,1])==0){
      chr3_dxy_slide[i,n]<-"nan"
    } else {
      chr3_dxy_slide[i,n]<-mean(x[,n+2])
    }
  }
}
chr4<-tab[tab[,1]=="Chromosome_4",]
chr4_slide<-seq(1, 67590000, by=10000)
chr4_dxy_slide<-data.frame()
for(i in 1:length(chr4_slide)){
  x<-chr4[chr4[,2]>chr4_slide[i] & chr4[,2]<(chr4_slide[i]+500000),]
  for(n in 1:50){
    if(length(x[,1])==0){
      chr4_dxy_slide[i,n]<-"nan"
    } else {
      chr4_dxy_slide[i,n]<-mean(x[,n+2])
    }
  }
}
chr5<-tab[tab[,1]=="Chromosome_5",]
chr5_slide<-seq(1, 88020000, by=10000)
chr5_dxy_slide<-data.frame()
for(i in 1:length(chr5_slide)){
  x<-chr5[chr5[,2]>chr5_slide[i] & chr5[,2]<(chr5_slide[i]+500000),]
  for(n in 1:50){
    if(length(x[,1])==0){
      chr5_dxy_slide[i,n]<-"nan"
    } else {
      chr5_dxy_slide[i,n]<-mean(x[,n+2])
    }
  }
}
chr6<-tab[tab[,1]=="Chromosome_6",]
chr6_slide<-seq(1, 69450000, by=10000)
chr6_dxy_slide<-data.frame()
for(i in 1:length(chr6_slide)){
  x<-chr6[chr6[,2]>chr6_slide[i] & chr6[,2]<(chr6_slide[i]+500000),]
  for(n in 1:50){
    if(length(x[,1])==0){
      chr6_dxy_slide[i,n]<-"nan"
    } else {
      chr6_dxy_slide[i,n]<-mean(x[,n+2])
    }
  }
}
chr7<-tab[tab[,1]=="Chromosome_7",]
chr7_slide<-seq(1, 66640000, by=10000)
chr7_dxy_slide<-data.frame()
for(i in 1:length(chr7_slide)){
  x<-chr7[chr7[,2]>chr7_slide[i] & chr7[,2]<(chr7_slide[i]+500000),]
  for(n in 1:50){
    if(length(x[,1])==0){
      chr7_dxy_slide[i,n]<-"nan"
    } else {
      chr7_dxy_slide[i,n]<-mean(x[,n+2])
    }
  }
}
chr8<-tab[tab[,1]=="Chromosome_8",]
chr8_slide<-seq(1, 82060000, by=10000)
chr8_dxy_slide<-data.frame()
for(i in 1:length(chr8_slide)){
  x<-chr8[chr8[,2]>chr8_slide[i] & chr8[,2]<(chr8_slide[i]+500000),]
  for(n in 1:50){
    if(length(x[,1])==0){
      chr8_dxy_slide[i,n]<-"nan"
    } else {
      chr8_dxy_slide[i,n]<-mean(x[,n+2])
    }
  }
}
chr9<-tab[tab[,1]=="Chromosome_9",]
chr9_slide<-seq(1, 99980000, by=10000)
chr9_dxy_slide<-data.frame()
for(i in 1:length(chr9_slide)){
  x<-chr9[chr9[,2]>chr9_slide[i] & chr9[,2]<(chr9_slide[i]+500000),]
  for(n in 1:50){
    if(length(x[,1])==0){
      chr9_dxy_slide[i,n]<-"nan"
    } else {
      chr9_dxy_slide[i,n]<-mean(x[,n+2])
    }
  }
}

slide<-NULL
for(i in 1:10){
slide<-c(slide,chr1_dxy_slide[,1+2])
}
for(i in 1:10){
slide<-c(slide,chr2_dxy_slide[,1+2])
}
for(i in 1:10){
slide<-c(slide,chr3_dxy_slide[,1+2])
}
for(i in 1:10){
slide<-c(slide,chr4_dxy_slide[,1+2])
}
for(i in 1:10){
slide<-c(slide,chr5_dxy_slide[,1+2])
}
for(i in 1:10){
slide<-c(slide,chr6_dxy_slide[,1+2])
}
for(i in 1:10){
slide<-c(slide,chr7_dxy_slide[,1+2])
}
for(i in 1:10){
slide<-c(slide,chr8_dxy_slide[,1+2])
}
for(i in 1:10){
slide<-c(slide,chr9_dxy_slide[,1+2])
}

slide<-gsub("nan","NA", slide)
slide<-as.numeric(slide)

quantile(slide,prob=c(0.01,0.99), na.rm=T)

	# Step 5: plot the dxy
		# this is done in python using script: dxy_plotntIIvsntIIIa_500kb10kb.ipynb
		# all necessary files are provided in the the folder genome_scans/fxy

	# Step 6: bean plots
		# get the real data in 100kb windows for bean plots in R
R
tab<-read.table("ntII-ntIIIa_filtered_dxy.frq", header =T)

chr1<-tab[tab[,1]=="Chromosome_1",]
chr1_slide<-seq(1, 76910000, by=10000)
chr1_dxy_slide<-NULL
for(i in 1:length(chr1_slide)){
  x<-chr1[chr1[,2]>chr1_slide[i] & chr1[,2]<(chr1_slide[i]+100000),]
    if(length(x[,1])==0){
      chr1_dxy_slide[i]<-"nan"
    } else {
      chr1_dxy_slide[i]<-mean(x[,24])
    }
}
chr1_2<-data.frame(chr=rep("Chromosome_1", length(chr1_dxy_slide)), pos=chr1_slide+50000, dxy=chr1_dxy_slide)
chr2<-tab[tab[,1]=="Chromosome_2",]
chr2_slide<-seq(1, 88290000, by=10000)
chr2_dxy_slide<-NULL
for(i in 1:length(chr2_slide)){
  x<-chr2[chr2[,2]>chr2_slide[i] & chr2[,2]<(chr2_slide[i]+100000),]
    if(length(x[,1])==0){
      chr2_dxy_slide[i]<-"nan"
    } else {
      chr2_dxy_slide[i]<-mean(x[,24])
    }
}
chr2_2<-data.frame(chr=rep("Chromosome_2", length(chr2_dxy_slide)), pos=chr2_slide+50000, dxy=chr2_dxy_slide)
chr3<-tab[tab[,1]=="Chromosome_3",]
chr3_slide<-seq(1, 89510000, by=10000)
chr3_dxy_slide<-NULL
for(i in 1:length(chr3_slide)){
  x<-chr3[chr3[,2]>chr3_slide[i] & chr3[,2]<(chr3_slide[i]+100000),]
    if(length(x[,1])==0){
      chr3_dxy_slide[i]<-"nan"
    } else {
      chr3_dxy_slide[i]<-mean(x[,24])
    }
}
chr3_2<-data.frame(chr=rep("Chromosome_3", length(chr3_dxy_slide)), pos=chr3_slide+50000, dxy=chr3_dxy_slide)
chr4<-tab[tab[,1]=="Chromosome_4",]
chr4_slide<-seq(1, 67590000, by=10000)
chr4_dxy_slide<-NULL
for(i in 1:length(chr4_slide)){
  x<-chr4[chr4[,2]>chr4_slide[i] & chr4[,2]<(chr4_slide[i]+100000),]
    if(length(x[,1])==0){
      chr4_dxy_slide[i]<-"nan"
    } else {
      chr4_dxy_slide[i]<-mean(x[,24])
    }
}
chr4_2<-data.frame(chr=rep("Chromosome_4", length(chr4_dxy_slide)), pos=chr4_slide+50000, dxy=chr4_dxy_slide)
chr5<-tab[tab[,1]=="Chromosome_5",]
chr5_slide<-seq(1, 88020000, by=10000)
chr5_dxy_slide<-NULL
for(i in 1:length(chr5_slide)){
  x<-chr5[chr5[,2]>chr5_slide[i] & chr5[,2]<(chr5_slide[i]+100000),]
    if(length(x[,1])==0){
      chr5_dxy_slide[i]<-"nan"
    } else {
      chr5_dxy_slide[i]<-mean(x[,24])
    }
}
chr5_2<-data.frame(chr=rep("Chromosome_5", length(chr5_dxy_slide)), pos=chr5_slide+50000, dxy=chr5_dxy_slide)
chr6<-tab[tab[,1]=="Chromosome_6",]
chr6_slide<-seq(1, 69450000, by=10000)
chr6_dxy_slide<-NULL
for(i in 1:length(chr6_slide)){
  x<-chr6[chr6[,2]>chr6_slide[i] & chr6[,2]<(chr6_slide[i]+100000),]
    if(length(x[,1])==0){
      chr6_dxy_slide[i]<-"nan"
    } else {
      chr6_dxy_slide[i]<-mean(x[,24])
    }
}
chr6_2<-data.frame(chr=rep("Chromosome_6", length(chr6_dxy_slide)), pos=chr6_slide+50000, dxy=chr6_dxy_slide)
chr7<-tab[tab[,1]=="Chromosome_7",]
chr7_slide<-seq(1, 66640000, by=10000)
chr7_dxy_slide<-NULL
for(i in 1:length(chr7_slide)){
  x<-chr7[chr7[,2]>chr7_slide[i] & chr7[,2]<(chr7_slide[i]+100000),]
    if(length(x[,1])==0){
      chr7_dxy_slide[i]<-"nan"
    } else {
      chr7_dxy_slide[i]<-mean(x[,24])
    }
}
chr7_2<-data.frame(chr=rep("Chromosome_7", length(chr7_dxy_slide)), pos=chr7_slide+50000, dxy=chr7_dxy_slide)
chr8<-tab[tab[,1]=="Chromosome_8",]
chr8_slide<-seq(1, 82060000, by=10000)
chr8_dxy_slide<-NULL
for(i in 1:length(chr8_slide)){
  x<-chr8[chr8[,2]>chr8_slide[i] & chr8[,2]<(chr8_slide[i]+100000),]
    if(length(x[,1])==0){
      chr8_dxy_slide[i]<-"nan"
    } else {
      chr8_dxy_slide[i]<-mean(x[,24])
    }
}
chr8_2<-data.frame(chr=rep("Chromosome_8", length(chr8_dxy_slide)), pos=chr8_slide+50000, dxy=chr8_dxy_slide)
chr9<-tab[tab[,1]=="Chromosome_9",]
chr9_slide<-seq(1, 99980000, by=10000)
chr9_dxy_slide<-NULL
for(i in 1:length(chr9_slide)){
  x<-chr9[chr9[,2]>chr9_slide[i] & chr9[,2]<(chr9_slide[i]+100000),]
    if(length(x[,1])==0){
      chr9_dxy_slide[i]<-"nan"
    } else {
      chr9_dxy_slide[i]<-mean(x[,24])
    }
}
chr9_2<-data.frame(chr=rep("Chromosome_9", length(chr9_dxy_slide)), pos=chr9_slide+50000, dxy=chr9_dxy_slide)

ntIIntIIIa_dxy_slide<-rbind(chr1_2,chr2_2,chr3_2,chr4_2,chr5_2,chr6_2,chr7_2,chr8_2,chr9_2)
write.table(ntIIntIIIa_dxy_slide,"ntII-ntIIIa_100kb10kb_dxy.txt", col.names=T, row.names=F, sep="\t", dec=".")

		# get the start and end positions of each chromosomes
sed 's/"//g' ntII-ntIIIa_100kb10kb_dxy.txt > temp && mv temp ntII-ntIIIa_100kb10kb_dxy.txt
cat ${files}/chr.txt | while read line ; do grep "$line" ntII-ntIIIa_100kb10kb_dxy.txt | head -1 >> chr_start.txt ; done
cat ${files}/chr.txt | while read line ; do grep "$line" ntII-ntIIIa_100kb10kb_dxy.txt | tail -1 >> chr_end.txt ; done

		# make files with the windows to extract in R
R
start<-read.table("chr_start.txt", dec=".", sep="\t")
end<-read.table("chr_end.txt", dec=".", sep="\t")
tab<-cbind(start[,1:2], end[,2])
colnames(tab)<-c("chr","start","end")
chr<-list()
for(i in 1:length(tab[,1])){
  chr[[i]]<-data.frame(seq(tab[i,2],tab[i,3],by=100000))
  chr[[i]][,2]<-paste("Chromosome_",i,sep="")
  chr[[i]]<-data.frame(chr[[i]][,2],chr[[i]][,1])
  chr[[i]][,3]<-chr[[i]][,2]+99999
  chr[[i]][,3]<-format(chr[[i]][,3], scientific = FALSE)
  colnames(chr[[i]])<-c("chr","start","end")
}
for(i in 1:length(tab[,1])){
  write.table(chr[[i]][,1:2], paste("chr",i,"_100kb_windows.txt", sep=""),col.names=T, row.names=F, sep="\t", dec=".")
}

		#extract the windows I like
ls *windows.txt | sed 's/_100kb_windows.txt//g' > chr2.txt
cat chr2.txt | while read line ; do sed 's/"//g' "$line"_100kb_windows.txt > temp && mv temp "$line"_100kb_windows.txt ; done
cat chr2.txt | while read line1 ; do cat "$line1"_100kb_windows.txt | while read line2 ; do grep "$line2" ntII-ntIIIa_100kb10kb_dxy.txt >> "$line1"_100kb_windows_dxy.txt ; done ; done
cp chr1_100kb_windows_dxy.txt 100kb_windows_dxy.txt
cat chr3.txt | while read line ; do cat 100kb_windows_dxy.txt "$line"_100kb_windows_dxy.txt > temp && mv temp 100kb_windows_dxy.txt ; done

		# Plot the results in R
R
library(beanplot)
			# import data
tab<-read.table("100kb_windows_dxy.txt", header=F, dec=".", sep="\t")
colnames(tab)<-c("CHROM","BIN_START","dxy")
			# define the values for the windows with DE genes
genes<-data.frame(genes=c("pck","ppdk","asp-at"),fst=c(0.0108180082 ,0.018175574,0.0119036925))
			# plot
beanplot(tab$dxy, col="grey80", beanlines="median", beanlinewd=0.5, horizontal = TRUE, xlab="dxy")
for(i in 1:length(genes[,1])){
  abline(v=genes[i,2], col="red")
  text(x=genes[i,2], y=0.75, as.factor(genes[i,1]))
}
abline(v=quantile(tab$dxy, probs = c(0.05,0.95), na.rm=T), col="blue")
abline(v=quantile(tab$dxy, probs = c(0.01,0.98),na.rm=T), col="green")


