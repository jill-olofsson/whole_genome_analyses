# This is the genetic analyses of all samples

# Analyses 1: Phylogeny of nuclear SNPs
	# step 0: define location of the program unless they are already in your path
RAXML=/home/$USER/path_to_raxml
	# step 1: execute raxml
${RAXML}/raxmlHPC-SSE3 -s merged_plate01-12_ASEM_C4_v1_unique_f2_covfilt_allsites_TEmask_chr_max99_SNP_mac56_missing90_thin1k.phy -n merged_plate01-12_ASEM_C4_v1_unique_f2_covfilt_allsites_TEmask_chr_max99_SNP_mac56_missing90_thin1k -m GTRGAMMA -p 12345 -f a -# 100 -x 12345
	# phy file provided in files/phy

# Analyses 2: Preform pca on the nuclear SNPs
	# step 0: covert vcf file to structure file using PDGSpider
		# define location of programs unless they are already in your path
PDGspider=/home/$USER/path_to_PDGSpider
coverter=/home/$USER/path_to_converter_file

sed 's/0\/0:0,.,./0\/0:0,20,20/g' merged_plate01-12_ASEM_C4_v1_unique_f2_covfilt_allsites_TEmask_chr_max99_SNP_mac56_missing90_thin1k.recode.vcf >  merged_plate01-12_ASEM_C4_v1_unique_f2_covfilt_allsites_TEmask_chr_max99_SNP_mac56_missing90_thin1k_fixPL.vcf

java -Xmx1024m -Xms512m -jar ${PDGSpider}/PGDSpider2-cli.jar -inputfile merged_plate01-12_ASEM_C4_v1_unique_f2_covfilt_allsites_TEmask_chr_max99_SNP_mac56_missing90_thin1k_fixPL.vcf inputformat VCF -outputfile merged_plate01-12_ASEM_C4_v1_unique_f2_covfilt_allsites_TEmask_chr_max99_SNP_mac56_missing90_thin1k_fixPL.str -outputformat STRUCTURE -spid ${coverter}/vcf-to-str.spid
	# str file provided in files/pca

	# Step 1: Preform the pca using adgenet in R
		# prepare a file with all the samples in the order they appear and meta data such as the nuclear clade assignment of each sample (provided)
		# open R and load adegenet
R
library("adegenet")
		# load the data
			# need to know number of individuals and loci in the structure file
data1<-read.structure("merged_plate01-12_ASEM_C4_v1_unique_f2_covfilt_allsites_TEmask_chr_max99_SNP_mac56_missing90_thin1k_fixPL.str", n.ind=561, n.loc=18513, onerowperind = F, col.lab=1, col.pop=0, ask=FALSE)
samp<-read.table("str_order_20200421.txt", header=T, sep="\t", dec=".")
		# preform the pca
			# attach the individual and the nt-clade assignment to the data
indNames(data1)<-as.vector(samp$samp)
pop(data1)<-as.vector(samp$nt)
			# replace missing data with mean
x.data1 <- scaleGen(data1, NA.method="mean", scale=FALSE)
			# do PCA keeping 4 axes
pca<-dudi.pca(x.data1, center=FALSE, scale=FALSE,scannf = FALSE, nf = 4)
			# example of how to plot the results
				#define colors
col<-c("gold", "forestgreen","orange", "blue", "red", "black")
				#plot PC 1 vs 2
s.class(pca$li,fac=pop(data1),axesel=FALSE, cstar=0, cpoint=1, col=col, xax=1, yax=2)
title("PCA 1-2")
				#plot PC 1 vs 3
s.class(pca$li,fac=pop(data1),axesel=FALSE, cstar=0, cpoint=1, col=col, xax=1, yax=3)
title("PCA 1-3")
				#plot PC 2 vs 3
s.class(pca$li,fac=pop(data1),axesel=FALSE, cstar=0, cpoint=1, col=col, xax=2, yax=3)
title("PCA 2-3")


# Analyses 3: admixture using NGSAdmix
	# Step 0: define location of programs  needed unless they are already in your path
BEAGLE=/home/$USER/path_to_beagle
NGSadmix=/home/$USER/path_to_ngsadmix

	# Step 1: covert the vcf file to a beagle file
cat merged_plate01-12_ASEM_C4_v1_unique_f2_covfilt_allsites_TEmask_chr_max99_SNP_mac56_missing90_thin1k.recode.vcf | java -jar ${Beagle}/vcf2beagle.jar 0 merged_plate01-12_ASEM_C4_v1_unique_f2_covfilt_allsites_TEmask_chr_max99_SNP_mac56_missing90_thin1k
awk '{print $1, $3, $4}' merged_plate01-12_ASEM_C4_v1_unique_f2_covfilt_allsites_TEmask_chr_max99_SNP_mac56_missing90_thin1k.markers > markers.txt
java -jar ${BEAGLE}/beagle2gprobs.jar markers.txt merged_plate01-12_ASEM_C4_v1_unique_f2_covfilt_allsites_TEmask_chr_max99_SNP_mac56_missing90_thin1k.bgl.gz 0  > merged_plate01-12_ASEM_C4_v1_unique_f2_covfilt_allsites_TEmask_chr_max99_SNP_mac56_missing90_thin1k.beagle.gz

	# Step 2: Run NGSadmix for K 1-10
${NGSadmix}/NGSadmix -likes merged_plate01-12_ASEM_C4_v1_unique_f2_covfilt_allsites_TEmask_chr_max99_SNP_mac56_missing90_thin1k.beagle.gz -K 1 -outfiles admix_indmax99_chr_SNPs_TEmask_mac56_thin1kb_90missing_01 -seed 12345

${NGSadmix}/NGSadmix -likes merged_plate01-12_ASEM_C4_v1_unique_f2_covfilt_allsites_TEmask_chr_max99_SNP_mac56_missing90_thin1k.beagle.gz -K 2 -outfiles admix_indmax99_chr_SNPs_TEmask_mac56_thin1kb_90missing_02 -seed 12345

${NGSadmix}/NGSadmix -likes merged_ASEM_C4_v1_unique_f2_covfilt_allsites_TEmask_chr_max99_SNP_mac56_missing90_thin1k.beagle.gz -K 3 -outfiles admix_indmax99_chr_SNPs_TEmask_mac56_thin1kb_90missing_03 -seed 12345

${NGSadmix}/NGSadmix -likes merged_plate01-12_ASEM_C4_v1_unique_f2_covfilt_allsites_TEmask_chr_max99_SNP_mac56_missing90_thin1k.beagle.gz -K 4 -outfiles admix_indmax99_chr_SNPs_TEmask_mac56_thin1kb_90missing_04 -seed 12345

${NGSadmix}/NGSadmix -likes merged_plate01-12_ASEM_C4_v1_unique_f2_covfilt_allsites_TEmask_chr_max99_SNP_mac56_missing90_thin1k.beagle.gz -K 5 -outfiles admix_indmax99_chr_SNPs_TEmask_mac56_thin1kb_90missing_05 -seed 12345

${NGSadmix}/NGSadmix -likes merged_plate01-12_ASEM_C4_v1_unique_f2_covfilt_allsites_TEmask_chr_max99_SNP_mac56_missing90_thin1k.beagle.gz -K 6 -outfiles admix_indmax99_chr_SNPs_TEmask_mac56_thin1kb_90missing_06 -seed 12345

${NGSadmix}/NGSadmix -likes merged_plate01-12_ASEM_C4_v1_unique_f2_covfilt_allsites_TEmask_chr_max99_SNP_mac56_missing90_thin1k.beagle.gz -K 7 -outfiles admix_indmax99_chr_SNPs_TEmask_mac56_thin1kb_90missing_07 -seed 12345

${NGSadmix}/NGSadmix -likes merged_plate01-12_ASEM_C4_v1_unique_f2_covfilt_allsites_TEmask_chr_max99_SNP_mac56_missing90_thin1k.beagle.gz -K 8 -outfiles admix_indmax99_chr_SNPs_TEmask_mac56_thin1kb_90missing_08 -seed 12345

${NGSadmix}/NGSadmix -likes merged_plate01-12_ASEM_C4_v1_unique_f2_covfilt_allsites_TEmask_chr_max99_SNP_mac56_missing90_thin1k.beagle.gz -K 9 -outfiles admix_indmax99_chr_SNPs_TEmask_mac56_thin1kb_90missing_09 -seed 12345

${NGSadmix}/NGSadmix -likes merged_plate01-12_ASEM_C4_v1_unique_f2_covfilt_allsites_TEmask_chr_max99_SNP_mac56_missing90_thin1k.beagle.gz -K 10 -outfiles admix_indmax99_chr_SNPs_TEmask_mac56_thin1kb_90missing_10 -seed 12345

	# Step 3: Example on how to plot the results in R
		# read in data
samp<-read.table("str_order_20200421.txt", header=T, sep="\t", dec=".")
admix_k3<-t(as.matrix(read.table("admix_indmax99_chr_SNPs_TEmask_mac56_thin1kb_90missing_03.qopt")))
		# sort the data so individuals appear in the desired order
admix2_k3<-admix_k3[,order(pop[,12])]
		# plot the results in a barplot
barplot(admix2_k3,col=c("red","gold", "blue"),space=0,border=NA,xlab="Individuals",ylab="admixture")

# Analyses 4: Identity by descent
	# this is just an example using a small set of populations (all populations assigned to the nuclear clade I in the manuscript).

	# Step 1: calculate the pairwise Fst between all populations
		# to prepare
			# make a list of all populations (ntI.txt)
			# for each population make a file with all the individuals in that population (pop.txt where pop is replaced with each of the population names) save in the foled ${comp} as defined below
			# for each population make a file with all the populations to which it should be compares (pop_vs.txt where pop is replaced with each of the population names) save in the foled ${comp} as defined below
ntI=/home/$USER/files/IBD/
comp=/home/$USER/files/IBD/ntI
		# use vcftools to calculate Fst in sliding windows of length 500kb slide 10kb
mkdir output
mkdir output/original
cat ${ntI}/ntI.txt | while read lineA ; do cat ${comp}/"$lineA"_vs.txt | while read lineB ; do vcftools --gzvcf merged_plate01-12_ASEM_C4_v1_unique_f2_covfilt_allsites_TEmask_chr_max99_SNP_mac56_missing90_thin1k.recode.vcf --max-missing 0.1 --weir-fst-pop ${comp}/"$lineA".txt --weir-fst-pop ${comp}/"$lineB".txt --fst-window-size 500000 --fst-window-step 10000 --mac 3 --out output/origina/500kb10kb_90missing_mac3_"$lineA"_"$lineB" ; done ; done
		# insert windows with missing data into the file to be able to average across all possible windows later
			#list all output files
ls output/original | sed 's/output\/original\///g' | sed 's/500kb10kb_90missing_mac3_//g' | sed 's/.windowed.weir.fst.gz//g' > files.txt
			# make column for merging
mkdir output/merge
cat files.txt | while read line ; do zcat output/original/500kb10kb_90missing_mac3_"$line".windowed.weir.fst.gz  | awk 'BEGIN{OFS="\t"}; {print $1":"$2, $0}'| sort -k2,2 -k3,3n | awk 'BEGIN{OFS="\t"}; {print $1, $2, $3, $4, $7}' > output/merge/500kb10kb_90missing_mac3_"$line"_tomerge.txt ; done
			# merge the fst files with the bed file containign all possible windows
				# merge the files (insert windows not in the original fst-file
				# replace the MEAN_FST with the comparision name
				# replace all negative Fst's with 0
cat files_allcomp.txt | while read line ; do awk 'BEGIN{OFS="\t"}; BEGIN{OFS="\t"}; NR==FNR{a[$1]=$0 ; next} {{ found=0; for(i=1;i<=NF;i++) { if($i in a) { print a[$i]; found=1; break; } } if (!found) {print $0}}}' output/merge/500kb10kb_90missing_mac3_"$line"_tomerge.txt wind_500Kb10Kb_tomerge.bed | sed 's/MEAN_FST/'$line'/' | awk 'BEGIN{OFS="\t"};$5<0 {$5=0} 1'> temp && mv temp output/merge/500kb10kb_90missing_mac3_"$line"_tomerge.txt ; done
		# calculate mean of each comparision
			# remove lines with NA in the 5th columnlength()
			# calculate mean and put in a file
cat files_allcomp.txt | while read line ; do awk 'BEGIN{OFS="\t"};$5 != "NA"' output/merge/500kb10kb_90missing_mac3_"$line"_tomerge.txt | awk '{ sum += $5 } END {print sum / NR }' >> mean_fst_90missing_mac3_allcomp.txt ; done
			# add information about the comparision
paste files_allcomp.txt mean_fst_90missing_mac3_allcomp.txt > temp && mv temp mean_fst_90missing_mac3_allcomp.txt

	# Step 2: Calculate the geographical distance in km between each population
		# this is done in R
R
		# define function to do the calculations: OBS! these functions are for java I have converted themt o R code
			# from https://stackoverflow.com/questions/365826/calculate-distance-between-2-gps-coordinates
degreesToRadians<-function(degrees) {
  return(degrees*pi/180)
}
distanceInKmBetweenEarthCoordinates<-function(lat1, lon1, lat2, lon2) {
  earthRadiusKm<-6371
  dLat<-degreesToRadians(lat2-lat1)
  dLon<-degreesToRadians(lon2-lon1)
  lat1<-degreesToRadians(lat1)
  lat2<-degreesToRadians(lat2)
  a<-sin(dLat/2)*sin(dLat/2)+sin(dLon/2)*sin(dLon/2)*cos(lat1)*cos(lat2)
  c<-2*atan2(sqrt(a),sqrt(1-a))
  return(earthRadiusKm*c)
}

		# read table with GPS coordinates for each population
tab<-read.table("gps_table_Africa_RAD_clean_for_IBD.txt", header=T, sep="\t", dec=".")
			# want to go down the table and compare the first row with all others below then the 2nd row with all others below etc. 
			#define matrix
geo<-data.frame(matrix(data=NA, nrow=length(tab[,1]), ncol=length(tab[,1])+1))
				#add pop names to first columns
					# change if you pop names ar not in the first column
geo[,1]<-tab[,1]
				#add pop names to column headers
					# change if you pop names are not in the first column
colnames(geo)<-c("NA",as.vector(tab[,1]))

			# fill in the matrix
				# this makes a distance matrix it only addes to the lower half of the matrix, the upper half remains NA
for(n in 1:length(tab[,1])){
vec<-rep("NA", length(tab[,1]))	# define an empty vector to put distances in
  for(i in n:length(tab[,1])){
  vec[i]<-distanceInKmBetweenEarthCoordinates(tab[n,6],tab[n,7],tab[i,6],tab[i,7])	# add to vector from pop where you compare to, change column numbers to where you have your lat long coordinates
  geo[,n+1]<-vec	# add vector to distance matrix
}}
colnames(geo)[1]<-"pops"
			# write out distance matrix
write.table(geo, "geo_dist_all_African_ms_20191119.txt", sep="\t", na="NA", dec=".", row.names=F, col.names=T)
			# want a table of values instead
geo_2<-list()
for(n in 1:(length(tab[,1])-1)){
geo_2[[n]]<-data.frame(pop1=rep(tab[n,1], length(tab[n:length(tab[,1]),1])-1), pop2=tab[(n+1):length(tab[,1]),1])
vec<-NULL
for(i in n:(length(tab[,1])-1)){
 dist<-distanceInKmBetweenEarthCoordinates(tab[n,6],tab[n,7],tab[i+1,6],tab[i+1,7])
 vec<-c(vec,dist)
}
geo_2[[n]][,3]<-vec
colnames(geo_2[[n]])[3]<-"dist"
geo_2[[n]]$comp<-paste(geo_2[[n]][,1],geo_2[[n]][,2], sep="_")
}
dist<-data.frame()
for(i in 1:length(geo_2)){
dist<-rbind(dist,geo_2[[i]])
}
dist<-dist[,c(4,1,2,3)]
write.table(dist,"dist_table_Africa_RAD_20191119.txt", col.names=T,row.names=F, sep="\t")

	# Step 3: Put all the information into one big table (this is provided here for all the comparisions)
		# file: allcomp_IDB_tab_90missing_mac3.txt

	# Step 4: Calculate the correlations using Mantel tests in R
R
		# identify the samples belonging to each clade
			# read in a table with the group codes
groupcode<-read.table('group_code',header=F)
			# define groups
C3<-names(table(as.vector(groupcode$V1[groupcode$V2=="C3" & groupcode$V3>2]))) # clade I
Int<-names(table(as.vector(groupcode$V1[groupcode$V2=="Int" & groupcode$V3>2]))) # clade II
C4iiia<-names(table(as.vector(groupcode$V1[groupcode$V2=="C4iiia" & groupcode$V3>2]))) # clade IIIa
C4iiib<-names(table(as.vector(groupcode$V1[groupcode$V2=="C4iiib" & groupcode$V3>2]))) # clade IIIb
C4iv<-names(table(as.vector(groupcode$V1[groupcode$V2=="C4iv" & groupcode$V3>2]))) # clade IV

		# upload the data containing pairwise genetic and geographic distances
data<-read.table('allcomp_IDB_tab_90missing_mac3.txt',header=T)

		# for populations assigned to different nuclear clades, example of C3 (=clade I) vs C4iiib (=clade IIIb)
C3C4iiib_dist<-matrix(nrow=length(C3),ncol=length(C4iiib))
C3C4iiib_fst<-matrix(nrow=length(C3),ncol=length(C4iiib))
for(i in 1:length(C3)){
for(j in 1:length(C4iiib)){
C3C4iiib_dist[i,j]<-data$dist[(data$pop1==C3[i] & data$pop2==C4iiib[j]) | (data$pop2==C3[i] & data$pop1==C4iiib[j])]
C3C4iiib_fst[i,j]<-data$fst[(data$pop1==C3[i] & data$pop2==C4iiib[j]) | (data$pop2==C3[i] & data$pop1==C4iiib[j])]
}
}
permu<-c()
for(i in 1:9999){
permu[i]<-cor.test(as.vector(as.matrix(as.data.frame(C3C4iiib_fst)[sample(1:dim(C3C4iiib_fst)[1]),sample(1:dim(C3C4iiib_fst)[2])])),as.vector(C3C4iiib_dist),method=c("spearman"))$estimate
}
(1+length(permu[permu>cor.test(as.vector(as.matrix(C3C4iiib_fst)),as.vector(as.matrix(C3C4iiib_dist)),method=c("spearman"))$estimate]))/(9999+1)
summary(lm(as.vector(as.matrix(C3C4iiib_fst))~as.vector(as.matrix(C3C4iiib_dist))))

		# for populations assigned to the name nuclear clade, example of C3 (=clade I)
C3C3_dist<-matrix(nrow=length(C3),ncol=length(C3))
C3C3_fst<-matrix(nrow=length(C3),ncol=length(C3))
for(i in 1:length(C3)){
for(j in 1:length(C3)){
if(i==j){
C3C3_dist[i,j]<-0
C3C3_fst[i,j]<-0
}
if(i!=j){
C3C3_dist[i,j]<-data$dist[(data$pop1==C3[i] & data$pop2==C3[j]) | (data$pop2==C3[i] & data$pop1==C3[j])]
C3C3_fst[i,j]<-data$fst[(data$pop1==C3[i] & data$pop2==C3[j]) | (data$pop2==C3[i] & data$pop1==C3[j])]
}
}
}
permu<-c()
for(i in 1:9999){
newvec<-sample(1:dim(C3C3_fst)[1])
newmat<-as.matrix(as.data.frame(C3C3_fst)[newvec,newvec])
permu[i]<-cor.test(newmat[row(C3C3_fst)>col(C3C3_fst)],C3C3_dist[row(C3C3_fst)>col(C3C3_fst)],method=c("spearman"))$estimate
}
(length(permu[permu>cor.test(C3C3_fst[row(C3C3_fst)>col(C3C3_fst)],C3C3_dist[row(C3C3_fst)>col(C3C3_fst)],method=c("spearman"))$estimate])+1)/(9999+1)
summary(lm(C3C3_fst[row(C3C3_fst)>col(C3C3_fst)]~C3C3_dist[row(C3C3_fst)>col(C3C3_fst)]))
