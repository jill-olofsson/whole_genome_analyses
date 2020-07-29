####Analyses of RADSeq
	# Publication: Olofsson et al 2020 

# Step 1: define location of data and programs
data=/home/$USER/path_to_raw_data
adapt=/home/$USER/path_to_fasta_file_with_adaptor_sequence
barcodes=/home/$USER/path_to_file_with_barcodes
TRIM=/home/$USER/path_to_Trimmomatic-0.36
stacks=/home/$USER/path_to_stacks

# Step 2: Trim and clean raw data with Trimmomatic
	# input is the two raw read files in fastq.gz format
	# need to specify 4 output files: 2 for each read, 1 is the reads which have a pair and the other are the unpaired reads
	# input trailing ILLUMINACLIP are settings for the program
	# need to specify the adaptor sequences in a fasta file
java -jar $TRIM/trimmomatic-0.36.jar PE -phred33 $data/data_R1.fastq.gz $data/data_R2.fastq.gz data_R1_paired.fq.gz data_R1_unpaired.fq.gz data_R2_paired.fq.gz data_R2_unpaired.fq.gz ILLUMINACLIP:${adapt}/adaptor_trimmomatic.fa:2:30:10:1:true LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

# Step 3: De-multiplex using the process_radtag in stacks
	# http://catchenlab.life.illinois.edu/stacks/comp/process_radtags.php
	# only use the paired reads from above
mkdir demultiplex
${stacks}/process_radtags -1 fata_R1_paired.fq.gz -2 fata_R2_paired.fq.gz -i gzfastq -b ${barcodes}/data_barcodes -o demultiplex -r -e ecoRI --inline_null

