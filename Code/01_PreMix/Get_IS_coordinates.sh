#HTLV
#Sequences of pre-integration genomic sequences were obtained directly from supplementary data of the publication.

#HIV
#Zhyvoloup et al. 2017
#IS coordinates were obtained from a publication supplementary table. Replicates of DMSO treated samples of wt and capsid N74D mutant were joined to create wt and N74D IS sets. Coordinates of both IS sets were joined to create a single file.

#Vansant et al. 2020
#Coordinates were obtained from GSE135295_ledgins1_integration_features table. Mock, LEDGIN+ and LEDGF-KD data sets were created joining bulk, GFP+ and GFP- samples from same treatment group marked by the key present in field 8 of the feature table as follows: wt (S168, S169, S178), LEDGIN+ (S170, S171, S179), LEDGF-KD (S174, S175, S181).

#Demeulemeester et al. 2016
#Coordinates of IS in the form of text file were obtained upon request from the authors of the study. Data obtained by transduction of different cell types with identical variant were mixed. The IS coordinates of the IN variants were transformed to BED format with custom script. 

#MLV
#Raw reads were obtained from Sequence Read Archive using:
fastq-dump --split-e command.

#Sequences were trimmed using cutadapt with following sequences to be trimmed:
LTR3nest=TGACTACCCGTCAGCGGGGGTC
LTR3rest=TTTCA
LTR5nest=CAAACCTACAGGTGGGGTCTTTC
LTR5rest=A
Adaptor=AGTCCCTTAAGCGGAGCCCT
#First 5’/3’ LTR sequences from primer (LTR*nest) were trimmed with
cutadapt -m 11 -O 20 -e 0.1 -j 0 --trim-n
#Then rest of each LTR (LTR*rest) was trimmed using 
cutadapt -m X -O 20 -e 0 -j 0
#where X in -m option equals to LTR*rest sequence length. Adapter was removed from sequences using
cutadapt -m 11 -O 10 -e 0.1 -j 0
#Finally sequences containing inner proviral sequences were removed using
cutadapt -m 11 -O 10 -e 0.1 -j 0
#The resulting LTR5_F_full.fastq and LTR3_F_full.fastq were mapped to hg38 human genome assembly using
bowtie2 with -p 20 -q --no-unal -x hg38
#Next, alignments of mapped reads were filtered to start at position 0 with
grep -v "MD:Z:0"
#Alignments with a single hit in the genome were then sorted with
grep -v "XS:i:" 
#Final BED file was created with samtools view, bedtools bamtobed and sort commands.


#MVV
#BED formatted IS coordinates retrieved from MVV vector infected HEK293T cells were downloaded from Gene Expression Omnibus (study accession GSE196042).

#PFV
#BED formatted IS coordinates were downloaded from Gene Expression Omnibus (study accession GSE97973). 

#ASLV
#Malhotra et al. 2017
#SAM formatted alignments were obtained from the authors of the study upon request. First, data from experiments with identical vectors and cells were joined (48h and 120h collection time points). For the compatibility with downstream tools, SAM headers from hg19 (HeLa)/galGal4 (CEF) were added to alignments. Next, alignments of mapped reads were filtered to start at position 0 with grep -v "MD:Z:0". Alignments with a single hit in the genome were then sorted with grep -v "XS:i:" command. Final BED file was created with samtools view, bedtools bamtobed and sort commands. BED ranges were converted to contain single LTR-proximal position with custom awk script:
awk '{if ($6 == "+") print $1"\t"$1":"$2"_"$6; else print $1"\t"$1":"$3"_"$6;}'
#If the mapped strand was “+”, the start coordinate was used, otherwise the end coordinate was used. The awk command was followed by:
sort | uniq -c | awk '{print $3"\t"$1}'
#to count occurrences of each IS. Only IS coordinates with more than one occurrence were selected for further analysis.

#Moiani et al. 2014 
#Raw reads were obtained from Sequence Read Archive using
fastq-dump --split-e command.
#Sequences were trimmed using cutadapt with following sequences to be trimmed:
LTR=TTGGTGTGCACCTGGGTTGATGGCCGGACCGTTGATTCCCTGACGACTACGAGCACCTGCATGAAGCAGAAGG
LTRend=^CTTCA
ADAPT=GTCCCTTAAC
MseADAPT=TTAGTCC
#First, the majority of LTR sequence was trimmed:
cutadapt -g $LTR -m 20 -O 30 -e 0.1 -j 0 --trim-n
#Next, the end of LTR was trimmed if the sequence was exactly matching the LTRend:
cutadapt -g $LTRend -m 15 -O 5 -e 0 -j 0
#Finally, the adaptor sequence was removed from the reads:
cutadapt -a $MseADAPT -m 15 -O 7 -e 0.1 -j 0
#In adapter-containing sequences, the end of the read was modified to contain MseI site instead of MseI-Adapter junction using sed:
sed 's/TTAGTCC$/TTAA/'
#The reads were mapped to hg38 genome assembly using
bowtie2 with -p 20
#parameter set.
#In the next step, alignments were filtered using custom code. First, the header of SAm file was removed with
grep -v "^@"
#and only alignments starting at position 0 were selected using
grep -v "MD:Z:0" 
#and only reads with single alignment were further selected with
grep -v "XS:i:"
#Final BED file was created with SAMtools view, BEDtools bamtobed and sort commands.
#To calculate the frequency of each IS in the IS set, first, the awk was used to convert IS into ISIDs:
awk '{if ($6 == "+") print $1"\t"$1":"$2+1"_"$6; else print $1"\t"$1":"$3"_"$6;}'
#and then counted with:
sort | uniq -c | awk '{print $3"\t"$1}' .
#IS with at least 5 occurrences were selected and the BED coordinates were transformed to represent the central position of the target site with awk:
awk -v tsd=6 'BEGIN{half=tsd/2} {if ($4 =="+") print $1"\t"$2+half-1"\t"$2+half"\t"$4"\t"$5"\t"$6; else print $1"\t"$2-half"\t"$2-half+1"\t"$4"\t"$5"\t"$6}' .
