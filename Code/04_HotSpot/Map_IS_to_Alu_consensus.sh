#Map intra Alu IS to consesnsus ALu sequence to analyze position of IS relative to consesnsus

#Tools:
#bowtie2, samtools (view), bedtools (bamtobed)

#Output:
#${SAMOUTNAME}.bed - coordinates of mapped reads relative to consensus as a reference
#${SAMOUTNAME}_PosFreq.txt - calculated frequencies of mapped reads per position in consesnus

#---------------------
#PART I
#---------------------
#Build bowtie index files from consesus fasta
#indexes are called "Alu"
bowtie-build Alu_consensus.fa Alu

#---------------------
#PART II
#---------------------
#Set variables

#WD - Working Directory
#GENOME - 	name of genome assembly
#       -	in publication: hg38
#GENOMFADIR -	path to the directory where genome assembly sequence in FASTA format is stored
#VIRUS -	base name of sample
#INPUTBEDDIR -	path to directory containing file with IS coordinates
#INPUTBEDNAME -name of IS coordinate file
#QUERYFANAME -	name of FASTA file with sequences that will be mapped
#SAMOUTNAME -	name for output SAM file

WD= #example: ./

GENOME=hg38
GENOMFADIR= #example: ${WD}/Data/Annotations/${GENOME}/FASTA

VIRUS=HIV_Zhyvoloup_M08_PPM07
INPUTBEDDIR=${WD}/Results/HotSpot
INPUTBEDNAME=${VIRUS}_in_hg38_rmsk_Alu

QUERYFANAME=${VIRUS}_in_hg38_rmsk_Alu

SAMOUTNAME=${QUERYFANAME}_to_consensus

#---------------------
#PART III
#---------------------
#Map sequences
bowtie2 -f -L 5 -N 1 -i S,1,0.2 --score-min L,0,-2 --all -x Alu -U ${QUERYFANAME}.fa -S ${SAMOUTNAME}.sam

#Calulate numbers of alignment per sequence
samtools view -F 0x4 ${SAMOUTNAME}.sam | cut -f1 | sort | uniq -c > ${QUERYFANAME}_AlignCount.txt

#Select alignments of single-mapped reads
grep "^@" ${SAMOUTNAME}.sam > ${SAMOUTNAME}_single.sam
grep --no-group-separator "$(awk '$1 == 1 {print $2}' ${QUERYFANAME}_AlignCount.txt)" ${SAMOUTNAME}.sam >> ${SAMOUTNAME}_single.sam

#Select FASTA reads of multimappers
grep -A1 --no-group-separator "$(awk '$1 > 1 {print $2}' ${QUERYFANAME}_AlignCount.txt)" ${QUERYFANAME}.fa > ${SAMOUTNAME}.fa

#Map multimappers from 1st mapping
bowtie2 -f -L 5 -N 1 -i S,1,0.2 --all -x Alu -U ${SAMOUTNAME}.fa -S ${SAMOUTNAME}_mII.sam

#Calulate numbers of alignment per sequence
samtools view -F 0x4 ${SAMOUTNAME}_mII.sam | cut -f1 | sort | uniq -c > ${QUERYFANAME}_AlignCount.txt

#Select alignments of single-mapped reads
grep --no-group-separator "$(awk '$1 == 1 {print $2}' ${QUERYFANAME}_AlignCount.txt)" ${SAMOUTNAME}.sam >> ${SAMOUTNAME}_single.sam

#Convert alignments in SAM to BED
samtools view ${SAMOUTNAME}_single.sam -b | bedtools bamtobed -i stdin | awk -v OFS="\t" '{print $0,"IS"}' > ${SAMOUTNAME}.bed

#Remove temporary files
rm ${SAMOUTNAME}.sam ${SAMOUTNAME}_mII.sam ${QUERYFANAME}_AlignCount.txt ${SAMOUTNAME}.fa

#Reduce ranges to middle base
#Calcuate frequency of positions
awk '{print $2"\t"$3}' ${SAMOUTNAME}.bed | sort | uniq -c | awk '{print $2"\t"$3"\t"$1}' | sort -k3nr > ${SAMOUTNAME}_PosFreq.txt


echo ""
echo "Created $(wc -l ${SAMOUTNAME}.bed | cut -d" " -f1) ranges from $(grep -c ">" ${QUERYFANAME}.fa) sequences."
echo ""
head ${SAMOUTNAME}.bed
echo "..."
tail ${SAMOUTNAME}.bed
echo ""
echo "Position counts in ${SAMOUTNAME}_PosFreq.txt file:"
head ${SAMOUTNAME}_PosFreq.txt
echo ""
head ""

