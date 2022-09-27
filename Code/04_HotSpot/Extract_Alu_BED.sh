#This code extracts coordinates of Alu repeats from RepeatMasker using awk code

#Annotations for the repeats is downloaded as a Repeatmasker (rmsk.txt) file from the ucsc goldenpath

#Set variables for path, directories and file names
#GENOME -	name of genome assembly
#RMSKDIR -	path to directory with RepeatMasker file
#ANNDIR -	path to the directory where BED file of extracted coordinates is stored

GENOME= #example: hg38
RMSKDIR= #example: ../Data/${GENOME}/Repeats
ANNDIR= #Data/Annotations/${GENOME}/

#Select only Alu repeats and export in BED format
awk 'BEGIN{OFS="\t"} $13 == "Alu" {print $6,$7,$8,$11,$2,$10}' ${RMSKDIR}/rmsk_${GENOME}.txt > ${ANNDIR}/${GENOME}_rmsk_Alu.bed

