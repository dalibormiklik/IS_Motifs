#Separate positions of patterns that are located inside given feature (for instance Alu repeats)

#Tools: bedtools (intersect)

#Output:	${PATPOSFILE}_in_${FEATNAME}.bed - coordinates of the intra-feature sequence motifs
#		${PATPOSFILE}_in_${FEATNAME}_len1.bed - ${PATPOSFILE}_in_${FEATNAME}.bed reduced to single (central) position

#WD - Working Directory
#GENOME - 	name of genome assembly
#       -	in publication: hg18, hg19 or hg38
#FEATDIR - 	path to the directory where BED file with feature coordinates (Alu) is stores
#FEATNAME -	name of a BED file with genomic features
#         -	in publication: ${GENOME}_rmsk_Alu created by Extract_Alu.sh
#PATPOSFILE -	name of the BED file containing coordinates of the sequence pattern; created by Locate_pattern_in_genome.sh
#

GENOME=hg38
FEATNAME=${GENOME}_rmsk_Alu
PATPOSFILE=${GENOME}_HIV_palMotif

WD= #example: ./
FEATDIR= #example: ${WD}/Data/Annotations/BED

#Run bedtools intersect

echo "Running bedtools intersect."
echo "Pattern file: ${PATPOSFILE}"
echo "Feature: ${FEATNAME}"

echo "..."

bedtools intersect -wo -a ${PATPOSFILE}.bed -b ${FEATDIR}/${FEATNAME}.bed > ${PATPOSFILE}_in_${FEATNAME}.txt

echo ""
echo "${PATPOSFILE}_in_${FEATNAME}.txt ready:"
echo ""

#Create BED file of pattern positions from Intersect file
#Strand is the strand of a feature
#Name is concatenation of feature name and feature genomic position
#Score field is filled by pattern 
echo "Creating BED file"
echo "..."

awk -v OFS="\t" '{print $1,$2,$3,$7":"$8"-"$9"_"$10,$4,$12}' ${PATPOSFILE}_in_${FEATNAME}.txt | sort -k1,1 -k2,2n > ${PATPOSFILE}_in_${FEATNAME}.bed

#Create 1bp-wide features
awk -v OFS="\t" '{print $1, $2+6, $3-6, $4, $5, $6}' ${PATPOSFILE}_in_${FEATNAME}.bed > ${PATPOSFILE}_in_${FEATNAME}_len1.bed

echo ""
echo "Files ${PATPOSFILE}_in_${FEATNAME}.bed"
echo "      ${PATPOSFILE}_in_${FEATNAME}_len1.bed"
echo "Created"
echo ""

head ${PATPOSFILE}_in_${FEATNAME}.bed
echo "..."
tail ${PATPOSFILE}_in_${FEATNAME}.bed

