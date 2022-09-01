#use seqkit tool to localize all ocurences of sequence pattern in a genome
#define pattern in regExp style ( "."  instead of "N")

PATTERN=CT..G...C..AG
GENOME=hg38
OUTNAME=${GENOME}_HIV_palMotif

#Path to the directory with genome assembly sequence stored in FASTA file
WD= #example: ./
GENOMFADIR= #example: ${WD}/Data/Annotations/${GENOME}/FASTA

#Run seqkit
echo "Running seqkit locate."
echo "Pattern: ${PATTERN}"
echo "Genome: ${GENOME}"

echo "..."

cat ${GENOMFADIR}/${GENOME}.fa | seqkit locate -i -r -p ${PATTERN} --bed > ${OUTNAME}.bed

echo ""
echo "${OUTNAME}.bed ready:"
head ${OUTNAME}.bed
echo "..."
tail ${OUTNAME}.bed
