#Set files and directory names

WD= #~/IS_Motifs

SEQDIR=${WD}/Data/IS/Alpha_Moiani
SEQFILENAME=Alpha_Moiani_is26

OUTDIR=${WD}/Results/HotSpot/Alpha_Moiani

rm ${OUTDIR}/${SEQFILENAME}_Motif.txt

MOTIF=TTCA
NUCSBEFORE=6

egrep -B1 "^[ACGTN]{${NUCSBEFORE}}${MOTIF}" ${SEQDIR}/${SEQFILENAME}.fa | grep -v "^--" | \
awk -F'\n' 'BEGIN{RS=">";OFS="\t"} {print $1,$2}' | tail -n+2 | \
awk -F'[:_]' -v motif=${MOTIF} -v nnum=${NUCSBEFORE} 'BEGIN{OFS="\t"} {print $5":"$6,"N"nnum""motif,$3}' >> ${OUTDIR}/${SEQFILENAME}_Motif.txt

MOTIF=TGAA
NUCSBEFORE=16

egrep -B1 "^[ACGTN]{${NUCSBEFORE}}${MOTIF}" ${SEQDIR}/${SEQFILENAME}.fa | grep -v "^--" | \
awk -F'\n' 'BEGIN{RS=">";OFS="\t"} {print $1,$2}' | tail -n+2 | \
awk -F'[:_]' -v motif=${MOTIF} -v nnum=${NUCSBEFORE} 'BEGIN{OFS="\t"} {print $5":"$6,"N"nnum""motif,$3}' >> ${OUTDIR}/${SEQFILENAME}_Motif.txt

