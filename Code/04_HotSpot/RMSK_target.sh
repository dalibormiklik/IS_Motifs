#Calculate frequency of IS in repetitive elements

#Output:
#uniqNum_${VIRUS}_${MIX_COMP}_RMSK.txt - table containing the frequences of targeted repeats 

#GENOME	name of genome assembly
#RMSKDIR	path to the directory with the repeatmasker data
#RMSK		path to the BED file with repeatmasker coordinates
#ISDIR		path to the directory with the IS sequences
#VIRUS		base name for the sample
#MIX_COMP	name of the component to which IS sequnces are associated
#FADIR/	full path to the FASTA file with IS sequences
#  FILENAME/
#  FILE
#SEQDIR/	path to the txt file containing sequences of IS associated with the component
#   FILENAME/
#   FILE

#---------------------
#PART I
#---------------------
#Set variables
WD= #example: ./

GENOME=hg38

RMSKDIR= #example: ${WD}/Data/Repeats
RMSK=${RMSKDIR}/rmsk.bed

#Set sequence variables
ISDIR= #example: ${WD}/Results/HotSpot/
VIRUS=HIV_Zhyvoloup
MIX_COMP=M08_PPM07

FADIR= #example: ${WD}/Data/${VIRUS}
FAFILENAME=${VIRUS}_is26_CAP
SEQDIR= #example: ${WD}/Results
SEQFILENAME=ISseq_${VIRUS}_${MIX_COMP}

FAFILE=${FADIR}/${FAFILENAME}.fa
SEQFILE=${SEQDIR}/${SEQFILENAME}.txt

#---------------------
#PART II
#---------------------
#Extract fasta entries of selected IS

if [ -f "${SEQFILE}" ]; then
 
 echo -e "${SEQFILENAME}.txt in\n${SEQDIR}\nexists.\n"
 
 if [ -f "${FAFILE}" ]; then
 
  echo -e "${FAFILENAME}.txt in\n${FADIR}\n already exists.\nSkipping FASTA creation.\n"
  
 else
  
  echo -e "Creating ${FAFILENAME}.txt to\n${FADIR}\n."
 
  for ISSEQ in $(cat $SEQFILE)
   do
   grep --no-group-separator -B 1 "${ISSEQ}" ${FAFILE} >> ISseq_${VIRUS}_${MIX_COMP}.fa
   done
 
  echo -e "\nDONE\n"
  
  fi

else
 
 if [ -f "${FAFILE}" ]; then
  
  echo -e "${FAFILENAME}.txt in\n${FADIR}\n already exists.\nSkipping FASTA creation.\n"
  
 else
  
  echo -e "${SEQFILENAME}.txt in\n${SEQDIR}\nand\n${FAFILENAME}.txt in\n${FADIR}\nDO NOT EXIST.\n"
  echo "Stoping the script"
  exit

fi

#Create BED file from FASTA headers
grep "^>" ${FAFILE} | awk 'BEGIN{FS="::"; OFS="\t"} {print $3,$4,$1,0,"*"}' | awk 'BEGIN{FS="-"; OFS="\t"} {print $1,$2}' | sed 's/Chr/chr/g' | sed 's/>//g' | sort -k1,1 -k2,2n | uniq > ${SEQFILENAME}.bed
grep "^>" ${FAFILE} | awk 'BEGIN{FS="Chr"} {print "chr"$2}' | awk 'BEGIN{FS=":"; OFS="\t"} {print $1,$2}' | awk 'BEGIN{FS="-"; OFS="\t"} {print $1,$2}' | sed 's/>//g' | sort -k1,1 -k2,2n | uniq > ${SEQFILENAME}.bed
ISTOTAL=$(wc -l ${SEQFILENAME}.bed | awk '{print $1}')

#---------------------
#PART III
#---------------------
#Run BEDtools intercept between selected IS and Repeatmasker
##output with targeted repeats
bedtools intersect -wo -a ${ISDIR}/${VIRUS}_is26_sort.bed -b ${RMSKDIR}/rmsk.bed > ${SEQFILENAME}_RMSK.txt
##output number of IS outside of any repeat
bedtools intersect -v -a ${ISDIR}/${VIRUS}_is26_sort.bed -b ${RMSKDIR}/rmsk.bed > ${SEQFILENAME}_outRMSK.txt

#Create table with counts and "full" repeats names (family etc)
RUNINAM=$(awk 'BEGIN{FS="\t"} {print $10}' ISseq_${VIRUS}_${MIX_COMP}_RMSK.txt | sort | uniq)

rm uniqNum_${VIRUS}_${MIX_COMP}_RMSK.txt
#loop on unique repeat names and count occurence in IS RMSK data
for RNAMi in $(echo $RUNINAM)
 do
 RNUM=$(awk -v Rnam=${RNAMi} 'BEGIN{FS="\t"} $10==Rnam {print $4}' ${SEQFILENAME}_RMSK.txt | sort | uniq | wc -l | awk '{print $1}')
 awk -v Rnum=$RNUM -v Rnam=${RNAMi} 'BEGIN{OFS="\t"} $1==Rnam {print Rnum,Rnam,$2,$3}' ${RMSKDIR}/uniq_rmsk_names.txt >> uniqNum_${VIRUS}_${MIX_COMP}_RMSK_tmp.txt
 done

#Add row with number of IS not overlapping any repeat
RNUM=$(wc -l ${SEQFILENAME}_outRMSK.txt | awk '{print $1}')
echo -e ${RNUM}"\tNone\tNone\tNone" >> uniqNum_${VIRUS}_${MIX_COMP}_RMSK_tmp.txt

sort -n -k1 -r uniqNum_${VIRUS}_${MIX_COMP}_RMSK_tmp.txt > uniqNum_${VIRUS}_${MIX_COMP}_RMSK.txt
rm uniqNum_${VIRUS}_${MIX_COMP}_RMSK_tmp.txt
head uniqNum_${VIRUS}_${MIX_COMP}_RMSK.txt

