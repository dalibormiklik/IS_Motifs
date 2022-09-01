#Run code from directory dedicated to sample

#In part I, variables are set.
#In part II, calculations are performed.

#tools used in the code: bedtools (shuffle, intersect, getfasta)

#output:
#${FEATNAME}_target_freq.txt -	frequency of IS ranges intersect with the feature ranges

#---------------------
#PART I
#---------------------

#Paths to the directories/files need to exist before the code is ran
#for example: IS directory (${WD}/Data/${VIRUS}) with the IS range coordinates BED file (${VIRUS}_is${SEQLEN}_sort) need to be already present

#Set sequence variables
#If more samples are entered, separate names by "," (colon)
#WD - Working Directory
#SMPLSET -	can be one of predefined sets of HIV IS: Zhyvoloup, Demeulemeester or Vansant
#SEQLEN -	lengt of IS sequence range
#       -	in publication SEQLEN=12 is used (needs to be created by Create_Ranges.sh)

WD= #example: ./
SMPLSET= #example: Zhyvoloup

SEQLEN=12

#Select Alu entries from repeatmasker, save as a BED file and extract fasta sequences
#Set file with repeatmasker

#GENOME - 	name of genome assembly
#       -	in publication: hg18, hg19 or hg38
#FEATNAME -	name of a BED file with genomic features
#         -	in publication: ${GENOME}_rmsk_Alu created by Extract_Alu.sh
#GENOMFADIR -	path to the directory where genome assembly sequence in FASTA format is stored
#ChromLen -	sorted txt file with chromoseome lengths (needed for shuffled controls to be cerated)
#FEATDIR - 	path to the directory where BED file with feature coordinates (Alu) is stores

GENOME= #example: hg38
FEATNAME= #example: ${GENOME}_rmsk_Alu

GENOMFADIR= #example: ${WD}/Data/Annotations/${GENOME}/FASTA
ChromLen= #example: ${WD}/Data/Annotations/${GENOME}/chromInfo_sort.txt
FEATDIR= #example: ${WD}/Data/Annotations/BED

#Set number of itteration for shuffled control used for motif targeting
SHUFFLE_ITTERS=100

#Set SMPLNAME according to name of selected set of IS given by SMPLSET variable
#The preset SMPLNAME represents sets of HIV-1 IS coordinates analyzed in the publication

if [ $(echo $SMPLSET) == "Zhyvoloup" ]; then
 SMPLNAME=HIV_Zhyvoloup,HIV_Zhyvoloup_N74D
 echo ${SMPLNAME}
 fi
 
if [ $(echo $SMPLSET) == "Vansant" ]; then
 SMPLNAME=HIV_Vansant_SupT1_wt_0_mix,HIV_Vansant_SupT1_wt_6_mix,HIV_Vansant_SupT1_KD_0_mix
 echo ${SMPLNAME}
 fi
 
if [ $(echo $SMPLSET) == "Demeulemeester" ]; then SMPLNAME=HIV_Demeulemeester_WT,HIV_Demeulemeester_S119A,HIV_Demeulemeester_S119T,HIV_Demeulemeester_S119T_R231G,HIV_Demeulemeester_R231G,HIV_Demeulemeester_R231K,HIV_Demeulemeester_R231L,HIV_Demeulemeester_R231Q,HIV_Demeulemeester_S119P,HIV_Demeulemeester_S119G,HIV_Demeulemeester_S119I,HIV_Demeulemeester_S119K,HIV_Demeulemeester_S119R,HIV_Demeulemeester_S119V
 echo ${SMPLNAME}
 fi


#---------------------
#PART II
#---------------------


#Number of samples per set
SMPLTOTAL=$(echo ${SMPLNAME} | sed 's/,/\n/g' | wc -l | cut -d" " -f1)
SMPLACT=1

#Run loop for each SMPLNAME
#Particular SMPLNAME is stored in VIRUS variable
 
for VIRUS in $(echo ${SMPLNAME} | sed 's/,/ /g')
 do
 
 echo ""
 echo "Analyze ${VIRUS}: ${SMPLACT}/${SMPLTOTAL}"
 echo ""
 
 #Set base name for IS-containing BED file + directory where the BED file is found
 ISNAME=${VIRUS}_is${SEQLEN}_sort
 ISDIR=${WD}/Data/${VIRUS}
 
 #Create randomly shuffled sequences for all IS
 echo "- Create shuffled control sites"
 echo "  ..."
 bedtools shuffle -i ${ISDIR}/${ISNAME}.bed -g ${ChromLen} | sort -k1,1 -k2,2n > ./${ISNAME}_shuffle.bed

 ###Targeting of repeat elements
 echo "- Targeting ${FEATNAME}"
 echo "  ..."

 #Run BEDtools intercept between selected IS and Repeatmasker
 ##output with targeted repeats
 echo " - IS intersect"
 echo "   ..."
 bedtools intersect -wo -a ${ISDIR}/${ISNAME}.bed -b ${FEATDIR}/${FEATNAME}.bed > ${VIRUS}_is${SEQLEN}_in_${FEATNAME}.txt
 cut -f 1-5,12 ${VIRUS}_is${SEQLEN}_in_${FEATNAME}.txt > ${VIRUS}_is${SEQLEN}_in_${FEATNAME}.bed
 cut -f 7-11,12 ${VIRUS}_is${SEQLEN}_in_${FEATNAME}.txt | uniq > ${FEATNAME}_with_${VIRUS}_is${SEQLEN}.bed

 #Get shuffle control IS in Alu
 echo " - Shuffle intersect"
 echo "   ..."
 bedtools intersect -wo -a ./${ISNAME}_shuffle.bed -b ${FEATDIR}/${FEATNAME}.bed > ${VIRUS}_is${SEQLEN}_in_${FEATNAME}_shuffle.txt
 cut -f 1-5,12 ${VIRUS}_is${SEQLEN}_in_${FEATNAME}_shuffle.txt > ${VIRUS}_is${SEQLEN}_in_${FEATNAME}_shuffle.bed
 cut -f 7-11,12 ${VIRUS}_is${SEQLEN}_in_${FEATNAME}_shuffle.txt | uniq > ${FEATNAME}_with_${VIRUS}_is${SEQLEN}_shuffle.bed
 
 echo " Intersection with ${FEATNAME} done."
 echo ""
 
 ##Create shuffled controls targeting repeat elements (shuffle intra element IS)
 echo "- Create feature-matched shuffled IS"
 echo "  ..."
 
 #Target random IS to defined (BED) regions
 #Target all repeats
 echo "  - Target all features"
 echo "    ..."
 
 TARGETBED=${FEATDIR}/${FEATNAME}.bed
 
 bedtools shuffle -i ${VIRUS}_is${SEQLEN}_in_${FEATNAME}.bed -g ${ChromLen} -incl ${TARGETBED} -f 0.6 | sort -k1,1 -k2,2n > ./${VIRUS}_is${SEQLEN}_in_${FEATNAME}_shuffle_all.bed

#Target exact repeats targeted by IS
 echo "  - Target IS-matched features"
 echo "    ..."
 
 TARGETBED=${FEATNAME}_with_${VIRUS}_is${SEQLEN}.bed
 
 bedtools shuffle -i ${VIRUS}_is${SEQLEN}_in_${FEATNAME}.bed -g ${ChromLen} -incl ${TARGETBED} -f 0.6 | sort -k1,1 -k2,2n > ./${VIRUS}_is${SEQLEN}_in_${FEATNAME}_shuffle_same.bed
 
 echo ""
 
 #
 #Get sequences of selected IS in repeat
 echo "- Get genomic sequences of intra-feature IS"
 echo "  -IS"
 echo "  ..."
 bedtools getfasta -name -s -fi ${GENOMFADIR}/${GENOME}.fa -bed ${VIRUS}_is${SEQLEN}_in_${FEATNAME}.bed -fo ${VIRUS}_is${SEQLEN}_in_${FEATNAME}.fa
 awk -F'\n' 'BEGIN{RS=">";OFS="\t"} {print $1,$2}' ${VIRUS}_is${SEQLEN}_in_${FEATNAME}.fa | tail -n+2 > ${VIRUS}_is${SEQLEN}_in_${FEATNAME}_fa_tab.txt
 
 echo "  -shuffle-all"
 echo "  ..." 
 bedtools getfasta -name -s -fi ${GENOMFADIR}/${GENOME}.fa -bed ${VIRUS}_is${SEQLEN}_in_${FEATNAME}_shuffle_all.bed -fo ${VIRUS}_is${SEQLEN}_in_${FEATNAME}_shuffle_all.fa
 awk -F'\n' 'BEGIN{RS=">";OFS="\t"} {print $1,$2}' ${VIRUS}_is${SEQLEN}_in_${FEATNAME}_shuffle_all.fa | tail -n+2 > ${VIRUS}_is${SEQLEN}_in_${FEATNAME}_shuffle_all_fa_tab.txt

 echo "  -shuffle-feature-matched"
 echo "  ..."
 bedtools getfasta -name -s -fi ${GENOMFADIR}/${GENOME}.fa -bed ${VIRUS}_is${SEQLEN}_in_${FEATNAME}.bed -fo ${VIRUS}_is${SEQLEN}_in_${FEATNAME}_shuffle_same.fa
 awk -F'\n' 'BEGIN{RS=">";OFS="\t"} {print $1,$2}' ${VIRUS}_is${SEQLEN}_in_${FEATNAME}_shuffle_same.fa | tail -n+2 > ${VIRUS}_is${SEQLEN}_in_${FEATNAME}_shuffle_same_fa_tab.txt
 
 echo " Done."
 echo ""
 
 #Create tables of IS sequences in rmsk from FASTA files
 INPUTFILE=${VIRUS}_is${SEQLEN}_in_${FEATNAME}
 awk -F'\n' 'BEGIN{RS=">";OFS="\t"} {print $1,$2}' ${INPUTFILE}.fa | tail -n+2 > ${INPUTFILE}_fa_tab.txt
 awk -F'\n' 'BEGIN{RS=">";OFS="\t"} {print $1,$2}' ${INPUTFILE}_shuffle_all.fa | tail -n+2 > ${INPUTFILE}_shuffle_all_fa_tab.txt
 awk -F'\n' 'BEGIN{RS=">";OFS="\t"} {print $1,$2}' ${INPUTFILE}_shuffle_same.fa | tail -n+2 > ${INPUTFILE}_shuffle_same_fa_tab.txt

 
 #If specified create shuffled control composed of multiple random sequences
 #number of itterations is defined by variable SHUFFLE_ITTERS
 
  if [ $(echo ${SHUFFLE_ITTERS}) -gt 1 ]
  then
   
   echo "Create ${SHUFFLE_ITTERS} columns of shuffled sequences."
   echo ""
   
   INPUTFILE=${VIRUS}_is${SEQLEN}_in_${FEATNAME}
   SHUFFLEOUTNAME=${INPUTFILE}_shuffle_all
   #Target random IS to defined (BED) regions
   ##Target all repeats
   TARGETBED=${FEATDIR}/${FEATNAME}.bed
    
   #Create first column of final file with itterations
   echo "Creating 1st column."
   echo "..."
   cut -f2 ${SHUFFLEOUTNAME}_fa_tab.txt > ${SHUFFLEOUTNAME}_itter1_tmp_tab.txt
     
   ACT_i=2
    
   while [ "${ACT_i}" -le "${SHUFFLE_ITTERS}" ]
   do
    echo "${VIRUS}: ${SMPLACT}/${SMPLTOTAL}"
    echo "Itteration: ${ACT_i}"
    echo "..."
    
    ##Create random BED regions
    bedtools shuffle -i ${VIRUS}_is${SEQLEN}_in_${FEATNAME}.bed -g ${ChromLen} -incl ${TARGETBED} -f 0.6 | sort -k1,1 -k2,2n > ./${SHUFFLEOUTNAME}_itter${ACT_i}.bed
    ##Get FASTA seq of target sequences
    bedtools getfasta -name -s -fi ${GENOMFADIR}/${GENOME}.fa -bed ${SHUFFLEOUTNAME}_itter${ACT_i}.bed -fo ${SHUFFLEOUTNAME}_itter${ACT_i}.fa
    ##Create table
    awk -F'\n' 'BEGIN{RS=">";OFS="\t"} {print $2}' ${SHUFFLEOUTNAME}_itter${ACT_i}.fa | tail -n+2 > ${SHUFFLEOUTNAME}_itter${ACT_i}_col.txt
    
    #rm all itteratipon files but FA col.txt
    rm ./${SHUFFLEOUTNAME}_itter${ACT_i}.bed ${SHUFFLEOUTNAME}_itter${ACT_i}.fa 
    
    paste ${SHUFFLEOUTNAME}_itter$(($ACT_i-1))_tmp_tab.txt ${SHUFFLEOUTNAME}_itter${ACT_i}_col.txt > ${SHUFFLEOUTNAME}_itter${ACT_i}_tmp_tab.txt
    
    rm ${SHUFFLEOUTNAME}_itter$(($ACT_i-1))_tmp_tab.txt ${SHUFFLEOUTNAME}_itter${ACT_i}_col.txt
    
    echo "Itteration: ${ACT_i} DONE"
    echo ""
   
    #if [ "$ACT_i" -lt  "${SHUFFLE_ITTERS}" ]; then
     ACT_i=$(($ACT_i+1))
     #fi
    
   done 
   
  fi
  
 mv ${SHUFFLEOUTNAME}_itter${SHUFFLE_ITTERS}_tmp_tab.txt ${SHUFFLEOUTNAME}_${SHUFFLE_ITTERS}x_tab.txt
 ls -la *itter*
 
 head ${SHUFFLEOUTNAME}_${SHUFFLE_ITTERS}x_tab.txt

#Done for sample (virus)
SMPLACT=$((SMPLACT+1))
done


#Remove frequency table if present
rm ${FEATNAME}_target_freq.txt

echo -e "Calculate IS frequencies:\n"
echo -e "\t-IS set: ${SMPLSET}\n"
echo -e "\t-IS samples in set: ${SMPLTOTAL}\n"
echo -e "\t-Feature: ${FEATNAME}\n"
echo -e "\t-Saving into ${FEATNAME}_target_freq.txt\n"

echo -e "\t-Save into ${FEATNAME}_target_freq.txt\n"

echo -e "RUN ->\n"

SMPLACT=1
for VIRUS in $(echo ${SMPLNAME} | sed 's/,/ /g')
 do
 
 echo -e "Set ${SMPLSET}\n"
 echo -e "Sample ${SMPLSET}: ${SMPLACT}/${SMPLACT}\n"
 
 ISNAME=${VIRUS}_is${SEQLEN}_sort
 ISDIR=${WD}/IS_seq/${VIRUS}
 
 #Calculate number of IS in repeat
 #Header line - do not include if stats of more samples are going to be present
 #echo -e "Group\tVirus\tAlu\tAll" >> Alu_target_freq.txt
 echo -e "Shuffle\t${VIRUS}\t$(wc -l ${VIRUS}_is${SEQLEN}_in_${FEATNAME}_shuffle.bed | cut -d " " -f1)\t$(wc -l ${ISNAME}_shuffle.bed | cut -d " " -f1)" >> ${FEATNAME}_target_freq.txt
 
 echo -e "IS\t${VIRUS}\t$(wc -l ${VIRUS}_is${SEQLEN}_in_${FEATNAME}.bed | cut -d " " -f1)\t$(wc -l ${ISDIR}/${ISNAME}.bed | cut -d " " -f1)" >> ${FEATNAME}_target_freq.txt
 
 echo -e "DONE ${SMPLACT}/${SMPLACT}\n"
 tail -2 ${FEATNAME}_target_freq.txt
 echo ""
  
 SMPLACT=$((SMPLACT+1))
 
 done

echo -e "Feature targeting for ${FEATNAME} DONE\n"

cat ${FEATNAME}_target_freq.txt

