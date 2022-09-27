#Calculate distance of IS found in Alu to nearest targeted Alu-derived Palindromic motif

#In part I, variables are set.
#In part II, calculations are performed.

#tools used in the code: bedtools (closest, shuffle)

#output: ${VIRUS}_in_${FEATNAME}_palMotif_distance.txt - table of IS distances to the nearest sequence motifs (only distances up to $MAXDIST are stored)
#		+ table containing distances of shuffled controls

#---------------------
#PART I
#---------------------

#Set sequence variables
#If more samples are entered, separate names by "," (colon)
#WD - Working Directory
#SMPLSET -	can be one of predefined sets of HIV IS: Zhyvoloup, Demeulemeester or Vansant
#SEQLEN -	lengt of IS sequence range
#       -	in publication SEQLEN=12 is used (needs to be created by Create_Ranges.sh)

WD= #example: ./
SMPLSET= #example: Zhyvoloup

SEQLEN=12

#GENOME - 	name of genome assembly
#       -	in publication: hg18, hg19 or hg38
#FEATNAME -	name of a BED file with genomic features
#         -	in publication: ${GENOME}_rmsk_Alu created by Extract_Alu.sh
#FEATDIR - 	path to the directory where BED file with feature coordinates (Alu) is stores
#ChromLen -	sorted txt file with chromoseome lengths (needed for shuffled controls to be cerated)
#MOTIFCOORD -	name of file with the coordinates of the sequences motif; created by Locate_pattern_in_genome.sh

GENOME= #example: hg38
ChromLen= #example: ${WD}/Data/Annotations/${GENOME}/chromInfo_sort.txt

FEATNAME= #example: ${GENOME}_rmsk_Alu
FEATDIR= #example: ${WD}/Data/Annotations/BED

MOTIFCOORD= #example: ${FEATDIR}/${GENOME}_HIV_palMotif_in_${FEATNAME}_len1

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
 
if [ $(echo $SMPLSET) == "Demeulemeester" ]; then
 SMPLNAME=HIV_Demeulemeester_WT,HIV_Demeulemeester_S119A,HIV_Demeulemeester_S119T,HIV_Demeulemeester_S119T_R231G,HIV_Demeulemeester_R231G,HIV_Demeulemeester_R231K,HIV_Demeulemeester_R231L,HIV_Demeulemeester_R231Q,HIV_Demeulemeester_S119P,HIV_Demeulemeester_S119G,HIV_Demeulemeester_S119I,HIV_Demeulemeester_S119K,HIV_Demeulemeester_S119R,HIV_Demeulemeester_S119V
 echo ${SMPLNAME}
 fi


#---------------------
#PART II
#---------------------

#Number of samples per set
SMPLTOTAL=$(echo ${SMPLNAME} | sed 's/,/\n/g' | wc -l | cut -d" " -f1)
SMPLACT=1

#Run loop
for VIRUS in $(echo ${SMPLNAME} | sed 's/,/ /g')
 do
 
 #Set base name for IS-containing BED file + directory where the BED file is found
 ISNAME=${VIRUS}_is${SEQLEN}_sort
 ISDIR=${WD}/IS_seq/${VIRUS}
 
 echo "Calculate distances of ${VIRUS} (${SMPLACT}/${SMPLTOTAL})"
 echo "to nearest ${MOTIFCOORD}"
 echo "..."
 
 
 #Prior distance calculation, shrink BED intervals of IS to single nucleotide position in the middle of interval
 ##All IS
 ###IS
 awk -v OFS="\t" '{print $1, $2+6, $3-6, "IS_"$4, $5, "*"}' ${ISDIR}/${ISNAME}.bed | sort -k1,1 -k2,2n | bedtools closest -a stdin -b ${MOTIFCOORD}.bed -D b | awk -F"\t" '!_[$4]++' | awk -v OFS="\t" '$8 != -1' > ${VIRUS}_palMotif_distance.txt
 ###Shuffle control
 awk -v OFS="\t" '{print $1, $2+6, $3-6, "Shuffle_"$4, $5, "*"}' ${ISNAME}_shuffle.bed | sort -k1,1 -k2,2n | bedtools closest -d -a stdin -b ${MOTIFCOORD}.bed | awk -F"\t" '!_[$4]++'| awk -v OFS="\t" '$8 != -1' >> ${VIRUS}_palMotif_distance.txt
 
 #Filter distance results
 #Select single line for each IS (unique column 4 with IS name)
 #MAXDIST sets tha maximum IS-FEATURE distance that will be kept 
 MAXDIST=1000
 awk -v OFS="\t" -v maxdist=${MAXDIST} '$8 != -1 && $13 >=(-1 * maxdist) && $13 <= maxdist' ${VIRUS}_palMotif_distance.txt | awk -F"\t" '!_[$4]++' > ${VIRUS}_palMotif_distance_up${MAXDIST}.txt
 
 ##Intra-Alu IS
 #IS
 awk -v OFS="\t" '{print $1, $2+6, $3-6, "IS_"$4, $5, $6}' ${VIRUS}_is${SEQLEN}_in_${FEATNAME}.bed | bedtools closest -a stdin -b ${MOTIFCOORD}.bed -D b | awk -F"\t" '!_[$4]++' > ${VIRUS}_in_${FEATNAME}_palMotif_distance.txt
 
 ##Shuffle control - same Alu
 if [ $(echo ${SHUFFLE_ITTERS}) -gt 1 ]
 then
  
  echo "Create ${SHUFFLE_ITTERS} columns of shuffled distances."
  echo ""
  
  INPUTFILE=${VIRUS}_is${SEQLEN}_in_${FEATNAME}
  SHUFFLEOUTNAME=${INPUTFILE}_shuffle_same_palMotif_distance
  #Target random IS to defined (BED) regions
  ##Target all repeats
  TARGETBED=${FEATNAME}_with_${VIRUS}_is${SEQLEN}
  
  #Create first column of final file with itterations
  echo "Creating 1st column."
  echo "..."
  bedtools shuffle -i ${INPUTFILE}.bed -g ${ChromLen} -incl ${TARGETBED}.bed -f 0.6 | sort -k1,1 -k2,2n | awk -v OFS="\t" '{print $1, $2+6, $3-6, $4, $5, $6}' | bedtools closest -a stdin -b ${MOTIFCOORD}.bed -D b | awk -F"\t" '!_[$4]++' | cut -f13 > ${SHUFFLEOUTNAME}_itter1_tmp_tab.txt
    
  ACT_i=2
  
  while [ "${ACT_i}" -le "${SHUFFLE_ITTERS}" ]
  do
   echo "${VIRUS} (${SMPLACT}/${SMPLTOTAL})"
   echo " Itteration: ${ACT_i} / ${SHUFFLE_ITTERS}"
   echo " ..."
   
   ##Create ACT_i th column of distances
   bedtools shuffle -i ${INPUTFILE}.bed -g ${ChromLen} -incl ${TARGETBED}.bed -f 0.6 | sort -k1,1 -k2,2n | awk -v OFS="\t" '{print $1, $2+6, $3-6, $4, $5, $6}' | bedtools closest -a stdin -b ${MOTIFCOORD}.bed -D b | awk -F"\t" '!_[$4]++' | cut -f13 > ${SHUFFLEOUTNAME}_itter${ACT_i}_col.txt
   
   #Join i-th column with previous column file   
   paste ${SHUFFLEOUTNAME}_itter$(($ACT_i-1))_tmp_tab.txt ${SHUFFLEOUTNAME}_itter${ACT_i}_col.txt > ${SHUFFLEOUTNAME}_itter${ACT_i}_tmp_tab.txt
   
   rm ${SHUFFLEOUTNAME}_itter$(($ACT_i-1))_tmp_tab.txt ${SHUFFLEOUTNAME}_itter${ACT_i}_col.txt
   
   echo " Itteration: ${ACT_i} / ${SHUFFLE_ITTERS}  DONE"
   echo ""
   
   #Next step...
   ACT_i=$(($ACT_i+1))
      
  done 
  
 fi
 
 mv ${SHUFFLEOUTNAME}_itter${SHUFFLE_ITTERS}_tmp_tab.txt ${SHUFFLEOUTNAME}_${SHUFFLE_ITTERS}x_tab.txt
 ls -la *itter*

 head ${SHUFFLEOUTNAME}_${SHUFFLE_ITTERS}x_tab.txt

 echo "${VIRUS} (${SMPLACT}/${SMPLTOTAL})"
 echo "calculations DONE."
 echo ""
 
 SMPLACT=$((SMPLACT+1))
 done

echo -e "\nIS Set ${SMPLSET}\nDONE\n"

