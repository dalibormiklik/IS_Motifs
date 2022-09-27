#The code for transformation of position-wise IS BED ranges to IS-surrounding ranges
#dist variables are set to produce ranges of length 26
#tsd variable set the retrovirus-specific target site duplication length
#xtr variable is 1 if tsd is odd

#ISFILE is the name of a BED file containing single-position IS coordinate

#HIV
##Zhyvoloup et al. 2017
#INVAR is either "wt" or "n74d". This variable selects the virus variant.
INVAR="wt"
awk -v tsd=2 -v dist=26 -v xtr=1 -v invar=$INVAR 'BEGIN{start=(dist/2)+xtr-tsd; end=tsd+(dist/2)} $5 == invar && $3 == "+" {print $1"\t"$2-start"\t"$2+end"\t"$4"\t""1""\t"$4}' ${ISFILE}.bed > ${OUTFILE}.bed
awk -v tsd=2 -v dist=26 -v xtr=1 -v invar=$INVAR 'BEGIN{start=(dist/2)+xtr+tsd; end=(dist/2)-tsd} $5 == invar && $3 == "-" {print $1"\t"$2-start"\t"$2+end"\t"$4"\t""2""\t"$4}' ${ISFILE}.bed >> ${OUTFILE}.bed

##Vansant et al. 2020
awk -v tsd=2 -v dist=26 -v xtr=1 'BEGIN{start=(dist/2)-xtr-tsd+1+1; end=(dist/2)+xtr+tsd-1} $6 == "-" {print $1"\t"$2-start"\t"$2+end"\t"$4"\t"$5"\t"$6}' ${ISFILE}.bed > ${OUTFILE}.bed
awk -v tsd=2 -v dist=26 -v xtr=1 'BEGIN{start=(dist/2)+xtr+tsd+1; end=(dist/2)-tsd-xtr} $6 == "+" {print $1"\t"$2-start"\t"$2+end"\t"$4"\t"$5"\t"$6}' ${ISFILE}.bed >> ${OUTFILE}.bed

##Demeulemeester et al. 2016
awk -v tsd=2 -v dist=26 -v xtr=1 'BEGIN{start=(dist/2)+xtr-tsd; end=tsd+(dist/2)} $6 == "+" {print $1"\t"$2-start+1"\t"$2+end+1"\t"$4"\t""1""\t"$4}' ${ISFILE}.bed > ${OUTFILE}.bed
 awk -v tsd=2 -v dist=26 -v xtr=1 'BEGIN{start=(dist/2)+xtr+tsd; end=(dist/2)-tsd} $6 == "-" {print $1"\t"$2-start"\t"$2+end"\t"$4"\t""2""\t"$4}' ${ISFILE}.bed >> ${OUTFILE}.bed

#MLV
awk -v tsd=2 -v dist=26 -v xtr=0 'BEGIN{start=(dist/2)+xtr-tsd; end=(dist/2)-tsd} {print $1"\t"$2-start"\t"$3+end"\t"$4"\t""1""\t"$4}' ${ISFILE}.bed > ${OUTFILE}.bed

#MVV
awk -v dist=26 'BEGIN{start=(dist/2)-1; end=(dist/2)+1} {print $1"\t"$2-start"\t"$2+end"\t"$4"\t""1""\t"$6}' ${ISFILE}.bed > ${OUTFILE}.bed

#PFV
awk -v tsd=4 -v dist=26 -v xtr=0 'BEGIN{start=(dist/2)+xtr; end=(dist/2)} {print $1"\t"$2-start"\t"$3+end-1"\t"$5"\t""1""\t"$4}' ${ISFILE}.bed > ${OUTFILE}.bed

#ASLV
##Malhotra et al. 2017
awk -v tsd=6 -v dist=26 -v xtr=0 'BEGIN{start=(dist/2)+xtr-tsd; end=tsd+(dist/2)} $4 == "+" {print $1"\t"$2-start"\t"$2+end"\t"$5"\t""1""\t"$5}' ${ISFILE}.bed > ${OUTFILE}.bed
awk -v tsd=6 -v dist=26 -v xtr=0  'BEGIN{start=(dist/2)+xtr+tsd; end=(dist/2)-tsd} $4 == "-" {print $1"\t"$2-start"\t"$2+end"\t"$5"\t""2""\t"$5}' ${ISFILE}.bed >> ${OUTFILE}.bed

##Moiani et al. 2014
awk -v tsd=6 -v dist=26 -v xtr=0 'BEGIN{start=(dist/2)+xtr; end=(dist/2)} {print $1"\t"$2-start"\t"$3+end-1"\t"$5"\t""1""\t"$4}' ${ISFILE}.bed > ${OUTFILE}.bed

