#!/bash

#Code should be run from the directory where directories containing the EM algorith recors are present
#(so they cen be listed by ls -d in DATASETS variable)

set -e

#Read names of viral data sets
#DATASETS should contain the full name of the directory in which the report is saved
#To find the full name, grep function looking for keys of the name is uses

DATASETS_KEY=MVV_HEK

DATASETS=$( ls -d */ | grep "${MVV_HEK}" | sed 's/\///g')

for DATASET in $(echo $DATASETS)
 do
 
 echo "Working on mixtures of ${DATASET}"
 echo ""
 
 #Check if directories of mixtures are not called PPM instead of M
 if compgen -G "${DATASET}/PPM*" > /dev/null
  then
   echo "Renaming PPM to M"
   echo ""
   for MIXNAME in ${DATASET}/PPM*
    do
    mv $MIXNAME ${MIXNAME/\/PPM/\/M}
    done
  fi
  
  echo "Mixture directories present:"
  ls -d ${DATASET}/M*
 
 for MIXPATH in $(ls -d ${DATASET}/M*)
  do
  
  #Save name of actual mixture
  MIXNAME=${MIXPATH##*/}
  
  #Check if final file PPM.txt already exists
  if compgen -G "${MIXPATH}/PPM.txt" > /dev/null; then
   echo "${MIXNAME} PPM.txt already exists."
   echo "Removing PPM.txt from ${MIXNAME}."
   rm ${MIXPATH}/PPM.txt
   echo "Creating new PPM.txt"
   echo ""
   fi
   
    echo "Extracting PPMs for $MIXPATH"
    echo ""
   
    FILENAME=${MIXPATH}/PARAM.GEN
  
    head ./${FILENAME}
  
    #Create table with all PPMs divided by one line with PPM number
    tail -n +3 ${FILENAME} | sed 's/ \+ /\t/g' | sed 's/^\t//g' > ./${MIXPATH}/PPMs_tmp.txt
    
    #Extract maximum number of components in mixture from PARAM file name
    MAXMIX=${MIXNAME#M*}
    
    for PPMNUM in `seq -w 00 $MAXMIX`
     do
     
     #Save PPM name as variable
     PPM="PPM"$PPMNUM
     
     echo "Extracting ${MIXNAME} ${PPM}..."
     echo ""
     
     #Extract PPM "header" - first row of "tmp"
     PPMHEAD=$(head -1 ./${MIXPATH}/PPMs_tmp.txt)
     
     #Extract weight of component
     echo $PPMHEAD > ./${MIXPATH}/${PPM}_Wm_tmp.txt
     
     #Set number of rows for PPM
     #number of rows is strored in second field of header
     # + 1 is added to extract all rows of PPM
     PPMROW=$(echo $PPMHEAD | awk '{print $2+1}')
     
     #Extract PPM
     #Add column with PPM name
     sed -n -e 2,${PPMROW}p ./${MIXPATH}/PPMs_tmp.txt | awk -v ppm=${PPM} 'BEGIN{OFS="\t"} {print ppm,$0}' > ./${MIXPATH}/${PPMNUM}_PPM_tmp.txt
     
     #Delete PPM prom temporary file
     sed -i 1,${PPMROW}d ./${MIXPATH}/PPMs_tmp.txt
     
     echo "${MIXNAME} ${PPM} extracted."
     echo ""
     
     done
    
    #Clear - remove PPM temporary files
    #rm ./${MIXPATH}/PPM*tmp.txt
    
    #Join all PPMs into single file
    echo "Joining ${MIXNAME} PPMs."
    echo ""
    cat ./${MIXPATH}/*_PPM_tmp.txt > ./${MIXPATH}/PPM.txt
    
    #Join all weights into single file
    cat ./${MIXPATH}/*_Wm_tmp.txt > ./${MIXPATH}/Wm.txt
    
    #Clear - remove all temporary files
    rm ./${MIXPATH}/*tmp.txt
  
    echo "Clearing ${MIXPATH} directory..."
    echo ""

    
    echo "${MIXNAME} mixture done."
    echo "Moving..."
    echo ""
    
   #fi
  
  #Finish "for MIXPATH"
  done
   
 echo "${DATASET} done."
 echo ""
  
done

echo "Data ready"
echo ""
echo "Lets look at the logos :-)"
echo ""



