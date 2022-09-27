#Scripts for selecting particular repeats from repeatmasker and transform outut to BED format
#Scripts are performed in the directory where rmsk.txt file downloaded from UCSC Goldenpath is saved

#all repeats to BED
awk 'BEGIN{OFS="\t"} {print $6,$7,$8,$11,$2,$10}' rmsk.txt > rmsk.bed

#Alu
awk 'BEGIN{OFS="\t"} $13 == "Alu" {print $6,$7,$8,$11,$2,$10}' rmsk.txt > rmsk_Alu.bed

