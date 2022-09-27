#BEDtools getfasta tool was used to retrieve genomic sequences. IS BED ranges and FASTA files with sequences of chromosomes served as input. In resulting FASTA files, all nucleotides were converted to capitals with



cat ${ISSEQ}.fa | sed 's/a/A/g' | sed 's/c/C/g' | sed 's/g/G/g' | sed 's/t/T/g' > ${ISSEQ}_CAP.fa

#Sequences in the FASTA file were first converted to table, where sequences are placed into single column and each sequence occupies a single row:
grep -v ">" ${ISSEQ}_CAP.fa > ${ISSEQ}_is26.txt

#Finally the sequences were converted to numeric strings:

sed 's/[Aa]/1/g' ${ISSEQ}_is26.txt | sed 's/[Cc]/2/g' | sed 's/[Gg]/3/g' | sed 's/[Tt]/4/g' | sed 's/[Nn]/5/g' ${ISSEQ}_is26_num.txt

