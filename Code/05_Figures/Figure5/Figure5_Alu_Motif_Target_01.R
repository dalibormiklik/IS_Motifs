#Part_I of Figure4
#Characterization of HIV-1 M08 PPM04-associated sequences

#Specify IS set and mixture component
virus <- "HIV_Zhyvoloup"
cm <- "PPM07"
mixture <- "M08"
cm.ppm <- paste0(mixture, "_", cm)
gnm <- "hg38"
repeatname <- "Alu"

#Panels A,B,C

#Panel A
#Frequency in repeats
source(paste0(scrd, "Figure5_Alu_Motif_Target_01_A_RepeatTarget.R"))

#Panel B
#Distribution in Alu consensus
source(paste0(scrd, "Figure5_Alu_Motif_Target_01_B_ConsPos.R"))

#Panel C
#Frequency of sequence motifs in M08_PPM07 component
source(paste0(scrd, "Figure5_Alu_Motif_Target_01_C_CompSeq.R"))

#Create Part_I plot
F4.pI <- ggarrange(p0, p.cmp.rmsk, p0, p.cons, p0, p.mfreq, p0,
                   nrow = 1, widths = c(p0.w, 1.5, p0.w, 3, p0.w, 2, p0.w),
                   labels = c("", "A", "", "B", "", "C"), hjust = 1)
