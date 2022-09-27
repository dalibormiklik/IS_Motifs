#Figure4 Part_III
#Plot HIV-1 targeting of Alu and intra-Alu motifs
#retargeted variants + IN mutants
#HIV-1 IS from 3 different studies

#Panels G, H

#Panel G
#Frequency in Alu repeats
source(paste0(scrd,"Figure4_Alu_Motif_Target_03_G_AluFreq.R"))

#Panel H
#Figure4 S1
#Frequency in intra-Alu palindromic motif
source(paste0(scrd,"Figure4_Alu_Motif_Target_03_H_MotifDist.R"))

#Frequency of motifs in intra-Alu IS
#Panel I - frequencies of top motifs in selected samples
#Supplementary Figure with full list of variants
#Panel I of the main Figure  is created as part of the "S2" code
source(paste0(scrd,"Figure4_Alu_Motif_Target_S2_MotifFreq.R"))

F4.pIII <- ggarrange(ggarrange(p0, p.freq.alu,
                    p0, p.freq.in.motif, p0,
                    nrow = 1, widths = c(p0.w , 1, p0.w / 2, 1, p0.w),
                    labels = c("", "G", "", "H", ""), hjust = 1),
                    p.m.freq.sel,
                    nrow = 2, heights = c(1, 1),
                    labels = c("", "I"))
