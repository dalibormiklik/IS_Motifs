#Figure4 Part_II
#Plot HIV-1 wt targeting of Alu and intra-Alu motifs

#First, load and calculate frequency data
#Motif sequence is defined here
source(paste0(scrd,"Figure5_Alu_Motif_Target_02_01_GetData.R"))

#Create Panels D, E, F

#Panel D
#Frequency of HIV-1 wt IS in Alu
source(paste0(scrd,"Figure5_Alu_Motif_Target_02_D_AluFreq.R"))

#Panel E
#Frequency of HIV-1 wt IS sequence motifs
source(paste0(scrd,"Figure5_Alu_Motif_Target_02_E_MotifFreq.R"))

#Panel F
#Frequency of HIV wt/N74D IS in Alu
source(paste0(scrd,"Figure5_Alu_Motif_Target_02_F_MotifDistance.R"))

F4.pII <- ggarrange(p0, ggarrange(p.freq.alu, p0,
                                  nrow = 2, heights = c(3,1)),
                    p0, p.motif.alu, p0, p.dist.motif,
                    nrow = 1, widths = c(p0.w, 1, p0.w, 2, p0.w, 2.8),
                    labels = c("", "D", "", "E", "", "F"), hjust = 1)
