#This code runs "Main" scripts creating individual Figures
#Prior this code, '00_General.R' code should be run to load necessary functions etc

#Figure1 - PDef & sequence logos
source(paste0(scr.wd, "Figure2/Figure2_00_SeqLogo_Main.R"))

#Figure2 - Sequence logos of full mixtures
source(paste0(scr.wd, "Figure3/Figure3_00_SeqLogo_Main.R"))
#FigureS2 - Sequence logos of all full mixtures
source(paste0(scr.wd, "Figure3/Figure3_S2_SeqLogo_All.R"))

#Figure3 - Positional nucleotide combinations
source(paste0(scr.wd, "Figure4/Figure4_00_PosComb_Main.R"))

#Figure4 - Characterization of HIV-1 hotspot
source(paste0(scr.wd, "Figure5/Figure5_Alu_Motif_Target_00_Main.R"))
