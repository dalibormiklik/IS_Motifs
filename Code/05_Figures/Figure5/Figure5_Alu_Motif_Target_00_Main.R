#Create plots describing identification of HIV-1 integration hotspot in Alu repeats:

#Part_I: identification of intra-Alu hotspot
#Part_II: Frequency of integration into hotspot in the whole HIV-1 IS set
#Part_III: Hotspot targeting in retargeted viruses and IN mutants

#Set directory for scripts
scrd <- paste0(scr.wd, "Figure4/")

#Empty plot
p0 <- ggplot() + theme_void() + theme(plot.background = element_rect(fill = "white", colour =  "white"))
p0.w <- 0.2

#Part_I
source(paste0(scrd, "Figure5_Alu_Motif_Target_01.R"))

#Part_II
source(paste0(scrd, "Figure5_Alu_Motif_Target_02.R"))

#Part_III
source(paste0(scrd, "Figure5_Alu_Motif_Target_03.R"))

F4 <- ggarrange(p0, F4.pI, F4.pII, F4.pIII,
                nrow = 4, heights = c(0.1, 1.2, 1, 2.5))

#Save figure
#A4 = 210 / 297 mm
a4w <- 210
a4h <- 297
ggsave("Figure4_HIV_Alu_Motif_Target.png", F4,
       path = paste0(result.wd, "Figures"),
       width = a4w, height = ((a4h * 2/3) / 5) * 4, units = "mm", dpi = "retina",
       device = "png", limitsize = FALSE)

#Supplementary:
#Create Figure S7 where KLID of nucleotide combinations are plotted
source(paste0(scrd, "Figure5_S7_PosCombs_KLID.R"))
