#Produce Supplementary Figure S8:
#Comparing ASLV data sets

#Compare nucleotide combination KLID values
source(paste0(scr.wd, "FigureS8/FigureS8_PosCombs_KLID.R"))

#Logos of Alpha_Moiani M08
source(paste0(scr.wd, "FigureS8/FigureS8_AlphaMoiani_Logo.R"))

#Frequences of TGAA and TTCA sequence motifs at tDNA strands
source(paste0(scr.wd, "FigureS8/FigureS8_AlphaMoiani_MotifFreq.R"))

#Empty plot
p0 <- ggplot() + theme_void() + theme(plot.background = element_rect(fill = "white", colour =  "white"))


p.s8.aslv <- ggarrange(ggarrange(p.aslv.a,
                                 ggarrange(p0, p0, p0,
                                           p0, p.alpha.mot, p0,
                                           p0, p0, p0,
                                           ncol = 3, nrow = 3,
                                           widths = c(0.5, 3, 1.5), heights = c(0.25, 5, 1)#,
                                           #labels = c(rep("", 4), "C", rep("", 4))
                                           ),
                                 nrow = 2, heights = c(2,1),
                                 labels = c("A", "C")
                                 ),
                       p8.logo,
                       ncol = 2, widths = c(2,1),
                       labels = c("", "B")
                       )

#Save figure:
#Save figure
#A4 = 210 / 297 mm
a4w <- 210
a4h <- 297
ggsave("FigureS8_ASLVs.png", p.s8.aslv,
       path = paste0(result.wd, "Figures/"),
       width = a4w, height = a4h * 4/5, units = "mm", dpi = "retina",
       device = "png", limitsize = FALSE, bg = "white")
