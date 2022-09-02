#Create supplementary figure S1
#All PPM00, most palindromic and closest RC-pairs
#This code is run after the code for main Fig.1 where other variables and objects are set 

#Single code setting full ranges for logos from Fig1C,D should be run

#Limits for plotting positions
#pos.rng can be "all" or vector of boundary positions
pos.rng <- c(-13,13)
source(paste0(scr.wd, "Figure1/Figure1_03_SeqLogo_Pal.R"))


#Create supplementary figure w/ all IS sets:
p.seq.logo.sup <- ggarrange(p1.a, p0,
                            p1.pdef.pannel, p0,
                            p1.pdef.max.rc.min.pannel,
                            nrow = 5,
                            heights = c(3, p0.w, 3, p0.w, 6),
                            labels = c("A", "", "B", "", "C"))

#Save figure
#A4 = 210 / 297 mm
a4w <- 210
a4h <- 297
ggsave("FigureS1_SeqLogo.png", p.seq.logo.sup,
       path = paste0(result.wd,"/Figures"),
       width = a4w, height = a4h * 4/5, units = "mm", dpi = "retina",
       device = "png", limitsize = FALSE)
