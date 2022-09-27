#Create KLID logos for components

#Set virus/sample folder name
# = which row of smpl.info table
select.smpl <- c(7,2,11,15,13,12,5,16)

#select.ppm = fill particular PPM ("PPM00", "PPM01", ...) or "mix" for all PPMs of mixture
#select.mix = "all" (all mixures) | name of mixture that will be selected
#include.ppm0 = Set if PPM0 or mixture PPMs will be shown
#logo.columns = divide graphical output into # of columns
select.ppm <- "PPM00"
#mix.name <- "M02"
logo.columns <- 1
smpl.title <- "virus"

#Limits for plotting positions
#pos.rng can be "all" or vector of boundary positions
pos.rng <- c(-8,8)
ylim <- c(-0.5, 1.6)

#Empty plot
p0 <- ggplot() + theme_void() + theme(plot.background = element_rect(fill = "white", colour =  "white"))
p0.w <- 0.2

select.mix <- "M02"
mix.name <- select.mix
include.ppm0 <- TRUE

source(paste0(scr.wd, "Figure1/Figure1_01_SeqLogo_PPM0.R"))

#Set theme for plots w/out y axis
theme_y_blank <- theme(axis.text.y = element_blank(),
                       axis.title.y = element_blank(),
                       axis.line.y = element_blank(),
                       axis.ticks.y = element_blank())


p1.a <- ggarrange(ggarrange(p1.list[[1]],
                            p1.list[[2]] + theme_y_blank,
                            p1.list[[3]] + theme_y_blank,
                            p1.list[[4]] + theme_y_blank,
                            ncol = 4, nrow = 1),
                  ggarrange(p1.list[[5]],
                            p1.list[[6]] + theme_y_blank,
                            p1.list[[7]] + theme_y_blank,
                            p1.list[[8]] + theme_y_blank,
                            ncol = 4, nrow = 1),
                  nrow = 2, heights = c(3,3))

#Create Panel B with PPM distances
source(paste0(scr.wd, "Figure1/Figure1_02_PPMdist.R"))
p1.b <- min.dist.p

#Set more narrow range of positions
pos.rng <- c(-8,8)
#Panel C: most palindromic components
source(paste0(scr.wd, "Figure1/Figure1_03_SeqLogo_Pal.R"))

p1.c <- ggarrange(p1.pdef.list[[1]],
                            p1.pdef.list[[2]] + theme_y_blank,
                            p1.pdef.list[[3]] + theme_y_blank,
                            p1.pdef.list[[4]] + theme_y_blank,
                            ncol = 4, nrow = 1)

pos.rng <- c(-8,8)
#Panel D: closest RC-pairs of components
p1.d <- ggarrange(p1.pdef.max.rc.min.list[[1]],
                  p1.pdef.max.rc.min.list[[2]] + theme_y_blank,
                  p1.pdef.max.rc.min.list[[3]] + theme_y_blank,
                  p1.pdef.max.rc.min.list[[4]] + theme_y_blank,
                  ncol = 4, nrow = 1)

#Compose figure 1
#Fig.1 WITH component sequence logos
#p.seq.logo <- ggarrange(p0, p1.a, p0, p1.b, p0, p1.c, p0, p1.d,
#                        nrow = 8,
#                        heights = c(0.75, 3, p0.w, 2, p0.w, 1.5, p0.w, 3),
#                        labels = c("A", "", "", "B", "", "C", "", "D"))

#Fig.1 WITHOUT component sequence logos
p.seq.logo <- ggarrange(p0, p1.a, p0, p1.b,
                        nrow = 4,
                        heights = c(0.75, 3, p0.w, 2),
                        labels = c("A", "", "", "B"))

#Save figure
#A4 = 210 / 297 mm
a4w <- 210
a4h <- 297
ggsave("Figure1_SeqLogo.png", p.seq.logo,
       path = paste0(result.wd,"/Figures"),
       width = a4w, height = a4h * 2/5, units = "mm", dpi = "retina",
       device = "png", limitsize = FALSE)

#Create full supplementary figure1
source(paste0(scr.wd, "Figure1/Figure1_S1_SeqLogo.R"))
