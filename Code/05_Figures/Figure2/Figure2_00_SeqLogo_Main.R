#Create KLID logos for components

#Set virus/sample folder name
# = which row of smpl.info table
select.smpl.name <- c("HIV_Zhyvoloup", "HTLV_Kirk",
                      "MLV_DeRavin", "MVV_HEK",
                      "PFV_WT_HT1080", "PFV_invitro",
                      "RCASC_CEF", "Random")

select.smpl <- sapply(select.smpl.name,
                      function(x) {
                        which(smpl.info$seqName == x)},
                      USE.NAMES = FALSE)

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

source(paste0(scr.wd, "Figure1/Figure2_01_SeqLogo_PPM0.R"))

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
source(paste0(scr.wd, "Figure2/Figure2_02_PPMdist.R"))
p1.b <- min.dist.p

#Fig.1 WITHOUT component sequence logos
p.seq.logo <- ggarrange(p0, p1.a, p0, p1.b,
                        nrow = 4,
                        heights = c(0.75, 3, p0.w, 2),
                        labels = c("A", "", "", "B"))

#Save figure
#A4 = 210 / 297 mm
a4w <- 210
a4h <- 297
ggsave("Figure2_SeqLogo.png", p.seq.logo,
       path = paste0(result.wd,"/Figures"),
       width = a4w, height = a4h * 2/5, units = "mm", dpi = "retina",
       device = "png", limitsize = FALSE)

#Create full supplementary figure1
source(paste0(scr.wd, "Figure2/Figure2_S1_SeqLogo.R"))
