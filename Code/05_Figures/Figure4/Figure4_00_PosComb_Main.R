#Compose Figure3:
#Frequencies of positional combinations

#Set virus/sample folder name
# = which row of smpl.info table
select.smpl.name <- c("HIV_Zhyvoloup", "HTLV_Kirk",
                      "MLV_DeRavin", "MVV_HEK")

select.smpl <- sapply(select.smpl.name,
                      function(x) {
                        which(smpl.info$seqName == x)},
                      USE.NAMES = FALSE)

#Create data.frame with positional combination frequencies
sp.range <- seq(-13,-1, by = 1)
smpl.title <- "virus" # whole | virus | virus_celltype
df <- as.data.frame(do.call(rbind,
                            lapply(select.smpl,
                                   function(ss) {
                                     #Load sample variables
                                     virus <- smpl.info$seqName[ss]
                                     vrs <- smpl.info$Virus[ss]
                                     ctrl <- smpl.info$Control[ss]
                                     sln <- smpl.info$seqLen[ss]
                                     celltype <- smpl.info$Cell[ss]
                                     
                                     #Target site duplication (length)
                                     tsd <- unname(retro.tsd[names(retro.tsd) == vrs])
                                     half.tsd <- floor(tsd/2)
                                     
                                     #Load sequences
                                     #Load IS
                                     is.seq <- as(read.table(paste0(data.wd, "IS/", virus, "/", virus, "_is", sln,".txt"),
                                                             sep = "\t", stringsAsFactor = FALSE, header = FALSE),
                                                  "DataFrame")
                                     
                                     is.nuc.mat <- NucMat(is.seq[,1])
                                     colnames(is.nuc.mat) <- RelPosNames(is.nuc.mat)
                                     
                                     ran.seq.name <- paste0(ctrl,"_seq_CAP")
                                     ran.seq <- as(read.table(paste0(data.wd, "IS/Ran/", ran.seq.name, ".txt"),
                                                              sep = "\t", stringsAsFactor = FALSE, header = FALSE),
                                                   "DataFrame")
                                     
                                     ran.nuc.mat <- NucMat(ran.seq[,1])
                                     colnames(ran.nuc.mat) <- RelPosNames(ran.nuc.mat)
                                     
                                     
                                     #Calculate dinucleotide frequency and enrichment (fold, KLID)
                                     sp.range <- as.numeric(colnames(is.nuc.mat)[colnames(is.nuc.mat) < 0])
                                     dinuc.freq <- ListToDF(DinucAtPosFreq(sp.range, sample.seq.mat = is.nuc.mat, control.seq.mat = ran.nuc.mat),
                                                            method = "dinuc")
                                     
                                     #Add virus-specific columns to data:
                                     ##Virus name
                                     if(smpl.title == "whole") {dinuc.freq$Virus <- smpl.info$seqName[ss]}
                                     if(smpl.title == "virus") {
                                       if(celltype == "invitro") {
                                         dinuc.freq$Virus <- paste0(vrs,"iv")
                                       } else {
                                         dinuc.freq$Virus <- vrs
                                       }
                                     }
                                     if(smpl.title == "virus_celltype") {dinuc.freq$Virus <- paste0(vrs, "_", celltype)}
                                     
                                     ##position relative to cleavage site
                                     cs.rel.pos <- as.numeric(levels(dinuc.freq$Pos)[dinuc.freq$Pos]) + half.tsd
                                     cs.rel.pos[cs.rel.pos >= 0] <- cs.rel.pos[cs.rel.pos >= 0] + 1
                                     dinuc.freq$tsdPos <- cs.rel.pos
                                     
                                     dinuc.freq
                                   })))

#Create panel 1: positional KLID
source(paste0(scr.wd, "Figure4/Figure4_01_KLID_dinucleotide.R"))

#Create panel 2: positional FREQuencies
source(paste0(scr.wd, "Figure4/Figure4_02_FREQ_dinucleotide.R"))

#Create panel 3: table of frequencies
library(png)
library(grid)
img.2 <- readPNG(paste0(result.wd, "Figures/Figure4_PositionalCombs_table.png"), native = TRUE)
img.2 <- rasterGrob(img.2, interpolate=TRUE)
p0 <- ggplot() + theme_void() + theme(plot.background = element_rect(fill = "white", colour =  "white"))
t1 <- p0 +
  annotation_custom(img.2, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf)
  
#Create panel 4: Most and least frequent nucleotides
source(paste0(scr.wd, "Figure4/Figure4_04_Most-Least_Freq.R"))


#Final plot
p.pos.comb <- ggarrange(ggarrange(p3.a, p0, p3.b,
                                  labels = c("A", "", "B"),
                                  ncol = 3, widths = c(1.5, 0.1, 2.5)),
                        ggarrange(t1, p3.d, p0,
                                  labels = c("C", "D", ""),
                                  widths = c(c(1, 4, 0.1)),
                                  nrow = 1),
                        nrow = 2, heights = c(2.2,1))

#Save figure
#A4 = 210 / 297 mm
a4w <- 210
a4h <- 297
ggsave("Figure4_PositionalCombs.png", p.pos.comb,
       path = paste0(result.wd, "Figures/"),
       width = a4w, height = a4h * 4/5, units = "mm", dpi = "retina",
       device = "png", limitsize = FALSE, bg = "white")
