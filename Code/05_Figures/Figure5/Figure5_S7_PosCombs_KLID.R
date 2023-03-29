#Create supplementary figure showing positional combinations of HIV IN variants

#Load sample info
vrs <- "HIV"
smpl.pattern <- paste0(vrs, "_Demeulemeester") 
ctrl <- "Ran9k_hg19_is27"
sln <- 26
tsd <- 5

#Select directory names containing sample pattern in name
select.smpl <- dir(path = sq.wd,
                   full.names = FALSE, recursive = FALSE,
                   pattern = smpl.pattern)

#Create data.frame with positional combination frequencies
sp.range <- seq(-13,-1, by = 1)
df.hivinvar <- as.data.frame(do.call(rbind,
                            lapply(1:length(select.smpl),
                                   function(ss) {
                                     #Load sample variables
                                     virus <- select.smpl[ss]
                                     smpl <- unlist(strsplit(virus, "_"))
                                     smpl <- paste(smpl[3:length(smpl)], collapse = "_")
                                     
                                     #Target site duplication (length)
                                     half.tsd <- floor(tsd/2)
                                     
                                     #Load sequences
                                     #Load IS
                                     is.seq <- as(read.table(paste0(sq.wd, virus, "/", virus, "_is", sln,".txt"),
                                                             sep = "\t", stringsAsFactor = FALSE, header = FALSE),
                                                  "DataFrame")
                                     
                                     is.nuc.mat <- NucMat(is.seq[,1])
                                     colnames(is.nuc.mat) <- RelPosNames(is.nuc.mat)
                                     
                                     ran.seq.name <- paste0(ctrl,"_seq_CAP")
                                     ran.seq <- as(read.table(paste0(sq.wd, "Ran/", ran.seq.name, ".txt"),
                                                              sep = "\t", stringsAsFactor = FALSE, header = FALSE),
                                                   "DataFrame")
                                     
                                     ran.nuc.mat <- NucMat(ran.seq[,1])
                                     colnames(ran.nuc.mat) <- RelPosNames(ran.nuc.mat)
                                     ran.nuc.mat <- ran.nuc.mat[,which(colnames(ran.nuc.mat) %in% colnames(is.nuc.mat))]
                                     
                                     #Calculate dinucleotide frequency and enrichment (fold, KLID)
                                     sp.range <- as.numeric(colnames(is.nuc.mat)[colnames(is.nuc.mat) < 0])
                                     dinuc.freq <- ListToDF(DinucAtPosFreq(sp.range, sample.seq.mat = is.nuc.mat, control.seq.mat = ran.nuc.mat),
                                                            method = "dinuc")
                                     
                                     #Add virus-specific columns to data:
                                     ##Virus name
                                     if(smpl.title == "whole") {dinuc.freq$Virus <- smpl.info$Name[ss]}
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
                                     
                                     dinuc.freq$Sample <- smpl
                                     
                                     dinuc.freq
                                     
                                   })))

#Set the levels for sample
#Set Virus column as column with sample names
lvls.smpl <- unique(df.hivinvar$Sample)
lvls.smpl <- c("WT",
               lvls.smpl[grep("^R231",lvls.smpl)],
               lvls.smpl[grep("S119",lvls.smpl)])
lvls.smpl <- lvls.smpl[-grep("_", lvls.smpl)]
lvls.virus <- c(lvls.smpl, unique(df.hivinvar$Sample)[grep("_", unique(df.hivinvar$Sample))])

#Create data.frame with data used for plotting
select.pos <- c(-5:2)
data.to.plot <- df.hivinvar[df.hivinvar$tsdPos %in% select.pos,]

data.to.plot$Virus <- factor(data.to.plot$Sample, levels = lvls.virus)
data.to.plot$tsdPos <- factor(data.to.plot$tsdPos, levels = select.pos)
data.to.plot$N1 <- factor(data.to.plot$N1, levels = nucs)
data.to.plot$N2 <- factor(data.to.plot$N2, levels = nucs)

klid <- as.data.frame(do.call(rbind,
                              lapply(levels(data.to.plot$Virus),
                                     function(lvls.v) {
                                       #Derive (per postion) KLID
                                       k <- sapply(lvls.pos,
                                                   function(p) {
                                                     s <- sum(data.to.plot$KLID[data.to.plot$tsdPos == p &
                                                                                  data.to.plot$Virus == lvls.v])
                                                     if(length(s) == 0) {s <- c()}
                                                     s
                                                   })
                                       
                                       df <- data.frame(tsdPos = lvls.pos,
                                                        KLID = k,
                                                        Virus = lvls.v)
                                       df[df == 0] <- NA
                                       df
                                     })))
klid <- klid[!is.na(klid$KLID),]
klid$tsdPos <- factor(klid$tsdPos, levels = select.pos)
klid$Virus <- factor(klid$Virus, levels = lvls.virus)

#PLOT
p.klid <- 
  ggplot(data = data.to.plot) +
  geom_col(data = klid, aes(x = tsdPos, y = KLID, group = Virus),
           fill = "gray90", color = "gray40", width = 0.75, na.rm = TRUE) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = abs(min(select.pos)) + 0.5, color = "black", lty = 2) +
  scale_shape_manual(values=c(21:25)) +
  scale_colour_manual(values = nuc.cols) +
  geom_quasirandom(aes(x = tsdPos, y = KLID, group = Virus, color = N2, shape = N1),
                   size = 2, stroke = 1.5, fill = alpha("gray",.5)) +
  geom_text(aes(label = Virus), x = 1.75, y = 55, fontface = "bold", size = 3) +
  xlab("Strand transfer site-relative position") +
  #scale_x_discrete(labels = lvls.pos) +
  guides(shape = guide_legend(order=1),
         color = guide_legend(order=2)) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0, colour = "black", face = "bold"),
        axis.line.x = element_blank(),
        axis.title.x = element_text(colour = "black", face = "bold", size=10, angle = 0),
        axis.title.y = element_text(colour = "black", face = "bold", size=10),
        axis.text.x = element_text(colour = "black", face = "bold", size=10, angle = 0),
        axis.text.y = element_text(colour = "black", face = "bold", size=10, angle = 0),
        legend.position=c("top"), legend.box = "vertical",
        legend.title = element_text(colour = "black", size=9, face="bold"),
        legend.text = element_text(colour = "black", size=8, face="bold"),
        legend.spacing = unit(0, units = "points"),
        strip.background = element_blank(),
        strip.text = element_blank(),
        rect = element_rect(fill = "transparent"),
        panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "white", color = "white")
  ) +
  facet_wrap("Virus", ncol = 3, scales = "fixed", strip.position = "right")

#Save figure
#A4 = 210 / 297 mm
a4w <- 210
a4h <- 297
ggsave("FigureS7_HIV-INvars_PositionalCombs_KLID.png", p.klid,
       path = paste0(result.wd, "Figures/"),
       width = a4w, height = a4h * 2/3, units = "mm", dpi = "retina",
       device = "png", limitsize = FALSE, bg = "white")
