#Create bar plots showing the most and the least frequent nucleotide per STR-relative position

select.pos <- c(-3:-1, 1:3)
sel.vir <- c("HIV", "HTLV", "MLV", "MVV")

lvls.virus <- sel.vir

#Define functions

#Select "most" or "least" frequent nucleotide per position
most.or.least <- "most"

if(most.or.least == "most") {FreqDF <- MostFreqDF}
if(most.or.least == "least") {FreqDF <- LeastFreqDF}

df.2 <- as.data.frame(do.call(rbind,
                              lapply(select.smpl,
                                     function(act.smpl) {
                                       #Load sample variables
                                       seqname <- smpl.info$seqName[act.smpl]
                                       vrs <- smpl.info$Virus[act.smpl]
                                       gnm <- smpl.info$Genome[act.smpl]
                                       sln <- smpl.info$seqLen[act.smpl]
                                       ctrl.name <- smpl.info$Control[act.smpl]
                                       celltype <- smpl.info$Cell[act.smpl]
                                       
                                       #Target site duplication (length)
                                       tsd <- unname(retro.tsd[names(retro.tsd) == vrs])
                                       half.tsd <- floor(tsd/2)
                                       
                                       #Load sequences
                                       #Load IS
                                       is.seq <- as(read.table(paste0(sq.wd, seqname, "/", seqname, "_is", sln,".txt"),
                                                               sep = "\t", stringsAsFactor = FALSE, header = FALSE),
                                                    "DataFrame")
                                       
                                       is.nuc.mat <- NucMat(is.seq[,1])
                                       colnames(is.nuc.mat) <- RelPosNames(is.nuc.mat)
                                       
                                       ran.seq.name <- paste0(ctrl.name, "_seq_CAP")
                                       ran.seq <- as(read.table(paste0(sq.wd, "Ran/", ran.seq.name, ".txt"),
                                                                sep = "\t", stringsAsFactor = FALSE, header = FALSE),
                                                     "DataFrame")
                                       
                                       ran.nuc.mat <- NucMat(ran.seq[,1])
                                       colnames(ran.nuc.mat) <- RelPosNames(ran.nuc.mat)
                                       
                                       #Calculate dinucleotide frequency and enrichment (fold, KLID)
                                       #Select maximum of "x.ranges" positions solist matrices have same number of rows
                                       x.ranges <- 13
                                       sp.range <- as.numeric(colnames(is.nuc.mat)[colnames(is.nuc.mat) <= 0])
                                       sp.range <- sp.range[(1 + abs(x.ranges - length(sp.range))):length(sp.range)]
                                       flist <- NucAtPosFreq(sp.range, nucs, sample.seq.mat = is.nuc.mat, control.seq.mat = ran.nuc.mat, output.number = "count")
                                       dinuc.freq <- FreqDF(flist)
                                       
                                       #Format table
                                       dinuc.freq$Pos <- factor(dinuc.freq$Pos, levels = sp.range)
                                       
                                       #Add virus-specific columns to data:
                                       ##Virus name
                                       if(celltype == "invitro") {
                                         dinuc.freq$Virus <- paste0(vrs,"iv")
                                       } else {
                                         dinuc.freq$Virus <- vrs
                                       }
                                       ##position relative to cleavage site
                                       cs.rel.pos <- as.numeric(levels(dinuc.freq$Pos)[dinuc.freq$Pos]) + half.tsd
                                       cs.rel.pos[cs.rel.pos >= 0] <- cs.rel.pos[cs.rel.pos >= 0] + 1
                                       dinuc.freq$tsdPos <- cs.rel.pos
                                       
                                       #set unique names for rows
                                       rownames(dinuc.freq) <- paste0(dinuc.freq$tsdPos, dinuc.freq$N1, dinuc.freq$N2,
                                                                      "_", dinuc.freq$Virus)
                                       
                                       dinuc.freq
                                     })))

#Create supplementary figure where all samples are plotted
source(paste0(scr.wd, "Figure3/Figure3_S4_Most-Least_Freq.R"))
p3.s1 <- p3.sup.freq

select.pos <- lvls.pos[between(lvls.pos, -3, 3)]
data.to.plot <- df.2[df.2$tsdPos %in% select.pos &df.2$Virus %in% sel.vir,]

data.to.plot$IS_prop <- unlist(data.to.plot$IS_prop)
data.to.plot$Ctrl_prop <- unlist(data.to.plot$Ctrl_prop)
data.to.plot$N2_prop <- unlist(data.to.plot$N2_prop)

data.to.plot$Virus <- factor(data.to.plot$Virus, levels = lvls.virus)
data.to.plot$tsdPos <- factor(data.to.plot$tsdPos, levels = select.pos)
data.to.plot$N1 <- factor(data.to.plot$N1, levels = names(nuc.cols))
data.to.plot$N2 <- factor(data.to.plot$N2, levels = names(nuc.cols))

uniq.data <- data.to.plot[!duplicated(paste0(data.to.plot$Virus, data.to.plot$tsdPos)),]
uniq.data$TextPos <- sapply(1:nrow(uniq.data),
                            function(x) {
                              prop.is <- uniq.data$IS_prop[x]
                              prop.ctrl <- uniq.data$Ctrl_prop[x]
                              if(prop.is > prop.ctrl) {
                                prop.is + 12
                              } else {
                                prop.ctrl + 12
                              }
                            })

data.to.plot.m <- data.to.plot
uniq.data.m <- uniq.data

p3.m <- 
  ggplot(data.to.plot.m, aes(x = tsdPos)) +
  geom_hline(yintercept = c(25, 50, 75, 100), colour = "gray60", linetype = 3) +
  geom_bar(data = uniq.data.m, aes(y = IS_prop, fill = N1), colour = NA, alpha = 0.75, stat = "identity") +
  geom_bar(aes(y = N2_prop, fill = N2), colour = NA, alpha = 1, stat = "identity") +
  geom_point(data = uniq.data.m, aes(y = Ctrl_prop)) +
  geom_vline(xintercept = abs(min(select.pos)) + 0.5, color = "black", lty = 2) +
  geom_text(data = uniq.data.m, aes(label = N1), y = uniq.data.m$TextPos, fontface = "bold", size = 3.5) +
  scale_colour_manual(values = nuc.cols) +
  scale_fill_manual(values = nuc.cols) +
  geom_hline(yintercept = 0) +
  coord_cartesian(ylim = c(0, 120)) +
  scale_y_continuous(breaks = seq(0, 100, 25), expand = c(0, 0)) +
  ggtitle("Most frequent") +
  xlab("STR site-relative position") +
  ylab("Nucleotide\nfrequency [%]") +
  guides(fill=guide_legend(title="RC\npartner")) +
  theme_classic() +
  theme(plot.title = element_blank(),
        axis.line.x = element_blank(),
        axis.title.x = element_text(colour = "black", face = "bold", size=8, angle = 0),
        axis.title.y = element_blank(),
        axis.text.x = element_text(colour = "black", face = "bold", size=8, angle = 0),
        axis.text.y = element_text(colour = "black", face = "bold", size=8, angle = 0),
        legend.position="right",
        legend.title = element_text(colour="black", size=8, face="bold"),
        legend.text = element_text(colour="black", size=8, face="bold"),
        legend.key.size = unit(8, "points"),
        panel.grid = element_line(colour = "gray"),
        strip.background = element_blank(),
        strip.text = element_text(colour = "black", size=10, face="bold"),
        strip.placement = "outside",
        #rect = element_rect(fill = "transparent"),
        #panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "white", color = "white")
  ) +
  facet_wrap(~Virus, nrow = 1)


#Select "most" or "least" frequent nucleotide per position
select.pos <- c(-3:-1, 1:3)
sel.vir <- c("HIV", "HTLV", "MLV", "MVV")

most.or.least <- "least"

if(most.or.least == "most") {FreqDF <- MostFreqDF}
if(most.or.least == "least") {FreqDF <- LeastFreqDF}

df.2 <- as.data.frame(do.call(rbind,
                              lapply(select.smpl,
                                     function(act.smpl) {
                                       #Load sample variables
                                       seqname <- smpl.info$seqName[act.smpl]
                                       vrs <- smpl.info$Virus[act.smpl]
                                       gnm <- smpl.info$Genome[act.smpl]
                                       sln <- smpl.info$seqLen[act.smpl]
                                       ctrl.name <- smpl.info$Control[act.smpl]
                                       celltype <- smpl.info$Cell[act.smpl]
                                       
                                       #Target site duplication (length)
                                       tsd <- unname(retro.tsd[names(retro.tsd) == vrs])
                                       half.tsd <- floor(tsd/2)
                                       
                                       #Load sequences
                                       #Load IS
                                       is.seq <- as(read.table(paste0(sq.wd, seqname, "/", seqname, "_is", sln,".txt"),
                                                               sep = "\t", stringsAsFactor = FALSE, header = FALSE),
                                                    "DataFrame")
                                       
                                       is.nuc.mat <- NucMat(is.seq[,1])
                                       colnames(is.nuc.mat) <- RelPosNames(is.nuc.mat)
                                       
                                       ran.seq.name <- paste0(ctrl.name, "_seq_CAP")
                                       ran.seq <- as(read.table(paste0(sq.wd, "Ran/", ran.seq.name, ".txt"),
                                                                sep = "\t", stringsAsFactor = FALSE, header = FALSE),
                                                     "DataFrame")
                                       
                                       ran.nuc.mat <- NucMat(ran.seq[,1])
                                       colnames(ran.nuc.mat) <- RelPosNames(ran.nuc.mat)
                                       
                                       #Calculate dinucleotide frequency and enrichment (fold, KLID)
                                       #Select maximum of "x.ranges" positions solist matrices have same number of rows
                                       x.ranges <- 13
                                       sp.range <- as.numeric(colnames(is.nuc.mat)[colnames(is.nuc.mat) <= 0])
                                       sp.range <- sp.range[(1 + abs(x.ranges - length(sp.range))):length(sp.range)]
                                       flist <- NucAtPosFreq(sp.range, nucs, sample.seq.mat = is.nuc.mat, control.seq.mat = ran.nuc.mat, output.number = "count")
                                       dinuc.freq <- FreqDF(flist)
                                       
                                       #Format table
                                       dinuc.freq$Pos <- factor(dinuc.freq$Pos, levels = sp.range)
                                       
                                       #Add virus-specific columns to data:
                                       ##Virus name
                                         if(celltype == "invitro") {
                                           dinuc.freq$Virus <- paste0(vrs,"iv")
                                         } else {
                                           dinuc.freq$Virus <- vrs
                                         }
                                       ##position relative to cleavage site
                                       cs.rel.pos <- as.numeric(levels(dinuc.freq$Pos)[dinuc.freq$Pos]) + half.tsd
                                       cs.rel.pos[cs.rel.pos >= 0] <- cs.rel.pos[cs.rel.pos >= 0] + 1
                                       dinuc.freq$tsdPos <- cs.rel.pos
                                       
                                       #set unique names for rows
                                       rownames(dinuc.freq) <- paste0(dinuc.freq$tsdPos, dinuc.freq$N1, dinuc.freq$N2,
                                                                      "_", dinuc.freq$Virus)
                                       
                                       dinuc.freq
                                     })))

#Create supplementary figure where all samples are plotted
source(paste0(scr.wd, "Figure3/Figure3_S4_Most-Least_Freq.R"))
p3.s2 <- p3.sup.freq

select.pos <- lvls.pos[between(lvls.pos, -3, 3)]
data.to.plot <- df.2[df.2$tsdPos %in% select.pos &df.2$Virus %in% sel.vir,]

data.to.plot$IS_prop <- unlist(data.to.plot$IS_prop)
data.to.plot$Ctrl_prop <- unlist(data.to.plot$Ctrl_prop)
data.to.plot$N2_prop <- unlist(data.to.plot$N2_prop)

data.to.plot$Virus <- factor(data.to.plot$Virus, levels = lvls.virus)
data.to.plot$tsdPos <- factor(data.to.plot$tsdPos, levels = select.pos)
data.to.plot$N1 <- factor(data.to.plot$N1, levels = names(nuc.cols))
data.to.plot$N2 <- factor(data.to.plot$N2, levels = names(nuc.cols))

uniq.data <- data.to.plot[!duplicated(paste0(data.to.plot$Virus, data.to.plot$tsdPos)),]
uniq.data$TextPos <- sapply(1:nrow(uniq.data),
                            function(x) {
                              prop.is <- uniq.data$IS_prop[x]
                              prop.ctrl <- uniq.data$Ctrl_prop[x]
                              if(prop.is > prop.ctrl) {
                                prop.is + 12
                              } else {
                                prop.ctrl + 12
                              }
                            })

data.to.plot.l <- data.to.plot
uniq.data.l <- uniq.data

p3.l <- 
  ggplot(data.to.plot.l, aes(x = tsdPos)) +
  geom_hline(yintercept = c(25, 50, 75, 100), colour = "gray60", linetype = 3) +
  geom_bar(data = uniq.data.l, aes(y = IS_prop, fill = N1), colour = NA, alpha = 0.75, stat = "identity") +
  geom_bar(aes(y = N2_prop, fill = N2), colour = NA, alpha = 1, stat = "identity") +
  geom_point(data = uniq.data.l, aes(y = Ctrl_prop)) +
  geom_vline(xintercept = abs(min(select.pos)) + 0.5, color = "black", lty = 2) +
  geom_text(data = uniq.data.l, aes(label = N1), y = uniq.data.l$TextPos, fontface = "bold", size = 3.5) +
  scale_colour_manual(values = nuc.cols) +
  scale_fill_manual(values = nuc.cols) +
  geom_hline(yintercept = 0) +
  coord_cartesian(ylim = c(0, 120)) +
  scale_y_continuous(breaks = seq(0, 100, 25), expand = c(0, 0)) +
  ggtitle("Most frequent") +
  xlab("STR site-relative position") +
  ylab("Nucleotide\nfrequency [%]") +
  guides(fill=guide_legend(title="RC\npartner")) +
  theme_classic() +
  theme(plot.title = element_blank(),
        axis.line.x = element_blank(),
        axis.title.x = element_text(colour = "black", face = "bold", size=8, angle = 0),
        axis.title.y = element_blank(),
        axis.text.x = element_text(colour = "black", face = "bold", size=8, angle = 0),
        axis.text.y = element_text(colour = "black", face = "bold", size=8, angle = 0),
        legend.position="right",
        legend.title = element_text(colour="black", size=8, face="bold"),
        legend.text = element_text(colour="black", size=8, face="bold"),
        legend.key.size = unit(8, "points"),
        panel.grid = element_line(colour = "gray"),
        strip.background = element_blank(),
        strip.text = element_text(colour = "black", size=10, face="bold"),
        strip.placement = "outside",
        plot.background = element_rect(fill = "white", color = "white")
  ) +
  facet_wrap(~Virus, nrow = 1)


p3.d <- ggarrange(p3.m, p3.l, nrow = 2,
                common.legend = TRUE, legend = "right")

p3.d <- annotate_figure(p3.d,
                      left = text_grob("Nucleotide frequency [%]", rot = 90, vjust = 1, size = 9, color = "black", face = "bold"))

#Save supplementary figure
p3.s.freq <- ggarrange(p3.s1, p3.s2,
                       nrow = 2, heights = c(1,1),
                       labels = c("A", "B"))

ggsave("FigureS5_PositionalCombs_Most-Least.png", p3.s.freq,
       path = paste0(result.wd, "Figures/"),
       width = a4w, height = a4h * 4/5, units = "mm", dpi = "retina",
       device = "png", limitsize = FALSE)


