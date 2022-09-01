#Create plots:
# i) panel H: frequency of intra-Alu IS in palindromic motif
# ii) Supplementary figure of distribution around the motif
# iii) Supplementary figure of sequence motifs of intra-Alu sequences

#Calculate frequencies of distances
motif.name <- "CT..G...C..AG"
dist.range <- dist.range
dist.bin <- dist.bin
y.range <- c(0, 10)
is.denom.factor <- 100

plot.legend <- FALSE

#Save distance distribution as supplementary figure Fig4_S1
source(paste0(scrd,"Figure4_Alu_Motif_Target_S1_MotifDistance.R"))

#Plot only frequency in the motif (Distance = 0)
data.to.plot <- d.freq[d.freq$Distance == 0,]
data.to.plot$Study <- factor(data.to.plot$Study)
data.to.plot$Sample <- factor(data.to.plot$Sample, levels = unique(data.to.plot$Sample)[smpl.order])

ctrl.col <- "#999999"

y.lab.2 <- "intra-Alu IS"

p.freq.in.motif <- ggplot() +
  scale_colour_manual(values = cols) +
  scale_fill_manual(values = alpha(cols, .4)) +
  geom_col(data = data.to.plot[data.to.plot$Group == "IS",],
           aes(x = Sample, y = Perc, colour = Sample, fill = Sample),
           show.legend = plot.legend) +
  geom_col(data = data.to.plot[data.to.plot$Group == "Shuffle",],
           aes(x = Sample, y = Perc), colour = ctrl.col, fill = alpha(ctrl.col, .3),
           show.legend = FALSE) +
  xlab("IS sample") +
  ylab(paste0(y.lab.1, y.lab.2)) +
  coord_cartesian(ylim = y.range) +
  geom_hline(yintercept = 0) +
  theme_classic() +
  theme(plot.title = element_text(colour = "black", face = "bold", size=10, hjust = 0.5),
        axis.line.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_text(colour = "black", face = "bold", size=8),
        axis.text.x = element_text(colour = "black", face = "bold", size=8, hjust = 1, angle = 45),
        axis.text.y = element_text(colour = "black", face = "bold", size=8, angle = 0),
        #axis.ticks.x = element_blank(),
        strip.background = element_blank(),
        strip.text = element_blank(),
        panel.spacing = unit(0.2, "lines"),
        text = element_text(colour = "black"),
        legend.title = element_blank(),
        legend.background = element_rect(fill="transparent"),
        legend.text = element_text(colour = "black", face="bold"),
        legend.direction = "vertical",
        legend.key.size = unit(0.3, "cm"),
        legend.position = c(.5,.9),
  ) + facet_grid(~Study, scales = "free_x", space = "free")


