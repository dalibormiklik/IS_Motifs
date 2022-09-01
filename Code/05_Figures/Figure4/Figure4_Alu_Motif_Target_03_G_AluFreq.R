#Create plot w/ frequency of HIV IS in Alu repeats
#Show sample from different sources

#Select samples
data.to.plot <- alu.f.tab

smpl.names <- levels(data.to.plot$Sample)

y.lab.2 <- "intra-Alu IS"
y.range <- c(0, 20)
plot.legend <- FALSE

p.freq.alu <- ggplot() +
  scale_colour_manual(values = cols) +
  scale_fill_manual(values = alpha(cols, .3)) +
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

