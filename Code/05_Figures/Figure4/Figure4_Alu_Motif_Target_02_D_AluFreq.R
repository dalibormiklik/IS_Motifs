#Plotting frequency of integration in Alu repeats

#Create table where file and sample names are connected in single object
smpl.tab.study <- smpl.tab[smpl.tab$Study == "Zhyvoloup" & smpl.tab$Sample == "WT",]
smpl.tab.study$Sample <- "HIV"

#Calculate frequencies of Alu targeting
#Set type of chart: point, histo
plot.title <- "Alu targeting"
y.lab.2 <- "IS in Alu"
y.range <- c(0, 20)
is.denom.factor <- 100

plot.legend <- TRUE

#Format plotting data
data.to.plot <- alu.f.tab[alu.f.tab$Virus == smpl.tab.study$Virus,]
data.to.plot$Virus <- factor(data.to.plot$Virus)
data.to.plot$Sample <- factor(sapply(as.character(data.to.plot$Group),
                                     function(gr) {
                                       if(gr == "Shuffle") {"Random"} else {"HIV"}
                                     }),
                              levels = c("HIV", "Random"))

#Set style of y axis labels
if(is.denom.factor == 100) {
  y.lab.1 <- "% of "
} else {
  if(is.denom.factor > 100) {
    y.lab.1 <- paste0("IS per 10^", log10(is.denom.factor), " ")
  } else {
    y.lab.1 <- paste0("IS per ", 10 * log10(is.denom.factor), " ")
  }
}

smpl.cols <- cols[names(cols) == "WT"]
names(smpl.cols) <- "HIV"
ctrl.col <- cols[length(cols)]

#Create plot
p.freq.alu <- ggplot() +
  scale_colour_manual(values = c(smpl.cols, ctrl.col)) +
  scale_fill_manual(values = alpha(c(smpl.cols, ctrl.col), .3)) +
  geom_col(data = data.to.plot[data.to.plot$Group == "IS",],
           aes(x = Virus, y = Perc, colour = Sample, fill = Sample),
           show.legend = plot.legend) +
  geom_col(data = data.to.plot[data.to.plot$Group == "Shuffle",],
           aes(x = Virus, y = Perc), colour = ctrl.col, fill = alpha(ctrl.col, .3),
           show.legend = FALSE) +
  xlab("IS sample") +
  ylab(paste0(y.lab.1, y.lab.2)) +
  coord_cartesian(ylim = y.range, xlim = c(0.5, 2)) +
  geom_hline(yintercept = 0) +
  theme_classic() +
  theme(plot.title = element_text(colour = "black", face = "bold", size=10, hjust = 0.5),
        axis.line.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_text(colour = "black", face = "bold", size=8),
        axis.text.x = element_blank(),
        axis.text.y = element_text(colour = "black", face = "bold", size=8, angle = 0),
        axis.ticks.x = element_blank(),
        strip.background = element_blank(),
        strip.text = element_blank(),
        panel.spacing = unit(0, "lines"),
        text = element_text(colour = "black"),
        legend.title = element_blank(),
        legend.background = element_rect(fill="transparent"),
        legend.text = element_text(colour = "black", face="bold"),
        legend.direction = "vertical",
        legend.key.size = unit(0.3, "cm"),
        legend.position = c(.5,.9)
  ) + facet_wrap(~Virus, nrow = 1, scales = "free_x")
