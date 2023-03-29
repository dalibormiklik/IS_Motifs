#Plotting frequency of integration around places where palindromic Alu motif occurs

#Calculate frequencies of distances
motif.name <- "CT..G...C..AG"
y.lab.2 <- "Alu IS"
dist.range <- dist.range
dist.bin <- dist.bin
y.range <- c(0, 6)
is.denom.factor <- 100

plot.legend <- FALSE

#Define distance windows
d.cat <- seq(dist.range[1], dist.range[2], by = dist.bin)

#Format plotting data
data.to.plot <- d.freq[d.freq$Study == "Zhyvoloup" & d.freq$Sample == "WT",]

#Convert Sample columns to factor 
data.to.plot$Sample <- factor(data.to.plot$Sample)

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
#Set breaks on x axis
x.brks <- c(round(min(dist.range)*3/4) , 0, round(max(dist.range)*3/4))

#Set colours
smpl.cols <- cols[names(cols) %in% levels(data.to.plot$Sample)]
ctrl.col <- cols[length(cols)]

#Create plot
p.dist.motif <- ggplot() +
  scale_colour_manual(values = c(smpl.cols, ctrl.col)) +
  scale_fill_manual(values = alpha(c(smpl.cols, ctrl.col), .3)) +
  geom_col(data = data.to.plot[data.to.plot$Group == "IS",],
           aes(x = Distance, y = Perc, colour = Sample, fill = Sample),
           show.legend = plot.legend) +
  geom_col(data = data.to.plot[data.to.plot$Group == "Shuffle",],
           aes(x = Distance, y = Perc),
           colour = ctrl.col, fill = alpha(ctrl.col, .3),
           show.legend = FALSE) +
  xlab(paste0("Distance to ", motif.name, " motif [bp]")) +
  ylab(paste0(y.lab.1, y.lab.2)) +
  coord_cartesian(xlim = dist.range, ylim = y.range) +
  geom_hline(yintercept = 0) +
  scale_x_continuous(breaks = x.brks) +
  theme_classic() +
  theme(plot.title = element_text(colour = "black", face = "bold", size=10, hjust = 0.5),
        axis.line.x = element_blank(),
        axis.title.x = element_text(colour = "black", face = "bold", size=8),
        axis.title.y = element_text(colour = "black", face = "bold", size=8),
        axis.text.x = element_text(colour = "black", face = "bold", size=8, hjust = 0.5),
        axis.text.y = element_text(colour = "black", family = "mono", face = "bold", size=8, angle = 0),
        strip.background = element_blank(),
        strip.text = element_blank(),
        text = element_text(colour = "black"),
        legend.title = element_blank(),
        legend.background = element_rect(fill="transparent"),
        legend.text = element_text(colour = "black", face="bold"),
        legend.direction = "vertical",
        legend.key.size = unit(0.3, "cm"),
        legend.position = c(.85,.9)
  ) + facet_wrap(~Sample, ncol = 1)


