#Plot whole distance range
#Format plotting data

data.to.plot <- d.freq
data.to.plot$Study <- factor(data.to.plot$Study)
data.to.plot$Sample <- factor(data.to.plot$Sample, levels = unique(data.to.plot$Sample)[smpl.order])

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
y.lab.2 <- "intra-Alu IS"
x.brks <- c(round(min(dist.range)*3/4) , 0, round(max(dist.range)*3/4))

smpl.cols <- cols[-length(cols)]
ctrl.col <- cols[length(cols)]

#Create plot
p.dist.list <- lapply(levels(d.freq$Study),
                      function(s.name) {
                        data.to.plot <- d.freq[d.freq$Study == s.name,]
                        
                        smpl.cols <- cols[names(cols) %in% unique(as.character(data.to.plot$Sample))]
                        ctrl.col <- cols[length(cols)]
                        
                        ggplot() +
                          scale_colour_manual(values = smpl.cols.all) +
                          scale_fill_manual(values = alpha(smpl.cols.all, .3)) +
                          geom_vline(xintercept = 0, colour = "gray", linetype="dashed") +
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
                                axis.title.x = element_text(colour = "black", face = "bold", size=10),
                                axis.title.y = element_text(colour = "black", face = "bold", size=10),
                                axis.text.x = element_text(colour = "black", face = "bold", size=10, hjust = 0.5),
                                axis.text.y = element_text(colour = "black", face = "bold", size=10, angle = 0),
                                strip.background = element_blank(),
                                strip.text = element_text(face = "bold"),
                                #strip.text = element_blank(),
                                text = element_text(colour = "black"),
                                legend.title = element_blank(),
                                legend.background = element_rect(fill="transparent"),
                                legend.text = element_text(colour = "black", face="bold"),
                                legend.direction = "vertical",
                                legend.key.size = unit(0.3, "cm"),
                                legend.position = c(.85,.9)
                          ) +  facet_wrap(~Sample, ncol = 3)
                        
                      })

p.dist.all <- ggarrange(ggarrange(p.dist.list[[1]] + theme(axis.title.y = element_text(colour = "white", face = "bold", size=10)),
                                  p0,
                                  ncol = 2, widths = c(2.3,1)),
                        p.dist.list[[2]] + theme(axis.title.y = element_text(colour = "white", face = "bold", size=10)),
                        p.dist.list[[3]],
                        nrow = 3, heights = c(1,1,4.5),
                        labels = c("A", "B", "C")
                        #labels = c("Zhyvoloup", "Vansant", "Demeulemeester")
)

#Save figure
#A4 = 210 / 297 mm
a4w <- 210
a4h <- 297
ggsave("FigureS6_HIV_Alu_Motif_Dist_All.png", p.dist.all,
       path = paste0(wd, "Text/Figures"),
       width = a4w, height = a4h, units = "mm", dpi = "retina",
       device = "png", limitsize = FALSE)
