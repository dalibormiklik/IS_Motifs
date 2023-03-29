#Plot frequency of component repeats in repetitive elements 

#Frequency of component IS in repeats

is.rmsk <- read.table(paste0(result.wd, "/HotSpot/", cm.ppm, "/uniqNum_", virus, "_", cm.ppm, "_RMSK.txt"),
                      header = FALSE, quote = "", sep = "\t", stringsAsFactors = TRUE,
                      col.names = c("NumIS", "PropIS","repName", "repClass", "repFamily"))

rmsk.target.freq <- data.frame(Repeats = c("Any", "None"),
                               numberIS = c(Rep = sum(is.rmsk$NumIS[is.rmsk$repFamily != "None"]),
                                            None = is.rmsk$NumIS[is.rmsk$repFamily == "None"]))
rmsk.target.freq$PercIS <- 100 * rmsk.target.freq$numberIS / sum(rmsk.target.freq$numberIS)

data.to.plot <- RepTargetFreq(is.rmsk, "repFamily")

#Filter data to plot
data.to.plot <- data.to.plot[data.to.plot$repFamily != "None",]
data.to.plot <- data.to.plot[1:5,]
data.to.plot$repFamily <- factor(data.to.plot$repFamily, levels = unique(data.to.plot$repFamily))

cols <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

#Create plot
p.cmp.rmsk <- ggplot(data.to.plot, aes(group = repFamily, fill = repFamily, color = repFamily)) +
  scale_fill_manual(values = alpha(cols, .3)) +
  scale_color_manual(values = cols) +
  geom_col(aes(x = repFamily, y = PercIS), size = 1, show.legend = FALSE) +
  geom_hline(yintercept = 0) +
  ylab("% HIV IS\nof M08 C07") +
  xlab("Repeat family") +
  coord_cartesian(ylim = c(0,100)) +
  theme_classic() +
  theme(axis.line.x = element_blank(),
        #axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_text(colour = "black", face = "bold", size=8),
        axis.text.x = element_text(colour = "black", face = "bold", size=8, angle = 45, hjust = 1, vjust = 1),
        axis.text.y = element_text(colour = "black", face = "bold", size=8, angle = 0),
        strip.background = element_blank(),
        strip.text = element_text(colour = "black", size=10, face="bold"),
        text = element_text(colour = "black"),
        legend.title = element_blank(),
        legend.text = element_text(colour = "black", face="bold"),
        legend.key.size = unit(0.3, "cm"),
        legend.position = c(.15,.9)
  )

