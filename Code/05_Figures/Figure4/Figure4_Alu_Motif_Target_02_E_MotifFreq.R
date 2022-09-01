#Get motif variants and their frequencies
#Sequences are extracted from sample-associated sequences

#Select control sample: number of itterations performed to create shuffled control
ctrl.itter <- 100

#Calculate frequencies of the motifs
#inside lapply function: reads "fa_tab.txt" files of IS and shuffled controls
#motif.pos = which positions in target sequences are considered to create motif

#motif.pos <- c(8,9,12,16,19,20)
motif.pos <- c(1,2,5,9,12,13)

#Select data from which the plot will be created
#Select top.x number of most frequent motifs
top.x <- 10
data.to.plot <- m.freq[m.freq$Study == "Zhyvoloup" & m.freq$Sample == "WT",]
data.to.plot <- head(data.to.plot[order(data.to.plot$Count, decreasing = TRUE),], top.x)
data.to.plot <- data.to.plot[order(data.to.plot$Perc, decreasing = FALSE),]

#Convert Motif and Sample columns to factor 
data.to.plot$Motif <- factor(data.to.plot$Motif)
data.to.plot$Sample <- factor(data.to.plot$Sample)

#Set colours
smpl.cols <- cols[names(cols) %in% levels(data.to.plot$Sample)]
ctrl.col <- cols[length(cols)]

#Plot
p.motif.alu <- ggplot(data = data.to.plot, aes(y = reorder_within(Motif, Perc, Sample))) +
  scale_colour_manual(values = unname(c(smpl.cols, ctrl.col)), breaks = names(c(smpl.cols, ctrl.col))) +
  scale_fill_manual(values = alpha(unname(c(smpl.cols, ctrl.col)), .3), breaks = names(c(smpl.cols, ctrl.col))) +
  geom_col(aes(x = Perc, colour = Sample, fill = Sample), show.legend = FALSE) +
  geom_col(aes(x = Perc_Shuffle), colour = ctrl.col, fill = alpha(ctrl.col, .3)) +
  geom_vline(xintercept = 0, size = 1) +
  xlab("% intra Alu sequences") +
  ylab(element_blank()) +
  coord_cartesian(xlim = c(0,6)) +
  scale_y_reordered() +
  theme_classic() +
  theme(axis.line.y = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_text(colour = "black", face = "bold", size=8),
        axis.text.y = element_text(colour = "black", family = "mono", face = "bold", size=8, angle = 0),
        axis.text.x = element_text(colour = "black", face = "bold", size=8, angle = 0),
        strip.background = element_blank(),
        strip.text = element_blank(),
        text = element_text(colour = "black"),
        legend.title = element_blank(),
        legend.background = element_rect(fill="transparent"),
        legend.text = element_text(colour = "black", face="bold"),
        legend.key.size = unit(0.3, "cm"),
        legend.position = c(.1,.9)
  ) + facet_wrap(~Sample, ncol = 1, scale = "free")
