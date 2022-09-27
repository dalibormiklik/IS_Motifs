#Create plots:
# i) Supplementary figure of distribution around the motif

#Select data from which the plot will be created
#Select top.x number of most frequent motifs
top.x <- 10
data.to.plot <-  do.call(rbind,
                         lapply(studies,
                                function(s.name) {
                                  smpl.names.study <- unique(as.character(m.freq$Sample[m.freq$Study == s.name]))
                                  #Select "top.x" motives of the study sample
                                  do.call(rbind,
                                          lapply(smpl.names.study,
                                                 function(smpl) {
                                                   smpl.df <- comp.motif.df[comp.motif.df$Sample == smpl,]
                                                   head(smpl.df[order(smpl.df$Count, decreasing = TRUE),], top.x)
                                                 }))
                                  
                                  
                                }))

data.to.plot <- data.to.plot[order(data.to.plot$Perc, decreasing = FALSE),]

#Convert Motif and Sample columns to factor 
data.to.plot$Motif <- factor(data.to.plot$Motif)
data.to.plot$Sample <- factor(data.to.plot$Sample, levels = smpl.names)
#Set colours
smpl.cols <- cols[-length(cols)]
ctrl.col <- cols[length(cols)]

#Create list of ggplots
p.m.freq.list <- lapply(levels(m.freq$Study),
                      function(s.name) {
                        
                        comp.motif.df <- m.freq[m.freq$Study == s.name,]
                        smpl.names.study <- unique(as.character(comp.motif.df$Sample))
                        
                        #Select "top.x" motives of the study sample
                        data.to.plot <- do.call(rbind,
                                                lapply(smpl.names.study,
                                               function(smpl) {
                                                 smpl.df <- comp.motif.df[comp.motif.df$Sample == smpl,]
                                                 head(smpl.df[order(smpl.df$Count, decreasing = TRUE),], top.x)
                                               }))

                        data.to.plot <- data.to.plot[order(data.to.plot$Perc, decreasing = FALSE),]
                        #Convert Motif and Sample columns to factor 
                        data.to.plot$Motif <- factor(data.to.plot$Motif)
                        data.to.plot$Sample <- factor(data.to.plot$Sample, levels = smpl.names)
                        
                        smpl.cols <- cols[names(cols) %in% unique(as.character(data.to.plot$Sample))]
                        ctrl.col <- cols[length(cols)]
                        
                        ggplot(data = data.to.plot, aes(y = reorder_within(Motif, Perc, Sample))) +
                          scale_colour_manual(values = unname(cols), breaks = names(cols)) +
                          scale_fill_manual(values = alpha(unname(cols), .3), breaks = names(cols)) +
                          geom_col(aes(x = Perc, colour = Sample, fill = Sample), show.legend = FALSE) +
                          geom_col(aes(x = Perc_Shuffle), colour = ctrl.col, fill = alpha(ctrl.col, .3)) +
                          geom_vline(xintercept = 0, size = 1) +
                          xlab("% intra Alu sequences") +
                          ylab(element_blank()) +
                          coord_cartesian(xlim = c(0,10)) +
                          scale_y_reordered() +
                          theme_classic() +
                          theme(axis.line.y = element_blank(),
                                axis.title.y = element_blank(),
                                axis.title.x = element_text(colour = "black", face = "bold", size=10),
                                axis.text.y = element_text(colour = "black", family = "mono", face = "bold", size=8, angle = 0),
                                axis.text.x = element_text(colour = "black", face = "bold", size=8, angle = 0),
                                strip.text = element_text(colour = "black", size=10, face="bold"),
                                strip.background = element_blank(),
                                text = element_text(colour = "black"),
                                legend.title = element_blank(),
                                legend.background = element_rect(fill="transparent"),
                                legend.text = element_text(colour = "black", face="bold"),
                                legend.key.size = unit(0.3, "cm"),
                                legend.position = c(.1,.9)
                          ) +  facet_wrap(~Sample, ncol = 3, scale = "free")
                        
                      })

#Arrange plots into pannels
p.m.freq.all <- ggarrange(ggarrange(p.m.freq.list[[1]], p0,
                                    ncol = 2, widths = c(2.3,1)),
                          p.m.freq.list[[2]],
                          p.m.freq.list[[3]],
                          nrow = 3, heights = c(1,1,4.5),
                          labels = c("A", "B", "C")
                        #labels = c("Zhyvoloup", "Vansant", "Demeulemeester")
)

#Save figure
#A4 = 210 / 297 mm
a4w <- 210
a4h <- 297
ggsave("FigureS7_HIV_Alu_Motif_Freq_All.png", p.m.freq.all,
       path = paste0(wd, "Text/Figures"),
       width = a4w, height = a4h, units = "mm", dpi = "retina",
       device = "png", limitsize = FALSE)

#Create panel with selected mutants
select.m.freq.df <- data.frame(Study = rep(levels(m.freq$Study)[3]),
                               Sample = c("S119G", "S119K", "S119R", "R231K"))

#Select "top.x" motives of the study sample
data.to.plot <- do.call(rbind,
                        lapply(1:nrow(select.m.freq.df),
                               function(w.r) {
                                 
                                 s.name <- select.m.freq.df$Study[w.r]
                                 
                                 comp.motif.df <- m.freq[m.freq$Study == s.name,]
                                 smpl.names.study <- select.m.freq.df$Sample[w.r]
                                 
                                 do.call(rbind,
                                         lapply(smpl.names.study,
                                                function(smpl) {
                                                  smpl.df <- comp.motif.df[comp.motif.df$Sample == smpl,]
                                                  head(smpl.df[order(smpl.df$Count, decreasing = TRUE),], top.x)
                                                }))
                                 
                               }))

data.to.plot <- data.to.plot[order(data.to.plot$Perc, decreasing = FALSE),]
#Convert Motif and Sample columns to factor 
data.to.plot$Motif <- factor(data.to.plot$Motif)
data.to.plot$Sample <- factor(data.to.plot$Sample, levels = smpl.names)

smpl.cols <- cols[names(cols) %in% unique(as.character(data.to.plot$Sample))]
ctrl.col <- cols[length(cols)]

p.m.freq.sel <- ggplot(data = data.to.plot, aes(y = reorder_within(Motif, Perc, Sample))) +
  scale_colour_manual(values = unname(cols), breaks = names(cols)) +
  scale_fill_manual(values = alpha(unname(cols), .3), breaks = names(cols)) +
  geom_col(aes(x = Perc, colour = Sample, fill = Sample), show.legend = FALSE) +
  geom_col(aes(x = Perc_Shuffle), colour = ctrl.col, fill = alpha(ctrl.col, .3)) +
  geom_vline(xintercept = 0, size = 1) +
  xlab("% intra Alu sequences") +
  ylab(element_blank()) +
  coord_cartesian(xlim = c(0,10)) +
  scale_y_reordered() +
  theme_classic() +
  theme(axis.line.y = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_text(colour = "black", face = "bold", size=8),
        axis.text.y = element_text(colour = "black", family = "mono", face = "bold", size=8, angle = 0),
        axis.text.x = element_text(colour = "black", face = "bold", size=6, angle = 45, vjust = 1, hjust = 1),
        strip.text = element_text(colour = "black", size=8, face="bold"),
        strip.background = element_blank(),
        text = element_text(colour = "black"),
        legend.title = element_blank(),
        legend.background = element_rect(fill="transparent"),
        legend.text = element_text(colour = "black", face="bold"),
        legend.key.size = unit(0.3, "cm"),
        legend.position = c(.1,.9)
  ) +  facet_wrap(~Sample, nrow = 1, scale = "free")

