#Select and print frequencies
#Set tsdPos positions, virus data set and limits
sel.pos <- c("-3", "-2", "-2", "-2")
sel.vir <- c("HIV", "HTLV", "MLV", "MVV")

lvls.virus <- sel.vir

ymax <- 63
ylines <- seq(10, ymax ,by = 10)

#Plot frequencies of dinucleotides on complementary positions of IS

#Create data table with data to plot
data.to.plot <- as.data.frame(do.call(rbind,
                                      lapply(1:length(sel.vir),
                                             function(w.vrs) {
                                               df[df$Virus == sel.vir[w.vrs] & df$tsdPos == sel.pos[w.vrs],]
                                             })))

data.to.plot$Virus <- factor(data.to.plot$Virus, levels = lvls.virus)
data.to.plot$tsdPos <- factor(data.to.plot$tsdPos, levels = select.pos)
data.to.plot$TextAt <- sapply(1:nrow(data.to.plot),
                              function(x) {
                                max(c(data.to.plot$IS[x], data.to.plot$Ran[x])) + 3
                              })
data.to.plot$PosLabel <- paste0("[", data.to.plot$tsdPos,"]")
data.to.plot$VirusPos <- paste0(data.to.plot$Virus, " ", data.to.plot$PosLabel)
data.to.plot <- data.to.plot[order(data.to.plot$Dinuc),]

brks <- seq(0, 10 * floor(ymax/10), by = 10)

#Create charts
p3.1 <- 
    ggplot(data.to.plot, aes(x = reorder(Dinuc, Ran))) +
    geom_hline(yintercept = ylines, colour = "gray") +
    geom_col(aes(y = IS, fill = Dinuc), alpha = .5, colour = NA) +
    geom_col(aes(y = Ran), fill = alpha("gray90", .6), color = "gray50") +
    geom_col(aes(y = IS), fill = NA, color = "black") +
    geom_label(aes(label = tsdPos, fontface = 2), x = 15, y = max(data.to.plot$TextAt), size = 3, fill = "gray90", colour = "black") +
    labs(y = "% of sequences", x = "Positional nucleotide combinations") +
    coord_flip() +
    scale_y_continuous(limits = c(0,ymax), breaks = brks, position = "right") +
    guides(size = "none", colour = "none", fill = "none") +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5, colour = "black", face = "bold", angle = 90),
          axis.line.x = element_line(colour = "gray"),
          axis.title.x = element_text(size = 8, colour = "black", face = "bold"),
          axis.title.y = element_text(size = 8, colour = "black", face = "bold"),
          axis.text.x = element_text(size = 8, vjust = 0.5, colour = "black", face = "bold", angle = 0),
          axis.line.y = element_blank(),
          axis.text.y = element_text(colour = "black", family = "mono", face = "bold", size=8, angle = 0),
          panel.grid = element_line(colour = "gray"),
          strip.background = element_blank(),
          strip.text = element_blank(),
          rect = element_rect(fill = "transparent"),
          panel.background = element_rect(fill = "transparent"),
          plot.background = element_rect(fill = "white", color = "white")
    ) +
    facet_wrap("Virus", ncol = 1, strip.position = "right")

#Select and print frequencies
#Set tsdPos positions, virus data set and limits
sel.pos <- c("1", "1", "1", "1", "1")
sel.vir <- sel.vir
ymax <- 23
ylines <- seq(10, ymax ,by = 10)

#Plot frequencies of dinucleotides on complementary positions of IS

#Create data table with data to plot
data.to.plot <- as.data.frame(do.call(rbind,
                                      lapply(1:length(sel.vir),
                                             function(w.vrs) {
                                               df[df$Virus == sel.vir[w.vrs] & df$tsdPos == sel.pos[w.vrs],]
                                             })))

data.to.plot$Virus <- factor(data.to.plot$Virus, levels = lvls.virus)
data.to.plot$tsdPos <- factor(data.to.plot$tsdPos, levels = select.pos)
data.to.plot$TextAt <- sapply(1:nrow(data.to.plot),
                              function(x) {
                                max(c(data.to.plot$IS[x], data.to.plot$Ran[x])) + 3
                              })
data.to.plot$PosLabel <- paste0("[", data.to.plot$tsdPos,"]")
data.to.plot$VirusPos <- paste0(data.to.plot$Virus, " ", data.to.plot$PosLabel)
data.to.plot <- data.to.plot[order(data.to.plot$Dinuc),]

brks <- seq(0, 10 * floor(ymax/10), by = 10)

#Create charts
p3.2 <- 
  ggplot(data.to.plot, aes(x = reorder(Dinuc, Ran))) +
  geom_hline(yintercept = ylines, colour = "gray") +
  geom_col(aes(y = IS, fill = Dinuc), alpha = .5, colour = NA) +
  geom_col(aes(y = Ran), fill = alpha("gray90", .6), color = "gray50") +
  geom_col(aes(y = IS), fill = NA, color = "black") +
  geom_label(aes(label = tsdPos, fontface = 2), x = 15, y = max(data.to.plot$TextAt), size = 3, fill = "gray90", colour = "black") +
  labs(y = "% of sequences", x = element_blank()) +
  coord_flip() +
  scale_y_continuous(limits = c(0,ymax), breaks = brks, position = "right") +
  guides(size = "none", colour = "none", fill = "none") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5, colour = "black", face = "bold", angle = 90),
        axis.title.x = element_text(size = 8, colour = "black", face = "bold"),
        axis.line.x = element_line(colour = "gray"),
        axis.text.x = element_text(size = 8, vjust = 0.5, colour = "black", face = "bold", angle = 0),
        axis.line.y = element_blank(),
        axis.text.y = element_text(colour = "black", family = "mono", face = "bold", size=8, angle = 0),
        panel.grid = element_line(colour = "gray"),
        strip.background = element_blank(),
        strip.text = element_text(colour = "black", size=10, face="bold"),
        rect = element_rect(fill = "transparent"),
        panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "white", color = "white")
  ) +
  #facet_grid(~ Virus)
  facet_wrap("Virus", ncol = 1, strip.position = "right")

#Arrange plots for both positional categories
p3.b <- ggarrange(p3.1, p3.2, ncol = 2, widths = c(3, 1.4))


#Create Supplementary Figure with all frequencies in single figure
source(paste0(scr.wd, "Figure3/Figure3_S2_FREQ_dinucleotide.R"))
