select.pos <- lvls.pos[between(lvls.pos, -6, 3)]

lvls.virus <- unique(df.2$Virus)

data.to.plot <- df.2[df.2$tsdPos %in% select.pos,]

data.to.plot$IS_prop <- unlist(data.to.plot$IS_prop)
data.to.plot$Ctrl_prop <- unlist(data.to.plot$Ctrl_prop)
data.to.plot$N2_prop <- unlist(data.to.plot$N2_prop)

data.to.plot$Virus <- factor(data.to.plot$Virus, levels = lvls.virus)
data.to.plot$tsdPos <- factor(data.to.plot$tsdPos, levels = select.pos)
data.to.plot$N1 <- factor(data.to.plot$N1, levels = nucs)
data.to.plot$N2 <- factor(data.to.plot$N2, levels = nucs)

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

p3.sup.freq <- 
  ggplot(data.to.plot, aes(x = tsdPos)) +
  geom_hline(yintercept = c(25, 50, 75), colour = "gray60", linetype = 3) +
  geom_bar(data = uniq.data, aes(y = IS_prop, fill = N1), colour = NA, alpha = 0.75, stat = "identity") +
  geom_bar(aes(y = N2_prop, fill = N2), colour = NA, alpha = 1, stat = "identity") +
  geom_point(data = uniq.data, aes(y = Ctrl_prop)) +
  geom_vline(xintercept = abs(min(select.pos)) + 0.5, color = "black", lty = 2) +
  geom_text(data = uniq.data, aes(label = N1), y = uniq.data$TextPos, fontface = "bold", size = 3.5) +
  scale_colour_manual(values = c("darkgreen", "blue", "red", "gold")) +
  scale_fill_manual(values = c("darkgreen", "blue", "red", "gold")) +
  geom_hline(yintercept = 0) +
  coord_cartesian(ylim = c(0, 120)) +
  scale_y_continuous(breaks = seq(0, 100, 25), expand = c(0, 0)) +
  ggtitle("Most frequent") +
  xlab("Cleavage site-relative position") +
  ylab("Nucleotide\nfrequency [%]") +
  guides(fill=guide_legend(title="RC\npartner")) +
  theme_classic() +
  theme(plot.title = element_blank(),
        axis.line.x = element_blank(),
        axis.title.x = element_text(colour = "black", face = "bold", size=8, angle = 0),
        axis.title.y = element_text(colour = "black", face = "bold", size=8, angle = 90),
        axis.text.x = element_text(colour = "black", face = "bold", size=8, angle = 0),
        axis.text.y = element_text(colour = "black", face = "bold", size=8, angle = 0),
        legend.position="right",
        legend.title = element_text(colour="black", size=10, face="bold"),
        legend.text = element_text(colour="black", size=8, face="bold"),
        legend.key.size = unit(8, "points"),
        panel.grid = element_line(colour = "gray"),
        strip.background = element_blank(),
        strip.text = element_text(colour = "black", size=10, face="bold"),
        strip.placement = "outside",
        rect = element_rect(fill = "transparent"),
        panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "white", color = "white")
  ) +
  facet_wrap(~Virus, ncol = 4)
