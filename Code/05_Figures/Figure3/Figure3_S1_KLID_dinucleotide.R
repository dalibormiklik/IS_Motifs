#Create combined bar plot + points showing dinucleotide combinations KLID
select.pos <- c(-6:-1, 1:3)

lvls.virus <- unique(df$Virus)

#Create data.frame with data used for plotting
data.to.plot <- df[df$tsdPos %in% select.pos,]

data.to.plot$Virus <- factor(data.to.plot$Virus, levels = lvls.virus)
data.to.plot$tsdPos <- factor(data.to.plot$tsdPos, levels = select.pos)
data.to.plot$N1 <- factor(data.to.plot$N1, levels = names(nuc.cols))
data.to.plot$N2 <- factor(data.to.plot$N2, levels = names(nuc.cols))

klid <- as.data.frame(do.call(rbind,
                              lapply(lvls.virus,
                                     function(lvls.v) {
                                       #Derive (per postion) KLID
                                       k <- sapply(lvls.pos,
                                                   function(p) {
                                                     s <- sum(data.to.plot$KLID[data.to.plot$tsdPos == p &
                                                                                  data.to.plot$Virus == lvls.v])
                                                     if(length(s) == 0) {s <- c()}
                                                     s
                                                   })
                                       
                                       df <- data.frame(tsdPos = lvls.pos,
                                                        KLID = k,
                                                        Virus = lvls.v)
                                       df[df == 0] <- NA
                                       df
                                     })))
klid <- klid[!is.na(klid$KLID),]
klid$tsdPos <- factor(klid$tsdPos, levels = select.pos)
klid$Virus <- factor(klid$Virus, levels = lvls.virus)

#PLOT
p3.s1 <- 
  ggplot(data = data.to.plot) +
  geom_col(data = klid, aes(x = tsdPos, y = KLID, group = Virus),
           fill = "gray90", color = "gray40", width = 0.75, na.rm = TRUE) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = abs(min(select.pos)) + 0.5, color = "black", lty = 2) +
  scale_shape_manual(values=c(21:25)) +
  scale_colour_manual(values = nuc.cols) +
  geom_quasirandom(aes(x = tsdPos, y = KLID, group = Virus, color = N2, shape = N1),
                   size = 2, stroke = 1.5, fill = alpha("gray",.5)) +
  geom_text(aes(label = Virus), x = length(select.pos) - 1, y = 90, fontface = "bold", size = 3) +
  xlab("STR site-relative position") +
  coord_cartesian(ylim = c(0, 120)) +
  guides(shape = guide_legend(order=1),
         color = guide_legend(order=2)) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0, colour = "black", face = "bold"),
        axis.line.x = element_blank(),
        axis.title.x = element_text(colour = "black", face = "bold", size=8, angle = 0),
        axis.title.y = element_text(colour = "black", face = "bold", size=8),
        axis.text.x = element_text(colour = "black", face = "bold", size=8, angle = 0),
        axis.text.y = element_text(colour = "black", face = "bold", size=8, angle = 0),
        legend.position=c("top"), legend.box = "vertical",
        legend.title = element_text(colour = "black", size=9, face="bold"),
        legend.text = element_text(colour = "black", size=8, face="bold"),
        legend.spacing = unit(0, units = "points"),
        strip.background = element_blank(),
        strip.text = element_blank(),
        rect = element_rect(fill = "transparent"),
        panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "white", color = "white")
  ) +
  facet_wrap(~Virus, ncol = 2, scales = "fixed", strip.position = "right")

#Save figure
#A4 = 210 / 297 mm
a4w <- 210
a4h <- 297
ggsave("FigureS3_PositionalCombs_KLID.png", p3.s1,
       path = paste0(result.wd, "Figures/"),
       width = a4w, height = a4h * 3/4, units = "mm", dpi = "retina",
       device = "png", limitsize = FALSE, bg = "white")

