#Create supplementary figure with all PDefs per mixture

#pdef.tab object created in Figure1_02_PPMdist.R script is used

#Transform character columns to factors
data.to.plot <- pdef.tab
data.to.plot$M[which(data.to.plot$PPM1 == "PPM00")] <- "PPM0"
data.to.plot <- unique(data.to.plot)
data.to.plot$M <- factor(data.to.plot$M, levels = unique(data.to.plot$M))
data.to.plot$Virus <- factor(data.to.plot$Virus, levels = unique(data.to.plot$Virus))
data.to.plot$PPM1 <- factor(data.to.plot$PPM1, levels = unique(data.to.plot$PPM1))

#Create plot
pdef.all.p <- ggplot() +
  geom_hline(yintercept = 0) +
  geom_hline(yintercept = 0.05, colour = "gray", linetype = "dashed") +
  geom_quasirandom(data = data.to.plot,
                   aes(x = M, y = PPMdist), colour = "gray50",
                   varwidth = TRUE, shape = 21, size = 2, stroke = 1.5, fill = NA, groupOnX=TRUE) +
  #geom_blank(data = empty.df, aes(x = M, y = PPMdist)) +
  coord_cartesian(ylim = c(0,2)) +
  #scale_colour_manual(values = cbPalette) +
  #coord_cartesian(ylim = c(0,2)) +
  labs(colour = "IS set:") +
  ylab("Palindromic defect (PDef)") +
  xlab("Mixture") +
  theme_classic() +
  theme(axis.line.x = element_blank(),
        axis.title.x = element_text(colour = "black", face = "bold", size = 10),
        axis.title.y = element_text(colour = "black", face = "bold", size = 10, margin = margin(r = 10)),
        axis.text.x = element_text(colour = "black", face = "bold", size = 10, angle = 45, hjust = 1, vjust = 1.25),
        axis.text.y = element_text(colour = "black", face = "bold", size = 10, angle = 0),
        axis.ticks.x =  element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(size = 12, face = "bold"),
        text = element_text(colour = "black"),
        legend.title = element_text(colour = "black", face="bold", size = 10),
        legend.background = element_rect(fill="transparent"),
        legend.text = element_text(colour = "black", face="bold"),
        legend.direction = "vertical",
        legend.key.size = unit(0.3, "cm"),
  ) +
  facet_wrap(~Virus, scales = "free")

#Save figure
#A4 = 210 / 297 mm
a4w <- 210
a4h <- 297
ggsave("FigureS1_PDef.png", pdef.all.p,
       path = paste0(result.wd,"/Figures"),
       width = a4w, height = a4h * 3/4, units = "mm", dpi = "retina",
       device = "png", limitsize = FALSE)
