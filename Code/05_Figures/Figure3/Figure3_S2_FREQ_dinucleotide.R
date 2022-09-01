#Select and print frequencies
#Set tsdPos positions, virus data set and limits
sel.pos <- c("-3", "-2", "-1", "1", "2", "3")
#sel.vir <- rep("HTLV", length(sel.pos))
ymax <- 65
ylines <- seq(10, ymax, by = 10)

lvls.virus <- unique(df$Virus)

#Plot frequencies of dinucleotides on complementary positions of IS

#Divide plots into 2 parts to better fit the page
s2.part.A <- lvls.virus[1:4]
s2.part.B <- lvls.virus[-which(lvls.virus %in% s2.part.A)]
s2.part.list <- list(A = s2.part.A,
                     B = s2.part.B)

lapply(names(s2.part.list),
       function(s2.part) {
         sel.vir <- unname(unlist(s2.part.list[names(s2.part.list) == s2.part]))
         
         #Create data table with data to plot
         data.to.plot <- as.data.frame(do.call(rbind,
                                               lapply(1:length(sel.vir),
                                                      function(w.vrs) {
                                                        df[df$Virus == sel.vir[w.vrs] & df$tsdPos %in% sel.pos,]
                                                      })))
         
         #Create data table with data to plot
         #data.to.plot <- data.to.plot[data.to.plot$tsdPos %in% sel.pos,]
         
         data.to.plot$Virus <- factor(data.to.plot$Virus, levels = sel.vir)
         data.to.plot$tsdPos <- factor(data.to.plot$tsdPos, levels = sel.pos)
         data.to.plot$TextAt <- sapply(1:nrow(data.to.plot),
                                       function(x) {
                                         max(c(data.to.plot$IS[x], data.to.plot$Ran[x])) + 3
                                       })
         data.to.plot$PosLabel <- paste0("[", data.to.plot$tsdPos,"]")
         data.to.plot$VirusPos <- paste0(data.to.plot$Virus, " ", data.to.plot$PosLabel)
         data.to.plot <- data.to.plot[order(data.to.plot$Dinuc),]
         
         brks <- seq(0, 10 * floor(ymax/10), by = 10)
         
         #Create charts
         p3.s2 <- 
           ggplot(data.to.plot, aes(x = reorder(Dinuc, Ran))) +
           geom_hline(yintercept = ylines, colour = "gray") +
           geom_col(aes(y = IS, fill = Dinuc), alpha = .5, colour = NA) +
           geom_col(aes(y = Ran), fill = alpha("gray90", .6), color = "gray50") +
           geom_col(aes(y = IS), fill = NA, color = "black") +
           geom_label(aes(label = tsdPos, fontface = 2), x = 15, y = max(data.to.plot$TextAt), size = 3, fill = "gray90", colour = "black") +
           labs(y = "% of sequences", x = "Positional nucleotide combinations") +
           coord_flip() +
           scale_y_continuous(limits = c(0,ymax), breaks = brks, position = "left") +
           guides(size = "none", colour = "none", fill = "none") +
           theme_classic() +
           theme(plot.title = element_text(hjust = 0.5, colour = "black", face = "bold", angle = 90),
                 axis.line.x = element_line(colour = "gray"),
                 axis.title.x = element_text(size = 10, colour = "black", face = "bold"),
                 axis.title.y = element_text(size = 10, colour = "black", face = "bold"),
                 axis.text.x = element_text(size = 8, vjust = 0.5, colour = "black", face = "bold", angle = 0),
                 axis.line.y = element_blank(),
                 axis.text.y = element_text(size=8, colour = "black", family = "mono", face = "bold", angle = 0),
                 panel.grid = element_line(colour = "gray"),
                 strip.background = element_blank(),
                 strip.text.x = element_text(colour = "black", face = "bold", size=12, angle = 0),
                 strip.text.y = element_blank(),
                 rect = element_rect(fill = "transparent"),
                 panel.background = element_rect(fill = "transparent"),
                 plot.background = element_rect(fill = "white", color = "white")
           ) +
           #facet_wrap(Virus ~ tsdPos, ncol = 3, strip.position = "top")
           facet_grid(tsdPos ~ Virus)
         
         #Save figure
         #A4 = 210 / 297 mm
         a4w <- 210
         a4h <- 297
         ggsave(paste0("FigureS4_PositionalCombs_Freq_", s2.part, ".png"), p3.s2,
                path = paste0(result.wd, "Figures/"),
                width = a4w / 4 * length(sel.vir), height = a4h * 0.8, units = "mm", dpi = "retina",
                device = "png", limitsize = FALSE, bg = "white")
         
         
       })

