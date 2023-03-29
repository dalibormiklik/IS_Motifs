sp.range <- seq(-13,-1, by = 1)
smpl.title <- "whole" # whole | virus | virus_celltype

#virus.smpl <- unique(smpl.info$Virus[smpl.info$Virus != "Random"])
virus.smpl <- c("HIV", "HTLV", "MLV", "PFV", "MVV", "ASLV")
#virus.smpl <- c("HIV", "HTLV")

#Create combined bar plot + points showing dinucleotide combinations KLID
select.pos <- c(-6:-1, 1:3)


df.list <- lapply(virus.smpl,
                  function(vs) {
                    select.smpl.vir <- which(smpl.info$Virus == vs)
                    
                    #Create data.frame with positional combination frequencies
                    df <- as.data.frame(do.call(rbind,
                                                lapply(select.smpl.vir,
                                                       function(ss) {
                                                         #Load sample variables
                                                         virus <- smpl.info$seqName[ss]
                                                         vrs <- smpl.info$Virus[ss]
                                                         ctrl <- smpl.info$Control[ss]
                                                         sln <- smpl.info$seqLen[ss]
                                                         celltype <- smpl.info$Cell[ss]
                                                         
                                                         #Target site duplication (length)
                                                         tsd <- unname(retro.tsd[names(retro.tsd) == vrs])
                                                         half.tsd <- floor(tsd/2)
                                                         
                                                         #Load sequences
                                                         #Load IS
                                                         is.seq <- as(read.table(paste0(data.wd, "IS/", virus, "/", virus, "_is", sln,".txt"),
                                                                                 sep = "\t", stringsAsFactor = FALSE, header = FALSE),
                                                                      "DataFrame")
                                                         
                                                         is.nuc.mat <- NucMat(is.seq[,1])
                                                         colnames(is.nuc.mat) <- RelPosNames(is.nuc.mat)
                                                         
                                                         ran.seq.name <- paste0(ctrl,"_seq_CAP")
                                                         ran.seq <- as(read.table(paste0(data.wd, "IS/Ran/", ran.seq.name, ".txt"),
                                                                                  sep = "\t", stringsAsFactor = FALSE, header = FALSE),
                                                                       "DataFrame")
                                                         
                                                         ran.nuc.mat <- NucMat(ran.seq[,1])
                                                         colnames(ran.nuc.mat) <- RelPosNames(ran.nuc.mat)
                                                         
                                                         #Set equal number of columns in the PPMs
                                                         # if random PPM has more columns than IS PPM
                                                         nc.ran <- ncol(ran.nuc.mat)
                                                         nc.is <- ncol(is.nuc.mat)
                                                         if(nc.ran > nc.is) {
                                                           c.dif <- nc.ran - nc.is
                                                           if(c.dif %% 2 == 0) {
                                                             c.start <- 1 + (c.dif / 2)
                                                             c.end <- nc.ran - (c.dif / 2)
                                                             ran.nuc.mat <- ran.nuc.mat[,c.start:c.end]
                                                           } else {
                                                             c.start <- 2 + floor(c.dif / 2)
                                                             c.end <- nc.ran - floor(c.dif / 2)
                                                             ran.nuc.mat <- ran.nuc.mat[,c.start:c.end]
                                                             colnames(ran.nuc.mat) <- RelPosNames(ran.nuc.mat)
                                                           }
                                                         }
                                                         
                                                         #Calculate dinucleotide frequency and enrichment (fold, KLID)
                                                         sp.range <- as.numeric(colnames(is.nuc.mat)[colnames(is.nuc.mat) < 0])
                                                         dinuc.freq <- ListToDF(DinucAtPosFreq(sp.range, sample.seq.mat = is.nuc.mat, control.seq.mat = ran.nuc.mat),
                                                                                method = "dinuc")
                                                         
                                                         #Add virus-specific columns to data:
                                                         ##Virus name
                                                         if(smpl.title == "whole") {dinuc.freq$Virus <- smpl.info$seqName[ss]}
                                                         if(smpl.title == "virus") {
                                                           if(celltype == "invitro") {
                                                             dinuc.freq$Virus <- paste0(vrs,"iv")
                                                           } else {
                                                             dinuc.freq$Virus <- vrs
                                                           }
                                                         }
                                                         if(smpl.title == "virus_celltype") {dinuc.freq$Virus <- paste0(vrs, "_", celltype)}
                                                         
                                                         ##position relative to cleavage site
                                                         cs.rel.pos <- as.numeric(levels(dinuc.freq$Pos)[dinuc.freq$Pos]) + half.tsd
                                                         cs.rel.pos[cs.rel.pos >= 0] <- cs.rel.pos[cs.rel.pos >= 0] + 1
                                                         dinuc.freq$tsdPos <- cs.rel.pos
                                                         
                                                         dinuc.freq
                                                       })))
                    
                    #Create data.frame with data used for plotting
                    data.to.plot <- df[df$tsdPos %in% select.pos,]
                    lvls.virus <- unique(df$Virus)
                    
                    data.to.plot$Virus <- factor(data.to.plot$Virus, levels = lvls.virus)
                    data.to.plot$tsdPos <- factor(data.to.plot$tsdPos, levels = select.pos)
                    data.to.plot$N1 <- factor(data.to.plot$N1, levels = names(nuc.cols))
                    data.to.plot$N2 <- factor(data.to.plot$N2, levels = names(nuc.cols))
                    
                    klid <- as.data.frame(do.call(rbind,
                                                  lapply(lvls.virus,
                                                         function(lvls.v) {
                                                           #Derive (per postion) KLID
                                                           k <- sapply(select.pos,
                                                                       function(p) {
                                                                         s <- sum(data.to.plot$KLID[data.to.plot$tsdPos == p &
                                                                                                      data.to.plot$Virus == lvls.v])
                                                                         if(length(s) == 0) {s <- c()}
                                                                         s
                                                                       })
                                                           
                                                           df <- data.frame(tsdPos = select.pos,
                                                                            KLID = k,
                                                                            Virus = lvls.v)
                                                           df[df == 0] <- NA
                                                           df
                                                         })))
                    klid <- klid[!is.na(klid$KLID),]
                    klid$tsdPos <- factor(klid$tsdPos, levels = select.pos)
                    klid$Virus <- factor(klid$Virus, levels = lvls.virus)
                    
                    ggplot(data = data.to.plot) +
                      geom_col(data = klid, aes(x = tsdPos, y = KLID, group = Virus),
                               fill = "gray90", color = "gray40", width = 0.75, na.rm = TRUE) +
                      geom_hline(yintercept = 0) +
                      geom_vline(xintercept = abs(min(select.pos)) + 0.5, color = "black", lty = 2) +
                      scale_shape_manual(values=c(21:25)) +
                      scale_colour_manual(values = nuc.cols) +
                      geom_quasirandom(aes(x = tsdPos, y = KLID, group = Virus, color = N2, shape = N1),
                                       size = 2, stroke = 1.5, fill = alpha("gray",.5)) +
                      #geom_text(aes(label = Virus), x = length(select.pos) - 1, y = 90, size = 3) +
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
                            strip.text = element_text(colour = "black", face = "bold"),
                            rect = element_rect(fill = "transparent"),
                            panel.background = element_rect(fill = "transparent"),
                            plot.background = element_rect(fill = "white", color = "white")
                      ) +
                      facet_wrap(~Virus, ncol = 4, scales = "fixed", drop = FALSE, strip.position = "top")
                    
                    
                    })
names(df.list) <- virus.smpl

p0 <- ggplot() + theme_void() + theme(plot.background = element_rect(fill = "white", colour =  "white"))

sapply(1:length(df.list),
       function(x) {
         df <- df.list[[x]]
         length(unique(df$Virus))
       })

p3.s1 <- ggarrange(p0,
                   ggarrange(df.list[[1]], p0, widths = c(4, 0), legend = "none"),
                   ggarrange(df.list[[2]], p0, widths = c(3.2, 0.8), legend = "none"),
                   ggarrange(df.list[[3]], p0, widths = c(1.5, 3.5), legend = "none"),
                   ggarrange(df.list[[4]], p0, widths = c(3.2, 0.8), legend = "none"),
                   ggarrange(df.list[[5]], p0, widths = c(1.5, 3.5), legend = "none"),
                   ggarrange(df.list[[6]], p0, widths = c(4, 0), legend = "none"),
          labels = c("",virus.smpl), heights = c(0.5, 1, 1, 1, 1, 1, 2),
          ncol = 1)

#To print only particular viral samples
#p3.s1 <- ggarrange(ggarrange(df.list[[1]], p0, widths = c(4,0)),
#                   ggarrange(df.list[[2]], p0, widths = c(3.2,0.8), legend = "none"),
#                   labels = virus.smpl, common.legend = TRUE,
#                   ncol = 1

#Save figure
#A4 = 210 / 297 mm
a4w <- 210
a4h <- 297
ggsave("FigureS3_PositionalCombs_KLID.png", p3.s1,
       path = paste0(result.wd, "Figures/"),
       width = a4w, height = a4h * 5/4, units = "mm", dpi = "retina",
       device = "png", limitsize = FALSE, bg = "white")
