#Calculate distances (RC and PDeficient) of PPMs
#PPMs are stored in ppm.tab object
#Distance is calculated by PPMdist function saved in Functions.R

#Calculate RC distances for all mixtures of the sample
rc.dist <- do.call(rbind,
                   lapply(select.smpl,
                          function(w.smpl) {
                            
                            #Load mixtures
                            virus.file.name <- smpl.info$Name[w.smpl]
                            virus <- smpl.info$Virus[w.smpl]
                            celltype <- smpl.info$Cell[w.smpl]
                            if(celltype == "invitro") {
                              virus <- paste0(virus,"iv")
                            }
                            
                            dir.list <- list.dirs(path = paste0(result.wd, "Mixtures/", virus.file.name, "/"), full.names = FALSE)
                            comp.mix.dir <- dir.list[grep("M[0-9]+", dir.list)]
                            
                            #...for each of the mixtures
                            do.call(rbind,
                                    lapply(comp.mix.dir,
                                           function(mix.name) {
                                             
                                             d <- comp.mix.dir[grep(mix.name, comp.mix.dir)]
                                             ppm.dir <- paste0(result.wd, "Mixtures/", virus.file.name, "/", d, "/")
                                             
                                             c.names <- c("PPM","Pos", "A", "C", "G", "T", "KLID")
                                             
                                             ppm.tab <- read.table(paste0(ppm.dir, "PPM.txt"),
                                                                   sep = "\t", stringsAsFactor = FALSE, header = FALSE,
                                                                   col.names = c.names)
                                             
                                             ppm.tab$PPM <- as.factor(ppm.tab$PPM)
                                             
                                             #Calulate distances
                                             rc.d <- PPMdist(ppm.tab, ppm.direct = "RC", select.pos = "all")
                                             
                                             #add columns for indetification of sample and mixture
                                             rc.d$M <- mix.name
                                             rc.d$Virus <- virus
                                             rc.d$wSmpl <- w.smpl
                                             
                                             #output
                                             rc.d

                                           }))
                            
                            
                            
                            
                          }))

#Select only PDef distances
pdef.tab <- rc.dist[rc.dist$DistType == "PDef",]
#Select single component with lowest PDef for each mixture
#!!!MVV M06 PPM01 is excluded, as it is nearly random PPM (which is highly palindromic)
pdef.min <- do.call(rbind,
                    lapply(unique(pdef.tab$Virus),
                           function(v) {
                             v.tab <- pdef.tab[pdef.tab$Virus == v,]
                             v.tab.0 <- v.tab[v.tab$PPM1 == "PPM00",][1,]
                             v.tab.0$M <- "PPM0"
                             v.tab.mix <- do.call(rbind,
                                                  lapply(unique(v.tab$M),
                                                         function(m) {
                                                           m.tab <- v.tab[v.tab$M == m & v.tab$PPM1 != "PPM00",]
                                                           m.tab <- m.tab[order(m.tab$PPMdist, decreasing = FALSE),]
                                                           if(v != "MVV" && m != "M06") {
                                                             m.tab[1,]
                                                           } else {
                                                             m.tab[m.tab$PPM1 != "PPM01" | m.tab$PPM2 != "PPM01",][1,]
                                                           }
                                                         }
                                                  ))
                             rbind(v.tab.0, v.tab.mix)
                             
                           }))

#Separate PPM0s from PDef.min object
pdef.ppm0 <- pdef.min[pdef.min$M == "PPM0",]
pdef.min <- pdef.min[pdef.min$M != "PPM0",]

#Set factors for separated PDef.min and PDef.PPM0 data.frames
pdef.ppm0$M <- factor(pdef.ppm0$M, levels = unique(pdef.ppm0$M))
pdef.ppm0$Virus <- factor(pdef.ppm0$Virus, levels = unique(pdef.ppm0$Virus))

pdef.min$M <- factor(pdef.min$M, levels = unique(pdef.min$M))
pdef.min$Virus <- factor(pdef.min$Virus, levels = unique(pdef.min$Virus))

#Select the closest RC-pair from the least palindromic sequences
#Select least palindromic sequences
pdef.max <- do.call(rbind,
                    lapply(unique(pdef.tab$Virus),
                           function(v) {
                             v.tab <- pdef.tab[pdef.tab$Virus == v,]
                             v.tab.0 <- v.tab[v.tab$PPM1 == "PPM00",][1,]
                             v.tab.0$M <- "PPM0"
                             v.tab.mix <- do.call(rbind,
                                                  lapply(unique(v.tab$M),
                                                         function(m) {
                                                           m.tab <- v.tab[v.tab$M == m & v.tab$PPM1 != "PPM00",]
                                                           m.tab <- m.tab[order(m.tab$PPMdist, decreasing = TRUE),]
                                                           if(v != "MVV" && m != "M06") {
                                                             m.tab[1,]
                                                           } else {
                                                             m.tab[m.tab$PPM1 != "PPM01" | m.tab$PPM2 != "PPM01",][1,]
                                                           }
                                                         }
                                                  ))
                             rbind(v.tab.0, v.tab.mix)
                             
                           }))
pdef.max <- pdef.max[pdef.max$M != "PPM0",]
pdef.max$M <- factor(pdef.max$M, levels = unique(pdef.max$M))
pdef.max$Virus <- factor(pdef.max$Virus, levels = unique(pdef.max$Virus))
pdef.max$DistType <- "Least palindromic"


#Select only RC distances (w/out PDef)
rc.dist.only <- rc.dist[rc.dist$DistType == "RC" &
                          rc.dist$PPM1 != "PPM00" & rc.dist$PPM2 != "PPM00",]

#Select minimum RC distance per sample per mixture
rc.dist.min <- do.call(rbind,
                       lapply(unique(rc.dist.only$Virus),
                              function(v) {
                                v.tab <- rc.dist.only[rc.dist.only$Virus == v,]
                                do.call(rbind,
                                        lapply(unique(v.tab$M),
                                               function(m) {
                                                 m.tab <- v.tab[v.tab$M == m,]
                                                 unique(m.tab[m.tab$PPMdist == min(m.tab$PPMdist), which(colnames(m.tab) %in% c("PPMdist", "M", "Virus"))])
                                               }))
                                
                                
                              }))
rc.dist.min$M <- factor(rc.dist.min$M, levels = unique(rc.dist.min$M))
rc.dist.min$Virus <- factor(rc.dist.min$Virus, levels = unique(rc.dist.min$Virus))

#Add discrimination column to each data.frame before they are joined
pdef.ppm0$DistType <- "PPM0"
pdef.min$DistType <- "PDef-low PPMs (palindromic)"
pdef.max$DistType <- "PDef-high PPMs (asymmetric)"
rc.dist.min$DistType <- "Closest RC-pairs"

#Join PDef and RC-pair data.frames
min.dist.tab <- rbind(pdef.ppm0, pdef.min, pdef.max)
min.dist.tab$DistType <- factor(min.dist.tab$DistType, levels = unique(min.dist.tab$DistType))

#Set color-blind-friendly pallette
cb.cols <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#999999")

#Create empty.df with limits for axes
empty.df <- data.frame(DistType = factor(unique(min.dist.tab$DistType), levels = unique(min.dist.tab$DistType)),
                       M = c("PPM0", rep("M02",2)),
                       PPMdist = c(0.1, 0.5, 2))

#Create plot
min.dist.p <- ggplot() +
  geom_hline(yintercept = 0) +
  geom_hline(yintercept = 0.05, colour = "gray", linetype = "dashed") +
  geom_quasirandom(data = min.dist.tab,
                   aes(x = M, y = PPMdist, colour = Virus),
                   varwidth = TRUE, shape = 21, size = 2, stroke = 1.5, fill = NA, groupOnX=TRUE) +
  geom_blank(data = empty.df, aes(x = M, y = PPMdist)) +
  scale_colour_manual(values = cb.cols) +
  #coord_cartesian(ylim = c(0,2)) +
  labs(colour = "IS set:") +
  ylab("Palindromic\ndefect (PDef)") +
  xlab("Mixture") +
  theme_classic() +
  theme(axis.line.x = element_blank(),
        axis.title.x = element_text(colour = "black", face = "bold", size = 8),
        axis.title.y = element_text(colour = "black", face = "bold", size = 8, margin = margin(r = 10)),
        axis.text.x = element_text(colour = "black", face = "bold", size = 8, hjust = 0.5),
        axis.text.y = element_text(colour = "black", face = "bold", size = 8, angle = 0),
        strip.background = element_blank(),
        strip.text = element_text(size = 9, face = "bold"),
        text = element_text(colour = "black"),
        legend.title = element_text(colour = "black", face="bold", size = 10),
        legend.background = element_rect(fill="transparent"),
        legend.text = element_text(colour = "black", face="bold"),
        legend.direction = "vertical",
        legend.key.size = unit(0.3, "cm"),
        ) +
        facet_row(~DistType, scales = "free", space = "free")

