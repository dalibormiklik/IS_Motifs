#Select samples that will be plotted
select.smpl.name <- c("HIV_Zhyvoloup", "HIV_Vansant_wt", "HIV_Melamed_cellcult", "HIV_Melamed_patient",
                    "HTLV_Kirk", "HTLV_Melamed_cellcult", "HTLV_Melamed_patient",
                    "MLV_DeRavin", "MVV_HEK",
                    "PFV_WT_HT1080", "PFV_invitro",
                    "RCASC_CEF", "Alpha_Moiani", "Random")

select.smpl <- sapply(select.smpl.name,
                    function(x) {
                      which(smpl.info$seqName == x)},
                    USE.NAMES = FALSE)


#select.smpl <- c(7:10, 17, 18, 2, 19, 20, 11, 15, 13, 12, 1, 5, 16)

#Calculate RC distances for all mixtures of the sample
name.by <- "seqname" # | "virus"

rc.dist <- do.call(rbind,
                   lapply(select.smpl,
                          function(w.smpl) {
                            
                            #Load mixtures
                            virus.file.name <- smpl.info$Name[w.smpl]
                            if(name.by == "seqname") {virus <- smpl.info$seqName[w.smpl]}
                            if(name.by == "virus") {
                              virus <- smpl.info$Virus[w.smpl]
                              if(celltype == "invitro") {
                                virus <- paste0(virus,"iv")
                                }
                              }
                            celltype <- smpl.info$Cell[w.smpl]
                            
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

cm <- paste0("M0", 2:8) 
#Select only PDef distances
pdef.tab <- rc.dist[rc.dist$DistType == "PDef" &
                      rc.dist$M %in% cm,]

#Select samples that will be plotted
pdef.smpl <- select.smpl


#Transform character columns to factors
data.to.plot <- pdef.tab[pdef.tab$wSmpl %in% pdef.smpl,]
vir.lvls <- unique(pdef.tab$Virus)

#data.to.plot$M[which(data.to.plot$PPM1 == "PPM00")] <- "PPM0"

pdef.ppm0 <- data.to.plot[data.to.plot$PPM1 == "PPM00",]
pdef.ppm0$Virus <- factor(pdef.ppm0$Virus, levels = vir.lvls)

data.to.plot <- data.to.plot[data.to.plot$PPM1 != "PPM00",]
data.to.plot <- unique(data.to.plot)
data.to.plot$M <- factor(data.to.plot$M, levels = unique(data.to.plot$M))
data.to.plot$Virus <- factor(data.to.plot$Virus, levels = vir.lvls)
data.to.plot$PPM1 <- factor(data.to.plot$PPM1, levels = unique(data.to.plot$PPM1))


#Create plot
pdef.all.p <- ggplot() +
  scale_fill_brewer(palette = "Dark2") +
  geom_hline(yintercept = 0) +
  geom_quasirandom(data = pdef.ppm0,
                   aes(x = M, y = PPMdist), fill = "black",
                   varwidth = TRUE, shape = 21, size = 1.9, stroke = 0.1, groupOnX=TRUE) +
  geom_quasirandom(data = data.to.plot,
                   aes(x = M, y = PPMdist, fill = M), alpha = 0.7,
                   varwidth = TRUE, shape = 21, size = 2, stroke = 0.5, groupOnX=TRUE) +
  #geom_blank(data = empty.df, aes(x = M, y = PPMdist)) +
  coord_cartesian(ylim = c(0, 2)) +
  ylab("Palindromic defect (PDef)") +
  xlab("Mixture") +
  theme_classic() +
  theme(axis.line.x = element_blank(),
        axis.title.x = element_text(colour = "black", face = "bold", size = 8),
        axis.title.y = element_text(colour = "black", face = "bold", size = 8, margin = margin(r = 10)),
        axis.text.x = element_text(colour = "black", face = "bold", size = 8, angle = 45, hjust = 1, vjust = 1.25),
        axis.text.y = element_text(colour = "black", face = "bold", size = 8, angle = 0),
        axis.ticks.x =  element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(size = 8, face = "bold"),
        text = element_text(colour = "black"),
        legend.position = "none",
        legend.title = element_text(colour = "black", face="bold", size = 10),
        legend.background = element_rect(fill="transparent"),
        legend.text = element_text(colour = "black", face="bold"),
        legend.direction = "vertical",
        legend.key.size = unit(0.3, "cm"),
  ) +
  facet_wrap(~factor(Virus, vir.lvls), scales = "fixed", ncol = 4)

#pdef.all.p

#Save figure
#A4 = 210 / 297 mm
a4w <- 210
a4h <- 297
ggsave("FigureS1_PDef.png", pdef.all.p,
       path = paste0(result.wd,"/Figures"),
       width = a4w, height = a4h * 3/5, units = "mm", dpi = "retina",
       device = "png", limitsize = FALSE)

