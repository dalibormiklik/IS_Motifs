#Select most palindromic components and plot the logos
#pdef.tab and pdef.min objects are created in Calculate_PPMdistance.R
head(pdef.tab)
head(pdef.min)

#Identify PPMs with minimal Self-RC distance of all components/mixtures
pdef.min.total <- do.call(rbind,
                          lapply(levels(pdef.min$Virus),
                                 function(v) {
                                   pdef.min.v <- pdef.min[pdef.min$M != "PPM0" & pdef.min$Virus == v,]
                                   d.min <- pdef.min.v[pdef.min.v$PPMdist == min(pdef.min.v$PPMdist),]
                                   pdef.tab[pdef.tab$Virus == v &
                                            pdef.tab$M == d.min$M &
                                            pdef.tab$PPMdist == d.min$PPMdist,]
                                 }))
select.smpl.sub <- unique(pdef.min.total$wSmpl)

p1.pdef.list <- lapply(select.smpl.sub,
                       function(w.smpl) {
                         
                         #Select row containing w.smpl number in wSmpl column
                         w.smpl.row <- pdef.min.total[pdef.min.total$wSmpl == w.smpl,]
                         
                         #Select Mixture, PPMs, etc
                         select.ppm <- unique(c(w.smpl.row$PPM1,w.smpl.row$PPM2))
                         select.mix <- w.smpl.row$M
                         mix.name <- select.mix
                         v.name <- w.smpl.row$Virus
                         include.ppm0 <- FALSE
                         logo.columns <- 1
                         smpl.title <- "virus"
                         
                         #Create list of ggplot logos
                         #Load info of sample
                         virus.file.name <- smpl.info$Name[w.smpl]
                         virus <- smpl.info$Virus[w.smpl]
                         celltype <- smpl.info$Cell[w.smpl]
                         
                         tsd <- unname(retro.tsd[names(retro.tsd) == virus])
                         
                         print("Creating PPMs:")
                         print(paste0("Virus: ", virus))
                         print(paste0("dataset: ", virus.file.name))
                         
                         #Load random PPM
                         #Control column in sample info table
                         ran.seq.file.name <- smpl.info$Control[w.smpl]
                         
                         ppm.ran <- as(read.table(paste0(s.wd, ran.seq.file.name, "_PPM.txt"),
                                                  sep = "\t", stringsAsFactor = FALSE, row.names = 1, header = TRUE, check.names = FALSE),
                                       "matrix")
                         
                         #Random frequency of nucleotides from ppm.ran
                         ran.nuc.freq <- apply(ppm.ran, 1, mean)
                         
                         #Load PPM
                         #Set possible component mixtures
                         #----
                         dir.list <- list.dirs(path = paste0(result.wd, "Mixtures/", virus.file.name, "/"), full.names = FALSE)
                         comp.mix.dir <- dir.list[grep("M[0-9]+", dir.list)]
                         print(paste("Mixtures defined:", paste(comp.mix.dir, collapse = ", ")))
                         #----
                         
                         if(exists("mix.name") && mix.name %in% comp.mix.dir) {
                           
                           print(paste0("Loading PPMs of mixture: ", mix.name))
                           
                           d <- comp.mix.dir[grep(mix.name, comp.mix.dir)]
                           ppm.dir <- paste0(result.wd, "Mixtures/", virus.file.name, "/", d, "/")
                           
                           c.names <- c("PPM","Pos", "A", "C", "G", "T", "KLID")
                           
                           ppm.tab <- read.table(paste0(ppm.dir, "PPM.txt"),
                                                 sep = "\t", stringsAsFactor = FALSE, header = FALSE,
                                                 col.names = c.names)
                           
                           ppm.tab$PPM <- as.factor(ppm.tab$PPM)
                           
                           ppm.names <- levels(ppm.tab$PPM)
                           
                           #Print info about mixture PPMs
                           print(paste("PPMs of mixture", mix.name, "loaded."))
                           print(paste(mix.name, "contains", length(levels(ppm.tab$PPM)), "PPMs"))
                           print(ppm.names)
                           
                           
                         } else {
                           
                           print("mix.name object not present or not defined in dataset")
                           stop(print("PPMs not loaded."))
                           
                         }
                         
                         if(any(select.mix == "all")) {
                           comp.mix.dir <- comp.mix.dir
                         } else {
                           if(any(select.mix %in% comp.mix.dir)) {
                             comp.mix.dir <- comp.mix.dir[comp.mix.dir %in% select.mix]
                           }
                         }
                         
                         #Run only if single value in comp.mix.dir
                         if(length(comp.mix.dir) != 1) {
                           if(length(comp.mix.dir) == 0) {
                             stop("No component mixture in comp.mix.dir object.")
                           }
                           if(length(comp.mix.dir) > 1) {
                             stop(paste0(length(comp.mix.dir), " component mixtures in comp.mix.dir object"))
                           }
                         } else {
                           mix.name <- comp.mix.dir
                         }
                         
                         #Load mixture PPMs
                         #source(paste0(scr.wd,'Load_mixPPM.R'))
                         
                         #Load mixture weights
                         #Extract weight for the PPM
                         wm <- read.table(paste0(result.wd, "Mixtures/", virus.file.name, "/",
                                                 mix.name, "/Wm.txt"),
                                          header = FALSE, stringsAsFactors = FALSE)
                         colnames(wm) <- c("cNum", "seqLen", "Wm")
                         
                         klid.tab <- ppm.tab
                         w.cols.nucs <- which(colnames(klid.tab) %in% nucs)
                         
                         #Calculate KLID values
                         #Random nucleotides frequencies are loaded in "Load_Mixtures.R" script
                         klid.tab[,w.cols.nucs] <- sapply(colnames(klid.tab)[w.cols.nucs],
                                                          function(n) {
                                                            pos.val <- klid.tab[,which(colnames(ppm.tab) == n)]
                                                            p.ran <- unname(ran.nuc.freq[names(ran.nuc.freq) == n])
                                                            pos.val * log(pos.val / p.ran)
                                                          })
                         #Create PPM-like matrix with KLID values
                         klid.ppm <- lapply(ppm.names,
                                            function(x) {
                                              klid <- t(klid.tab[klid.tab$PPM == x, which(colnames(klid.tab) %in% c("Pos", nucs))])
                                              klid <- PPMformat(klid)
                                              
                                              klid <- klid[order(rownames(klid)),]
                                              if(!all(pos.rng == "all")) {
                                                klid <- klid[,colnames(klid) %in% as.character(seq(pos.rng[1],pos.rng[2], by = 1))]
                                              }
                                              klid
                                            })
                         names(klid.ppm) <- ppm.names
                         
                         #Set positions and labels for x axis
                         #Expected that all PPMs have same number of columns
                         brks <- 1:ncol(klid.ppm[[1]])
                         lbls <- colnames(klid.ppm[[1]])
                         
                         #Set the boundaries for TSD rectangle
                         #saved in tsd.middle, half.tsd, is.middle
                         if(tsd %% 2 == 0) {
                           half.tsd <- tsd/2
                           is.middle <- max(brks)/2
                           tsd.range <- c(is.middle - half.tsd + 0.5,
                                          is.middle + half.tsd + 0.5
                           )
                         } else {
                           half.tsd <- tsd/2
                           is.middle <- ceiling(max(brks)/2)
                           tsd.range <- c(is.middle - half.tsd,
                                          is.middle + half.tsd
                           )
                         }
                         
                         #Select data for PPMs that the logo will be created for
                         if(select.ppm == "mix") {
                           logo.file.name <- paste0(virus.file.name, "_", mix.name)
                           if(include.ppm0) {
                             data.to.plot <- klid.ppm
                           } else {
                             data.to.plot <- klid.ppm[names(klid.ppm) != "PPM00"]
                           }
                           plot.width <- 330 * logo.columns
                         } else {
                           logo.file.name <- paste0(virus.file.name, "_", mix.name, "_", select.ppm)
                           data.to.plot <- klid.ppm[names(klid.ppm) == select.ppm]
                           plot.width <- 300
                         }
                         
                         #Modify names of List KLID PPMs to include mixture name
                         #names(data.to.plot) <- paste0(mix.name, ": ",names(data.to.plot))
                         
                         #Use geom_logo function to create data.frame for logo construction
                         logo.tab <- geom_logo(data.to.plot, method = "custom", seq_type = "dna", plot = FALSE)
                         
                         #Create table with annotations:
                         #component name (seq_group), weight (Wm), max y value per component
                         #Set name for PPM00 and change names for component
                         #component.name = preposition for component name (before number)
                         ppm0.name <- "PPM0"
                         component.name <- "C"
                         
                         ann <- as.data.frame(t(sapply(as.character(levels(logo.tab$seq_group)),
                                                       function(p) {
                                                         comp.num <- as.numeric(unlist(strsplit(p, "PPM"))[2])
                                                         if(comp.num == 0) {
                                                           wm.p <- 100
                                                         } else {
                                                           wm.p <- round(100 * wm$Wm[wm$cNum == comp.num])
                                                         }
                                                         ymax <- round(max(logo.tab$y[logo.tab$seq_group == p]), digits = 2)
                                                         out <- c(mix.name, p, wm.p, ymax)
                                                         out
                                                       })))
                         colnames(ann) <- c("Mix", "seq_group", "Wm", "yMax")
                         if(ppm0.name != "PPM00") {ann$seq_group <- sub("PPM00", ppm0.name, ann$seq_group)}
                         if(component.name != "PPM") {ann$seq_group <- sub("PPM", component.name, ann$seq_group)}
                         
                         #reorder annotations by component weights
                         ann <- ann[order(as.numeric(ann$Wm), decreasing = TRUE),]
                         ann$seq_group <- factor(ann$seq_group)
                         
                         #Set order of seq_group levels by weights in logo.tab
                         if(ppm0.name != "PPM00") {logo.tab$seq_group <- sub("PPM00", ppm0.name, logo.tab$seq_group)}
                         if(component.name != "PPM") {logo.tab$seq_group <- sub("PPM", component.name, logo.tab$seq_group)}
                         logo.tab$seq_group <- factor(logo.tab$seq_group,
                                                      levels = ann$seq_group)
                         
                         if(smpl.title == "whole") {gtitle <- paste0(mix.name, " of ", virus, " (", virus.file.name, ")")}
                         if(smpl.title == "virus") {
                           if(celltype == "invitro") {
                             gtitle <- paste0(virus,"iv", ": ", mix.name)
                           } else {
                             gtitle <- paste0(virus, ": ", mix.name)
                           }
                         }
                         if(smpl.title == "virus_celltype") {gtitle <- paste0(virus, "_", celltype)}
                         
                         ggplot() +
                           klid_logo(data = logo.tab) +
                           geom_label(data = ann, aes(x = 3, y = 1.6 - 0.1, label = seq_group, colour = NULL), label.size = 0, fontface = "bold", size = 3, fill = "gray", alpha = 0.3) +
                           geom_label(data = ann, aes(x = max(brks)-1, y = 1.6 - 0.1, label = Wm, colour = "white"), colour = "white", label.size = 0, fontface = "bold", size = 3, fill = "black", alpha = 0.9) +
                           theme_logo() +
                           scale_x_continuous(breaks = brks, labels = as.numeric(lbls)) +
                           geom_vline(xintercept = tsd.range, color = alpha("black", .5), linetype = "dashed") +
                           geom_hline(yintercept=0) +
                           ggtitle(gtitle) +
                           theme(plot.title = element_text(hjust = 0.5, size = 10, face = "bold", colour = "black"),
                                 #plot.subtitle = element_text(hjust = 0.5, size = 8, face = "bold", colour = "black"),
                                 #axis.text.x = element_text(colour = "black", face = "bold", size=5, angle = 0),
                                 axis.text.x = element_blank(),
                                 axis.text.y = element_text(colour = "black", face = "bold", size=8, angle = 0),
                                 axis.title.x = element_blank(),
                                 axis.title.y = element_text(colour = "black", face = "bold", size=8, angle = 90),
                                 axis.line.y = element_line(colour = "black"),
                                 axis.ticks.y = element_line(colour = "black"),
                                 strip.background = element_blank(),
                                 strip.text = element_blank(),
                                 plot.background = element_rect(fill = "white", color = "white")
                           ) +
                           facet_wrap(~factor(seq_group, levels = ann$seq_group), ncol = logo.columns)
                       })

p1.pdef.pannel <- ggarrange(ggarrange(p1.pdef.list[[1]],
                                      p1.pdef.list[[2]] + theme_y_blank,
                                      p1.pdef.list[[3]] + theme_y_blank,
                                      p1.pdef.list[[4]] + theme_y_blank,
                                      ncol = 4, nrow = 1),
                            ggarrange(p1.pdef.list[[5]],
                                      p1.pdef.list[[6]] + theme_y_blank,
                                      p1.pdef.list[[7]] + theme_y_blank,
                                      p1.pdef.list[[8]] + theme_y_blank,
                                      ncol = 4, nrow = 1),
                            nrow = 2)



#Identify PPMs with closest RC-pairs of all components/mixtures
rc.min.total <- do.call(rbind,
                          lapply(levels(rc.dist.min$Virus),
                                 function(v) {
                                   rc.dist.min.v <- rc.dist.min[rc.dist.min$M != "PPM0" & rc.dist.min$Virus == v,]
                                   d.min <- rc.dist.min.v[rc.dist.min.v$PPMdist == min(rc.dist.min.v$PPMdist),]
                                   rc <- rc.dist.only[rc.dist.only$Virus == v &
                                                        rc.dist.only$M == d.min$M &
                                                        rc.dist.only$PPMdist == d.min$PPMdist,]
                                   rc[!duplicated(paste0(rc$PPMdist, rc$M, rc$Virus)),]
                                 }))

p1.rc.min.list <- lapply(select.smpl.sub,
                       function(w.smpl) {
                         
                         #Select row containing w.smpl number in wSmpl column
                         w.smpl.row <- rc.min.total[rc.min.total$wSmpl == w.smpl,]
                         
                         #Select Mixture, PPMs, etc
                         select.ppm <- unique(c(w.smpl.row$PPM1, w.smpl.row$PPM2))
                         select.mix <- w.smpl.row$M
                         mix.name <- select.mix
                         v.name <- w.smpl.row$Virus
                         include.ppm0 <- FALSE
                         logo.columns <- 1
                         smpl.title <- "virus"
                         
                         #Create list of ggplot logos
                         #Load info of sample
                         virus.file.name <- smpl.info$Name[w.smpl]
                         virus <- smpl.info$Virus[w.smpl]
                         celltype <- smpl.info$Cell[w.smpl]
                         
                         tsd <- unname(retro.tsd[names(retro.tsd) == virus])
                         
                         print("Creating PPMs:")
                         print(paste0("Virus: ", virus))
                         print(paste0("dataset: ", virus.file.name))
                         
                         #Load random PPM
                         #Control column in sample info table
                         ran.seq.file.name <- smpl.info$Control[w.smpl]
                         
                         ppm.ran <- as(read.table(paste0(s.wd, ran.seq.file.name, "_PPM.txt"),
                                                  sep = "\t", stringsAsFactor = FALSE, row.names = 1, header = TRUE, check.names = FALSE),
                                       "matrix")
                         
                         #Random frequency of nucleotides from ppm.ran
                         ran.nuc.freq <- apply(ppm.ran, 1, mean)
                         
                         #Load PPM
                         #Set possible component mixtures
                         #----
                         dir.list <- list.dirs(path = paste0(result.wd, "Mixtures/", virus.file.name, "/"), full.names = FALSE)
                         comp.mix.dir <- dir.list[grep("M[0-9]+", dir.list)]
                         print(paste("Mixtures defined:", paste(comp.mix.dir, collapse = ", ")))
                         #----
                         
                         if(exists("mix.name") && mix.name %in% comp.mix.dir) {
                           
                           print(paste0("Loading PPMs of mixture: ", mix.name))
                           
                           d <- comp.mix.dir[grep(mix.name, comp.mix.dir)]
                           ppm.dir <- paste0(result.wd, "Mixtures/", virus.file.name, "/", d, "/")
                           
                           c.names <- c("PPM","Pos", "A", "C", "G", "T", "KLID")
                           
                           ppm.tab <- read.table(paste0(ppm.dir, "PPM.txt"),
                                                 sep = "\t", stringsAsFactor = FALSE, header = FALSE,
                                                 col.names = c.names)
                           
                           ppm.tab$PPM <- as.factor(ppm.tab$PPM)
                           
                           ppm.names <- levels(ppm.tab$PPM)
                           
                           #Print info about mixture PPMs
                           print(paste("PPMs of mixture", mix.name, "loaded."))
                           print(paste(mix.name, "contains", length(levels(ppm.tab$PPM)), "PPMs"))
                           print(ppm.names)
                           
                           
                         } else {
                           
                           print("mix.name object not present or not defined in dataset")
                           stop(print("PPMs not loaded."))
                           
                         }
                         
                         if(any(select.mix == "all")) {
                           comp.mix.dir <- comp.mix.dir
                         } else {
                           if(any(select.mix %in% comp.mix.dir)) {
                             comp.mix.dir <- comp.mix.dir[comp.mix.dir %in% select.mix]
                           }
                         }
                         
                         #Run only if single value in comp.mix.dir
                         if(length(comp.mix.dir) != 1) {
                           if(length(comp.mix.dir) == 0) {
                             stop("No component mixture in comp.mix.dir object.")
                           }
                           if(length(comp.mix.dir) > 1) {
                             stop(paste0(length(comp.mix.dir), " component mixtures in comp.mix.dir object"))
                           }
                         } else {
                           mix.name <- comp.mix.dir
                         }
                         
                         #Load mixture PPMs
                         #source(paste0(scr.wd,'Load_mixPPM.R'))
                         
                         #Load mixture weights
                         #Extract weight for the PPM
                         wm <- read.table(paste0(result.wd, "Mixtures/", virus.file.name, "/",
                                                 mix.name, "/Wm.txt"),
                                          header = FALSE, stringsAsFactors = FALSE)
                         colnames(wm) <- c("cNum", "seqLen", "Wm")
                         
                         klid.tab <- ppm.tab
                         w.cols.nucs <- which(colnames(klid.tab) %in% nucs)
                         
                         #Calculate KLID values
                         #Random nucleotides frequencies are loaded in "Load_Mixtures.R" script
                         klid.tab[,w.cols.nucs] <- sapply(colnames(klid.tab)[w.cols.nucs],
                                                          function(n) {
                                                            pos.val <- klid.tab[,which(colnames(ppm.tab) == n)]
                                                            p.ran <- unname(ran.nuc.freq[names(ran.nuc.freq) == n])
                                                            pos.val * log(pos.val / p.ran)
                                                          })
                         #Create PPM-like matrix with KLID values
                         klid.ppm <- lapply(ppm.names,
                                            function(x) {
                                              klid <- t(klid.tab[klid.tab$PPM == x, which(colnames(klid.tab) %in% c("Pos", nucs))])
                                              klid <- PPMformat(klid)
                                              
                                              klid <- klid[order(rownames(klid)),]
                                              if(!all(pos.rng == "all")) {
                                                klid <- klid[,colnames(klid) %in% as.character(seq(pos.rng[1],pos.rng[2], by = 1))]
                                              }
                                              klid
                                            })
                         names(klid.ppm) <- ppm.names
                         
                         #Set positions and labels for x axis
                         #Expected that all PPMs have same number of columns
                         brks <- 1:ncol(klid.ppm[[1]])
                         lbls <- colnames(klid.ppm[[1]])
                         
                         #Set the boundaries for TSD rectangle
                         #saved in tsd.middle, half.tsd, is.middle
                         if(tsd %% 2 == 0) {
                           half.tsd <- tsd/2
                           is.middle <- max(brks)/2
                           tsd.range <- c(is.middle - half.tsd + 0.5,
                                          is.middle + half.tsd + 0.5
                           )
                         } else {
                           half.tsd <- tsd/2
                           is.middle <- ceiling(max(brks)/2)
                           tsd.range <- c(is.middle - half.tsd,
                                          is.middle + half.tsd
                           )
                         }
                         
                         #Select data for PPMs that the logo will be created for
                         if(any(select.ppm == "mix")) {
                           logo.file.name <- paste0(virus.file.name, "_", mix.name)
                           if(include.ppm0) {
                             data.to.plot <- klid.ppm
                           } else {
                             data.to.plot <- klid.ppm[names(klid.ppm) != "PPM00"]
                           }
                           plot.width <- 330 * logo.columns
                         } else {
                           logo.file.name <- paste0(virus.file.name, "_", mix.name, "_", select.ppm)
                           data.to.plot <- klid.ppm[names(klid.ppm) %in% select.ppm]
                           plot.width <- 300
                         }
                         
                         #Modify names of List KLID PPMs to include mixture name
                         #names(data.to.plot) <- paste0(mix.name, ": ",names(data.to.plot))
                         
                         #Use geom_logo function to create data.frame for logo construction
                         logo.tab <- geom_logo(data.to.plot, method = "custom", seq_type = "dna", plot = FALSE)
                         
                         #Create table with annotations:
                         #component name (seq_group), weight (Wm), max y value per component
                         #Set name for PPM00 and change names for component
                         #component.name = preposition for component name (before number)
                         ppm0.name <- "PPM0"
                         component.name <- "C"
                         
                         ann <- as.data.frame(t(sapply(as.character(levels(logo.tab$seq_group)),
                                                       function(p) {
                                                         comp.num <- as.numeric(unlist(strsplit(p, "PPM"))[2])
                                                         if(comp.num == 0) {
                                                           wm.p <- 100
                                                         } else {
                                                           wm.p <- round(100 * wm$Wm[wm$cNum == comp.num])
                                                         }
                                                         ymax <- round(max(logo.tab$y[logo.tab$seq_group == p]), digits = 2)
                                                         out <- c(mix.name, p, wm.p, ymax)
                                                         out
                                                       })))
                         colnames(ann) <- c("Mix", "seq_group", "Wm", "yMax")
                         if(ppm0.name != "PPM00") {ann$seq_group <- sub("PPM00", ppm0.name, ann$seq_group)}
                         if(component.name != "PPM") {ann$seq_group <- sub("PPM", component.name, ann$seq_group)}
                         
                         #reorder annotations by component weights
                         ann <- ann[order(as.numeric(ann$Wm), decreasing = TRUE),]
                         ann$seq_group <- factor(ann$seq_group)
                         
                         #Set order of seq_group levels by weights in logo.tab
                         if(ppm0.name != "PPM00") {logo.tab$seq_group <- sub("PPM00", ppm0.name, logo.tab$seq_group)}
                         if(component.name != "PPM") {logo.tab$seq_group <- sub("PPM", component.name, logo.tab$seq_group)}
                         logo.tab$seq_group <- factor(logo.tab$seq_group,
                                                      levels = ann$seq_group)
                         
                         if(smpl.title == "whole") {gtitle <- paste0(mix.name, " of ", virus, " (", virus.file.name, ")")}
                         if(smpl.title == "virus") {
                           if(celltype == "invitro") {
                             gtitle <- paste0(virus,"iv", ": ", mix.name)
                           } else {
                             gtitle <- paste0(virus, ": ", mix.name)
                           }
                         }
                         if(smpl.title == "virus_celltype") {gtitle <- paste0(virus, "_", celltype)}
                         
                         ggplot() +
                           klid_logo(data = logo.tab) +
                           geom_label(data = ann, aes(x = 3, y = 1.6 - 0.1, label = seq_group, colour = NULL), label.size = 0, fontface = "bold", size = 3, fill = "gray", alpha = 0.3) +
                           geom_label(data = ann, aes(x = max(brks)-1, y = 1.6 - 0.1, label = Wm, colour = "white"), colour = "white", label.size = 0, fontface = "bold", size = 3, fill = "black", alpha = 0.9) +
                           theme_logo() +
                           scale_x_continuous(breaks = brks, labels = as.numeric(lbls)) +
                           geom_vline(xintercept = tsd.range, color = alpha("black", .5), linetype = "dashed") +
                           geom_hline(yintercept=0) +
                           ggtitle(gtitle) +
                           theme(plot.title = element_text(hjust = 0.5, size = 10, face = "bold", colour = "black"),
                                 #plot.subtitle = element_text(hjust = 0.5, size = 8, face = "bold", colour = "black"),
                                 #axis.text.x = element_text(colour = "black", face = "bold", size=5, angle = 0),
                                 axis.text.x = element_blank(),
                                 axis.text.y = element_text(colour = "black", face = "bold", size=8, angle = 0),
                                 axis.title.x = element_blank(),
                                 axis.title.y = element_text(colour = "black", face = "bold", size=8, angle = 90),
                                 axis.line.y = element_line(colour = "black"),
                                 axis.ticks.y = element_line(colour = "black"),
                                 strip.background = element_blank(),
                                 strip.text = element_blank(),
                                 plot.background = element_rect(fill = "white", color = "white")
                           ) +
                           facet_wrap(~factor(seq_group, levels = ann$seq_group), ncol = logo.columns)
                       })

#Find closest PPM to max(PDef) PPMs
pdef.max.rc.min <- do.call(rbind,
                           lapply(1:nrow(pdef.max),
                                  function(w.row) {
                                    #Select actual row
                                    r <- pdef.max[w.row,]
                                    #Select RC nearest PPM from the mixture
                                    rc <- rc.dist.only[rc.dist.only$PPM1 == r$PPM1 & rc.dist.only$M == r$M & rc.dist.only$Virus == r$Virus,]
                                    if(nrow(rc) > 1) {
                                      rc <- rc[order(rc$PPMdist, decreasing = FALSE),][1,]
                                    }
                                    #Select Pdef of nearest PPM to PPM1
                                    r2 <- pdef.tab[pdef.tab$PPM1 == rc$PPM2 & pdef.tab$M == rc$M & pdef.tab$Virus == rc$Virus,]
                                    #Create data.frame with both PDef values and RC distances between PPM1 & PPM2
                                    data.frame(PPM1 = rc$PPM1,
                                               PPM2 = rc$PPM2,
                                               PPM1_PDef = r$PPMdist,
                                               PPM2_Pdef = r2$PPMdist,
                                               RC_dist = rc$PPMdist,
                                               M = r$M,
                                               Virus = r$Virus,
                                               wSmpl = r$wSmpl)
                                  }))
pdef.max.rc.min$M <- factor(pdef.max.rc.min$M, levels = unique(pdef.max.rc.min$M))
pdef.max.rc.min$Virus <- factor(pdef.max.rc.min$Virus, levels = unique(pdef.max.rc.min$Virus))

pdef.max.rc.min.total <- do.call(rbind,
                        lapply(levels(pdef.max.rc.min$Virus),
                               function(v) {
                                 pdef.max.rc.min.v <- pdef.max.rc.min[pdef.max.rc.min$M != "PPM0" & pdef.max.rc.min$Virus == v,]
                                 d.min <- pdef.max.rc.min.v[pdef.max.rc.min.v$PPM1_PDef == max(pdef.max.rc.min.v$PPM1_PDef),]
                                 pdef.max.rc.min.v[rownames(pdef.max.rc.min.v) == rownames(d.min),]
                               }))


p1.pdef.max.rc.min.list <- lapply(select.smpl.sub,
                         function(w.smpl) {
                           
                           #Select row containing w.smpl number in wSmpl column
                           w.smpl.row <- pdef.max.rc.min.total[pdef.max.rc.min.total$wSmpl == w.smpl,]
                           
                           #Select Mixture, PPMs, etc
                           select.ppm <- unique(c(w.smpl.row$PPM1, w.smpl.row$PPM2))
                           select.mix <- w.smpl.row$M
                           mix.name <- select.mix
                           v.name <- w.smpl.row$Virus
                           include.ppm0 <- FALSE
                           logo.columns <- 1
                           smpl.title <- "virus"
                           
                           #Create list of ggplot logos
                           #Load info of sample
                           virus.file.name <- smpl.info$Name[w.smpl]
                           virus <- smpl.info$Virus[w.smpl]
                           celltype <- smpl.info$Cell[w.smpl]
                           
                           tsd <- unname(retro.tsd[names(retro.tsd) == virus])
                           
                           print("Creating PPMs:")
                           print(paste0("Virus: ", virus))
                           print(paste0("dataset: ", virus.file.name))
                           
                           #Load random PPM
                           #Control column in sample info table
                           ran.seq.file.name <- smpl.info$Control[w.smpl]
                           
                           ppm.ran <- as(read.table(paste0(s.wd, ran.seq.file.name, "_PPM.txt"),
                                                    sep = "\t", stringsAsFactor = FALSE, row.names = 1, header = TRUE, check.names = FALSE),
                                         "matrix")
                           
                           #Random frequency of nucleotides from ppm.ran
                           ran.nuc.freq <- apply(ppm.ran, 1, mean)
                           
                           #Load PPM
                           #Set possible component mixtures
                           #----
                           dir.list <- list.dirs(path = paste0(result.wd, "Mixtures/", virus.file.name, "/"), full.names = FALSE)
                           comp.mix.dir <- dir.list[grep("M[0-9]+", dir.list)]
                           print(paste("Mixtures defined:", paste(comp.mix.dir, collapse = ", ")))
                           #----
                           
                           if(exists("mix.name") && mix.name %in% comp.mix.dir) {
                             
                             print(paste0("Loading PPMs of mixture: ", mix.name))
                             
                             d <- comp.mix.dir[grep(mix.name, comp.mix.dir)]
                             ppm.dir <- paste0(result.wd, "Mixtures/", virus.file.name, "/", d, "/")
                             
                             c.names <- c("PPM","Pos", "A", "C", "G", "T", "KLID")
                             
                             ppm.tab <- read.table(paste0(ppm.dir, "PPM.txt"),
                                                   sep = "\t", stringsAsFactor = FALSE, header = FALSE,
                                                   col.names = c.names)
                             
                             ppm.tab$PPM <- as.factor(ppm.tab$PPM)
                             
                             ppm.names <- levels(ppm.tab$PPM)
                             
                             #Print info about mixture PPMs
                             print(paste("PPMs of mixture", mix.name, "loaded."))
                             print(paste(mix.name, "contains", length(levels(ppm.tab$PPM)), "PPMs"))
                             print(ppm.names)
                             
                             
                           } else {
                             
                             print("mix.name object not present or not defined in dataset")
                             stop(print("PPMs not loaded."))
                             
                           }
                           
                           if(any(select.mix == "all")) {
                             comp.mix.dir <- comp.mix.dir
                           } else {
                             if(any(select.mix %in% comp.mix.dir)) {
                               comp.mix.dir <- comp.mix.dir[comp.mix.dir %in% select.mix]
                             }
                           }
                           
                           #Run only if single value in comp.mix.dir
                           if(length(comp.mix.dir) != 1) {
                             if(length(comp.mix.dir) == 0) {
                               stop("No component mixture in comp.mix.dir object.")
                             }
                             if(length(comp.mix.dir) > 1) {
                               stop(paste0(length(comp.mix.dir), " component mixtures in comp.mix.dir object"))
                             }
                           } else {
                             mix.name <- comp.mix.dir
                           }
                           
                           #Load mixture PPMs
                           #source(paste0(scr.wd,'Load_mixPPM.R'))
                           
                           #Load mixture weights
                           #Extract weight for the PPM
                           wm <- read.table(paste0(result.wd, "Mixtures/", virus.file.name, "/",
                                                   mix.name, "/Wm.txt"),
                                            header = FALSE, stringsAsFactors = FALSE)
                           colnames(wm) <- c("cNum", "seqLen", "Wm")
                           
                           klid.tab <- ppm.tab
                           w.cols.nucs <- which(colnames(klid.tab) %in% nucs)
                           
                           #Calculate KLID values
                           #Random nucleotides frequencies are loaded in "Load_Mixtures.R" script
                           klid.tab[,w.cols.nucs] <- sapply(colnames(klid.tab)[w.cols.nucs],
                                                            function(n) {
                                                              pos.val <- klid.tab[,which(colnames(ppm.tab) == n)]
                                                              p.ran <- unname(ran.nuc.freq[names(ran.nuc.freq) == n])
                                                              pos.val * log(pos.val / p.ran)
                                                            })
                           #Create PPM-like matrix with KLID values
                           klid.ppm <- lapply(ppm.names,
                                              function(x) {
                                                klid <- t(klid.tab[klid.tab$PPM == x, which(colnames(klid.tab) %in% c("Pos", nucs))])
                                                klid <- PPMformat(klid)
                                                
                                                klid <- klid[order(rownames(klid)),]
                                                if(!all(pos.rng == "all")) {
                                                  klid <- klid[,colnames(klid) %in% as.character(seq(pos.rng[1],pos.rng[2], by = 1))]
                                                }
                                                klid
                                              })
                           names(klid.ppm) <- ppm.names
                           
                           #Set positions and labels for x axis
                           #Expected that all PPMs have same number of columns
                           brks <- 1:ncol(klid.ppm[[1]])
                           lbls <- colnames(klid.ppm[[1]])
                           
                           #Set the boundaries for TSD rectangle
                           #saved in tsd.middle, half.tsd, is.middle
                           if(tsd %% 2 == 0) {
                             half.tsd <- tsd/2
                             is.middle <- max(brks)/2
                             tsd.range <- c(is.middle - half.tsd + 0.5,
                                            is.middle + half.tsd + 0.5
                             )
                           } else {
                             half.tsd <- tsd/2
                             is.middle <- ceiling(max(brks)/2)
                             tsd.range <- c(is.middle - half.tsd,
                                            is.middle + half.tsd
                             )
                           }
                           
                           #Select data for PPMs that the logo will be created for
                           if(any(select.ppm == "mix")) {
                             logo.file.name <- paste0(virus.file.name, "_", mix.name)
                             if(include.ppm0) {
                               data.to.plot <- klid.ppm
                             } else {
                               data.to.plot <- klid.ppm[names(klid.ppm) != "PPM00"]
                             }
                             plot.width <- 330 * logo.columns
                           } else {
                             logo.file.name <- paste0(virus.file.name, "_", mix.name, "_", select.ppm)
                             data.to.plot <- klid.ppm[names(klid.ppm) %in% select.ppm]
                             plot.width <- 300
                           }
                           
                           #Modify names of List KLID PPMs to include mixture name
                           #names(data.to.plot) <- paste0(mix.name, ": ",names(data.to.plot))
                           
                           #Use geom_logo function to create data.frame for logo construction
                           logo.tab <- geom_logo(data.to.plot, method = "custom", seq_type = "dna", plot = FALSE)
                           
                           #Create table with annotations:
                           #component name (seq_group), weight (Wm), max y value per component
                           #Set name for PPM00 and change names for component
                           #component.name = preposition for component name (before number)
                           ppm0.name <- "PPM0"
                           component.name <- "C"
                           
                           ann <- as.data.frame(t(sapply(as.character(levels(logo.tab$seq_group)),
                                                         function(p) {
                                                           comp.num <- as.numeric(unlist(strsplit(p, "PPM"))[2])
                                                           if(comp.num == 0) {
                                                             wm.p <- 100
                                                           } else {
                                                             wm.p <- round(100 * wm$Wm[wm$cNum == comp.num])
                                                           }
                                                           ymax <- round(max(logo.tab$y[logo.tab$seq_group == p]), digits = 2)
                                                           out <- c(mix.name, p, wm.p, ymax)
                                                           out
                                                         })))
                           colnames(ann) <- c("Mix", "seq_group", "Wm", "yMax")
                           if(ppm0.name != "PPM00") {ann$seq_group <- sub("PPM00", ppm0.name, ann$seq_group)}
                           if(component.name != "PPM") {ann$seq_group <- sub("PPM", component.name, ann$seq_group)}
                           
                           #reorder annotations by component weights
                           ann <- ann[order(as.numeric(ann$Wm), decreasing = TRUE),]
                           ann$seq_group <- factor(ann$seq_group)
                           
                           #Set order of seq_group levels by weights in logo.tab
                           if(ppm0.name != "PPM00") {logo.tab$seq_group <- sub("PPM00", ppm0.name, logo.tab$seq_group)}
                           if(component.name != "PPM") {logo.tab$seq_group <- sub("PPM", component.name, logo.tab$seq_group)}
                           logo.tab$seq_group <- factor(logo.tab$seq_group,
                                                        levels = ann$seq_group)
                           
                           if(smpl.title == "whole") {gtitle <- paste0(mix.name, " of ", virus, " (", virus.file.name, ")")}
                           if(smpl.title == "virus") {
                             if(celltype == "invitro") {
                               gtitle <- paste0(virus,"iv", ": ", mix.name)
                             } else {
                               gtitle <- paste0(virus, ": ", mix.name)
                             }
                           }
                           if(smpl.title == "virus_celltype") {gtitle <- paste0(virus, "_", celltype)}
                           
                           ggplot() +
                             klid_logo(data = logo.tab) +
                             geom_label(data = ann, aes(x = 3, y = 1.6 - 0.1, label = seq_group, colour = NULL), label.size = 0, fontface = "bold", size = 3, fill = "gray", alpha = 0.3) +
                             geom_label(data = ann, aes(x = max(brks)-1, y = 1.6 - 0.1, label = Wm, colour = "white"), colour = "white", label.size = 0, fontface = "bold", size = 3, fill = "black", alpha = 0.9) +
                             theme_logo() +
                             scale_x_continuous(breaks = brks, labels = as.numeric(lbls)) +
                             geom_vline(xintercept = tsd.range, color = alpha("black", .5), linetype = "dashed") +
                             geom_hline(yintercept=0) +
                             ggtitle(gtitle) +
                             theme(plot.title = element_text(hjust = 0.5, size = 10, face = "bold", colour = "black"),
                                   #plot.subtitle = element_text(hjust = 0.5, size = 8, face = "bold", colour = "black"),
                                   #axis.text.x = element_text(colour = "black", face = "bold", size=5, angle = 0),
                                   axis.text.x = element_blank(),
                                   axis.text.y = element_text(colour = "black", face = "bold", size=8, angle = 0),
                                   axis.title.x = element_blank(),
                                   axis.title.y = element_text(colour = "black", face = "bold", size=8, angle = 90),
                                   axis.line.y = element_line(colour = "black"),
                                   axis.ticks.y = element_line(colour = "black"),
                                   strip.background = element_blank(),
                                   strip.text = element_blank(),
                                   plot.background = element_rect(fill = "white", color = "white")
                             ) +
                             facet_wrap(~factor(seq_group, levels = ann$seq_group), ncol = logo.columns)
                         })

#Create figure panel
p1.rc.min.pannel <- ggarrange(ggarrange(p1.rc.min.list[[1]],
                                        p1.rc.min.list[[2]] + theme_y_blank,
                                        p1.rc.min.list[[3]] + theme_y_blank,
                                        p1.rc.min.list[[4]] + theme_y_blank,
                                      ncol = 4, nrow = 1),
                            ggarrange(p1.rc.min.list[[5]],
                                      p1.rc.min.list[[6]] + theme_y_blank,
                                      p1.rc.min.list[[7]] + theme_y_blank,
                                      p1.rc.min.list[[8]] + theme_y_blank,
                                      ncol = 4, nrow = 1),
                            nrow = 2)

p1.pdef.max.rc.min.pannel <- ggarrange(ggarrange(p1.pdef.max.rc.min.list[[1]],
                                                 p1.pdef.max.rc.min.list[[2]] + theme_y_blank,
                                                 p1.pdef.max.rc.min.list[[3]] + theme_y_blank,
                                                 p1.pdef.max.rc.min.list[[4]] + theme_y_blank,
                                        ncol = 4, nrow = 1),
                              ggarrange(p1.pdef.max.rc.min.list[[5]],
                                        p1.pdef.max.rc.min.list[[6]] + theme_y_blank,
                                        p1.pdef.max.rc.min.list[[7]] + theme_y_blank,
                                        p1.pdef.max.rc.min.list[[8]] + theme_y_blank,
                                        ncol = 4, nrow = 1),
                              nrow = 2)

