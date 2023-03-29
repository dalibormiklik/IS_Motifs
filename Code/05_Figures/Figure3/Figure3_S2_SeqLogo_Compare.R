#Create KLID logos for components

#Set virus/sample folder name
# = which row of smpl.info table

#Select samples that will be plotted
select.smpl.name <- c("HIV_Zhyvoloup", "HIV_Vansant_wt", "HIV_Melamed_cellcult", "HIV_Melamed_patient",
                      "HTLV_Kirk", "HTLV_Melamed_cellcult", "HTLV_Melamed_patient"
                      )

select.smpl <- sapply(select.smpl.name,
                      function(x) {
                        which(smpl.info$seqName == x)},
                      USE.NAMES = FALSE)

smpl.title <- "whole"

pos.rng <- c(-13,13)
ylim <- c(-0.5, 1.6)

#Set theme for plots w/out y axis
theme_y_blank <- theme(axis.text.y = element_blank(),
                       axis.title.y = element_blank(),
                       axis.line.y = element_blank(),
                       axis.ticks.y = element_blank())

#Which mixtures to plot (number of components)
cm.all <- 2:8

dir.create(paste0(result.wd,"/Figures/S2_SeqLogo_Compare"))
for(cm.num in cm.all) {
  
  cm <- paste0("M0", cm.num)
  #Constrict data.frame with info about samples
  select.mix <- rep(cm, length(select.smpl))
  
  select.smpl.df <- data.frame(PPM = rep("mix", length(select.smpl)),
                               M = select.mix,
                               wSmpl = select.smpl)
  
  #Create list of ggplots
  p2.logo.list <- lapply(select.smpl,
                         function(w.smpl) {
                           
                           #Select row containing w.smpl number in wSmpl column
                           info.df <- select.smpl.df
                           w.smpl.row <- info.df[info.df$wSmpl == w.smpl,]
                           
                           #Select Mixture, PPMs, etc
                           select.ppm <- unique(w.smpl.row[,grep("PPM([0-9]|$)", colnames(w.smpl.row))])
                           select.mix <- w.smpl.row$M
                           mix.name <- select.mix
                           v.name <- w.smpl.row$Virus
                           include.ppm0 <- FALSE
                           logo.columns <- 1
                           smpl.title <- "seqname"
                           
                           #Create list of ggplot logos
                           #Load info of sample
                           virus.file.name <- smpl.info$Name[w.smpl]
                           virus <- smpl.info$Virus[w.smpl]
                           seqname <- smpl.info$seqName[w.smpl]
                           celltype <- smpl.info$Cell[w.smpl]
                           
                           tsd <- unname(retro.tsd[names(retro.tsd) == virus])
                           
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
                           source(paste0(scr.wd,'Load_mixPPM.R'))
                           
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
                                                if(all(pos.rng != "all")) {
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
                           if(smpl.title == "seqname") {gtitle <- seqname}
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
  
  #Build a pannels
  p2.logo.1 <- ggarrange(p2.logo.list[[1]],
                         p2.logo.list[[2]] + theme_y_blank,
                         p2.logo.list[[3]] + theme_y_blank,
                         p2.logo.list[[4]] + theme_y_blank,
                         ncol = 4, nrow = 1,
                         widths = c(1, rep(0.83, 3)))
  
  p2.logo.2 <- ggarrange(p2.logo.list[[5]],
                         p2.logo.list[[6]] + theme_y_blank,
                         p2.logo.list[[7]] + theme_y_blank,
                         p0 + theme_y_blank,
                         ncol = 4, nrow = 1,
                         widths = c(1, rep(0.83, 3)))
  
  #A4 = 210 / 297 mm
  a4w <- 210
  a4h <- 297
  #Calculate the height of single panel
  p.high <- a4h * 3/4 * 1/8 * cm.num
  
  #Save pannels in single image or to multiple images if size does not fit A4
  if((p.high * 2) <= a4h) {
    p2.logo.pannel <- ggarrange(p2.logo.1,
                                p2.logo.2,
                                nrow = 2)
    #Save figure with combined plots
    ggsave(paste0("FigureS2_Compare_", cm, ".png"), p2.logo.pannel,
           path = paste0(result.wd,"/Figures/S2_SeqLogo_Compare"),
           width = a4w, height = p.high * 2, units = "mm", dpi = "retina",
           device = "png", limitsize = FALSE)
    
  } else {
    #Save figure with combined plots
    ggsave(paste0("FigureS2_Compare_", cm, "_A.png"), p2.logo.1,
           path = paste0(result.wd,"/Figures/S2_SeqLogo_Compare"),
           width = a4w, height = p.high, units = "mm", dpi = "retina",
           device = "png", limitsize = FALSE)
    
    ggsave(paste0("FigureS2_Compare_", cm, "_B.png"), p2.logo.2,
           path = paste0(result.wd,"/Figures/S2_SeqLogo_Compare"),
           width = a4w, height = p.high, units = "mm", dpi = "retina",
           device = "png", limitsize = FALSE)
    
    
  }
  
  
  
}

