#Create "empty" plot for y axis plotting

klid_logo <- function(data = NULL) {
  
  seq_type <- attr(data, "seq_type")
  col_scheme <- 'nucleotide'
  cs <- data.frame(
    letter = names(nuc.cols),
    col = cb.nuc.3,
    stringsAsFactors = FALSE
  )
  cs$group <- cs$letter
  attr(cs, 'cs_label') <- col_scheme
  class(cs) = c('data.frame','ggseqlogo_cs')
  legend_title <- attr(cs, "cs_label")
  data <- merge(data, cs, by = "letter", all.x = T)
  data <- data[order(data$order), ]
  #set color scheme
  colscale_gradient = is.numeric(cs$group)
  colscale_opts = NULL
  tmp = cs[!duplicated(cs$group) & !is.na(cs$group), ]
  col_map = unlist(split(tmp$col, tmp$group))
  colscale_opts = scale_fill_manual(values = col_map, name = legend_title, 
                                    na.value = "gray")
  #set other options
  guides_opts = guides(fill = "none")
  y_lim = c(-0.5, 1.6)
  extra_opts = NULL
  y_lab = "KLID"
  x_lab = "Position"
  
  data$group_by = with(data, interaction(seq_group, letter, 
                                         position))
  
  #Create ggplot layer for logo
  logo_layer <- layer(stat = "identity", data = data,
                      mapping = aes_string(x = "x", y = "y", fill = "group", group = "group_by"),
                      geom = "polygon", 
                      position = "identity", show.legend = NA, inherit.aes = F, 
                      params = list(na.rm = T)
  )
  
  breaks_fun = function(lim) {
    1:floor(lim[2]/1.05)
  }
  
  logo_list <- list(logo_layer, scale_x_continuous(breaks = breaks_fun,
                                                   labels = identity),
                    ylab(y_lab), xlab(x_lab), colscale_opts, 
                    guides_opts, coord_cartesian(ylim = y_lim), extra_opts)
  return(logo_list)
}

#Create list of ggplot logos
p1.list <- lapply(select.smpl,
                  function(w.smpl) {
                    
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
                    dir.list <- list.dirs(path = paste0(data.wd, "Data/", virus.file.name, "/"), full.names = FALSE)
                    comp.mix.dir <- dir.list[grep("M[0-9]+", dir.list)]
                    print(paste("Mixtures defined:", paste(comp.mix.dir, collapse = ", ")))
                    #----
                    
                    if(exists("mix.name") && mix.name %in% comp.mix.dir) {
                      
                      print(paste0("Loading PPMs of mixture: ", mix.name))
                      
                      d <- comp.mix.dir[grep(mix.name, comp.mix.dir)]
                      ppm.dir <- paste0(data.wd, "Data/", virus.file.name, "/", d, "/")
                      
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
                    
                    #Use geom_logo function to create data.frame for logo construction
                    logo.tab <- geom_logo(data.to.plot, method = "custom", seq_type = "dna", plot = FALSE)

                    if(smpl.title == "whole") {gtitle <- paste0(mix.name, " of ", virus, " (", virus.file.name, ")")}
                    if(smpl.title == "virus") {
                      if(celltype == "invitro") {
                        gtitle <- paste0(virus,"iv")
                      } else {
                        gtitle <- virus
                        }
                      }
                    if(smpl.title == "virus_celltype") {gtitle <- paste0(virus, "_", celltype)}
                    
                    ggplot() +
                      klid_logo(data = logo.tab) +
                      theme_logo() +
                      scale_x_continuous(breaks = brks, labels = as.numeric(lbls)) +
                      geom_vline(xintercept = tsd.range, color = alpha("black", .5), linetype = "dashed") +
                      geom_hline(yintercept=0) +
                      ggtitle(gtitle) +
                      theme(plot.title = element_text(hjust = 0.5, size = 10, face = "bold", colour = "black"),
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

#Empty plot for axis plot
p.ax <- ggplot() +
  klid_logo(data = logo.tab) +
  theme_logo() +
  #scale_x_continuous(breaks = brks, labels = as.numeric(lbls)) +
  geom_vline(xintercept = tsd.range, color = alpha("black", .5), linetype = "dashed") +
  geom_hline(yintercept=0) +
  ggtitle(gtitle) +
  theme(plot.title = element_text(hjust = 0.5, size = 10, face = "bold", colour = "white"),
        axis.text.x = element_text(colour = "white", face = "bold", size=5, angle = 0),
        axis.text.y = element_text(colour = "black", face = "bold", size=8, angle = 0),
        axis.title.x = element_blank(),
        axis.title.y = element_text(colour = "black", face = "bold", size=8, angle = 90),
        axis.line.y = element_line(colour = "black"),
        axis.ticks.y = element_line(colour = "black"),
        strip.background = element_blank(),
        strip.text = element_blank(),
        plot.background = element_rect(fill = "white", color = "white")
  ) +
  facet_wrap(~factor(seq_group), ncol = logo.columns)

p.ax.w<- 0.2


