#Create sequence logo with KLID values
#Modified function from ggseqlogo package

klid_logo <- function(data = NULL) {
  
  seq_type <- attr(data, "seq_type")
  col_scheme <- 'nucleotide'
  cs <- data.frame(
    letter = nucs,
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

#Functions for calculation of positional combinations
#Together they create 'df' object:
# nucleotide combination frequencies at STR site-relative positions

DinucAtPosFreq <- function(select.pos, sample.seq.mat = is.nuc.mat, control.seq.mat = ran.nuc.mat) {
  
  #check if colnames of is and random matrices are the same
  if(all(colnames(sample.seq.mat) == colnames(control.seq.mat))) {
    
    #Count number of all sequences
    #"nuc.mat" matrices shoul be predefined
    all.is.count <- nrow(sample.seq.mat)
    all.ran.count <- nrow(control.seq.mat)
    
    #Define nucleotides and complement counterparts
    nucs <- c("A", "C", "T", "G")
    rc.nucs <- c("T", "G", "A", "C")
    
    #Table with dinucleotides
    n.tab <- as.data.frame(do.call(rbind,lapply(nucs,
                                                function(n1) {
                                                  n1.rc <- rc.nucs[nucs == n1]
                                                  di.f <- unname(sapply(rc.nucs,
                                                                        function(rc.n2) {
                                                                          paste0(n1,rc.n2, collapse = "")
                                                                        }))
                                                  di.r <- unname(sapply(rc.nucs,
                                                                        function(rc.n2) {
                                                                          paste0(rc.n2,n1, collapse = "")
                                                                        }))
                                                  
                                                  a <- cbind(rep(n1, length(nucs)),
                                                             nucs,
                                                             di.f, di.r
                                                  )
                                                  colnames(a) <- c("5", "3", "FW", "RV")
                                                  a
                                                })
    ))
    rownames(n.tab) <- n.tab[,which(colnames(n.tab) == "FW")]
    
    #Create list of position probabilities
    pos.list <- lapply(select.pos,
                       function(pos) {
                         
                         n.freq.mat <- t(sapply(rownames(n.tab),
                                                function(nuc) {
                                                  #Select nucleotide and position
                                                  #Set left and right values
                                                  w.row <- which(rownames(n.tab) == nuc)
                                                  neg.nuc <- n.tab[w.row, which(colnames(n.tab) == "5")]
                                                  pos.nuc <- n.tab[w.row, which(colnames(n.tab) == "3")]
                                                  if(pos < 0) {
                                                    neg.pos <- pos
                                                    pos.pos <- abs(pos)
                                                  } else {
                                                    neg.pos <- -1 * pos
                                                    pos.pos <- pos
                                                  }
                                                  
                                                  #Select which positions are equal to relative positions
                                                  w.col.neg <- which(colnames(sample.seq.mat) == neg.pos)
                                                  w.col.pos <- which(colnames(sample.seq.mat) == pos.pos)
                                                  
                                                  #CALCULATE frequencies of sequences in % of all sequences
                                                  n.freq <- c((100 * nrow(rbind(sample.seq.mat[sample.seq.mat[,w.col.neg] == neg.nuc & sample.seq.mat[,w.col.pos] == pos.nuc,])) /
                                                                 all.is.count),
                                                              (100 * nrow(rbind(control.seq.mat[control.seq.mat[,w.col.neg] == neg.nuc & control.seq.mat[,w.col.pos] == pos.nuc,])) /
                                                                 all.ran.count)
                                                  )
                                                  names(n.freq) <- c("IS", "Ctrl")
                                                  n.freq
                                                  #end of sapply(select.nuc)
                                                })
                                         #rbind
                         )
                         
                         rownames(n.freq.mat) <- rownames(n.tab)
                         #Return frequency matrix
                         t(n.freq.mat)
                         
                         #end of lapply(select.pos)
                       })
    names(pos.list) <- select.pos
    
    return(pos.list)
    
    print(paste0("Frequency of ", nuc, " at position ", pos,":"))
    print(n.freq)
    
  } else {
    print("IS and Random colnames are NOT eaqual!")
  }
}
ListToDF <- function(list.of.matrices, method = "dinuc") {
  df <- as.data.frame(
    do.call(rbind,
            lapply(1:length(list.of.matrices),
                   function(x) {
                     lnam <- names(list.of.matrices)[x]
                     l <- list.of.matrices[[x]]
                     l.is <- unname(l[which(rownames(l) == "IS"),])
                     l.ran <- unname(l[which(rownames(l) == "Ctrl"),])
                     #Fold increase in frequency to random (IS / Ran)
                     f <- l.is / l.ran
                     #KLID
                     k <- l.is * log(l.is / l.ran)
                     c.nam <- colnames(l)
                     
                     if(method == "dinuc") {
                       c.nam2 <- do.call(rbind,
                                         unname(sapply(c.nam,
                                                       function(n) {
                                                         strsplit(n,"")
                                                       })))
                       c.nam.all <- cbind(lnam,
                                          l.is,
                                          l.ran,
                                          f,
                                          k,
                                          c.nam,
                                          c.nam2
                       )
                     }
                     if(method == "single") {
                       c.nam.all <- cbind(lnam,
                                          l.is,
                                          l.ran,
                                          f,
                                          k,
                                          c.nam
                       )
                     }
                     
                     c.nam.all
                     
                   })
    )
  )
  
  if(method == "dinuc") {
    colnames(df) <- c("Pos", "IS", "Ran", "Enrich", "KLID", "Dinuc", "N1", "N2")
    df$Pos <- factor(df$Pos, levels = sp.range)
    df$IS <- as.numeric(df$IS)
    df$Ran <- as.numeric(df$Ran)
    df$Enrich <- as.numeric(df$Enrich)
    df$KLID <- as.numeric(df$KLID)
    df$N1 <- factor(df$N1, levels = nucs)
    df$N2 <- factor(df$N2, levels = nucs)
    
  }
  if(method == "single") {
    colnames(df) <- c("Pos", "IS", "Ran", "Enrich", "KLID", "Nuc")
    df$Pos <- factor(df$Pos, levels = sp.range)
    df$IS <- as.numeric(df$IS)
    df$Ran <- as.numeric(df$Ran)
    df$Enrich <- as.numeric(df$Enrich)
    df$KLID <- as.numeric(df$KLID)
    df$Nuc <- factor(df$Nuc, levels = nucs)
  }
  
  return(df)
  
}
NucAtPosFreq <- function(select.pos, select.nuc, sample.seq.mat = is.nuc.mat, control.seq.mat = ran.nuc.mat, output.number = "count") {
  
  #check if colnames of is and random matrices are the same
  if(all(colnames(sample.seq.mat) == colnames(control.seq.mat))) {
    
    #Count number of all sequences
    #"nuc.mat" matrices shoul be predefined
    all.is.count <- nrow(sample.seq.mat)
    all.ran.count <- nrow(control.seq.mat)
    
    #Set if frequency or sequence numbers will be output
    if(output.number == "percent") {
      multip <- 100
      denom.is <- all.is.count
      denom.ran <- all.ran.count 
    }
    if(output.number == "count") {
      multip <- 1
      denom.is <- 1
      denom.ran <- 1 
    }
    
    #Define nucleotides and complement counterparts
    nucs <- c("A", "C", "T", "G")
    rc.nucs <- c("T", "G", "A", "C")
    
    pos.list <- lapply(select.pos,
                       function(pos) {
                         
                         n.freq.mat <- t(sapply(select.nuc,
                                                function(nuc) {
                                                  #Select nucleotide and position
                                                  #Set left and right values
                                                  if(pos == 0) {
                                                    #CALCULATE frequencies of sequences in % of all sequences
                                                    w.col.pos <- which(colnames(sample.seq.mat) == pos)
                                                    n.freq <- c((multip * nrow(sample.seq.mat[sample.seq.mat[,w.col.pos] == nuc,]) /
                                                                   denom.is),
                                                                (multip * nrow(control.seq.mat[control.seq.mat[,w.col.pos] == nuc,]) /
                                                                   denom.ran)
                                                    )
                                                    n.pair <- rep(0, length(nucs))
                                                    names(n.pair) <- nucs
                                                    
                                                  } else {
                                                    if(pos < 0) {
                                                      neg.pos <- pos
                                                      pos.pos <- abs(pos)
                                                      neg.nuc <- nuc
                                                      pos.nuc <- rc.nucs[nucs == nuc]
                                                    } else {
                                                      neg.pos <- -1 * pos
                                                      pos.pos <- pos
                                                      neg.nuc <- rc.nucs[nucs == nuc]
                                                      pos.nuc <- nuc
                                                    }
                                                    
                                                    #Select which positions are equal to relative positions
                                                    w.col.neg <- which(colnames(sample.seq.mat) == neg.pos)
                                                    w.col.pos <- which(colnames(sample.seq.mat) == pos.pos)
                                                    
                                                    #Select sequences with nucleotide at position
                                                    #.both object stores sequences with nucleotides at both positions
                                                    ##IS
                                                    n.seq.is.neg <- sample.seq.mat[sample.seq.mat[,w.col.neg] == neg.nuc,]
                                                    n.seq.is.pos <- sample.seq.mat[sample.seq.mat[,w.col.neg] != neg.nuc & sample.seq.mat[,w.col.pos] == pos.nuc,]
                                                    n.seq.is.both <- n.seq.is.neg[n.seq.is.neg[,w.col.pos] == pos.nuc,]
                                                    ##Control group
                                                    n.seq.ctrl.neg <- control.seq.mat[control.seq.mat[,w.col.neg] == neg.nuc,]
                                                    n.seq.ctrl.pos <- control.seq.mat[control.seq.mat[,w.col.neg] != neg.nuc & control.seq.mat[,w.col.pos] == pos.nuc,]
                                                    n.seq.ctrl.both <- n.seq.ctrl.neg[n.seq.ctrl.neg[,w.col.pos] == pos.nuc,]
                                                    
                                                    #Calculate number of sequences
                                                    n.freq <- c((multip * sum(nrow(n.seq.is.neg), nrow(n.seq.is.pos)) / denom.is),
                                                                (multip * sum(nrow(n.seq.ctrl.neg), nrow(n.seq.ctrl.pos)) / denom.ran)
                                                    )
                                                    
                                                    #Calculate frequency of other nucleotides at position when nucleotide is "n"
                                                    #Calculated for RC - thus complementary nucleotides need to be selected
                                                    n.pair <- sapply(nucs,
                                                                     function(nuc.pair) {
                                                                       rc.nuc <- rc.nucs[nuc.pair == nucs]
                                                                       if(nuc.pair == nuc) {
                                                                         nrow(n.seq.is.both)
                                                                       } else {
                                                                         sum(nrow(n.seq.is.neg[n.seq.is.neg[,w.col.pos] == rc.nuc,]),
                                                                             nrow(n.seq.is.pos[n.seq.is.pos[,w.col.neg] == nuc.pair,])
                                                                         )
                                                                       }
                                                                     })
                                                    #percent of all sequences
                                                    n.pair <- 100 * n.pair / all.is.count
                                                    #Percent of sequences with "nuc"
                                                    #n.pair <- 100 * n.pair / sum(nrow(n.seq.is.neg), nrow(n.seq.is.pos))
                                                    
                                                  }
                                                  if(output.number == "percent") {
                                                    names(n.freq) <- c("IS", "Ctrl")
                                                    n.freq
                                                  }
                                                  #Calculate frequencies and significance from total numbers
                                                  if(output.number == "count") {
                                                    names(n.freq) <- c("IS_num", "Ctrl_num")
                                                    num.is <- n.freq[grep("IS", names(n.freq))]
                                                    num.ctrl <- n.freq[grep("Ctrl", names(n.freq))]
                                                    
                                                    is.prop <- num.is / all.is.count 
                                                    ctrl.prop <- num.ctrl / all.ran.count
                                                    n.perc <- 100 * c(is.prop, ctrl.prop)
                                                    names(n.perc) <- c("IS_prop", "Ctrl_prop")
                                                    
                                                    #Statistics of fold change (enrich) & KLID & proportion test (z-statistics) & chi-square
                                                    #z.stat <- prop.test(n.freq[1], all.is.count, n.freq[2], correct = FALSE)
                                                    #chi.stat <- chisq.test(rbind(c(num.is, all.is.count - num.is),
                                                    #                             c(num.ctrl, all.ran.count - num.ctrl))
                                                    #                       )
                                                    #stat.val <- c(is.prop / ctrl.prop,
                                                    #              is.prop * log(is.prop / ctrl.prop),
                                                    #              round(z.stat$p.value, 3),
                                                    #              round(chi.stat$p.value, 3))
                                                    #names(stat.val) <- c("Enrich", "KLID", "Pz", "Pchi")
                                                    
                                                    
                                                    #Create vector of calculated values for nucleotide "n"
                                                    n.freq <- c(n.freq,
                                                                n.perc,
                                                                #stat.val,
                                                                n.pair)
                                                    n.freq
                                                  }
                                                  
                                                  #end of sapply(select.nuc)
                                                })
                         )
                         
                         rownames(n.freq.mat) <- select.nuc
                         #Return frequency matrix
                         t(n.freq.mat)
                         
                         #end of lapply(select.pos)
                       })
    names(pos.list) <- select.pos
    
    return(pos.list)
    
    print(paste0("Frequency of ", nuc, " at position ", pos,":"))
    print(n.freq)
    
  } else {
    print("IS and Random colnames are NOT eaqual!")
  }
}

#Caclulate the most or the least frequent nucleotides at STR site-relative positions

MostFreqDF <- function(list.of.matrices) {
  df <- as.data.frame(
    do.call(rbind,
            lapply(1:length(list.of.matrices),
                   function(x) {
                     lnam <- names(list.of.matrices)[x]
                     lmat <- list.of.matrices[[x]]
                     ldf <- as.data.frame(t(lmat))
                     
                     r.nam <- rownames(lmat)
                     
                     w.col <- which(ldf$IS_prop == max(ldf$IS_prop))
                     nuc.name <- rownames(ldf)[w.col]
                     
                     new.cols <- c(lnam, nuc.name)
                     names(new.cols) <- c("Pos", "N1")
                     
                     nuc.stats <- ldf[w.col,]
                     sel.nuc.stats <- unlist(c(nuc.stats[grep("prop", names(nuc.stats))]))
                     new.mat.1 <- matrix(rep(sel.nuc.stats,
                                             length(nucs)),
                                         ncol = length(sel.nuc.stats), byrow = TRUE
                     )
                     colnames(new.mat.1) <- names(sel.nuc.stats)
                     
                     new.mat.2 <- matrix(rep(new.cols, length(nucs)),
                                         ncol = 2, byrow = TRUE)
                     colnames(new.mat.2) <- c("Pos", "N1")
                     
                     new.mat.3 <- t(sapply(nucs,
                                           function(n) {
                                             a <- c(n, nuc.stats[names(nuc.stats) == n])
                                             names(a) <- c("N2", "N2_prop")
                                             a
                                           }))
                     
                     new.mat <- as.data.frame(cbind(new.mat.1, new.mat.2, new.mat.3))
                     rownames(new.mat) <- paste0(new.mat$tsdPos, new.mat$N1, new.mat$N2)
                     new.mat
                     
                   }
            ))
  )
  
  return(df)
  
}
LeastFreqDF <- function(list.of.matrices) {
  df <- as.data.frame(
    do.call(rbind,
            lapply(1:length(list.of.matrices),
                   function(x) {
                     lnam <- names(list.of.matrices)[x]
                     lmat <- list.of.matrices[[x]]
                     ldf <- as.data.frame(t(lmat))
                     
                     r.nam <- rownames(lmat)
                     
                     w.col <- which(ldf$IS_prop == min(ldf$IS_prop))
                     nuc.name <- rownames(ldf)[w.col]
                     
                     new.cols <- c(lnam, nuc.name)
                     names(new.cols) <- c("Pos", "N1")
                     
                     nuc.stats <- ldf[w.col,]
                     sel.nuc.stats <- unlist(c(nuc.stats[grep("prop", names(nuc.stats))]))
                     new.mat.1 <- matrix(rep(sel.nuc.stats,
                                             length(nucs)),
                                         ncol = length(sel.nuc.stats), byrow = TRUE
                     )
                     colnames(new.mat.1) <- names(sel.nuc.stats)
                     
                     new.mat.2 <- matrix(rep(new.cols, length(nucs)),
                                         ncol = 2, byrow = TRUE)
                     colnames(new.mat.2) <- c("Pos", "N1")
                     
                     new.mat.3 <- t(sapply(nucs,
                                           function(n) {
                                             a <- c(n, nuc.stats[names(nuc.stats) == n])
                                             names(a) <- c("N2", "N2_prop")
                                             a
                                           }))
                     
                     new.mat <- as.data.frame(cbind(new.mat.1, new.mat.2, new.mat.3))
                     rownames(new.mat) <- paste0(new.mat$tsdPos, new.mat$N1, new.mat$N2)
                     new.mat
                     
                   }
            ))
  )
  
  return(df)
  
}

#Function used to produce Figure 4
#Calculations of frequencies in repeat elements

RepTargetFreq <- function(rmsk.tab, select.column.name) {
  #select analyzed column
  w.col <- which(colnames(rmsk.tab) == select.column.name)
  #Create data.frame with number of hits
  df <- as.data.frame(do.call(rbind,
                              lapply(levels(rmsk.tab[, select.column.name]),
                                     function(x) {
                                       s <- sum(rmsk.tab$NumIS[rmsk.tab[, select.column.name] == x])
                                       p <- sum(rmsk.tab$PropIS[rmsk.tab[, select.column.name] == x])
                                       v <- c(x, s, round(100 * p, 1))
                                       names(v) <- c(select.column.name, "NumIS", "PercIS")
                                       v
                                     })
  ))
  df[,2] <- as.numeric(df[,2])
  df[,3] <- as.numeric(df[,3])
  df <- df[order(df[,2], decreasing = TRUE),]
  df[,1] <- factor(df[,1], levels = unique(df[,1]))
  return(df)
}

#Define function creating table of motifs
MotifTab <- function(sequence.table, motif.pos) {
  do.call(rbind,
          lapply(1:nrow(sequence.table),
                 function(s) {
                   seq.entry <- toupper(as.character(unlist(sequence.table[s,])))
                   n.mat <- NucMat(seq.entry)
                   do.call(cbind, lapply(1:length(seq.entry),
                                         function(w.mat) {
                                           if(length(seq.entry) == 1) {n.mat.x <- n.mat}
                                           else {n.mat.x <- n.mat[w.mat,]}
                                           
                                           data.frame(Motif = paste0(n.mat.x[motif.pos[1]],
                                                                     n.mat.x[motif.pos[2]], "..", 
                                                                     n.mat.x[motif.pos[3]], "...",
                                                                     n.mat.x[motif.pos[4]], "..",
                                                                     n.mat.x[motif.pos[5]],
                                                                     n.mat.x[motif.pos[6]]))
                                           
                                         }))
                 })
          
  )
}


#calculate membership function of sequence to all PPMs in mixture
#nuc.m = matrix of nucleotide sequences as creted by NucMat function
#ppm.list = list of PPMs (components)
#wm.m = named vector of weights of PPMs

MembershipFunction <- function(nuc.m, ppm.list, wm.m) {
  
  #Check formats of inputs
  if(any(class(is.nuc.mat) == "matrix")) {
  } else {
    stop("nuc.mat input is not a matrix!")
  }
  if(class(ppm.list) != "list") {
    if(class(ppm.list) == "data.frame") {
      ppm.df <- ppm.list
      ppm.list <- lapply(unique(ppm.df$PPM),
                         function(ppm.name) {
                           p.l <- ppm.df[ppm.df$PPM == ppm.name,
                                         which(colnames(ppm.df) != "PPM")]
                           p.l
                         })
      names(ppm.list) <- unique(ppm.df$PPM)
    } else {
      stop("ppm.list input is not a list nor data.frame!")
    }
  }
  if(class(wm.m) != "numeric") {
    stop("wm.m input is not numeric!")
  }
  
  #Check if all PPMs have weights in wm.m
  if(length(ppm.list) != length(wm.m)) {
    stop("PPM list and weights are not of same length!")
  }
  
  #Check if all names of PPMs are in weigths
  if(all(names(ppm.list) %in% names(wm.m)) == FALSE) {
    stop("Some of PPM names are not in weigth vecor!")
  }
  
  #Calculate
  
  #P(x|m)
  pseq <- SeqProbByPPM(nuc.m, ppm.list)
  
  #calculate weighted probability
  w.p <- sapply(1:ncol(pseq),
                function(x) {
                  p.nam <- colnames(pseq)[x]
                  p <- pseq[,colnames(pseq) == p.nam]
                  w <- wm.m[names(wm.m) == p.nam]
                  p * w
                })
  
  sum.w.p <- apply(w.p, 1, sum)
  
  m.ship <- t(sapply(1:nrow(w.p),
                     function(x) {
                       w.p[x,] / sum.w.p[x] 
                     }))
  
  colnames(m.ship) <- colnames(pseq)
  
  return(m.ship)
}

#Relative position according to axis of symetry
#Left side is negative
RelPosNames <- function(input) {
  if(length(input) > 0) {
    if(length(ncol(input)) == 0) {
      n.col <- length(input)
    } else {
      n.col <- ncol(input)
    }
    
    #calculate values of relative positions
    if(n.col %% 2 == 0) {
      n.col.half <- n.col / 2
      n.pos <- c( seq(-1 * n.col.half, -1),
                  seq(1, n.col.half)
      )
    } else {
      
      n.col.half <- floor(n.col / 2)
      n.pos <- c( seq(-1 * n.col.half, -1),
                  0,
                  seq(1, n.col.half)
      )
    }
    
    return(n.pos)
    
    
  } else {
    print("Input has length 0. No column names can be created.")
  }
}


#Format PPM
#Resulting in nucleotides per rows and position per column
PPMformat <- function(ppm) {
  c.nam <- colnames(ppm)
  if(any(c.nam %in% nucs)) {
    
    #Search for nucletide names in culumn names
    if(all(nucs %in% c.nam)) {
      ppm <- t(ppm[,which(c.nam %in% nucs)])
      ppm <- ppm[order(rownames(ppm)),]

    } else {
      message("Only partial match between nucleotides and column names of PPM")
    }
    
  } else {
    #Search for nucletide names in rows
    r.nam <- rownames(ppm)
    if(any(r.nam %in% nucs)) {
      
      if(all(nucs %in% r.nam)) {
        ppm <- ppm[which(r.nam %in% nucs),]
        ppm <- ppm[order(rownames(ppm)),]

      } else {
        message("Only partial match between nucleotides and row names of PPM")
      }
      
    } else {
      message("No match between nucleotides and row or column names of PPM")
    }
  }
  
  #Set column names for PPM
  #Relative position according to axis of symetry
  colnames(ppm) <- RelPosNames(ppm)
  return(ppm)
}

#Break the sequence down to nucleotide matrix
NucMat <- function(seq.vec) {
  
  if(class(seq.vec) != "character") {
    stop(print(paste0("Sequence input in class ", class(seq.vec), " ! Must be character!")))
  }

  is.nuc.l <- strsplit(seq.vec, split = "")
  is.nuc.mat <- matrix(unlist(is.nuc.l, use.names=FALSE),
                       ncol = length(is.nuc.l[[1]]), byrow = TRUE)
  
  return(is.nuc.mat)
}

#Create PPM from nucleotide matrix (NucMat function)
#Create PPM
PPMfromSeqMat <- function(nucleotide.matrix) {
  if(is.matrix(nucleotide.matrix)) {
    nucs <- c("A", "C", "G", "T")
    ppm <- cbind(apply(nucleotide.matrix, 2,
                       function(x) {
                         sapply(nucs,
                                function(y) {
                                  length(x[x == y]) / length(x)
                                })
                       }))
    #ppm <- PPMformat(ppm)
    return(ppm)
  } else {
    stop("Input is not matrix!")
  }
}

#CALCULATE PERCENTAGE OF BASES IN SEQUENCES
NucFreq <- function(seq.vec) {
  
  if(class(seq.vec) != "character") {
    stop(print(paste0("Sequence input in class ", class(seq.vec), " ! Must be character!")))
  }

  #Create matrix of sequences where columns represent position in sequence
  
  is.nuc.mat <- NucMat(seq.vec)

  #For every row, count ratio of A,C,G,T
  nucs <- c("A", "C", "G", "T")
  print("Calculating frequencies for: ")
  for(i in 1:length(nucs)) {
    print(nucs[i])
    f <- apply(is.nuc.mat, 1, function(x) {
      length(x[x == nucs[i]]) / length(x)
    })
    if(i == 1) {
      nuc.freq <- f
    } else {
      nuc.freq <- cbind(nuc.freq,
                        f)
    }
  }
  
  colnames(nuc.freq) <- nucs
  rownames(nuc.freq) <- seq.vec
  
  print("NucFreq clculated")
  print(head(nuc.freq))
  
  return(nuc.freq)
  
}

#CALCULATE SEQUENCE PROBABILITY FROM PPM
#nuc.mat is culeotide matrix created by NucMat function
SeqProbByPPM <- function(nuc.mat, ppm.list) {
  
  print(paste0("Started at ", Sys.time()))

  for(i in 1:length(ppm.list)) {
    
    ppm <- t(ppm.list[[i]][,2:5])
    
    px.prod <- sapply(1:nrow(nuc.mat), function(x.seq) {
      
      x.vec <- nuc.mat[x.seq,]
      x.prob <- rep(NA, length(x.vec))
      ppm.nuc <- 1:nrow(ppm)
      nucs <- rownames(ppm)
      
      pos1 <- which(x.vec == nucs[1])
      pos2 <- which(x.vec == nucs[2])
      pos3 <- which(x.vec == nucs[3])
      pos4 <- which(x.vec == nucs[4])
      
      x.prob[pos1] <- ppm[1, pos1]
      x.prob[pos2] <- ppm[2, pos2]
      x.prob[pos3] <- ppm[3, pos3]
      x.prob[pos4] <- ppm[4, pos4]
      
      prod(x.prob)
    })
    #pseq.part <- data.frame(Pseq = px.prod,
    #                        PPM = names(ppm.list)[i])
    
    if(i == 1) {
      pseq <- px.prod
    } else {
      pseq <- cbind(pseq,
                    px.prod)
    }
  }
  
  colnames(pseq) <- names(ppm.list)
  
  print(paste0("Done    at ", Sys.time()))
  
  return(pseq)
}

#Calculate KLID and SKLID
#KLID
KLID <- function(test.ppm, bground.ppm) {
  klid.ppm <- test.ppm * log(test.ppm / bground.ppm)

  return(apply(klid.ppm, 2, function(x) {
    sum(x, na.rm = TRUE)
  })
  )
}

#SKLID
SKLID <- function(test.ppm, bground.ppm) {
  p <- test.ppm * log(test.ppm / bground.ppm)
  q <- bground.ppm * log(bground.ppm / test.ppm)
  q[!is.finite(q)] <- NA
  sklid.ppm <- p + q
  
  return(apply(sklid.ppm, 2, function(x) {
           sum(x, na.rm = TRUE)
           })
  )
}

#Transform data to reverse complement
#run.on can be "PPM" or "SEQ" according to input data
#SEQ should be in nuc.mat format
FWtoRC <- function(input.data, input.format = "PPM") {
  nucs <- c("A","C","G","T")
  #For SEQ
  if(input.format == "SEQ") {
    c.nam <- colnames(input.data)
    
    w.a <- (which(input.data == "A"))
    w.c <- (which(input.data == "C"))
    w.g <- (which(input.data == "G"))
    w.t <- (which(input.data == "T"))
    input.data[w.a] <- "T"
    input.data[w.c] <- "G"
    input.data[w.g] <- "C"
    input.data[w.t] <- "A"
    colnames(input.data) <- c.nam
    return(input.data)
  }
  
  #For PPM
  if(input.format == "PPM") {
    r.nam <- rownames(input.data)
    if(all(nucs %in% r.nam )) {
    } else {
      if(any(nucs %in% r.nam)) {
        message(paste0("Only ",
                       paste(nucs[which(nucs %in% r.nam)], collapse = ", "),
                       " found in row names!"))
      } else {
        message("No nucletodies found in row names!")
        print(class(input.data))
        print(head(input.data))
        stop()
      }
    }
    w.rows <- which(r.nam %in% nucs)
    to.rev <- input.data[w.rows,]
    c.nam <- colnames(to.rev)
    to.rev <- to.rev[,rev(1:ncol(to.rev))]
    r.nam <- rownames(to.rev)
    w.a <- (which(r.nam == "A"))
    w.c <- (which(r.nam == "C"))
    w.g <- (which(r.nam == "G"))
    w.t <- (which(r.nam == "T"))
    rownames(to.rev)[w.a] <- "T"
    rownames(to.rev)[w.c] <- "G"
    rownames(to.rev)[w.g] <- "C"
    rownames(to.rev)[w.t] <- "A"
    colnames(to.rev) <- c.nam
    return(to.rev)
  }
}

#PPMdist - calculate distances of probability matrices
#dist.input: is a data.frame with PPMs
#ppm1.name, ppm2.name: names of PPM for distance calculation
#ppm2 can be transformed to reverse-complement if ppm.direct == "RC"
#ppm.direct: "FW" for forward "RC" for Reverse-complement
#dist.pos: "all" for calculation of all positions or
#select.pos: range of distances from axis of symmetry (-6:6)
PPMdist <- function(dist.input, ppm.direct = "FW", select.pos = "all") {
  if(class(dist.input) != "data.frame") {
    message("Input is not of class data.frame")
    stop()
  }
  
  dist.tab <-   do.call(rbind,
                        lapply(levels(dist.input$PPM),
                               function(p.nam.1) {
                                 
                                 d.tab <- do.call(rbind,
                                                  lapply(levels(dist.input$PPM),
                                                         function(p.nam.2) {
                                                           #Select ppm1 and ppm2 by PPM names in PPM column
                                                           ##PPM1 - always in FW orientation
                                                           ppm1 <- PPMformat(dist.input[dist.input$PPM == p.nam.1,])
                                                           ##PPM2 - FW or transform to RC
                                                           ppm2 <- PPMformat(dist.input[dist.input$PPM == p.nam.2,])
                                                           if(ppm.direct %in% c("FW", "RC")) {
                                                             if(ppm.direct == "FW") {
                                                               ppm2 <- PPMformat(ppm2)
                                                             }
                                                             if(ppm.direct == "RC") {
                                                               ppm2 <- FWtoRC(ppm2, input.format = "PPM")
                                                               ppm2 <- PPMformat(ppm2)
                                                             }
                                                           } else {
                                                             stop("No valid value for ppm.direct - should be  FW or RC")
                                                           }
                                                           #Calculate distance
                                                           #Calculate for all positions or only for ones specified in selct.pos object
                                                           if(all(select.pos == "all")) {
                                                             ppm.dist <- sum(apply(abs(ppm1 - ppm2), 1, mean))
                                                             
                                                           } else {
                                                             ppm1 <- ppm1[,which(colnames(ppm1) %in% as.character(select.pos))]
                                                             ppm2 <- ppm2[,which(colnames(ppm2) %in% as.character(select.pos))]
                                                             ppm.dist <- sum(apply(abs(ppm1 - ppm2), 1, mean))
                                                           }
                                                           
                                                           #Include distance type
                                                           #If ppm.direct = RC and PPM names are equal, mark distance as PDef (palindrome deficit)
                                                           d.type <- ppm.direct
                                                           if(d.type == "RC" & p.nam.1 == p.nam.2) {
                                                             d.type <- "PDef"
                                                           }
                                                           
                                                           #Transform distance data to data.frame
                                                           ppm.dist.tab <- data.frame(PPM1 = p.nam.1,
                                                                                      PPM2 = p.nam.2,
                                                                                      PPMdist = ppm.dist,
                                                                                      DistType = d.type)
                                                           #output
                                                           ppm.dist.tab
                                                         }))
                               }))
  
  return(dist.tab)
  
}

#Create ReverseComplement of sequences
#input.seq should be of characer string vector class
SeqRC <- function(input.seq) {
  nuc.mat <- NucMat(input.seq)
  
  w.a <- which(nuc.mat == "A")
  w.c <- which(nuc.mat == "C")
  w.g <- which(nuc.mat == "G")
  w.t <- which(nuc.mat == "T")
  
  nuc.mat[w.a] <- "T"
  nuc.mat[w.t] <- "A"
  nuc.mat[w.c] <- "G"
  nuc.mat[w.g] <- "C"
  nuc.mat <- nuc.mat[,rev(1:ncol(nuc.mat))]
  input.seq.rc <- apply(nuc.mat, 1, function(x) {
    paste(x, collapse = "")
  })
}

#Function for creating KLID logo ggplot
#ppm.name = name of PPM in the format that is present in input ppm.list
#ppm.direct = logo will be plotted for forward "f" or reverse-complement "rc" orientation
#ppm.ran = PPM of random genomic sequences
#ranmat = if "mean" is entered KLID is calculated against mean nucleotide frequency

plotKLIDlogo <- function(ppm.name, ppm.direct = "f", ppm.ran, ranmat = "mean") {
  #Calculate Logo values
  nucs <- c("A", "C", "G", "T")
  
  #If mean of nucleotide frequencies is desired, enter "mean"
  ranmat <- "mean"
  
  #----
  ppm <- as.matrix(ppm.list[[which(names(ppm.list) == ppm.name)]])
  ppm <- t(ppm[,2:5])
  ppm <- PPMformat(ppm)
  ppm.ran <- ppm.ran[order(rownames(ppm.ran)),]
  
  if(ranmat == "mean") {
    #This creates the matrix where all columns are the same
    #each column represents means of nucletides frequency
    ran.nucs <- rownames(ppm.ran)
    ppm.ran <- matrix(apply(matrix(apply(ppm.ran, 1, mean),
                                   ncol = 1),
                            2, rep, ncol(ppm)),
                      ncol = ncol(ppm)
    )
    rownames(ppm.ran) <- ran.nucs
  }
  
  ppm <- ppm[order(rownames(ppm)),]
  class(ppm) <- "numeric"
  class(ppm.ran) <- "numeric"
  
  #Check the desired orientation of PPM
  #Change to RC if specified
  if(ppm.direct == "f") {
    
  } else {
    if(ppm.direct == "rc") {
      ppm <- ppm[,rev(1:ncol(ppm))]
      rownames(ppm) <- rev(rownames(ppm))
      
    } else {
      stop("Orientation of PPM needs to be specified: f or rc.")
    }
  }
  
  #Calculate KLID velues
  klid.ppm <- ppm * log(ppm / ppm.ran)
  
  y.min <- min(
    apply(klid.ppm[,11:16], 2, function(x) {
      x <- x[x < 0]
      sum(x, na.rm = TRUE)
    })
    , na.rm = TRUE)
  
  y.max <- max(
    apply(klid.ppm, 2, function(x) {
      x <- x[x >= 0]
      sum(x, na.rm = TRUE)
    })
    , na.rm = TRUE)
  if(y.max < 1) {
    y.max <- 1
  }
  
  #PLOT
  
  #Set the boundaries for TSD rectangle
  if(tsd %% 2 == 0) {
    half.tsd <- tsd/2
    is.middle <- ncol(ppm)/2
    tsd.range <- c(is.middle - half.tsd + 0.5,
                   is.middle + half.tsd + 0.5
    )
  } else {
    half.tsd <- tsd/2
    is.middle <- ceiling(ncol(ppm)/2)
    tsd.range <- c(is.middle - half.tsd,
                   is.middle + half.tsd
    )
  }
  
  y.max.rect <- max(
    apply(klid.ppm[,(tsd.range[1]+0.5):(tsd.range[2]-0.5)], 2, function(x) {
      x <- x[x >= 0]
      sum(x, na.rm = TRUE)
    })
    , na.rm = TRUE)
  
  #Create logo
  p <-
    ggplot() +
    geom_hline(xintercept = is.middle) +
    annotate('rect', xmin = tsd.range[1], xmax = tsd.range[2], ymin = y.min, ymax = y.max.rect, alpha = .5, col='gray', fill='gray') +
    geom_logo(klid.ppm, method = "custom") +
    ylab("KLID") +
    #scale_x_continuous(breaks = 1:ncol(ppm), labels = as.numeric(colnames(ppm))) +
    coord_cartesian(ylim = c(-0.25, y.max)) +
    geom_hline(yintercept=0) +
    ggtitle(ppm.name) +
    scale_y_continuous(breaks = seq(-0.2, y.max, by = 0.2), labels = seq(-0.2, y.max, by = 0.2)) +
    theme_logo() +
    theme(axis.text.x = element_blank())
  
  as_grob(p)
}
