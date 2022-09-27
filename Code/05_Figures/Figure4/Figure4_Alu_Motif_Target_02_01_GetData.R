#Quantify targeting of intra-Alu hotspot
#Load or Calculate:
# i) number of IS in Alu repeats
# ii) IS distances to the palindromic motif
# iii) Frequences of all IS intra-Alu motifs

wd.is.dist <- paste0(result.wd, "HotSpot/")

#Set names of studies used
studies <- c("Zhyvoloup", "Vansant", "Demeulemeester")

smpl.order <- c(1:6,12:17,8:11,7)

#Load and format table with Alu target frequency
#Derive data.frame smpl.tab from frequency table
##Alu freq
alu.f.tab <- do.call(rbind,
                     lapply(studies,
                            function(s.name) {
                              #Set path to study disrectory
                              s.dir <- paste0(wd.is.dist, "HIV_", s.name)
                              #Find file with frequencies
                              freq.file <- list.files(path = s.dir,
                                                      pattern = "*_rmsk_Alu_target_freq.txt")
                              #Load table
                              df <- read.table(paste0(s.dir, "/", freq.file),
                                               header = FALSE, quote = "", sep = "\t", stringsAsFactors = TRUE,
                                               col.names = c("Group", "Virus", "Alu", "All"))
                              df$Perc <- round(100 * df$Alu / df$All, digits = 2)
                              df$Study <- factor(s.name)
                              
                              df$Sample <- sapply(as.character(df$Virus),
                                                  function(v) {
                                                    smpl1 <- unlist(strsplit(v, paste0("HIV_", s.name,"[_]?")))[2]
                                                    smpl2 <- unlist(strsplit(smpl1, "_"))
                                                    if(is.na(smpl1)) {smpl1 <- "WT"}
                                                    if(length(smpl2) > 2 & smpl2[4] == "mix") {
                                                      if(smpl2[2] == "KD") {
                                                        smpl1 <- "LEDGF-KD"
                                                      } else {
                                                        if(smpl2[3] == 0) {
                                                          smpl1 <- "WT"
                                                        } else {
                                                          smpl1 <- "LEDGIN"
                                                        }
                                                      }
                                                    }
                                                    smpl1
                                                  }, USE.NAMES = FALSE)
                              #output
                              df
                            }))
alu.f.tab$Sample <- factor(alu.f.tab$Sample, levels = unique(alu.f.tab$Sample)[c(1:6,12:17,8:11,7)])
##smpl.tab
smpl.tab <- unique(alu.f.tab[,c("Virus", "Study", "Sample")])

smpl.names <- levels(smpl.tab$Sample)

#Set colors for samples
#colorblind-friendly palette is used
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
names(cbPalette) <- c("gray", "lightorange", "lightblue", "green", "yellow", "blue", "orange", "magenta")

smpl.cols.all <- c("#D55E00", rep("#0072B2", 11), rep("#009E73", 4), "#F0E442")
names(smpl.cols.all) <- smpl.names
cols <- c(smpl.cols.all, "#999999")
names(cols) <- c(smpl.names, "Random")

ctrl.col <- "#999999"


#Calculate frequencies of distances
motif.name <- "CT..G...C..AG"
y.lab.2 <- "intra-Alu IS in motif"
dist.range <- c(-33, 33)
dist.bin <- 1
y.range <- c(0, 10)
is.denom.factor <- 100

plot.legend <- FALSE

#Define distance windows
d.cat <- seq(dist.range[1], dist.range[2], by = dist.bin)

#Distribution of distances to the motif
d.freq <- do.call(rbind,
                  lapply(1:nrow(smpl.tab),
                         function(w.r) {
                           
                           v.name <- smpl.tab$Virus[w.r]
                           s.name <- smpl.tab$Study[w.r]
                           smpl <- smpl.tab$Sample[w.r]
                           
                           #Set path to study disrectory
                           s.dir <- paste0(wd.is.dist, "HIV_", s.name)
                           
                           #Find files with distances to motif
                           is.set.name <- list.files(path = s.dir,
                                                     pattern = paste0(v.name, "_in_.*_rmsk_Alu_palMotif_distance.txt"))
                           ctrl.set.name <- list.files(path = s.dir,
                                                       pattern = paste0(v.name, "_is12_.*_rmsk_Alu_shuffle_same_palMotif_distance_100x_tab.txt"))
                           
                           d.tab.c.names <- c(paste0("IS_", c("chrom", "start", "end", "name", "rep", "strand")),
                                              paste0("Motif_", c("chrom", "start", "end", "name", "AluStDist", "strand")),
                                              "Distance")
                           d.tab <- read.table(paste0(s.dir, "/", is.set.name),
                                               header = FALSE, quote = "", sep = "\t", stringsAsFactors = TRUE,
                                               col.names = d.tab.c.names)
                           d.tab.ctrl <- read.table(paste0(s.dir, "/",ctrl.set.name),
                                                    header = FALSE, quote = "", sep = "\t", stringsAsFactors = TRUE)
                           
                           
                           total.is <- nrow(d.tab)
                           total.ctrl <- nrow(d.tab.ctrl)
                           
                           all.is.dist <- d.tab$Distance
                           
                           d.freq.v <- as.data.frame(do.call(rbind, lapply(d.cat,
                                                                           function(d) {
                                                                             #Calculate number of IS in distance window
                                                                             #dist.bin needs to be >= 1
                                                                             if(dist.bin >= 1) {
                                                                               if(dist.bin == 1) {
                                                                                 #If distance is single value
                                                                                 count.is <- length(all.is.dist[all.is.dist == d])
                                                                                 count.ran <-apply(as.matrix(d.tab.ctrl), 2,
                                                                                                   function(dist.i) {
                                                                                                     length(dist.i[dist.i == d])
                                                                                                   })
                                                                               } else {
                                                                                 #If distance is window of length > 1
                                                                                 d.start <- d - (dist.bin / 2)
                                                                                 d.end <- d + (dist.bin / 2)
                                                                                 count.is <- length(all.is.dist[all.is.dist >= d.start & all.is.dist < d.end])
                                                                                 count.ran <- apply(as.matrix(d.tab.ctrl), 2,
                                                                                                    function(dist.i) {
                                                                                                      length(dist.i[dist.i >= d.start & dist.i < d.end])
                                                                                                    })
                                                                               }
                                                                             } else {
                                                                               stop("dist.bin object needs to be at least 1")
                                                                             }
                                                                             #Create data.frame with frequency values
                                                                             #Perc column contains percentages if is.denom.factor == 100
                                                                             data.frame(Group = c("Shuffle", "IS"),
                                                                                        Distance = c(d,d),
                                                                                        Count = c(mean(count.ran),
                                                                                                  count.is),
                                                                                        SD = c(sd(count.ran),
                                                                                               0),
                                                                                        Perc = c(is.denom.factor * mean(count.ran) / total.ctrl,
                                                                                                 is.denom.factor * count.is / total.is),
                                                                                        SD_Perc = c(is.denom.factor * sd(count.ran) / total.ctrl,
                                                                                                    0)
                                                                             )
                                                                           })))
                           
                           d.freq.v$Study <- s.name
                           d.freq.v$Sample <- smpl
                           d.freq.v
                           
                         }))

#Calculate motif frequencies
#Get motif variants and their frequencies
#Sequences are extracted from sample-associated sequences
#inside lapply function: reads "fa_tab.txt" files of IS and shuffled controls
#motif.pos = which positions in target sequences are considered to create motif
motif.pos <- c(1,2,5,9,12,13)
m.freq <- do.call(rbind,
                  lapply(studies,
                         function(s.name) {
                           #Set path to study disrectory
                           s.dir <- paste0(wd.is.dist, "HIV_", s.name)
                           #Find file with frequencies
                           freq.file <- list.files(path = s.dir,
                                                   pattern = "*_MotifSeq_target.txt")
                           
                           if(length(freq.file) != 0) {
                             
                             comp.motif.df <- read.table(paste0(s.dir, "/", freq.file),
                                                         header = TRUE, quote = "", sep = " ", stringsAsFactors = FALSE)
                             
                           } else {
                             print("Calculate motif frequencies across samples")
                             print(paste0("Study: ", s.name))
                             smpl.tab.study <- smpl.tab[smpl.tab$Study == s.name,]
                             comp.motif.df <- do.call(rbind,
                                                      lapply(1:nrow(smpl.tab.study),
                                                             function(smpl.n) {
                                                               #Load data for given sample
                                                               seq.tab <- read.table(paste0(s.dir, "/",
                                                                                            list.files(path = s.dir,
                                                                                                       pattern = paste0(smpl.tab.study$Virus[smpl.n],"_is.*rmsk_Alu_fa_tab.txt"))),
                                                                                     header = FALSE, quote = "", sep = "\t", stringsAsFactors = TRUE,
                                                                                     col.names = c("Name", "Seq"))
                                                               #Load shuffled control data
                                                               seq.tab.ran <- read.table(paste0(s.dir, "/",
                                                                                                list.files(path = s.dir,
                                                                                                           pattern = paste0(smpl.tab.study$Virus[smpl.n],"_is.*rmsk_Alu_shuffle_all_fa_tab.txt"))),
                                                                                         header = FALSE, quote = "", sep = "\t", stringsAsFactors = TRUE
                                                               )
                                                               
                                                               
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
                                                               
                                                               #Extract motifs from sequences
                                                               if(any(colnames(seq.tab) == "Seq")) {seq.tab <- data.frame(Seq = as.character(seq.tab[,which(colnames(seq.tab) == "Seq")]))}
                                                               if(any(colnames(seq.tab.ran) == "Seq")) {seq.tab.ran <- data.frame(Seq = as.character(seq.tab.ran[,which(colnames(seq.tab.ran) == "Seq")]))}
                                                               
                                                               motif.df <- MotifTab(seq.tab, motif.pos)
                                                               
                                                               #Calculate occurence of each motif
                                                               motif.uniq.df <- data.frame(Motif = motif.df$Motif[!duplicated(motif.df$Motif)])
                                                               motif.uniq.df$Count <- sapply(motif.uniq.df$Motif, function(x) {length(motif.df[motif.df$Motif == x,])})
                                                               motif.uniq.df$Motif <- factor(motif.uniq.df$Motif, levels = motif.uniq.df$Motif[order(motif.uniq.df$Count, decreasing = TRUE)])
                                                               
                                                               comp.motif.uniq.df <- motif.uniq.df
                                                               
                                                               total <- sum(comp.motif.uniq.df$Count)
                                                               
                                                               #Add data about random frequency of found motifs
                                                               #Create object with number of all sequences/shuffled control
                                                               total.ran <- nrow(seq.tab.ran)
                                                               
                                                               #Transform sequences to random motifs
                                                               motif.ran.df <- MotifTab(seq.tab.ran, motif.pos)
                                                               
                                                               #Create vector of all random motifs
                                                               motif.uniq.ran <- sort(unique(as.character(unlist(motif.ran.df))))
                                                               #Caclulate frequencies of each motif in each repeat(column) of shuffled control
                                                               motif.uniq.ran.freq <- do.call(rbind,
                                                                                              lapply(motif.uniq.ran,
                                                                                                     function(mtf) {
                                                                                                       sapply(1:ncol(motif.ran.df),
                                                                                                              function(w.col) {
                                                                                                                col.mtfs <- as.character(motif.ran.df[,w.col])
                                                                                                                length(which(col.mtfs == mtf))
                                                                                                              })
                                                                                                     }))
                                                               #Calculate mean, min and max frequency of every motif occuring in shuffled control
                                                               motif.uniq.ran.df <- data.frame(Motif = motif.uniq.ran,
                                                                                               Mean = apply(motif.uniq.ran.freq, 1, mean),
                                                                                               Min = apply(motif.uniq.ran.freq, 1, min),
                                                                                               Max = apply(motif.uniq.ran.freq, 1, max))
                                                               
                                                               comp.motif.uniq.df$Count_Shuffle <- sapply(as.character(motif.uniq.df$Motif),
                                                                                                          function(mtf) {
                                                                                                            w.row <- which(motif.uniq.ran.df$Motif == mtf)
                                                                                                            if(length(w.row) == 1) {
                                                                                                              motif.uniq.ran.df$Mean[w.row]
                                                                                                            } else {
                                                                                                              0
                                                                                                            }})
                                                               
                                                               comp.motif.uniq.df$Perc <- round(100 * comp.motif.uniq.df$Count / total, digits = 2)
                                                               comp.motif.uniq.df$Perc_Shuffle <- round(100 * comp.motif.uniq.df$Count_Shuffle / total, digits = 2)
                                                               
                                                               comp.motif.uniq.df$Sample <- smpl.tab.study$Sample[smpl.n]
                                                               
                                                               comp.motif.uniq.df
                                                               
                                                             }))
                             print("Save motif frequencies across samples")
                             write.table(comp.motif.df,
                                         paste0(s.dir, "/", s.name, "_MotifSeq_target.txt"),
                                         quote = FALSE, col.names = TRUE)
                             
                           }
                           
                           #add column w/ study name
                           comp.motif.df$Study <- factor(s.name)
                           #output
                           comp.motif.df
                           
                           
                         }))

