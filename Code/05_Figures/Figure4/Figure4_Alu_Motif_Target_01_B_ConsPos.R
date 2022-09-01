#Plot position of IS in Alu consensus sequence

#Load Table with positions in consensus sequences
#is.set.name <- paste0(virus, "_is20_in_rmsk_Alu_consensus.bed")

pos.freq.c.names <- c("start", "end", "count")

pos.freq <- read.table(paste0(result.wd, "HotSpot/",
                                 paste0(virus, "_M08_PPM07_in_hg38_rmsk_Alu_to_consensus_PosFreq.txt")),
                          header = FALSE, quote = "", sep = "\t", stringsAsFactors = TRUE,
                          col.names = pos.freq.c.names)

pos.freq$RangeLength <- pos.freq$end - pos.freq$start
pos.freq$RangeMiddle <- round(pos.freq$start + (pos.freq$end - pos.freq$start)/2)

head(pos.freq)

total.is <- sum(pos.freq$count)

#Calculate frequencies in set windows
#wnd.width - width of window (bumber of positions)
#rng.start - start of the analyzed range
#rng.end - end of the analyzed  range

wnd.width <- 1
rng.start <- 0
rng.end <- 290

rng <- seq(rng.start, rng.end, by = wnd.width)

#calculate frequencies in ranges
freq.in.rng <-   sapply(rng,
                        function(r) {
                          if(wnd.width == 1) {
                            100 * sum(pos.freq$count[pos.freq$RangeMiddle == r]) / total.is
                          } else {
                            if(r < max(rng)) {
                              100 * sum(pos.freq$count[pos.freq$RangeMiddle < (r + wnd.width) &
                                                    pos.freq$RangeMiddle >= r]) / total.is
                            } else {
                              #When r is the last value of the range, include that distance
                              100 * sum(pos.freq$count[pos.freq$RangeMiddle <= (r + wnd.width) &
                                                    pos.freq$RangeMiddle >= r]) / total.is
                              
                            }
                          }
                        })

head(freq.in.rng)

f.tab <- data.frame(Pos = factor(rng),
                    Freq = freq.in.rng)

#Select top frequent entries
top.f.tab <- head(f.tab[order(f.tab$Freq, decreasing = TRUE),],2)

cols <- c("#D55E00")

#Load sequences of the most frequent positions
top.seq <- read.table(paste0(result.wd, "HotSpot/", cm.ppm, "/",
                              paste0(virus, "_M08_PPM07_in_hg38_rmsk_Alu_to_consensus_TopConsSeq.txt")),
                       header = FALSE, quote = "", sep = "\t", stringsAsFactors = FALSE,
                       col.names = c("Pos", "Seq"))
seq.len <- 12
top.seq$Seq_Short <- sapply(top.seq$Seq,
                      function(s) {
                        nuc.vec <- unlist(strsplit(s, ""))
                        nuc.vec.mid <- ceiling(length(nuc.vec)/2)
                        paste0(nuc.vec[(nuc.vec.mid - (seq.len/2)):(nuc.vec.mid + (seq.len/2))], collapse = "")
                      })

#Plot IS freq in Alu percentage
p.cons <- ggplot(data = f.tab, aes(x = Pos, y = Freq)) +
  geom_col(colour = cols[1], fill = alpha(cols[1], .3)) +
  coord_cartesian(ylim = c(0, 80)) +
  scale_y_continuous(expand = c(0,0)) +
  scale_x_discrete(breaks = seq(50, 300, by = 50)) +
  ylab("% Mapped\nSequences") +
  xlab("Alu consensus position") +
  annotate("text" ,x=240, y = top.f.tab$Freq[1] + 10, label=top.seq$Seq_Short[1],
           size = 3, family = "mono", fontface = "bold") +
  annotate("text", x = 80, y = top.f.tab$Freq[2] + 10, label=top.seq$Seq_Short[2],
           size = 3, family = "mono", fontface = "bold") +
  theme_classic() +
  theme(axis.line.x = element_blank(),
        #axis.ticks.x = element_blank(),
        axis.title.x = element_text(colour = "black", face = "bold", size=8),
        axis.title.y = element_text(colour = "black", face = "bold", size=8),
        axis.text.x = element_text(colour = "black", face = "bold", size=8, angle = 0),
        axis.text.y = element_text(colour = "black", face = "bold", size=8, angle = 0),
        text = element_text(colour = "black"),
        legend.title = element_blank(),
  )
