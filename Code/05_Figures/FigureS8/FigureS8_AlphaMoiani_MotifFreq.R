#Calculate frequencies of Motif containing sequences in set

mot.freq.df <- read.table(paste0(result.wd, "HotSpot/Alpha_Moiani/Alpha_Moiani_is26_Motif.txt"),
                       sep = "\t", col.names = c("IS", "Seq", "Motif", "Strand"), stringsAsFactors = TRUE)

#calculate frequency of the motifs
mot.freq <- do.call(rbind,
                    lapply(levels(mot.freq.df$Motif),
                         function(mot) {
                           mot.total <- nrow(mot.freq.df[mot.freq.df$Motif == mot,])
                           do.call(rbind,
                                   lapply(levels(mot.freq.df$Strand),
                                        function(strnd) {
                                          data.frame(Motif = mot,
                                                     Strand = strnd,
                                                     Freq = 100 * nrow(mot.freq.df[mot.freq.df$Motif == mot & mot.freq.df$Strand == strnd,])/ mot.total
                                                     )
                                        }))
                         }))

mot.freq$Motif <- as.factor(mot.freq$Motif)
mot.freq$Strand <- as.factor(mot.freq$Strand)

#Create barplot
p.alpha.mot <- ggplot(data = mot.freq, aes(x = Motif, y = Freq, group = Strand, fill = Strand)) +
  geom_col(color = "black", position="dodge") +
  geom_hline(yintercept=0) +
  coord_cartesian(ylim = c(0,100), xlim = c(0.5,2.5), expand = FALSE) +
  ylab("% of motif-containing IS") +
  theme_classic() +
  theme(axis.line.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_text(size = 10, colour = "black"),
        axis.text.y = element_text(size = 10, colour = "black")
  )
