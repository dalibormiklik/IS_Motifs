#Get motif variants and their frequencies
#Sequences are extracted from 'virus_cm.ppm' associated sequences

#Set colors for samples
cols <- c("#D55E00")

smpl.tab <- data.frame(Virus = virus,
                       File = paste0(virus, "/", virus, "_" , cm.ppm,"_in_", gnm, "_rmsk_", repeatname),
                       Sample = cm.ppm)
smpl.names <- smpl.tab$Sample

#Calculate frequencies of the motifs
#inside lapply function: reads "fa_tab.txt" files of IS and shuffled controls
#motif.pos = which positions in target sequences are considered to create motif

motif.pos <- c(8,9,12,16,19,20)
#motif.pos <- 8:20

comp.motif.df <- do.call(rbind,
                         lapply(1:nrow(smpl.tab),
                                function(smpl.n) {
                                  #Load data for given sample
                                  seq.tab <- read.table(paste0(result.wd, "HotSpot/",
                                                              smpl.tab$File[smpl.n], "_fa_tab.txt"),
                                                        header = FALSE, quote = "", sep = "\t", stringsAsFactors = TRUE,
                                                        col.names = c("Name", "Seq"))
                                  
                                  #Extract motifs from sequences
                                  if(any(colnames(seq.tab) == "Seq")) {seq.tab <- data.frame(Seq = as.character(seq.tab[,which(colnames(seq.tab) == "Seq")]))}
                                  #if(any(colnames(seq.tab.ran) == "Seq")) {seq.tab.ran <- data.frame(Seq = as.character(seq.tab.ran[,which(colnames(seq.tab.ran) == "Seq")]))}
                                  
                                  motif.df <- MotifTab(seq.tab, motif.pos)
                                  
                                  #Calculate occurence of each motif
                                  motif.uniq.df <- data.frame(Motif = motif.df$Motif[!duplicated(motif.df$Motif)])
                                  motif.uniq.df$Count <- sapply(motif.uniq.df$Motif, function(x) {length(motif.df[motif.df$Motif == x,])})
                                  motif.uniq.df$Motif <- factor(motif.uniq.df$Motif, levels = motif.uniq.df$Motif[order(motif.uniq.df$Count, decreasing = TRUE)])
                                  
                                  comp.motif.uniq.df <- motif.uniq.df
                                  
                                  total <- sum(comp.motif.uniq.df$Count)
                                  
                                  comp.motif.uniq.df$Perc <- round(100 * comp.motif.uniq.df$Count / total, digits = 2)


                                  comp.motif.uniq.df
                                  
                                }))

#Select data from which the plot will be created
#Select top.x number of most frequent motifs
top.x <- 10
data.to.plot <-head(comp.motif.df[order(comp.motif.df$Count, decreasing = TRUE),], top.x)
data.to.plot <- data.to.plot[order(data.to.plot$Perc, decreasing = FALSE),]
#Convert Motif and Sample columns to factor 
data.to.plot$Motif <- factor(data.to.plot$Motif,
                             levels = data.to.plot$Motif[order(data.to.plot$Perc, decreasing = FALSE)])

#Plot
p.mfreq <- ggplot(data = data.to.plot, aes(x = Perc, y = Motif)) +
  scale_colour_manual(values = unname(cols)) +
  scale_fill_manual(values = alpha(unname(cols), .3)) +
  geom_col(color = cols, fill = alpha(cols, .3), show.legend = FALSE) +
  geom_vline(xintercept = 0, size = 1) +
  xlab("% component\nsequences") +
  ylab(element_blank()) +
  coord_cartesian(xlim = c(0,60)) +
  scale_y_reordered() +
  theme_classic() +
  theme(axis.line.y = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_text(colour = "black", face = "bold", size=8),
        axis.text.y = element_text(colour = "black", family = "mono", face = "bold", size=8, angle = 0),
        axis.text.x = element_text(colour = "black", face = "bold", size=8, angle = 0),
        text = element_text(colour = "black"),
        legend.title = element_blank(),
        legend.text = element_text(colour = "black", face="bold"),
  )

