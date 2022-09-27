#Calculate membership function (MFM) of the sequences to the mixture components
#Isolate sequences with MFM to selected component >= than threshold value

#Load selected sequences (Part I)
#-//- PPMs of mixture components (Part II)
#Calculate membership functions (Part III)
#Select sequences according to membership function (IV)

#Two files are created as an output:
#[...]_SeqMembFun.txt - membership function-containing table with all sequences in rows and components in columns
#"ISseq_[...].txt" - sequences with the

#select.smpl - row number of sample in smpl.info data.frame
#mix.name - name of the selected mixture
#cm - name of the selected component
#min.mf - treshold membership function for sequence selection
select.smpl <- 8
mix.name <- "M08"
cm <- "PPM07"
min.mf <- 0.9

cm.ppm <- paste0(mix.name, "_", cm)

#Set libraries nad objects with path names
library(GenomicRanges)

wd <- #'~/IS_mixtures/'
data.wd <- paste0(wd, 'Data/')
result.wd <- paste0(wd, 'Results/')
sq.wd <- paste0(data.wd, 'IS/')
scr.wd <- paste0(wd, 'Code/04_HotSpot/CompSeq/')

setwd(wd)

#Load functions
source(paste0(scr.wd, "Functions.R"))

#Load sample info
smpl.info <- read.table(paste0(data.wd, "Data/Sample_info.txt"),
                        header = TRUE, quote = "", sep = "\t", stringsAsFactors = FALSE)
smpl.info

nucs <- c("A", "T", "G", "C")


#PART I
#----
#Load sequences and create nucleotide matrix from the sequences

#Select names associated with sample
is.set.name <- smpl.info$Name[select.smpl]
virus <- smpl.info$seqName[select.smpl]
gnm <- smpl.info$Control[select.smpl]
sln <- smpl.info$seqLen[select.smpl]

#Run PARTs II-III only if table already not exists
if(file.exists(paste0(result.wd, "HotSpot/", cm.ppm, "/", virus, "_" , cm, "_SeqMembFun.txt"))) {
  
  print("MFM file exists:")
  print("Loading...")
  qseq <- read.table(paste0(result.wd, "HotSpot/", cm.ppm, "/", virus, "_" , cm, "_SeqMembFun.txt"),
                     sep = "\t", stringsAsFactor = FALSE, header = TRUE)
  print("qseq object ready.")
  
} else {
  
  print("MFM file DOES NOT exists:")
  print("Loading data for computation...")
  
  #Load IS sequences
  is.seq <- as(read.table(paste0(sq.wd, virus, "/", virus, "_is", sln,".txt"),
                          sep = "\t", stringsAsFactor = FALSE, header = FALSE),
               "DataFrame")
  
  #Create nucleotide matrix from sequences
  is.nuc.mat <- NucMat(is.seq[,1])
  colnames(is.nuc.mat) <- RelPosNames(is.nuc.mat)
  
  #Part II
  
  #Load names of the components in the selected mixture
  
  comp.mix.dir <- dir(path = paste0(result.wd, "Mixtures/", is.set.name), pattern = "[MB].*")
  source(paste0(scr.wd,'Load_mixPPM.R'))
  
  #PART III
  
  #Calculate membership functions of the sequences to mixture components
  
  print(paste0("Computing MFM for ", mix.name))
  
  #Load PPMs for mixture cm
  #Create ppm.list object
  ppm.tab <- ppm.tab[ppm.tab$PPM != "PPM00",]
  ppm.list <- lapply(unique(ppm.tab$PPM),
                     function(p) {
                       ppm.tab[ppm.tab$PPM == p, which(colnames(ppm.tab) != "PPM")]
                     })
  names(ppm.list) <- unique(ppm.tab$PPM)
  
  #Load component weights
  wm <- read.table(paste0(result.wd, "Mixtures/", is.set.name, "/",
                          mix.name, "/Wm.txt"),
                   header = FALSE, stringsAsFactors = FALSE)
  colnames(wm) <- c("cNum", "seqLen", "Wm")
  wm <- wm[wm$cNum != 0,]
  wm.m <- wm$Wm
  names(wm.m) <- paste0("PPM0", wm$cNum)
  
  #Calculate membership functions for the whole mixture
  print("Calculating membership functions")
  qseq <- MembershipFunction(is.nuc.mat, ppm.list, wm.m)
  
  #Transform qseq to data.frame (is.matrix)
  qseq <- as.data.frame(qseq)
  #Add column with sequences
  qseq$Seq <- is.seq[,1]
  
  #Save table with membership function values and sequences
  print(paste0("Saving ", cm, "_SeqMembFun.txt"))  
  write.table(qseq, paste0(result.wd, "HotSpot/", cm.ppm, "/", virus, "_" , cm, "_SeqMembFun.txt"), 
              quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
  }

#PART IV
#----
#Select sequences associated with given component
w.col <- which(colnames(qseq) == cm)
mfm.seq <- data.frame(Seq = qseq$Seq[qseq[,w.col] >= min.mf])

#Save table with membership function values and sequences
print(paste0("Saving ", cm, "_SeqMembFun.txt"))  
write.table(mfm.seq, paste0(result.wd, "HotSpot/", cm.ppm, "/", "ISseq_", virus, "_" , mix.name, "_", cm, ".txt"), 
            quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

print(paste0("ISseq_", virus, "_" , mix.name, "_", cm, ".txt Saved"))

print("...")
print("Mischief Managed.")
