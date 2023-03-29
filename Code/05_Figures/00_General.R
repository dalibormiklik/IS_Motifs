#Load libraries
#Load generaly important objects and functions

#Load libraries
library(GenomicRanges)
library(ggplot2)
library(ggseqlogo)
library(ggpubr)
library(dplyr)
library(ggbeeswarm)
library(ggforce)
library(tidytext)


#Set objects with path names
#Set the home path - where the IS_Motifs directory is present
#setwd(...)
wd <- paste0(getwd(),'/IS_Motifs/')
data.wd <- paste0(wd, 'Data/')
result.wd <- paste0(wd, 'Results/')
s.wd <- paste0(data.wd, 'IS/Ran/')
scr.wd <- paste0(wd, 'Code/05_Figures/')

setwd(wd)

#Load functions
source(paste0(scr.wd, "Functions.R"))

#Load lengths of target site dulication (TSD)
retro.tsd <- c(6, 5, 4, 6, 6, 4, 0)
names(retro.tsd) <- c("HTLV", "HIV", "MLV", "MVV",
                      "ASLV", "PFV", "Random")

#Load sample info
smpl.info <- read.table(paste0(data.wd, "Sample_info.txt"),
                        header = TRUE, quote = "", sep = "\t", stringsAsFactors = FALSE)
smpl.info

nucs <- c("A", "T", "G", "C")

#General color blind-friendly palette
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
#Custom color blind-friendly colours for A,T,G,C
cb.nuc.1 <- c("#009E4B", "#0072B2", "#D55E00", "#E69F00")
cb.nuc.2 <- c("#00AD58", "#026CA7", "#CE5E02", "#EAC700")
cb.nuc.3 <- c("#00A046", "#005AB5", "#DC3220", "#E6B200")

nuc.cols <- cb.nuc.3
names(nuc.cols) <- nucs
