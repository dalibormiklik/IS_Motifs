#Load Mixture PPMs
#mix.name= mixture name
#mix.name object needs to be set externally of this script
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

