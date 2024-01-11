#---------------------------------------#
#--- Pipeline for metabarcoding data ---#
#---------------------------------------#


# Date: 2023.05.29
# Primer pair: ill-1380F and ill-1510R
# Sequencing completed by Sequench (https://www.sequench.co.nz/)
# Run number CAW-22-11.1 (repeat/change details for other sequencing runs)


#--------------------------------#
##### LIBRARY INSTALL & LOAD #####
#--------------------------------#
#--- Library installation and setup ---#

# Installs BiocManager if don't already have it
# DO NOT INSTALL AGAIN IF YOU HAVE IT...STUFF BREAKS

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.16")

# Install phyloseq & dada2
# DO NOT INSTALL AGAIN IF YOU HAVE IT...STUFF BREAKS
BiocManager::install("phyloseq")
BiocManager::install("dada2", force = TRUE)
BiocManager::install("ShortRead")

# Loading libraries
library("phyloseq")
library("ggplot2")
library("ShortRead")
library("Biostrings")
library("dada2")

# Set working directory
### CHANGE ME ### This will be unique to your project, make sure it is correct
setwd("C:/Documents/CAW-22-11.1")


#--------------------#
##### DATA SETUP #####
#--------------------#

# If you have received your sequence data nested in multiple levels of folders
# the next step saves you a lot of copy and pasting and puts all the '.fastq' files
# into a single folder. 


#--- Get files out of sub folders ---#

###CHANGE ME###
# STEP 1: Set source destination - this needs to be the folder all the files and folders are nested in
source_folder <- "C:/Documents/CAW-22-11.1/v9"

# STEP 2: Make and Set destination folder - This is the folder you want all the '.fastq' files to end up in
seq_folder_name <- "sequence-files"
dir.create(seq_folder_name)
destination_folder <- "C:/Documents/CAW-22-11.1/sequence-files"

# STEP 3: Get list of files in source folder and subfolders:
file_list <- list.files(source_folder, recursive = TRUE, full.names = TRUE)

# STEP 4: Loop through each file - This finds and moves all the nested files to their new home
for (file in file_list) {
  # Construct the destination path by removing the source folder path
  destination_path <- file.path(destination_folder, basename(file))
  # Move the file to the destination folder
  file.rename(file, destination_path)
}
# File sorting COMPLETE


#--- Delete now empty sub folders ---#
# You don't have to do this step, but as a result of the last one you will now have a whole series of empty
# folders and sub-folder. For the humans out there that like a wee tidy, this goes out to you.

# STEP 1: Get the list of sub-directories in the source folder
# After the file moving loop, the code retrieves a list of sub-directories in the source folder
# Using list.dirs() with recursive = FALSE.
subdirs <- list.dirs(source_folder, recursive = TRUE, full.names = TRUE)

# STEP 2: Sort the sub-directories 
# This sorts remaining sub-directories in reverse order using 'sort()' 
# This is to ensure that the deepest ones are deleted first
subdirs <- sort(subdirs, decreasing = TRUE)

# STEP 3: This is a loop through each subdirectory
for (subdir in subdirs) {
  # Check if the subdirectory is empty
  if (length(list.files(subdir, recursive = TRUE)) == 0) {
    # Delete the empty subdirectory
    unlink(subdir, recursive = TRUE, force = TRUE)
  }
}


# Folder deleting COMPLETE, all neat and tidy



#--- Check files ---#
### CHANGE ME ### This will direct R to your sequence data in your new 'destination_folder', make sure it is correct
path <- "C:/Documents/CAW-22-11.1/sequence-files"
list.files(path)

# List names of forward and reverse reads
fnFs <- sort(list.files(path, pattern = "R1_001.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(path, pattern = "R2_001.fastq.gz", full.names = TRUE))

# View the list of names
fnFs
fnRs



#--- Primer pair setup and checking ---#

# This code specifies the primer pair used for your sequencing run. 
### CHANGE ME ### Make sure these are the primers you used. 
# CHANGE & REPEAT WHOLE SECTION FOR V4 region sequences.

FWD <- "CCCTGCCHTTTGTACACAC" #FORWARD: 1380F v9 
REV <- "CCTTCYGCAGGTTCACCTAC" #REVERSE: 1510R v9


#Here we find the forward and reverse compliments of primers in all sequences
allOrients <- function(primer) {
  # Then create all orientations of the input sequence
  require(Biostrings)
  dna <- DNAString(primer)  # The Biostrings works w/ DNAString objects rather than character vectors
  orients <- c(Forward = dna, Complement = complement(dna), Reverse = reverse(dna), 
               RevComp = reverseComplement(dna))
  return(sapply(orients, toString))  # Convert back to character vector
}

FWD.orients <- allOrients(FWD)
REV.orients <- allOrients(REV)
FWD.orients


# This says how many times the primers are found, and in what orientation
primerHits <- function(primer, fn) {
  # Counts number of reads in which the primer is found
  nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
  return(sum(nhits > 0))
}

rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs[[1]]), 
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs[[1]]))

# At this point you are hoping for all reads, or the vast majority, to be in the forward direction.
# If this is not the case, find Johnny P and grill him about why...it's most likely a you problem
# But he will know where to direct you to solve it.


#--------------------------------------#
##### PROCESSING STEP 1 - Cutadapt #####
#--------------------------------------#
#--- finds and sets up cutadapt --- #

### CHANGE ME ### 
# Please change the following file path to the cutadapt path on your machine
cutadapt <- "C:/Miniconda3/Scripts/cutadapt" 

# Run shell commands from R
system2(cutadapt, args = "--version") 


# This code determines the new file path for the cutadapt folder/files you are creating
path.cut <- file.path(path, "cutadapt18S")
if(!dir.exists(path.cut)) dir.create(path.cut) #does it exist, if it doesn't then make it
fnFs.cut <- file.path(path.cut, basename(fnFs))
fnRs.cut <- file.path(path.cut, basename(fnRs))



#--- Trimming primer sequences off all reads ---#

# Trim FWD off of R1 (forward reads) - 
R1.flags <- paste0("-g", " ^", FWD) 

# Trim REV off of R2 (reverse reads)
R2.flags <- paste0("-G", " ^", REV) 


#--- Run Cutadapt ---#
# This step can take some time

for(i in seq_along(fnFs)) {
  system2(cutadapt, args = c("-e 0.08 --discard-untrimmed", R1.flags, R2.flags,
                             "-o", fnFs.cut[i], "-p", fnRs.cut[i], # output files
                             fnFs[i], fnRs[i])) # input files
}

# You will now have a new 'cutadapt18s' folder within the folder your sequences are in.
# These have had the primers trimmed off each sequence.
# Next you are making sure the primer sequences have been removed and that you still have the same amount of files.

# Re-Check how many files still contain the primers after running cutadapt: 
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.cut[[1]]), 
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.cut[[1]]))


# Identify the forward and reverse fastq filenames have the format:
cutFs <- sort(list.files(path.cut, pattern = "R1_001.fastq.gz", full.names = TRUE))
cutRs <- sort(list.files(path.cut, pattern = "R2_001.fastq.gz", full.names = TRUE))

# Check that the number of forward and reverse files are still the same:
if(length(cutRs) == length(cutRs)) print("Forward and reverse files match. Go forth and explore")
if (length(cutRs) != length(cutRs)) stop("Forward and reverse files do not match. Better go back and have a check")

# Extract sample names, assuming filenames have format:
get.sample.name <- function(fname) strsplit(basename(fname), "_")[[1]][1]
sample.names <- unname(sapply(cutFs, get.sample.name))
head(sample.names)

#Check sample names
sample.names



#----------------------#
##### STOP POINT 1 #####
#----------------------#
# This stop point is designed to save the work-space environment so you cab pick up from here later if needed.
# It means you do not have to run the above code again and also just worth saving sometimes in case stuff goes awry


# generate folder in working directory for saved environments
folder_name <- "workspace_environments"
dir.create(folder_name)

# generate dynamic file name to track iterations
current_datetime <- format(Sys.time(), "%Y-%m-%d_%H-%M")
file_name_stop1 <- paste0("C:/Documents/CAW-22-11.1/workspace_environments/",current_datetime,"_CAW-22-11.1_stop1", ".RData")

# save workspace environment
save.image(file = file_name_stop1)



#-------------------------------------------------#
##### PROCESSING STEP 2 - Dada2 quality check #####
#-------------------------------------------------#


#--- LOAD WORKSPACE ENVIRONMENT ---#
### CHANGE ME ###
# Make sure you update this file link as you save out your stop point files.
#load("C:/Documents/CAW-22-11.1/workspace_environments/2023-05-29_14-38_CAW-22-11.1_stop1.RData")



#--- Dada2 - Quality Checking ---#

require(dada2)
require(ggplot2)

if(length(cutFs) <= 20) {
  fwd_qual_plots <- plotQualityProfile(cutFs) +
    scale_x_continuous(breaks=seq(0,250,10)) +
    scale_y_continuous(breaks=seq(0,40,2)) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
 
   rev_qual_plots <- plotQualityProfile(cutRs) +
    scale_x_continuous(breaks=seq(0,250,10)) +
    scale_y_continuous(breaks=seq(0,40,2)) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
} else {
  rand_samples <- sample(size = 20, 1:length(cutFs)) # grab 20 random
fwd_qual_plots <- plotQualityProfile(cutFs[rand_samples]) +
  scale_x_continuous(breaks=seq(0,250,10)) +
  scale_y_continuous(breaks=seq(0,40,2)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

rev_qual_plots <- plotQualityProfile(cutRs[rand_samples]) +
  scale_x_continuous(breaks=seq(0,250,10)) +
  scale_y_continuous(breaks=seq(0,40,2)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
}


# Shows the plotted quality profiles
fwd_qual_plots
rev_qual_plots  

# Good idea to export plots as images for later, you never know when you might need to whip them out.


#----------------------#
##### STOP POINT 2 #####
#----------------------#

# This stop point is designed to save the work-space environment so you cab pick up from here later if needed.
# It means you do not have to run the above code again and also just worth saving sometimes in case stuff goes awry


#generate dynamic file name to track iterations
current_datetime <- format(Sys.time(), "%Y-%m-%d_%H-%M") # This tags the date and time into the generated file

#This creates a file name with the date and time tag for this stop point
file_name_stop2 <- paste0("C:/Documents/CAW-22-11.1/workspace_environments/",current_datetime,"_CAW-22-11.1_stop2", ".RData")

#saves work space environment
save.image(file = file_name_stop2)


#--------------------------------------------------#
##### PROCESSING STEP 3 - Filtering & trimming #####
#--------------------------------------------------#

#--- LOAD WORKSPACE ENVIRONMENT ---#
### CHANGE ME ###
# Make sure you update this file link as you save out your stop point files.

#load("C:/Documents/CAW-22-11.1/workspace_environments/2023-05-29_15-28_CAW-22-11.1_stop2")



#--- Filtering and trimming ---#

filtpathF <- file.path(path.cut, "filtered", basename(cutFs))
filtpathR <- file.path(path.cut, "filtered", basename(cutRs))

out <- filterAndTrim(cutFs, filtpathF, cutRs, filtpathR,
                     ### CHANGE ME ### 'truncLen = ' This is the expected length of your sequences/what you are trimming them to.
                     # v9 = 150 max length
                     truncLen=c(130,130), maxEE=c(2,4), truncQ=2, maxN=0, rm.phix=TRUE,
                     compress=TRUE, verbose=TRUE, multithread=TRUE)

out

sample.names <- sapply(strsplit(basename(filtpathF), "_"), `[`, 1) # Assumes filename = samplename_XXX.fastq.gz
sample.namesR <- sapply(strsplit(basename(filtpathR), "_"), `[`, 1) # Assumes filename = samplename_XXX.fastq.gz
if(identical(sample.names, sample.namesR)) {print("Files are still matching.....congratulations")
} else {stop("Forward and reverse files do not match.")}
names(filtpathF) <- sample.names
names(filtpathR) <- sample.namesR


set.seed(100) # set seed to ensure that randomized steps are replicable


#--- Learn forward error rates ---#
# This can take time, go get a coffee, get out for a run or something
errF <- learnErrors(filtpathF, nbases=1e8, multithread=TRUE, verbose = TRUE)
## 100285500 total bases in 668570 reads from 24 samples will be used for learning the error rates.


#--- Learn reverse error rates ---#
# This will also take time, so maybe crack into some writing since you're all caffeinated or exercised.
errR <- learnErrors(filtpathR, nbases=1e8, multithread=TRUE, verbose = TRUE)

#--- Plotting forward error rates ---#
#Forward
errF_plot <- plotErrors(errF, nominalQ=TRUE)
errF_plot

#--- Plotting reverse error rates ---#
#Reverse
errR_plot <- plotErrors(errR, nominalQ=TRUE)
errR_plot


##### STOP POINT 3 #####
#--- Here you can save the work-space environment and pick up from here later if needed ---#
#--- This means you do not have to run the above code again
#--- Also just worth saving sometimes incase computer crashes


#generate dynamic file name to track iterations
current_datetime <- format(Sys.time(), "%Y-%m-%d_%H-%M")
file_name_stop3 <- paste0("C:/Documents/CAW-22-11.1/workspace_environments/",current_datetime,"_CAW-22-11.1_stop3", ".RData")

#save workspace environment
save.image(file = file_name_stop3)



#----------------------------------------------#
##### PROCESSING STEP 4 - Chimera removal ######
#----------------------------------------------#

#--- LOAD WORKSPACE ENVIRONMENT ---#
### CHANGE ME ###
# Make sure you update this file link as you save out your stop point files.

#load("C:/Documents/CAW-22-11.1/workspace_environments/2023-05-29_17-51_CAW-22-11.1_stop3.RData")



#----- Chimera removal -----#
# This code performs the de-replication of sequences in a FASTQ file specified by 'filtpathF' using the 'derepFastq()' function. 
# The resulting de-replicated data is stored in the variable 'derepF' or 'derepR'. 
# The 'verbose=TRUE' argument is used to display detailed progress information during the de-replication process.
derepF <- derepFastq(filtpathF, verbose=TRUE)
derepR <- derepFastq(filtpathR, verbose=TRUE)


#--- Checkpoint 1 ---#
#--- This is just to ensure saving of those longer processes, as sometimes the computer crashes after this point.
#generate dynamic file name to track iterations
current_datetime <- format(Sys.time(), "%Y-%m-%d_%H-%M")
file_name_stop3_checkpoint <- paste0("C:/Documents/CAW-22-11.1/workspace_environments/",current_datetime,"_CAW-22-11.1_stop3-checkpoint", ".RData")

#save workspace environment
save.image(file = file_name_stop3_checkpoint)

#--------------------#


# This code does the de-noising and chimera removal of amplicon sequence data using the 'dada()' function from the dada2 package. 
# The de-replicated data ('derepF' or 'derepR'), an error model ('errF' or 'errR'), and other parameters are provided as input to the function. 
# The de-noised sequences and related information are stored in the 'dadaF.pseudo' or 'dadaR.pseudo' object.
dadaF.pseudo <- dada(derepF, err=errF, multithread=TRUE, pool="pseudo")
# THIS BIT TAKES A LOOOOOONG TIME ^


dadaR.pseudo <- dada(derepR, err=errR, multithread=TRUE, pool="pseudo")
# THIS BIT TAKES A LOOOOOONG TIME ^


# Merging of paired-end amplicon sequence reads using the 'mergePairs()' function from the dada2 package. 
# The de-noised and de-replicated forward and reverse reads, along with specified parameters, are provided as input. 
# The merged sequences and related information are stored in the mergers object.
mergers <- mergePairs(dadaF.pseudo, derepF, dadaR.pseudo, derepR, maxMismatch = 1, minOverlap = 10, verbose=TRUE)
# THIS TAKES A WHILE TOO  ^


# line of code creates a sequence table (seqtab) from the merged amplicon sequences stored in the mergers object. 
# The sequence table provides information about the abundance of each unique sequence in the merged data, 
# allowing further analysis and downstream processing.
seqtab <- makeSequenceTable(mergers)

# Splits the name of a file or directory specified by path using the "-" character as the delimiter. 
# It then extracts the first two elements from each split part and stores them in the split.dir variable. 
split.dir <- sapply(strsplit(basename(path), "-"), `[`, 1:2)

# Combines the first and second elements from the `split.dir` variable into a single string, 
# separated by a period, and stores it in the `split.dir.name` variable.
split.dir.name  <- paste(split.dir[1], split.dir[2], sep=".")

# saves the R object `seqtab` as an RDS file, using the constructed file path and name that
# includes the directory path, base name, and a combination of split parts from `split.dir.name`.
saveRDS(seqtab, paste0(path, "seqtab", split.dir.name,".rds"))

# generates a frequency table that shows the distribution of sequence lengths in the `seqtab` object. 
# It counts the number of sequences with each length and provides a summary of their occurrence
table(nchar(getSequences(seqtab)))

# generates a histogram that visualizes the distribution of sequence lengths in the `seqtab` object. 
# It plots the frequency of different sequence lengths on the x-axis and the count or density on the y-axis, 
# providing an overview of the distribution pattern
hist(nchar(getSequences(seqtab)), main="Distribution of sequence lengths")

# filters the columns of the `seqtab` object based on the character length of their column names.
# It selects and retains only the columns with column names whose length is within your designated range.
### CHANGE ME ### the numbers in `seq(000,000)` to represent the expected character range for your sequences. 
# The resulting subset is assigned back to the `seqtab` object, effectively updating it with the filtered columns.
seqtab <- seqtab[,nchar(colnames(seqtab)) %in% seq(120,142)]

#  applies the `removeBimeraDenovo()` function to the `seqtab` object to remove chimeric sequences. 
# It utilizes multiple threads for faster computation and provides detailed output during the process. 
# The resulting sequence table without chimeras is assigned to the `seqtab.nochim` variable.
seqtab.nochim <- removeBimeraDenovo(seqtab, multithread=TRUE, verbose=TRUE)

### CHANGE ME ###
# make sure I match your file path to your working folder for this project
# You are saving a copy of the chimera free sequence table to your working directory
saveRDS(seqtab.nochim, "C:/Documents/CAW-22-11.1/CAW-22-11.1_j_v9_seqtab.nochim.rds")


### At this point, if you are not doing your own taxonomy files, send the .rds file to John.

##### STOP POINT 4 #####
#--- Here you can save the work-space environment and pick up from here later if needed ---#
#--- This means you do not have to run the above code again
#--- Also just worth saving sometimes incase computer crashes


#generate dynamic file name to track iterations
current_datetime <- format(Sys.time(), "%Y-%m-%d_%H-%M")
file_name_stop4 <- paste0("C:/Documents/CAW-22-11.1/workspace_environments/",current_datetime,"_CAW-22-11.1_stop4", ".RData")

#save workspace environment
save.image(file = file_name_stop4)

# At this point, repeat the above steps with all sequence runs. Change parameters for V4 such as primers etc.
# Only once they are all complete move to the next step!!!

#----------------------------------#
#--- Seqence table merge script ---#
#----------------------------------#

# This is to be run when you have multiple sequencing runs that were split to run through the your pipeline
# Only merge seqtab files of the same gene region
# Next step after this script  has been run is taxonomy

#----------------#
##### Set up #####
#----------------#

#--- Set working directory ---#
setwd("C:/Documents/R")

#--- Load libraries ---#
library("dada2")
#library("phyloseq")
#library("ggplot2")
#library("ShortRead")
#library("Biostrings")

#--- Read in seqtab files ---#
# Here I am defining the name of each data table `seqtab.nochim.file.name` 
# The reading in the existing sequence table for each sequence run using `readRDS()`
# make sure to change the file path to match where your files are!

# v9 region sequence tables
seqtab.nochim.CAW.22.11.1 <- readRDS("C:/Documents/R/CAW-22-11.1_j_v9_seqtab.nochim.rds")
seqtab.nochim.CAW.22.11.2 <- readRDS("C:/Documents/R/CAW-22-11.2_j_v9_seqtab.nochim.rds")
seqtab.nochim.CAW.22.11.3 <- readRDS("C:/Documents/R/CAW-22-11.3_j_v9_seqtab.nochim.rds")

seqtab.nochim.CAW.23.21.1 <- readRDS("C:/Documents/R/CAW-23-21_j_v9_seqtab.nochim.rds")

# v4 region sequence tables
seqtab.nochim.CAW.22.12.1 <- readRDS("C:/Documents/R/CAW-22-12.1_j_v4_seqtab.nochim.rds")
seqtab.nochim.CAW.22.12.2 <- readRDS("C:/Documents/R/CAW-22-12.2_j_v4_seqtab.nochim.rds")
seqtab.nochim.CAW.22.12.3 <- readRDS("C:/Documents/R/CAW-22-12.3_j_v4_seqtab.nochim.rds")

seqtab.nochim.CAW.23.22.1 <- readRDS("C:/Documents/R/CAW-23-22.1_j_v4_seqtab.nochim.rds")
seqtab.nochim.CAW.23.22.2 <- readRDS("C:/Documents/R/CAW-23-22.2_j_v4_seqtab.nochim.rds")

#----------------------------#
##### Merge seqtab files #####
#----------------------------#

# define input table groups
v9_tables <- list(seqtab.nochim.CAW.22.11.1, seqtab.nochim.CAW.22.11.2, seqtab.nochim.CAW.22.11.3, seqtab.nochim.CAW.23.21.1)
v4_tables <- list(seqtab.nochim.CAW.22.12.1, seqtab.nochim.CAW.22.12.2, seqtab.nochim.CAW.22.12.3, seqtab.nochim.CAW.23.22.1, seqtab.nochim.CAW.23.22.2)

# Merge all v9 region sequence tables
v9.seqtab.merge <- mergeSequenceTables(tables = v9_tables)
#check new file rownames order
rownames(v9.seqtab.merge)


# Merge all v4 region sequence tables
v4.seqtab.merge <- mergeSequenceTables(tables = v4_tables)
#check new file rownames order
rownames(v4.seqtab.merge)


# Save `.rds` file of each merged table to working directory
saveRDS(v9.seqtab.merge, "C:/Documents/R/2023-06-01_j_v9_merged-seqtab.nochim.rds")
saveRDS(v4.seqtab.merge, "C:/Documents/R/2023-06-01_j_v4_merged-seqtab.nochim.rds")


# Now you should have a new merged sequence table in your working directory, or the filepath you defined above
# Woohoo, all done.

#-------------------------------------#
#--- Metabarcoding Taxonomy script ---#
#-------------------------------------#
# Here we combine the taxonomy and sample metadata files and
# remove negative control sample data from the data set. Then you will have a delightful clean data-set
# You need to have a metadata file at this point that will have all the base information associated with your sequenced samples. 

#----------------#
##### Set up #####
#----------------#

#--- Load library ---#
#install.packages("dplyr")
library("dplyr")
library("phyloseq")
library("dada2")


#--- Read in data ---#
### CHANGE ME ### 
# Make sure file path and name is correct!!!

#--- v9 ---#
# read in `seqtab.nochim` file from your working directory
v9.seqtab.merge <- readRDS("C:/Documents/R/2023-06-01_j_v9_merged-seqtab.nochim.rds")
#read in metadata .csv file
v9.metadata <- read.csv("C:/Documents/R/JSPhd_sample-metadata.csv", header = TRUE, row.names = 1)
# read in taxonomy .rds file
v9.taxonomy <- readRDS("C:/Documents/R/v9.taxonomy.v5.rds")

#--- v4 ---#
# read in `seqtab.nochim` file from your working directory
v4.seqtab.merge <- readRDS("C:/Documents/R/2023-06-01_j_v4_merged-seqtab.nochim.rds")
#read in metadata .csv file
v4.metadata <- read.csv("C:/Documents/R/JSPhd_sample-metadata.csv", header = TRUE, row.names = 1)
# read in taxonomy .rds file
v4.taxonomy <- readRDS("C:/Documents/R/v4.taxonomy.v5.rds")

#--- Rename row names ---#
#rename `seqtab.nochim` rows - NOTE there is a way to just get row-names from the metadata file, commented out below.
# make sure these match the sample order in your sequence table file.


#--- v9 ---#
#re-check sample rowname order
rownames(v9.seqtab.merge)
#define new row names here
v9.new.rownames <- c(
  #CAW.22.11
  "N1a1", "N1a2", "N1a3", "N1b1", "N1b2", "N1b3", "N1c1", "N1c2", "N1c3",
  "N1d1", "N1d2", "N1d3", "N1e1", "N1e2", "N1e3", "N1f1", "N1f2", "N1f3",
  "N1g1", "N1g2", "N1g3", "N1h1", "N1h2", "N1h3", "N1i1", "N1i2", "N1i3",
  "N2a1", "N2a2", "N2a3", "N2b1", "N2b2", "N2b3", "N2c1", "N2c2", "N2c3",
  "N2d1", "N2d2", "N2d3", "N2e1", "N2e2", "N2e3", "N2f1", "N2f2", "N2f3",
  "N2g1", "N2g2", "N2g3", "N2h1", "N2h2", "N2h3", "N2i1", "N2i2", "N2i3",
  "N3a1", "N3a2", "N3a3", "N3b1", "N3b2", "N3b3", "N3c1", "N3c2", "N3c3",
  "N3d1", "N3d2", "N3d3", "N3e1", "N3e2", "N3e3", "N3f1", "N3f2", "N3f3",
  "N3g1", "N3g2", "N3g3", "N3h1", "N3h2", "N3h3", "N3i1", "N3i2", "N3i3",
  "N4a1", "N4a2", "N4a3", "N4b1", "N4b2", "N4b3", "N4c1", "N4c2", "N4c3",
  "N4d1", "N4d2", "N4d3", "ExtBL1", "PCRBL1", "H20BL1" ,
  "N4e1", "N4e2", "N4e3", "N4f1", "N4f2", "N4f3", "N4g1", "N4g2", "N4g3", 
  "N4h1", "N4h2", "N4h3", "N4i1", "N4i2", "N4i3", "N5a1", "N5a2", "N5a3", 
  "N5b1", "N5b2", "N5b3", "N5c1", "N5c2", "N5c3", "N5d1", "N5d2", "N5d3", 
  "N5e1", "N5e2", "N5e3", "N5f1", "N5f2", "N5f3", "N5g1", "N5g2", "N5g3", 
  "N5h1", "N5h2", "N5h3", "N5i1", "N5i2", "N5i3", "N6a1", "N6a2", "N6a3", 
  "N6b1", "N6b2", "N6b3", "N6c1", "N6c2", "N6c3", "N6d1", "N6d2", "N6d3",
  "N6e1", "N6e2", "N6e3", "N6f1", "N6f2", "N6f3", "N6g1", "N6g2", "N6g3",
  "N6h1", "N6h2", "N6h3", "N6i1", "N6i2", "N6i3", "H01", "H03", "H25",
  "H26", "H32", "H37", "H50", "H53", "H84", "H95", "H96", "H101",
  "H102", "H105", "H108", "H115", "H116", "H119", "H120", "H121", "H124",
  "H159", "H165", "H171", "ExtBL2", "PCRBL2", "H20BL2",
  "H172", "H178", "H216", "H252", "21_AP1", "21_AP2", "21_AP3", "21_AP8", "21_AP9",
  "21_AP10", "21_AP11", "H167", "ExtBL3", "PCRBL3", "H20BL3",
  #CAW.23.21
  "P2a-00", "P2a-25", "P2a-50", "P2a-75", "P2a-100", "P2a-125", "P2a-150", "P2a-175", "P2a-200",
  "P2a-225", "P2b-00", "P2b-25", "P2b-50", "P2b-75", "P2b-100", "P2b-125", "P2b-150", "P2b-175",
  "P2c-00", "P2a-CT", "P2a-CM", "P2a-CB", "P2b-CT", "P2b-CM", "P2b-CB", "P2c-CT", "P2c-CM",
  "P2c-CB", "P2a1", "P2a2", "P2a3", "P2b1", "P2b2", "P2b3", "P1a-00", "P1a-25",
  "P1a-50", "P1a-75", "P1a-100", "P1b-00", "P1b-25", "P1b-50", "P1b-75", "P1b-100", "P1a-CT",
  "P1a-CM", "P1a-CB", "P1b-CT", "P1b-CM", "P1b-CB", "P1c-CT", "P1c-CM", "P1c-CB", "P4a1",
  "P4a2", "P4a3", "P4a4", "P3a-CT", "P3a-CM", "P3a-CB", "P3b-CT", "P3b-CM", "P3b-CB",
  "P3c-CT", "P3c-CM", "P3c-CB", "P1a1", "P1a2", "P1a3", "P1b1", "P1b2", "P1b3",
  "P1c1", "P1c2", "P1c3", "P1d1", "P1d2", "P1d3", "P1e1", "P1e2", "P1e3", 
  "C1a1", "C1a2", "C1a3", "C1b1", "C1b2", "C1b3", "C1c1", "C1c2", "C1c3", 
  "EXT_BL","PCR_BL1","H20_BL1",
  "T1a1", "T1a2", "T1a3", "T1b1", "T1b2", "T1b3", "T1c1", "T1c2", "T1c3",
  "T1d1", "T1d2", "T1d3", "T1e1", "T1e2", "T1e3", "T1f1", "T1f2", "T1f3",
  "T1g1", "T1g2", "T1g3", "T1h1", "T1h2", "T1h3", "T1i1", "T1i2", "T1i3",
  "T2a1", "T2a2", "T2a3", "T2b1", "T2b2", "T2b3", "T2c1", "T2c2", "T2c3",
  "T2d1", "T2d2", "T2d3", "T2e1", "T2e2", "T2e3", "T2f1", "T2f2", "T2f3",
  "T2g1", "T2g2", "T2g3", "T2h1", "T2h2", "T2h3", "T2i1", "T2i2", "T2i3",
  "T3a1", "T3a2", "T3a3", "T3b1", "T3b2", "T3b3", "T3c1", "T3c2", "T3c3",
  "T3d1", "T3d2", "T3d3", "T3e1", "T3e2", "T3e3", "T3f1", "T3f2", "T3f3",
  "T3g1", "T3g2", "T3g3", "T3h1", "T3h2", "T3h3", "T3i1", "T3i2", "T3i3",
  "EXT_BL6","PCR_BL3", "H20_BL3"
)

#change rownames here
rownames(v9.seqtab.merge) <- v9.new.rownames


#--- v4 ---#
#re-check sample rowname order
rownames(v4.seqtab.merge)
#define new row names here
v4.new.rownames <- c(
  #CAW.22.12
  "N1a1", "N1a2", "N1a3", "N1b1", "N1b2", "N1b3", "N1c1", "N1c2", "N1c3",
  "N1d1", "N1d2", "N1d3", "N1e1", "N1e2", "N1e3", "N1f1", "N1f2", "N1f3",
  "N1g1", "N1g2", "N1g3", "N1h1", "N1h2", "N1h3", "N1i1", "N1i2", "N1i3",
  "N2a1", "N2a2", "N2a3", "N2b1", "N2b2", "N2b3", "N2c1", "N2c2", "N2c3",
  "N2d1", "N2d2", "N2d3", "N2e1", "N2e2", "N2e3", "N2f1", "N2f2", "N2f3",
  "N2g1", "N2g2", "N2g3", "N2h1", "N2h2", "N2h3", "N2i1", "N2i2", "N2i3",
  "N3a1", "N3a2", "N3a3", "N3b1", "N3b2", "N3b3", "N3c1", "N3c2", "N3c3",
  "N3d1", "N3d2", "N3d3", "N3e1", "N3e2", "N3e3", "N3f1", "N3f2", "N3f3",
  "N3g1", "N3g2", "N3g3", "N3h1", "N3h2", "N3h3", "N3i1", "N3i2", "N3i3",
  "N4a1", "N4a2", "N4a3", "N4b1", "N4b2", "N4b3", "N4c1", "N4c2", "N4c3",
  "N4d1", "N4d2", "N4d3", "v4_ExtBL1", "v4_PCRBL1", "v4_H20BL1" ,
  "N4e1", "N4e2", "N4e3", "N4f1", "N4f2", "N4f3", "N4g1", "N4g2", "N4g3", 
  "N4h1", "N4h2", "N4h3", "N4i1", "N4i2", "N4i3", "N5a1", "N5a2", "N5a3", 
  "N5b1", "N5b2", "N5b3", "N5c1", "N5c2", "N5c3", "N5d1", "N5d2", "N5d3", 
  "N5e1", "N5e2", "N5e3", "N5f1", "N5f2", "N5f3", "N5g1", "N5g2", "N5g3", 
  "N5h1", "N5h2", "N5h3", "N5i1", "N5i2", "N5i3", "N6a1", "N6a2", "N6a3", 
  "N6b1", "N6b2", "N6b3", "N6c1", "N6c2", "N6c3", "N6d1", "N6d2", "N6d3",
  "N6e1", "N6e2", "N6e3", "N6f1", "N6f2", "N6f3", "N6g1", "N6g2", "N6g3",
  "N6h1", "N6h2", "N6h3", "N6i1", "N6i2", "N6i3", "H01", "H03", "H25",
  "H26", "H32", "H37", "H50", "H53", "H84", "H95", "H96", "H101",
  "H102", "H105", "H108", "H115", "H116", "H119", "H120", "H121", "H124",
  "H159", "H165", "H171", "ExtBL2", "PCRBL2", "H20BL2",
  "H172", "H178", "H216", "H252", "21_AP1", "21_AP2", "21_AP3", "21_AP8", "21_AP9",
  "21_AP10", "21_AP11", "H167", "ExtBL3", "PCRBL3", "H20BL3",
  #CAW.23.22
  "P2a-00", "P2a-25", "P2a-50", "P2a-75", "P2a-100", "P2a-125", "P2a-150", "P2a-175", "P2a-200",
  "P2a-225", "P2b-00", "P2b-25", "P2b-50", "P2b-75", "P2b-100", "P2b-125", "P2b-150", "P2b-175",
  "P2c-00", "P2a-CT", "P2a-CM", "P2a-CB", "P2b-CT", "P2b-CM", "P2b-CB", "P2c-CT", "P2c-CM",
  "P2c-CB", "P2a1", "P2a2", "P2a3", "P2b1", "P2b2", "P2b3", "P1a-00", "P1a-25",
  "P1a-50", "P1a-75", "P1a-100", "P1b-00", "P1b-25", "P1b-50", "P1b-75", "P1b-100", "P1a-CT",
  "P1a-CM", "P1a-CB", "P1b-CT", "P1b-CM", "P1b-CB", "P1c-CT", "P1c-CM", "P1c-CB", "P4a1",
  "P4a2", "P4a3", "P4a4", "P3a-CT", "P3a-CM", "P3a-CB", "P3b-CT", "P3b-CM", "P3b-CB",
  "P3c-CT", "P3c-CM", "P3c-CB", "P1a1", "P1a2", "P1a3", "P1b1", "P1b2", "P1b3",
  "P1c1", "P1c2", "P1c3", "P1d1", "P1d2", "P1d3", "P1e1", "P1e2", "P1e3", 
  "C1a1", "C1a2", "C1a3", "C1b1", "C1b2", "C1b3", "C1c1", "C1c2", "C1c3", 
  "EXT_BL","PCR_BL1","H20_BL1",
  "T1c1", "T1c2", "T1c3",
  "T1d1", "T1d2", "T1d3", "T1e1", "T1e2", "T1e3", "T1f1", "T1f2", "T1f3",
  "T1g1", "T1g2", "T1g3", "T1h1", "T1h2", "T1h3", "T1i1", "T1i2", "T1i3",
  "T2a1", "T2a2", "T2a3", "T2b1", "T2b2", "T2b3", "T2c1", "T2c2", "T2c3",
  "T2d1", "T2d2", "T2d3", "T2e1", "T2e2", "T2e3", "T2f1", "T2f2", "T2f3",
  "T2g1", "T2g2", "T2g3", "T2h1", "T2h2", "T2h3", "T2i1", "T2i2", "T2i3",
  "T3a1", "T3a2", "T3a3", "T3b1", "T3b2", "T3b3", "T3c1", "T3c2", "T3c3",
  "T3d1", "T3d2", "T3d3", "T3e1", "T3e2", "T3e3", "T3f1", "T3f2", "T3f3",
  "T3g1", "T3g2", "T3g3", "T3h1", "T3h2", "T3h3", "T3i1", "T3i2", "T3i3",
  "EXT_BL6","PCR_BL3", "H20_BL3", "T1a1", "T1a2", "T1a3", "T1b1", "T1b2", "T1b3"
)
#change rownames here
rownames(v4.seqtab.merge) <- v4.new.rownames



#--- Merging data ---#

#combine metadata, sample information and metabarcoding files
v9 <- phyloseq(otu_table(v9.seqtab.merge, taxa_are_rows = FALSE), sample_data(v9.metadata), tax_table(v9.taxonomy))
v4 <- phyloseq(otu_table(v4.seqtab.nochim, taxa_are_rows = FALSE), sample_data(v4.metadata), tax_table(v4.taxonomy))



#---------------------------------------#
##### Subtracting negative controls #####
#---------------------------------------#

# This code segment subtracts sequences found in the negative control samples
# from all samples to account for contamination.


#--- V9 region ---#
# subsets neg control samples in merged data file
v9.negatives = subset_samples(v9, type =="control")
v9.negatives

# calculates the column sums for each variable in the `v9.negatives` OTU table and stores them in the `v9.negSums` vector.
v9.negSums <- colSums(otu_table(v9.negatives))

# Converts `v9.negSums` to a vector
v9.negSums_vec <- as.vector(v9.negSums)

# converts the object `otu_table(v9.negatives)` into a matrix and assigns it to the variable `v9.table`
v9.table = as(otu_table(v9.negatives), "matrix")

# converts the matrix object `v9.table` into a data frame and assigns it to the variable `v9.tabledf`
v9.tabledf = as.data.frame(v9.table)
v9.tabledf[,1:length(v9.tabledf)] <- sweep(v9.tabledf[,1:length(v9.tabledf)], 2,v9.negSums_vec)

# changes neg controls to zero
v9.tabledf <- replace(v9.tabledf, v9.tabledf < 0,0)

# removes zero samples (or negative controls)
v9.tabledf_noneg <- v9.tabledf[rowSums(v9.tabledf)!=0,]

# Creating a new phyloseq object with negative controls removed from all samples
# Then combining it with the appropriate metadata and taxonomy information
v9.cleaned <- phyloseq(otu_table(v9.tabledf_noneg, taxa_are_rows = FALSE),
                       sample_data(v9.metadata),
                       tax_table(v9.taxonomy))

v9.cleaned


#--------------------------------------------------#
##### Comparing original and clean sample sums #####
#--------------------------------------------------#

# calculates the total count sums for each sample in two different phyloseq objects (v9 and v9.cleaned) and 
# stores them in separate data frames. It then performs a left join operation to merge the sample sums from 
# both objects into a single data frame (control.rem.sums) for comparison.


org.ss <- as.data.frame(sample_sums(v9)) # calculate total count sums for each sample
org.ss$Names <- rownames(org.ss) # adds new column to data frame called names containing the sample names

new.ss <- as.data.frame(sample_sums(v9.cleaned)) # calculate total count sums for each sample
new.ss$Names <- rownames(new.ss) # adds new column to data frame called names containing the sample names

control.rem.sums <- dplyr::left_join(org.ss, new.ss, by="Names") # merging data frames based on matching `Name` values

colnames(control.rem.sums) <- c("Original", "Names", "New") # defining column names

control.rem.sums <- control.rem.sums %>%
  mutate(perc = New/Original*100) # finding the percent change between original and cleaned

control.rem.sums



#--- Save files ---#
#Saving clean data set to working directory
saveRDS(v9.cleaned,"C:/Documents/R/v9_cleaned.rds")

########################
#### ANALYSIS START ####
########################

library(grid)
library(tidyr)
library(tidyverse)
library(ggplot2)
library(phyloseq)

setwd("C:/Documents/R")

#--- Region efficiency ---#
# "Gene region comparison: V4 and V9 Small subunit rDNA in characterising Eukaryotic Microalgal communities"

v9 <- readRDS("C:/Documents/R/v9.cleaned2.rds")
v4 <- readRDS("C:/Documents/R/v4.cleaned2.rds")


v9.sums <- data.frame(sample_sums(v9))
write.csv(v9.sums, "v9sums.csv")

v4.sums <- data.frame(sample_sums(v4))
write.csv(v4.sums, "v4sums.csv")


# line of code filters the taxa in the `cleanv9` object and assigns the filtered result back to the `cleanv9` variable.
# It removes any taxa that have a sum of values less than or equal to 0, effectively removing taxa with no or negligible abundance in the dataset.
v9 <- filter_taxa(v9, function(x) sum(x) > 0, TRUE)
v4 <- filter_taxa(v4, function(x) sum(x) > 0, TRUE)

v9.taxa.df <- data.frame(tax_table(v9))
v4.taxa.df <- data.frame(tax_table(v4))


##### FUNCTION: Microalgae rare curves #####

# SUbset to EMC classes
taxa_filter <- c("Dinophyceae", # Dino
                 "Bacillariophyceae", "Coscinodiscophyceae", "Mediophyceae", # Dia
                 
                 "Bolidophyceae", "Syndiniophyceae", "Chrysophyceae", 
                 "Dictyochophyceae", "Eustigmatophyceae", "Raphidophyceae",
                 
                 "Picocystophyceae", "Pleurastrophyceae", "Prasinophyceae", # Chlorophyta
                 "Chlorodendrophyceae", "Chlorophyceae", "Pedinophyceae", "Trebouxiophyceae", "Ulvophyceae", # Chlorophyte - phylum; Chlorophytina - subphylum
                 "Mamiellophyceae", "Nephroselmidophyceae", "Pyramimonadophyceae", # Chlorophyte - phylum; Prasinophytina - subphylum
                 
                 "Pinguiophyceae", "Chloropicophyceae", 
                 "Prymnesiophyceae","Pavlovophyceae", "Rappephyceae", # Hapto
                 "Euglenida", "Xanthophyceae"
)

#make smaller supsets
dino <- c("Dinophyceae")
dia <- c("Bacillariophyceae", "Mediophyceae", "Coscinodiscophyceae")
dia.p <- c("Bacillariophyceae")
dia.rc <- c("Coscinodiscophyceae")
dia.pc <- c("Mediophyceae")
hap <- c("Pavlovophyceae", "Rappephyceae", "Prymnesiophyceae")
chl <- c("Chlorodendrophyceae", "Chlorophyceae", "Chloropicophyceae", "Mamiellophyceae", # Green algae
         "Nephroselmidophyceae","Pedinophyceae", "Trebouxiophyceae", "Pyramimonadophyceae")
oth <-c("Raphidophyceae", "Eustigmatophyceae","Chrysophyceae", 
        "Dictyochophyceae", "Pinguiophyceae", "Bolidophyceae", "Euglenida", "Xanthophyceae"
)

#--- TEMPERATE SITE ---#
n1.rare_comparison <- function(data, gene) {
  subset.1 <- subset_samples(data, siteID == "N1")
  subset.2 <- subset_taxa(subset.1, R5 %in% taxa_filter)
  merge.1 <- merge_samples(subset.2, "siteID")
  rare <- ggrare(merge.1, step = 500) +
    theme_minimal() +
    ylim(0, 450)
  
  data.2 <- ggplot_build(rare)$data[[1]]
  data.2 <- data.2 %>% mutate(region = gene) 
  
}

# --- TROPICAL SITE --- #
t3.rare_comparison <- function(data, gene) {
  subset.1 <- subset_samples(data, siteID == "T3")
  subset.2 <- subset_taxa(subset.1, R5 %in% taxa_filter)
  merge.1 <- merge_samples(subset.2, "siteID")
  rare <- ggrare(merge.1, step = 500) +
    theme_minimal() +
    ylim(0, 450)
  
  data.2 <- ggplot_build(rare)$data[[1]]
  data.2 <- data.2 %>% mutate(region = gene) 
  
}

# --- Combine plots --- #
combine_data_plot <- function(data1, data2, site) {
  combine <- rbind(data1, data2)
  plot <- ggplot(combine, aes(x = x, y = y, color = region)) +
    geom_line(aes(linetype = region), size = 1) +
    xlab("Sequence sample size") +
    ylab("ASV Abundance") +
    scale_color_manual(values = c("#a9d8d5", "#4b8cca"))+
    theme_minimal()+
    ylim(0, 2500)
  
}

# --- Run Functions ---#

N1.v9.rare <- n1.rare_comparison(v9, "v9")
N1.v4.rare <- n1.rare_comparison(v4, "v4")
T3.v9.rare <- t3.rare_comparison(v9, "v9")
T3.v4.rare <- t3.rare_comparison(v4, "v4")

N1.plot <- combine_data_plot(N1.v9.rare, N1.v4.rare, "N1")
T3.plot <- combine_data_plot(T3.v9.rare, T3.v4.rare, "T3")

combined_plots <- grid.arrange(N1.plot, T3.plot, ncol = 2)
combined_plots

##### FUNCTION: Eukaryote rare curves #####
# --- Temperate --- #
n1.euk_comparison <- function(data, gene) {
  subset.1 <- subset_samples(data, siteID == "N1")
  merge.1 <- merge_samples(subset.1, "siteID")
  rare <- ggrare(merge.1, step = 500) + # Plot
    theme_minimal() #+
  #ylim(0, 450)
  
  data.2 <- ggplot_build(rare)$data[[1]]
  data.2 <- data.2 %>% mutate(region = gene) 
  
}

# --- Tropical --- #
t3.euk_comparison <- function(data, gene) {
  subset.1 <- subset_samples(data, siteID == "T3")
  merge.1 <- merge_samples(subset.1, "siteID")
  rare <- ggrare(merge.1, step = 500) + # Plot
    theme_minimal() #+
  #ylim(0, 450)
  
  data.2 <- ggplot_build(rare)$data[[1]]
  data.2 <- data.2 %>% mutate(region = gene) 
  
}

# --- Run function --- #
N1.v9.euk <- n1.euk_comparison(v9, "v9")
N1.v4.euk <- n1.euk_comparison(v4, "v4")
t3.v9.euk <- t3.euk_comparison(v9, "v9")
t3.v4.euk <- t3.euk_comparison(v4, "v4")

N1.plot <- combine_data_plot(N1.v9.euk, N1.v4.euk, "N1")
T3.plot <- combine_data_plot(t3.v9.euk, t3.v4.euk, "T3")

euk.combined_plots <- grid.arrange(N1.plot, T3.plot, ncol = 2)
euk.combined_plots


##### SUMMARY STATISTICS #####
# --- Functions: create dataframe with summary stats --- #
std.error <- function(x) sd(x)/sqrt(length(x))

calculate_summary <- function(data, region, community) {
  total_ASVs <- sum(sample_sums(data))
  min_val <- min(sample_sums(data))
  max_val <- max(sample_sums(data))
  std_error <- std.error(sample_sums(data))
  mean_val <- mean(sample_sums(data))
  
  df <- data.frame(Site = unique(sample_data(data)$siteID),
                   Region = region,
                   Total_ASVs = total_ASVs,
                   Minimum = min_val,
                   Maximum = max_val,
                   Mean = mean_val,
                   Std_Error = std_error,
                   Community = community
  )
  
  return(df)
}

# --- Sub set data --- #
# EMC community #
v9.N1.emc <- subset_data <- v9 %>%
  subset_taxa(R5 %in% taxa_filter) %>%
  subset_samples(siteID == "N1")

v9.T3.emc <- subset_data <- v9 %>%
  subset_taxa(R5 %in% taxa_filter) %>%
  subset_samples(siteID == "T3")

v4.N1.emc <- subset_data <- v4 %>%
  subset_taxa(R5 %in% taxa_filter) %>%
  subset_samples(siteID == "N1")

v4.T3.emc <- subset_data <- v4 %>%
  subset_taxa(R5 %in% taxa_filter) %>%
  subset_samples(siteID == "T3")

# EUKARYOTIC COMMUNITY #
v9.N1 <- subset_samples(v9, siteID == "N1")
v9.T3 <- subset_samples(v9, siteID == "T3")
v4.N1 <- subset_samples(v4, siteID == "N1")
v4.T3 <- subset_samples(v4, siteID == "T3")

# DINOFLAGELLATES #
v9.N1.dino <- subset_data <- v9 %>%
  subset_taxa(R5 %in% dino) %>%
  subset_samples(siteID == "N1")

v4.N1.dino <- subset_data <- v4 %>%
  subset_taxa(R5 %in% dino) %>%
  subset_samples(siteID == "N1")

v9.t3.dino <- subset_data <- v9 %>%
  subset_taxa(R5 %in% dino) %>%
  subset_samples(siteID == "T3")

v4.t3.dino <- subset_data <- v4 %>%
  subset_taxa(R5 %in% dino) %>%
  subset_samples(siteID == "T3")

# DIATOMS #
v9.n1.dia <- subset_data <- v9 %>%
  subset_taxa(R5 %in% dia) %>%
  subset_samples(siteID == "N1")
v9.n1.diap <- subset_data <- v9 %>%
  subset_taxa(R5 %in% dia.p) %>%
  subset_samples(siteID == "N1")
v9.n1.diarc <- subset_data <- v9 %>%
  subset_taxa(R5 %in% dia.rc) %>%
  subset_samples(siteID == "N1")
v9.n1.diapc <- subset_data <- v9 %>%
  subset_taxa(R5 %in% dia.pc) %>%
  subset_samples(siteID == "N1")

v4.n1.dia <- subset_data <- v4 %>%
  subset_taxa(R5 %in% dia) %>%
  subset_samples(siteID == "N1")
v4.n1.diap <- subset_data <- v4 %>%
  subset_taxa(R5 %in% dia.p) %>%
  subset_samples(siteID == "N1")
v4.n1.diarc <- subset_data <- v4 %>%
  subset_taxa(R5 %in% dia.rc) %>%
  subset_samples(siteID == "N1")
v4.n1.diapc <- subset_data <- v4 %>%
  subset_taxa(R5 %in% dia.pc) %>%
  subset_samples(siteID == "N1")


v9.t3.dia <- subset_data <- v9 %>%
  subset_taxa(R5 %in% dia) %>%
  subset_samples(siteID == "T3")
v9.t3.diap <- subset_data <- v9 %>%
  subset_taxa(R5 %in% dia.p) %>%
  subset_samples(siteID == "T3")
v9.t3.diarc <- subset_data <- v9 %>%
  subset_taxa(R5 %in% dia.rc) %>%
  subset_samples(siteID == "T3")
v9.t3.diapc <- subset_data <- v9 %>%
  subset_taxa(R5 %in% dia.pc) %>%
  subset_samples(siteID == "T3")


v4.t3.dia <- subset_data <- v4 %>%
  subset_taxa(R5 %in% dia) %>%
  subset_samples(siteID == "T3")
v4.t3.diap <- subset_data <- v4 %>%
  subset_taxa(R5 %in% dia.p) %>%
  subset_samples(siteID == "T3")
v4.t3.diarc <- subset_data <- v4 %>%
  subset_taxa(R5 %in% dia.rc) %>%
  subset_samples(siteID == "T3")
v4.t3.diapc <- subset_data <- v4 %>%
  subset_taxa(R5 %in% dia.pc) %>%
  subset_samples(siteID == "T3")


# Chlorophtes #
v9.n1.chl <- subset_data <- v9 %>%
  subset_taxa(R5 %in% chl) %>%
  subset_samples(siteID == "N1")

v4.n1.chl <- subset_data <- v4 %>%
  subset_taxa(R5 %in% chl) %>%
  subset_samples(siteID == "N1")

v9.t3.chl <- subset_data <- v9 %>%
  subset_taxa(R5 %in% chl) %>%
  subset_samples(siteID == "T3")

v4.t3.chl <- subset_data <- v4 %>%
  subset_taxa(R5 %in% chl) %>%
  subset_samples(siteID == "T3")

# Hatoophtes #
v9.n1.hap <- subset_data <- v9 %>%
  subset_taxa(R5 %in% hap) %>%
  subset_samples(siteID == "N1")

v4.n1.hap <- subset_data <- v4 %>%
  subset_taxa(R5 %in% hap) %>%
  subset_samples(siteID == "N1")

v9.t3.hap <- subset_data <- v9 %>%
  subset_taxa(R5 %in% hap) %>%
  subset_samples(siteID == "T3")

v4.t3.hap <- subset_data <- v4 %>%
  subset_taxa(R5 %in% hap) %>%
  subset_samples(siteID == "T3")

# OTHER #
v9.n1.oth <- subset_data <- v9 %>%
  subset_taxa(R5 %in% oth) %>%
  subset_samples(siteID == "N1")

v4.n1.oth <- subset_data <- v4 %>%
  subset_taxa(R5 %in% oth) %>%
  subset_samples(siteID == "N1")

v9.t3.oth <- subset_data <- v9 %>%
  subset_taxa(R5 %in% oth) %>%
  subset_samples(siteID == "T3")

v4.t3.oth <- subset_data <- v4 %>%
  subset_taxa(R5 %in% oth) %>%
  subset_samples(siteID == "T3")


# Calculate summary statistics for each community
v4.N1.emc_sum <- calculate_summary(v4.N1.emc, "v4", "Emc")
v4.N1.euk_sum <- calculate_summary(v4.N1, "v4", "Euk")
v4.n1.dino_sum <- calculate_summary(v4.N1.dino, "v4", "Dino")
v4.n1.dia_sum <- calculate_summary(v4.n1.dia, "v4", "Dia")
v4.n1.diap_sum <- calculate_summary(v4.n1.diap, "v4", "Dia pennate")
v4.n1.diarc_sum <- calculate_summary(v4.n1.diarc, "v4", "Dia radial")
v4.n1.diapc_sum <- calculate_summary(v4.n1.diapc, "v4", "Dia polar")
v4.n1.chl_sum <- calculate_summary(v4.n1.chl, "v4", "chl")
v4.n1.hap_sum <- calculate_summary(v4.n1.hap, "v4", "hap")
v4.n1.oth_sum <- calculate_summary(v4.n1.oth, "v4", "oth")

v9.N1.emc_sum <- calculate_summary(v9.N1.emc, "v9", "Emc")
v9.N1.euk_sum <- calculate_summary(v9.N1, "v9", "Euk")
v9.n1.dino_sum <- calculate_summary(v9.N1.dino, "v9", "Dino")
v9.n1.dia_sum <- calculate_summary(v9.n1.dia, "v9", "Dia")
v9.n1.diap_sum <- calculate_summary(v9.n1.diap, "v9", "Dia pennate")
v9.n1.diarc_sum <- calculate_summary(v9.n1.diarc, "v9", "Dia radial")
v9.n1.diapc_sum <- calculate_summary(v9.n1.diapc, "v9", "Dia polar")
v9.n1.chl_sum <- calculate_summary(v9.n1.chl, "v9", "chl")
v9.n1.hap_sum <- calculate_summary(v9.n1.hap, "v9", "hap")
v9.n1.oth_sum <- calculate_summary(v9.n1.oth, "v9", "oth")

v4.T3.emc_sum <- calculate_summary(v4.T3.emc, "v4", "Emc")
v4.T3.euk_sum <- calculate_summary(v4.T3, "v4", "Euk")
v4.t3.dino_sum <- calculate_summary(v4.t3.dino, "v4", "Dino")
v4.t3.dia_sum <- calculate_summary(v4.t3.dia, "v4", "Dia")
v4.t3.diap_sum <- calculate_summary(v4.t3.diap, "v4", "Dia pennate")
v4.t3.diarc_sum <- calculate_summary(v4.t3.diarc, "v4", "Dia radial")
v4.t3.diapc_sum <- calculate_summary(v4.t3.diapc, "v4", "Dia polar")
v4.t3.chl_sum <- calculate_summary(v4.t3.chl, "v4", "chl")
v4.t3.hap_sum <- calculate_summary(v4.t3.hap, "v4", "hap")
v4.t3.oth_sum <- calculate_summary(v4.t3.oth, "v4", "oth")

v9.T3.emc_sum <- calculate_summary(v9.T3.emc, "v9", "Emc")
v9.T3.euk_sum <- calculate_summary(v9.T3, "v9", "Euk")
v9.t3.dino_sum <- calculate_summary(v9.t3.dino, "v9", "Dino")
v9.t3.dia_sum <- calculate_summary(v9.t3.dia, "v9", "Dia")
v9.t3.diap_sum <- calculate_summary(v9.t3.diap, "v9", "Dia pennate")
v9.t3.diarc_sum <- calculate_summary(v9.t3.diarc, "v9", "Dia radial")
v9.t3.diapc_sum <- calculate_summary(v9.t3.diapc, "v9", "Dia polar")
v9.t3.chl_sum <- calculate_summary(v9.t3.chl, "v9", "chl")
v9.t3.hap_sum <- calculate_summary(v9.t3.hap, "v9", "hap")
v9.t3.oth_sum <- calculate_summary(v9.t3.oth, "v9", "oth")

# Combine all summary data frames
summary_data <- rbind(v4.N1.euk_sum, v4.N1.emc_sum, v4.n1.dino_sum, v4.n1.diap_sum, v4.n1.diarc_sum, v4.n1.diapc_sum, v4.n1.chl_sum, v4.n1.hap_sum, v4.n1.oth_sum,
                      v9.N1.euk_sum, v9.N1.emc_sum, v9.n1.dino_sum, v9.n1.diap_sum, v9.n1.diarc_sum, v9.n1.diapc_sum, v9.n1.chl_sum, v9.n1.hap_sum, v9.n1.oth_sum,
                      v4.T3.euk_sum, v4.T3.emc_sum,  v4.t3.dino_sum, v4.t3.diap_sum, v4.t3.diarc_sum, v4.t3.diapc_sum, v4.t3.chl_sum, v4.t3.hap_sum, v4.t3.oth_sum,
                      v9.T3.euk_sum, v9.T3.emc_sum, v9.t3.dino_sum, v9.t3.diap_sum, v9.t3.diarc_sum, v9.t3.diapc_sum, v9.t3.chl_sum,v9.t3.hap_sum, v9.t3.oth_sum
)

write.csv(summary_data, file = "C:/Documents/R/23-07-12_chp1_sample_sums.csv", row.names = TRUE)

# You now have a CSV file in your working directory with all the summary statistics


##### Alpha diversity #####

# --- making data frame with all the alpha diversity info and plotting it
v4.subset <- subset_samples(v4, siteID %in% c("N1", "T3"))
v9.subset <- subset_samples(v9, siteID %in% c("N1", "T3"))

v9.emc <- subset_data <- v9 %>%
  subset_taxa(R5 %in% taxa_filter) %>%
  subset_samples(siteID %in% c("N1", "T3"))

v4.emc <- subset_data <- v4 %>%
  subset_taxa(R5 %in% taxa_filter) %>%
  subset_samples(siteID %in% c("N1", "T3"))

sample_data(v9.emc)
sample_sums(v9.emc)

# trim & rarefy
trim_rare_rich <- function(data, depth, ssu) {
  sample_data(data)$site.point <- paste(sample_data(data)$siteID, sample_data(data)$point, sep = "-") # make new column combining site & point info
  sample_data(data)$siteID <- factor(sample_data(data)$siteID, levels = c("N1", "T3")) # Make siteID factor
  merge <- merge_samples(data, "site.point") # merge by site.point
  trim <- prune_samples(sample_sums(merge)>= depth, merge) # trim to pre-defined depth, based off rare curve & sequencing depth
  rare <- rarefy_even_depth(trim, min(sample_sums(trim)), replace = F) # rarefy
  
  print(sample_sums(rare))
  
  rich <- estimate_richness(rare) # Get all alpha diversity measures
  rich$sampleID <- rownames(rich)
  sd <- data.frame(sample_data(rare))
  sd$sampleID <- rownames(sd)
  sd$sampleID <- gsub("-", "\\.", sd$sampleID)
  
  rich <- left_join(rich, sd, by = "sampleID")
  rich <- rich %>% mutate(region = ssu)
  
  print(rich)
}

#install.packages("gridExtra")
library(gridExtra)

# Within region comparison
v4.rare <- trim_rare_rich(v4.emc, 11000, "v4") # rarefy for v4/v4 comparison
v9.rare.comp <- trim_rare_rich(v9.emc, 11000, "v9") # rarefy for v4/v9 comparison

write.csv(v4.rare, file = "C:/Documents/R/23-08-09_chp1_v4-AlphaIndicies.csv", row.names = TRUE)
write.csv(v9.rare.comp, file = "C:/Documents/R/23-08-09_chp1_v9-AlphaIndicies.csv", row.names = TRUE)

rich.rare <- rbind(v4.rare, v9.rare.comp) # Combined data frame for v4/v9 comparison - all rarefied to 11000
region_colors <- c(v4 = "#a9d8d5", v9 = "#4b8cca")

# plot choa, shannon, invs simp, goods
alpha_plots <- function(data1) {
  obs <- ggplot() +
    geom_boxplot(data = data1, aes(x = interaction(region, siteID, sep = "-"), y = Observed, fill = region)) +
    theme_minimal() + scale_fill_manual(values = region_colors) +
    guides(fill = "none") + xlab("Site") + theme(legend.position = "none")
  chao <- ggplot() + geom_boxplot(data = data1, aes(x = interaction(region, siteID), y = Chao1, fill = region)) +
    theme_minimal() + scale_fill_manual(values = region_colors) +
    guides(fill = "none") + xlab("Site") + theme(legend.position = "none")
  sha <- ggplot() + geom_boxplot(data = data1, aes(x = interaction(region, siteID), y = Shannon, fill = region))+
    theme_minimal() + scale_fill_manual(values = region_colors) +
    guides(fill = "none") + xlab("Site") + theme(legend.position = "none")
  isim <- ggplot() + geom_boxplot(data = data1, aes(x = interaction(region, siteID), y = InvSimpson, fill = region)) +
    theme_minimal() + scale_fill_manual(values = region_colors) +
    guides(fill = "none") + xlab("Site") + theme(legend.position = "none")
  
  combined_plots <- grid.arrange(obs, chao, sha, isim, ncol = 2)
  combined_plots
}

site.labs <- c("N1", "T3")
labs <- c("v4:N1", "v9:N1", "v4:T3", "v9:T3")

alpha_plots <- function(data1,lab) {
  region_colors <- c(v4 = "#a9d8d5", v9 = "#4b8cca")
  measures <- c("Chao1", "Shannon", "InvSimpson")
  
  plots <- list()
  
  for (measure in measures) {
    plot <- ggplot(data1) +
      geom_boxplot(aes(x = interaction(region, siteID, sep = "-"), y = !!sym(measure), fill = region)) +
      theme_minimal() +
      scale_fill_manual(values = region_colors) +
      guides(fill = "none") +
      xlab("Site") +
      theme(legend.position = "none") +
      scale_x_discrete(labels = lab)
    
    plots[[measure]] <- plot
  }
  
  combined_plots <- grid.arrange(grobs = plots, ncol = 2)
  combined_plots
}

comp.plots <- alpha_plots(rich.rare, labs)

### --- Significance test - Wilcox --- ###

alpha <- read.csv("C:/Documents/R/23-10-10_chp1_AlphaIndicies.csv")
#rich.rare = estimate_richness(iceSamples_c3_trim.rare)
alpha.temp <- alpha[alpha$siteID == 1, ]
alpha.trop <- alpha[alpha$siteID == 2, ]

# CHao1 -
pairwise.wilcox.test(alpha.temp$Chao1, sample_data(alpha.temp)$region, p.adjust.method = "BH", paired = FALSE, exact = FALSE)
pairwise.wilcox.test(alpha.trop$Chao1, sample_data(alpha.trop)$region, p.adjust.method = "BH", paired = FALSE, exact = FALSE)

# Shannon -
pairwise.wilcox.test(alpha.temp$Shannon, sample_data(alpha.temp)$region, p.adjust.method = "BH", paired = FALSE, exact = FALSE)
pairwise.wilcox.test(alpha.trop$Shannon, sample_data(alpha.trop)$region, p.adjust.method = "BH", paired = FALSE, exact = FALSE)

# inv.simps -
pairwise.wilcox.test(alpha.temp$InvSimpson, sample_data(alpha.temp)$region, p.adjust.method = "BH", paired = FALSE, exact = FALSE)
pairwise.wilcox.test(alpha.trop$InvSimpson, sample_data(alpha.trop)$region, p.adjust.method = "BH", paired = FALSE, exact = FALSE)

#----------------------------------------------------------------------------------------------------#

##### Composition Plots #####

taxonomy_breaks <- c("Raphidophyceae", # brown algae
                     "Eustigmatophyceae", # Yellow-green
                     "Chrysophyceae", "Dictyochophyceae", # Golden algae
                     "Pinguiophyceae", #other algae
                     "Bolidophyceae", 
                     "Euglenida", "Xanthophyceae",
                     
                     "Pavlovophyceae", "Rappephyceae", "Prymnesiophyceae", # Haptophyte
                     
                     "Chlorodendrophyceae", "Chlorophyceae", "Chloropicophyceae", "Mamiellophyceae", # Green algae
                     "Nephroselmidophyceae","Pedinophyceae", "Trebouxiophyceae", "Pyramimonadophyceae",# Green algae
                     
                     "Dinophyceae", # dinoflagellates
                     
                     "Bacillariophyceae", "Coscinodiscophyceae", "Mediophyceae") # diatoms 


class_colours <- c(
  
  Eustigmatophyceae = "#574D90", Pinguiophyceae = "#574D90", Chrysophyceae = "#574D90",  Dictyochophyceae = "#574D90",Raphidophyceae =  "#574D90", 
  Euglenida = "#574D90", Xanthophyceae = "#574D90",
  
  Pavlovophyceae = "#9DB0CF", Rappephyceae = "#9DB0CF", Prymnesiophyceae = "#9DB0CF",
  
  Chlorodendrophyceae = "#CFE4E1", Chlorophyceae = "#CFE4E1", Chloropicophyceae = "#CFE4E1", Mamiellophyceae  =  "#CFE4E1", 
  Nephroselmidophyceae = "#CFE4E1", Pedinophyceae = "#CFE4E1",Trebouxiophyceae = "#CFE4E1",  Pyramimonadophyceae = "#CFE4E1",
  
  Dinophyceae  = "#6FB5CF",
  
  Bacillariophyceae = "#f26957", Coscinodiscophyceae = "#fcdeb4", Mediophyceae = "#fbb768" )

site_order <- c("N1", "T3")  # Replace with your desired site order


# Creates new object `phytoplankton_class` that combines the taxa in `phytoplankton` based on their Class level division.
plot_prep <- function(data, taxlvl, ssu) {
  mrg <- merge_samples(data, "siteID") 
  sample_data(mrg)$siteID <- c("N1", "T3")
  tax = tax_glom(mrg, taxlvl) # Groups by designated tax lvl
  perc = transform_sample_counts(tax, function(OTU) OTU/sum(OTU)*100) # calculates the percentages of each class
  # can change threshold - (> 1) currently at 1 % can go to zero if less abundance
  filtr = filter_taxa(perc, function(x) mean(x) > .0, TRUE)
  filtr_df <- psmelt(filtr)
  
  # Convert the Sample column to a factor with the desired order
  filtr_df$Sample <- factor(filtr_df$Sample, levels = site_order)
  filtr_df$R5 <- factor(filtr_df$R5, levels = taxonomy_breaks)
  
  filtr_df <- filtr_df %>%
    mutate(region = ssu)
} 

sample_data (v4.emc)


v9_df <- plot_prep(v9.emc, "R5", "v9")
v4_df <- plot_prep(v4.emc, "R5", "v4")

v9_df <- v9_df %>% mutate(site.point = paste(siteID, point, sep = "."))
v4_df <- v4_df %>% mutate(site.point = paste(siteID, point, sep = "."))


proportions_df <- rbind(v4_df, v9_df)
proportions_df



# Plot
composition <- ggplot(proportions_df, aes(x = Abundance, y = interaction(region, siteID, sep = "_"), fill = R5)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(name = "Taxonomy", values = class_colours) + 
  theme_minimal() +
  ylab("Site ID")
composition


####################
#### DENDROGRAM ####

# Load libraries #

library(tidyr)
library(tidyverse)
library(ggplot2)
#devtools::install_github("joey711/phyloseq")
library("phyloseq")

setwd("C:/Documents/R") # set Working directory

# Load data
v4 <- readRDS("C:/Documents/R/v4.cleaned2.rds")
v9 <- readRDS("C:/Documents/R/v9.cleaned2.rds")

# Select site to look at (change as desired)
EMC.site <- subset_samples(v9, siteID == "T3") # CHANGE SITE AS REQUIRED
sample_sums(EMC.site)

# Clean up subset
EMC.site <- filter_taxa(EMC.site, function(x) sum(x) > 0, TRUE)



#####---- DATA WRANGLING ----#####
#SUbset to Phytoplankton classes
# filter to R5
EMC.site <- subset_taxa(EMC.site, R5 %in% taxa_filter)
EMC.site <- filter_taxa(EMC.site, function(x) sum(x) > 0, TRUE)
sample_sums(EMC.site)

# Tally ASVs to species level
taxa <- data.frame(tax_table(EMC.site))
tally.dino <- taxa %>% group_by(R4, R5, R6, R7, R8, R9) %>% tally()

taxa.df <- data.frame(tax_table(EMC.site))
# removes all `NA's` in taxonomy.
tax.clean <- data.frame(tax_table(EMC.site))

for (i in 1:9){ tax.clean[,i] <- as.character(tax.clean[,i])}
tax.clean[is.na(tax.clean)] <- ""

for (i in 1:nrow(tax.clean)){
  if (tax.clean[i,2] == ""){
    kingdom <- paste("Unclassified_Domain_", tax.clean[i,1], sep = "") #R1
    tax.clean[i, 2:9] <- kingdom
  } else if (tax.clean[i,3] == ""){
    phylum <- paste("Unclassified_Supergroup_", tax.clean[i,2], sep = "") #R2
    tax.clean[i, 3:9] <- supergroup
  } else if (tax.clean[i,4] == ""){
    phylum <- paste("Unclassified_Division_", tax.clean[i,3], sep = "") #R3
    tax.clean[i, 4:9] <- division
  } else if (tax.clean[i,5] == ""){
    phylum <- paste("Unclassified_SubDivision_", tax.clean[i,4], sep = "") #R4
    tax.clean[i, 5:9] <- subdivision
  } else if (tax.clean[i,6] == ""){
    class <- paste("Unclassified_Class_", tax.clean[i,5], sep = "") #R5
    tax.clean[i, 6:9] <- class
  } else if (tax.clean[i,7] == ""){
    order <- paste("Unclassified_Order_", tax.clean[i,6], sep = "") #R6
    tax.clean[i, 7:9] <- order
  } else if (tax.clean[i,8] == ""){
    family <- paste("Unclassified_Family_", tax.clean[i,7], sep = "") #R7
    tax.clean[i, 8:9] <- family
  } else if (tax.clean[i,9] == ""){
    tax.clean$R9[i] <- paste("Unclassified_Genus_",tax.clean$R8[i], sep = "_") #R8
  } 
}


# Define new column names
taxo.divisions <- c("Domain", "Supergroup", "Division", "Sub-division", "Class", "Order", "Family", "Genus", "Species")
colnames(tax.clean) <- taxo.divisions

EMC.site.trim <- prune_samples(sample_sums(EMC.site)>= 10000, EMC.site) # trim to pre-defined depth, based off rare curve & sequencing depth
EMC.site.rare <- rarefy_even_depth(EMC.site.trim, min(sample_sums(EMC.site.trim)), replace = F)

sample_sums(EMC.site.rare)

EMC.site.rare.m <- merge_samples(EMC.site, "siteID")
EMC.site.rare.m.phyla <- tax_glom(EMC.site.rare.m, "R9")
EMC.site.rare.m.phyla_perc = transform_sample_counts(EMC.site.rare.m.phyla, function(OTU) OTU/sum(OTU)*100)
EMC.site.rare.m.phyla_perc.log <- transform_sample_counts(EMC.site.rare.m.phyla_perc, function(x) log(x+1))
phyla.asvtab <- as.data.frame(as.matrix(t(otu_table(EMC.site.rare.m.phyla_perc.log))))


#install.packages("ape")
#library(ape)
#remove.packages(phyloseq)
#detach(package:tidytree, unload = TRUE)

# Define the dendrogram as a character string
# Taxonomic information: To write Dendro list

dias_taxa_df <- tax.clean[tax.clean$Class %in% c("Bacillariophyceae"), ] # Change to look at each group
dias_counts <- table(dias_taxa_df$Family)
fam_names <- names(dias_counts)
fam_names
print(v9.dias_counts)

v9.df <- data.frame(v9.dias_counts)

dia_taxa_df <- tax.clean[tax.clean$`Sub-division` %in% c("Chlorophyta_X"), ]
Order_counts <- table(dia_taxa_df$Class, dia_taxa_df$Family)
fam_counts <- table(dia_taxa_df$Family)
fam_name <- names(fam_counts)
fam_name
Order_counts_with_totals <- addmargins(Order_counts)
print(Order_counts_with_totals)


#####---- v4 N1 UPDATED Dendro Setup----#####
#dendro_text <- "(((Chlorodendrales), (Mamiellales), (Pyramimonadales), (Chlorellales, Microthamniales, Watanabea-Clade)),
#((Dinophyceae_X, Dinophysiales, Gonyaulacales, Gymnodiniales, Peridiniales, Suessiales, Torodiniales, Unclass.Dinophyceae)),
#((Unclass.Bacillariophyceae),
# (Chaetocerotales, Thalassiosirales, Unclass.Mediophyceae)),
#((Coccolithales, Isochrysidales, Phaeocystales, Prymnesiales, Prymnesiophyceae_X, Unclass.Prymnesiophyceae)));"

dendro_text <- "(((Chlorodendrales), (Mamiellales), (Pyramimonadales), (Chlorellales, Microthamniales, Watanabea-Clade)),
((Dinophyceae_X, Dinophysiales, Gonyaulacales, Gymnodiniales, Peridiniales, Suessiales, Torodiniales)),
  ((Chaetocerotales, Thalassiosirales)),
((Coccolithales, Isochrysidales, Phaeocystales, Prymnesiales, Prymnesiophyceae_X)));"

tip_cols <- c( Chlorodendrales = "#2f9181", Mamiellales = "#2f9181", Pyramimonadales = "#2f9181", `Watanabea-Clade` = "#2f9181",
               Chlorellales = "#2f9181", Microthamniales = "#2f9181",
               
               Dinophyceae_X = "#517caf", Dinophysiales = "#517caf", Gonyaulacales = "#517caf", Gymnodiniales = "#517caf", 
               Peridiniales = "#517caf", Suessiales = "#517caf", Torodiniales = "#517caf", 
               #Unclass.Dinophyceae = "#517caf", 
               
               #Unclass.Bacillariophyceae = "#f26957",
               
               Chaetocerotales = "#ed8a54", Thalassiosirales = "#ed8a54", #Unclass.Mediophyceae = "#ed8a54",
               
               Coccolithales = "#7b90af", Isochrysidales = "#7b90af", Phaeocystales = "#7b90af", 
               Prymnesiales = "#7b90af", Prymnesiophyceae_X = "#7b90af" 
               #, Unclass.Prymnesiophyceae = "#7b90af"
)


#####---- v4 T3 Dendro setup ----#####
dendro_text <- "(((Chlorodendrales), (Trebouxiophyceae_X)),
((Dinophyceae_X, Dinophysiales, Gonyaulacales, Gymnodiniales, Peridiniales, Prorocentrales, Suessiales, Torodiniales, Unclass.Dinophyceae)),
((Rhabdonematales, Rhopalodiales)),
((Prymnesiales, Prymnesiophyceae_X, Unclass.Prymnesiophyceae)));"


tip_cols <- c(Chlorodendrales = "#2f9181", Trebouxiophyceae_X = "#2f9181",
              
              Dinophyceae_X = "#517caf", Dinophysiales = "#517caf", Gonyaulacales = "#517caf", Gymnodiniales = "#517caf", 
              Peridiniales = "#517caf", Prorocentrales = "#517caf", Suessiales = "#517caf", Torodiniales = "#517caf", 
              #Unclass.Dinophyceae = "#517caf",
              
              Rhabdonematales = "#f26957", Rhopalodiales = "#f26957",
              
              Prymnesiales = "#7b90af", Prymnesiophyceae_X = "#7b90af"
              #, Unclass.Prymnesiophyceae = "#7b90af" 
              
)

#####---- v9 N1 Dendro setupt ----#####

#v9 N1 dendro text
dendro_text <- "(((Chlorodendrales), (Chlamydomonadales, Sphaeropleales), (Chloropicales), (Dolichomastigales, Mamiellales, Unclass.Mamiellophyceae), 
(Nephroselmidales), (Pseudoscourfieldiales, Pyramimonadales), (Microthamniales, Watanabea-Clade)), 
((Dinophyceae_X, Dinophysiales, Gonyaulacales, Gymnodiniales, Peridiniales, Prorocentrales, Suessiales, Torodiniales, Unclass.Dinophyceae)), 
((Bacillariales, Fragilariales, Naviculales, Rhaphoneidales, Rhopalodiales, Surirellales, Unclass.Bacillariophyceae),
  (Corethrales, Coscinodiscales, Paraliales, Rhizosoleniales, Unclass.Coscinodiscophyceae),
  (Anaulales, Chaetocerotales, Cymatosirales, Eupodiscales, Hemiaulales, Lithodesmiales, Thalassiosirales, Unclass.Mediophyceae ),
  (Chrysophyceae_Clade-EC2H, Ochromonadales, Paraphysomonadales, Unclass.Chrysophyceae), (Dictyochophyceae_X),(Raphidophyceae_X)),
((Isochrysidales, Phaeocystales, Prymnesiales, Prymnesiophyceae_X, Unclass.Prymnesiophyceae)));"

#v9 N1
tip_cols <- c( Chlorodendrales = "#2f9181", Chlamydomonadales = "#2f9181", Sphaeropleales = "#2f9181", Chloropicales = "#2f9181", 
               Dolichomastigales = "#2f9181", Mamiellales = "#2f9181", Unclass.Mamiellophyceae = "#2f9181", Nephroselmidales = "#2f9181", 
               Pseudoscourfieldiales = "#2f9181", Pyramimonadales = "#2f9181", Microthamniales = "#2f9181", `Watanabea-Clade` = "#2f9181",
               
               Dinophyceae_X = "#517caf", Dinophysiales = "#517caf", Gonyaulacales = "#517caf", Gymnodiniales = "#517caf", 
               Peridiniales = "#517caf", Prorocentrales = "#517caf", Suessiales = "#517caf", Torodiniales = "#517caf", 
               Unclass.Dinophyceae = "#517caf",
               
               Bacillariales = "#f26957", Fragilariales = "#f26957", Naviculales = "#f26957",  
               Rhaphoneidales = "#f26957", Rhopalodiales = "#f26957", Surirellales = "#f26957", 
               Unclass.Bacillariophyceae = "#f26957",
               
               Corethrales = "#fcdeb4", Coscinodiscales = "#fcdeb4", Paraliales = "#fcdeb4", Rhizosoleniales = "#fcdeb4", 
               Unclass.Coscinodiscophyceae = "#fcdeb4", 
               
               Anaulales = "#fbb768", Chaetocerotales = "#fbb768", Cymatosirales = "#fbb768", Eupodiscales = "#fbb768", 
               Hemiaulales = "#fbb768", Lithodesmiales = "#fbb768", Thalassiosirales = "#fbb768",
               Unclass.Mediophyceae = "#fbb768",
               
               `Chrysophyceae_Clade-EC2H` = "#574D90", Ochromonadales = "#574D90", Paraphysomonadales = "#574D90", #Unclass.Chrysophyceae = "#574D90",
               Dictyochophyceae_X = "#574D90", Raphidophyceae_X = "#574D90", 
               
               Isochrysidales = "#9DB0CF", Phaeocystales = "#9DB0CF", Prymnesiales = "#9DB0CF", Prymnesiophyceae_X = "#9DB0CF", 
               Unclass.Prymnesiophyceae = "#9DB0CF"
)


#####---- v9 T3 Dendro setupt ----#####
#v9 T3 dendro text
dendro_text <- "(((Chlorodendrales), (Sphaeropleales), (Chloropicales), (Mamiellales), (Marsupiomonadales), (Pyramimonadales)), 
((Dinophysiales, Gonyaulacales, Gymnodiniales, Peridiniales, Prorocentrales, Suessiales, Torodiniales)), 
((Bacillariales, Cymbellales, Fragilariales, Licmophorales, Naviculales, Plagiogrammales, Rhabdonematales, Rhopalodiales),
  (Coscinodiscales, Melosirales, Stellarimales),
  (Ardissoneales, Chaetocerotales, Cymatosirales, Hemiaulales, Lithodesmiales, Probosciales, Thalassiosirales),
  (Parmales), (Chrysophyceae_Clade-EC2H, Hydrurales, Paraphysomonadales), (Dictyochophyceae_X), (Pinguiochrysidales), (Raphidophyceae_X)), 
((Pavlovales), (Prymnesiales, Prymnesiophyceae_X), 
  (Pavlomulinales)));"

# v9 T3
tip_cols <- c(
  Chlorodendrales = "#2f9181", Sphaeropleales = "#2f9181", #Unclass.Chlorophyceae = "#2f9181", 
  Chloropicales = "#2f9181",   Mamiellales = "#2f9181", Marsupiomonadales = "#2f9181", Pyramimonadales = "#2f9181", 
  
  Dinophysiales = "#517caf", Gonyaulacale = "#517caf", Gymnodiniales = "#517caf", Peridiniales = "#517caf", 
  Prorocentrales = "#517caf", Suessiales = "#517caf", Torodiniales = "#517caf", #Unclass.Dinophyceae = "#517caf",
  
  Bacillariales = "#f26957", Cymbellales = "#f26957", Fragilariales = "#f26957", Licmophorales = "#f26957", 
  Naviculales = "#f26957", Plagiogrammales = "#f26957", Rhabdonematales = "#f26957", Rhopalodiales = "#f26957" 
  , #Unclass.Bacillariophyceae = "#f26957",
  
  Coscinodiscales = "#fcdeb4", Melosirales = "#fcdeb4", Stellarimales = "#fcdeb4",
  
  Ardissoneales = "#fbb768", Chaetocerotales = "#fbb768", Cymatosirales = "#fbb768", Hemiaulales = "#fbb768", 
  Lithodesmiales = "#fbb768", Probosciales = "#fbb768", Thalassiosirales = "#fbb768"
  ,  #Unclass.Mediophyceae  = "#fbb768",
  
  Parmales = "#574D90", `Chrysophyceae_Clade-EC2H` = "#574D90", Hydrurales = "#574D90", Paraphysomonadales = "#574D90", 
  Dictyochophyceae_X = "#574D90", Pinguiochrysidales = "#574D90", Raphidophyceae_X = "#574D90",
  
  Pavlovales = "#9DB0CF", Prymnesiales = "#9DB0CF", Prymnesiophyceae_X = "#9DB0CF", #Unclass.Prymnesiophyceae = "#9DB0CF",
  Pavlomulinales = "#9DB0CF")



#####---- Dendrogram Plot ----#####

# Proceed with your code
dendro <- ape::read.tree(text = dendro_text)
# Specify text size
text_size <- 0.7

# Plot the dendrogram with customized settings
plot(dendro, type = "fan", cex = text_size, tip.color = tip_cols)



# Taxonomic information, make sure to do for all site/taxa group combinations:
taxa_df <- tax.clean[tax.clean$Domain %in% c("Eukaryota"), ]
class_counts <- table(taxa_df$Class)
order_counts <- table(taxa_df$Order)
family_counts <- table(taxa_df$Family)
genus_counts <- table(taxa_df$Genus)
species_counts <- table(taxa_df$Species)

print(order_counts)
sum(order_counts)
df <- data.frame(family_counts)

write.csv(class_counts, file = "C:/Documents/R/23-07-17_v9N1_classCounts.csv", row.names = TRUE)
write.csv(order_counts, file = "C:/Documents/R/23-07-17_v9N1_orderCounts.csv", row.names = TRUE)
write.csv(family_counts, file = "C:/Documents/R/23-07-17_v9N1_familyCounts.csv", row.names = TRUE)
write.csv(genus_counts, file = "C:/Documents/R/23-07-17_v9N1_genusCounts.csv", row.names = TRUE)
write.csv(species_counts, file = "C:/Documents/R/23-07-17_v9N1_speciesCounts.csv", row.names = TRUE)



# Print the counts of different families
print(family_counts)

# Print the counts of different genera
print(genus_counts)

# Print the counts of different species
print(species_counts)

# Access the names of different families
dino_taxa_df <- tax.clean[tax.clean$Class == "Chrysophyceae", ]
family_counts <- table(dino_taxa_df$Family)
family_names <- names(family_counts)
print(family_names)

# Access the names of different genera
class_names <- names(class_counts)
class_names

order_names <- names(Order_counts)
order_names

family_names <- names(family_counts)
family_names

genus_names <- names(genus_counts)
genus_names
# Access the names of different species
species_names <- names(species_counts)
species_names

#####---- Site Dendogram ----#####
# N1, both regions dendro text
dendro_text <- "(((Chlorodendrales), (Chlamydomonadales, Sphaeropleales), (Chloropicales), (Dolichomastigales, Mamiellales), 
  (Nephroselmidales), (Pseudoscourfieldiales, Pyramimonadales), (Chlorellales, Microthamniales, Watanabea-Clade)), 
((Dinophyceae_X, Dinophysiales, Gonyaulacales, Gymnodiniales, Peridiniales, Prorocentrales, Suessiales, Torodiniales)), 
((Bacillariales, Fragilariales, Naviculales, Rhaphoneidales, Rhopalodiales, Surirellales),
  (Corethrales, Coscinodiscales, Paraliales, Rhizosoleniales),
  (Anaulales, Chaetocerotales, Cymatosirales, Eupodiscales, Hemiaulales, Lithodesmiales, Thalassiosirales),
  (Chrysophyceae_Clade-EC2H, Ochromonadales, Paraphysomonadales), (Dictyochophyceae_X),(Raphidophyceae_X)),
((Coccolithales, Isochrysidales, Phaeocystales, Prymnesiales, Prymnesiophyceae_X)));"

# N1, both regions
tip_cols <- c( Chlorodendrales = "#2f9181", Chlamydomonadales = "#2f9181", Sphaeropleales = "#2f9181", Chloropicales = "#2f9181", 
               Dolichomastigales = "#2f9181", Mamiellales = "#2f9181", Nephroselmidales = "#2f9181", 
               Pseudoscourfieldiales = "#2f9181", Pyramimonadales = "#2f9181", Chlorellales = "#2f9181",  Microthamniales = "#2f9181", `Watanabea-Clade` = "#2f9181",
               
               Dinophyceae_X = "#517caf", Dinophysiales = "#517caf", Gonyaulacales = "#517caf", Gymnodiniales = "#517caf", 
               Peridiniales = "#517caf", Prorocentrales = "#517caf", Suessiales = "#517caf", Torodiniales = "#517caf",
               
               Bacillariales = "#f26957", Fragilariales = "#f26957", Naviculales = "#f26957",  
               Rhaphoneidales = "#f26957", Rhopalodiales = "#f26957", Surirellales = "#f26957",
               
               Corethrales = "#fcdeb4", Coscinodiscales = "#fcdeb4", Paraliales = "#fcdeb4", Rhizosoleniales = "#fcdeb4", 
               
               Anaulales = "#fbb768", Chaetocerotales = "#fbb768", Cymatosirales = "#fbb768", Eupodiscales = "#fbb768", 
               Hemiaulales = "#fbb768", Lithodesmiales = "#fbb768", Thalassiosirales = "#fbb768",
               
               `Chrysophyceae_Clade-EC2H` = "#574D90", Ochromonadales = "#574D90", Paraphysomonadales = "#574D90", #Unclass.Chrysophyceae = "#574D90",
               Dictyochophyceae_X = "#574D90", Raphidophyceae_X = "#574D90", 
               
               
               Coccolithales = "#9DB0CF", Isochrysidales = "#9DB0CF", Phaeocystales = "#9DB0CF", 
               Prymnesiales = "#9DB0CF", Prymnesiophyceae_X = "#9DB0CF"
)


#--------- T3 Both regions ----------#
#T3 dendro text
dendro_text <- "(((Chlorodendrales), (Sphaeropleales), (Chloropicales), (Mamiellales), (Marsupiomonadales), (Pyramimonadales), (Trebouxiophyceae_X)), 
((Dinophyceae_X, Dinophysiales, Gonyaulacales, Gymnodiniales, Peridiniales, Prorocentrales, Suessiales, Torodiniales)), 
((Bacillariales, Cymbellales, Fragilariales, Licmophorales, Naviculales, Plagiogrammales, Rhabdonematales, Rhopalodiales),
  (Coscinodiscales, Melosirales, Stellarimales),
  (Ardissoneales, Chaetocerotales, Cymatosirales, Hemiaulales, Lithodesmiales, Probosciales, Thalassiosirales),
  (Parmales), (Chrysophyceae_Clade-EC2H, Hydrurales, Paraphysomonadales), (Dictyochophyceae_X), (Pinguiochrysidales), (Raphidophyceae_X)), 
((Pavlovales), (Prymnesiales, Prymnesiophyceae_X), 
  (Pavlomulinales)));"

#T3
tip_cols <- c(
  Chlorodendrales = "#2f9181", Sphaeropleales = "#2f9181", 
  Chloropicales = "#2f9181",   Mamiellales = "#2f9181", Marsupiomonadales  = "#2f9181", Pyramimonadales = "#2f9181", Trebouxiophyceae_X = "#2f9181",
  
  Dinophysiales = "#517caf", Gonyaulacale = "#517caf", Gymnodiniales = "#517caf", Peridiniales = "#517caf", 
  Prorocentrales = "#517caf", Suessiales = "#517caf", Torodiniales = "#517caf",
  
  Bacillariales = "#f26957", Cymbellales = "#f26957", Fragilariales = "#f26957", Licmophorales = "#f26957", 
  Naviculales = "#f26957", Plagiogrammales = "#f26957", Rhabdonematales = "#f26957", Rhopalodiales = "#f26957", 
  
  Coscinodiscales = "#fcdeb4", Melosirales = "#fcdeb4", Stellarimales = "#fcdeb4",
  
  Ardissoneales = "#fbb768", Chaetocerotales = "#fbb768", Cymatosirales = "#fbb768", Hemiaulales = "#fbb768", 
  Lithodesmiales = "#fbb768", Probosciales = "#fbb768", Thalassiosirales = "#fbb768"
  ,  #Unclass.Mediophyceae  = "#fbb768",
  
  Parmales = "#574D90", `Chrysophyceae_Clade-EC2H` = "#574D90", Hydrurales = "#574D90", Paraphysomonadales = "#574D90", 
  Dictyochophyceae_X = "#574D90", Pinguiochrysidales = "#574D90", Raphidophyceae_X = "#574D90",
  
  Pavlovales = "#7b90af", Prymnesiales = "#7b90af", Prymnesiophyceae_X = "#7b90af",
  Pavlomulinales = "#7b90af")
