setwd("C:\\Users\\Tessa.Goerish\\OneDrive - Bio-Techne Corporation\\Desktop\\AVS554\\Example_Data_Files")

#Library Loadings ---
library(dada2); packageVersion('dada2')
library(Rqc)
library(Biostrings)
library(beepr) # if you want it to beep at you when it is done running
library(phyloseq); packageVersion("phyloseq")
library(vegan)
library(ape); packageVersion("ape")
library(plyr)
library(dplyr)
library(lme4)
library(lmerTest)
library(emmeans)
library(DESeq2)
library(permute)
library(randomForest)
library(rfPermute)
library(ggplot2)
library(RColorBrewer)
library(corrplot)
library(DECIPHER)
library(phangorn)
library(GenomeInfoDBData) # often R says it can't find package it just downloaded.  You can ether 1) manually check the box to enable this in the list of packages in R studio, or 2) manually enter the file location it gave you installed the package:
library("GenomeInfoDbData", lib.loc="~/R/win-library/4.0") #CHANGE ME

# Load raw data (Lab 2) ---------------------------
path_raw_F <- "C:\\Users\\Tessa.Goerish\\OneDrive - Bio-Techne Corporation\\Desktop\\AVS554\\Example_Data_Files\\Data\\scallop_Fb1_raw" 
path_raw_R <- "C:\\Users\\Tessa.Goerish\\OneDrive - Bio-Techne Corporation\\Desktop\\AVS554\\Example_Data_Files\\Data\\scallop_Fb2_raw"
fnsF <- list.files(path_raw_F)
fastqsF <- fnsF[grepl('.gz$', fnsF)] # CHANGE ME inside the '' if different file extensions, for example, might be '.fastq$'

fnsR <- list.files(path_raw_R)
fastqsR <- fnsR[grepl('.gz$', fnsR)] # CHANGE ME inside the '' if you have different file extensions

fnsF <- file.path(path_raw_F, fastqsF)
fnsR <- file.path(path_raw_R, fastqsR)

sample.names <- sapply(strsplit(basename(fastqsF), '_'), `[`, 1) 
names(fastqsF) <- sample.names 

sample.names2 <- sapply(strsplit(basename(fastqsR), '_'), `[`, 1) 
names(fastqsR) <- sample.names2


plotQualityProfile(fnsF[[1]]); beep()
plotQualityProfile(fnsR[[10]]); beep()

plotQualityProfile(fnsF, aggregate=T); beep() 
plotQualityProfile(fnsR, aggregate=T); beep()

## DO NOT RUN - Causes R to crash and terminate session :(
rqc(path = "C:\\Users\\Tessa.Goerish\\OneDrive - Bio-Techne Corporation\\Desktop\\AVS554\\Example_Data_Files\\Data\\scallop_Fb1_raw",
    sample = TRUE, n = 1e+06, group = NULL, top = 10, pair = NULL, pattern = ".fastq.gz",
    outdir = tempdir(), file = "rqc_report", openBrowser = TRUE, workers = multicoreWorkers())


### Summary of quality ----
# This dataset contains 58 samples Forward samples with 10114149 reads. There are 64 reverse samples with 5995530 reads. The reads are around 270 bases long. 
# The R1 quality goes below 30 around cycle 220 but going to trim at 250. The quality is good enough to use all R1.
# The R2 quality goes below 30 around cycle 140. The quality is okay for R2, trim sequence at 140. If not needed, drop the reverse reads.  

# You must include a Forward aggregate quality plot, and a Reverse aggregate quality plot if you have Revese reads, plus the Summary of Quality in an assignment on Brightspace and in your Manuscript. All other quality stats from this lecture are optional to include.

dir.create("C:\\Users\\Tessa.Goerish\\OneDrive - Bio-Techne Corporation\\Desktop\\AVS554\\Example_Data_Files\\Data\\Filtered_Data") #CHANGE ME
dir.create("C:\\Users\\Tessa.Goerish\\OneDrive - Bio-Techne Corporation\\Desktop\\AVS554\\Example_Data_Files\\Data\\Filtered_Data\\Filtered_forward")
dir.create("C:\\Users\\Tessa.Goerish\\OneDrive - Bio-Techne Corporation\\Desktop\\AVS554\\Example_Data_Files\\Data\\Filtered_Data\\Filtered_reverse")

path_filt_F <- "C:\\Users\\Tessa.Goerish\\OneDrive - Bio-Techne Corporation\\Desktop\\AVS554\\Example_Data_Files\\Data\\Filtered_Data\\Filtered_forward"
path_filt_R <- "C:\\Users\\Tessa.Goerish\\OneDrive - Bio-Techne Corporation\\Desktop\\AVS554\\Example_Data_Files\\Data\\Filtered_Data\\Filtered_reverse"


# run this whole chunk of code from filtoutput through verbose = TRUE)
filtoutput <- filterAndTrim(
  file.path(path_raw_F,fastqsF), # path to and file names of the raw reads
  file.path(path_filt_F, paste0(sample.names, "_F_filt.fastq.gz")), # set filtered file paths and filtered file names to create
  trimLeft=10, # cuts off the first XX bases from the F reads. Trim 10 for Illumina, 15 for 454 pyrosequencing.
  #trimRight=25, # cuts of last XX bases, or hash this out and use truncLen
  truncLen=250, # optional: cuts off end of F reads by trimming all to same length, use instead of trimRight
  #minLen=10, # optional: remove any sequences longer than max length, for use with 454 if don't want to truncate
  maxEE=2, # the maximum number of expected errors allowed in a read, always >1
  maxN= 0, # max number of ambiguous bases (N) allowed, DADA2 doesn't allow any!!
  rm.phix=TRUE, # remove any PhiX DNA (used as positive control),
  verbose=TRUE, matchIDs = TRUE); beep() # verbose = print to screen

# save the filtered variable for sample totals.  Your filtered files are already saved, this saves the summary R item to work with later without filtering again so you can track read numbers at each step.
saveRDS(filtoutput, "filtered_output_scallop_bac.rds") 

# check the dimensions of the variable you created, outputs two numbers: rows (# samples), columns (# of info it added)
dim(filtoutput) 

# take a look at the counts
filtoutput

# order the info by the first column (reads.in)
filtoutput[order(filtoutput[,1],decreasing=FALSE),]

# order the info by the second column (reads.out) which are filtered reads
filtoutput[order(filtoutput[,2],decreasing=FALSE),]

# get a sum total for raw reads in and filtered reads out
colSums(filtoutput)

# look at trends
filteroutputdf <- as.data.frame(filtoutput)

library(ggplot2)
ggplot(filteroutputdf) + 
  geom_point(aes(x = row.names(filteroutputdf), y = reads.in, color = "Reads In")) + 
  geom_point(aes(x = row.names(filteroutputdf), y = reads.out, color = "Reads Out")) + 
  scale_color_manual(values = c("Reads In" = "blue", "Reads Out" = "orange")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "Sample", y = "Reads", color = "Legend")
  
### Summary of filtering ----
# Filtering parameters: Trimleft=10 bp (illumina reads), and used TruncLen=250 since I knew the length I wanted to trim, Dropped Reverse reads due to discrepancy in sample labeling
# Largest raw reads = 757,154; Smallest raw reads = 53,246
# Largest filtered reads = 570,396; Smallest raw reads = 41,082
# total raw reads = 10,114,149
# total filtered reads = 8,022,887

set.seed(1234)
path_filt_F <- "C:\\Users\\Tessa.Goerish\\OneDrive - Bio-Techne Corporation\\Desktop\\AVS554\\Example_Data_Files\\Data\\Filtered_Data\\Filtered_forward"
#3. Set the names for Forward R1 filtered files. If your file names have a different pattern to the ones used in the example, you may need to update the strsplit portion.
fnsF <- list.files(path_filt_F, full.names = TRUE)
filtsF <- fnsF[grepl('.gz$', fnsF)] # CHANGE if different file extensions
sample.namesF <- sapply(strsplit(basename(filtsF), '_F_filt'), `[`, 1) # Assumes filename = samplename_XXX.fastq
names(filtsF) <- sample.namesF 

errF <- learnErrors(filtsF, nbases = 2e6, multithread=TRUE, randomize=TRUE); beep()
saveRDS(errF, 'error_profile_F.RDS')
plotErrors(errF, nominalQ=TRUE)
dadaFs <- dada(filtsF, err=errF, multithread=FALSE); beep()
seqtab <- makeSequenceTable(dadaFs); beep()
dim(seqtab)
svs_per_sample <- rowSums(seqtab > 0)
cat("SVs per sample range:", min(svs_per_sample), "to", max(svs_per_sample))
saveRDS(seqtab, "seqtab.rds")


### Summary of filtering ---- scallop tank bacteria
# Filtering parameters: Trimleft=10 bp (illumina reads), and used TruncLen=250 since I knew the length I wanted to trim, Dropped Reverse reads due to discrepancy in sample labeling
# Largest raw reads = 757,154; Smallest raw reads = 53,246
# Largest filtered reads = 570,396; Smallest raw reads = 41,082
# total raw reads = 10,114,149
# total filtered reads = 8,022,887

#Adjusted based filtering to 2e6 and multithread=TRUE due to the size of the dataset
#58 samples and 13,064 of SVs
#minimum SVs: 74, Maximum SVs: 1527

seqtab <- readRDS("C:\\Users\\Tessa.Goerish\\OneDrive - Bio-Techne Corporation\\Desktop\\AVS554\\Example_Data_Files\\seqtab.rds")
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=FALSE, verbose=TRUE) 
# CHANGE ME: My samples have 4330 bimeras out of 13064 total sequences

dim(seqtab.nochim) 
# CHANGE ME: 58 samples and 8734 SVs

sum(seqtab.nochim)/sum(seqtab)
saveRDS(seqtab.nochim, "C:\\Users\\Tessa.Goerish\\OneDrive - Bio-Techne Corporation\\Desktop\\AVS554\\Example_Data_Files\\seqtab.nochim.rds")
require(tidyr)
library(ggplot2)

seqtab.nochim <- readRDS("C:\\Users\\Tessa.Goerish\\OneDrive - Bio-Techne Corporation\\Desktop\\AVS554\\Example_Data_Files\\seqtab.nochim.rds")
filtoutput <- readRDS("C:\\Users\\Tessa.Goerish\\OneDrive - Bio-Techne Corporation\\Desktop\\AVS554\\Example_Data_Files\\filtered_output_scallop_bac.rds")

getN <- function(x) sum(getUniques(x))
meta <- read.csv("C:\\Users\\Tessa.Goerish\\OneDrive - Bio-Techne Corporation\\Desktop\\AVS554\\Example_Data_Files\\Scallop_metadata.csv", # CHANGE file name
                 header = TRUE, # specifies that the top row is column names
                 row.names = 1)



dim(filtoutput) # 58 rows
dim(seqtab.nochim) # 58 rows
dim(meta) # 58 rows

meta_cleaned <- meta[row.names(meta) %in% row.names(seqtab.nochim), ]
track <- cbind(filtoutput, rowSums(seqtab.nochim), meta$DNA_extract_tube_IDs)
track <- data.frame(DNA_extract_tube_IDs = rownames(track), track, row.names = NULL)
colnames(track) <- c("DNA_extract_tube_IDs", "reads.in","filtered", "nonchimeras")
rownames(track) <- rownames(meta)
head(track)
saveRDS(track, 'tracked_seqs.RDS')

track <- readRDS("C:\\Users\\Tessa.Goerish\\OneDrive - Bio-Techne Corporation\\Desktop\\AVS554\\Example_Data_Files\\tracked_seqs.RDS")

plotData <- as.data.frame(track) %>% gather(type,totals, reads.in, filtered, nonchimeras)
plotData$type <- factor(plotData$type, levels = c("reads.in", "filtered", "nonchimeras"))
ggplot(plotData,aes(x=DNA_extract_tube_IDs,y=as.numeric(totals))) + geom_jitter(aes(color=type)) + #CHANGE Sample_type to be the factor you previously chose
  ylab("Sequences") + xlab("DNA_extract_tube_IDs") 

##figure out the negative control (this is wrong) -> metadata: negative controls are 8b, 12b, 105, 129 -> best one looks like 12b
unqs.NC <- seqtab.nochim["12b",]
unqs.NC <- sort(unqs.NC[unqs.NC>0], decreasing=TRUE)
cat("DADA2 inferred", length(unqs.NC), "sequence variants present in the selected sample.\n")
# CHANGE ME: There are 227 SVs in the negative control
plot(unqs.NC, ylab="Number of seqs/SV, Neg Control", xlab="SVs", main="Number of sequences/SV in Negative Control") #CHANGE ME

#positive control -> WT is sample 10
unqs.PC <- seqtab.nochim["10",]
unqs.PC <- sort(unqs.PC[unqs.PC>0], decreasing=TRUE)
cat("DADA2 inferred", length(unqs.PC), "sequence variants present in the selected sample.\n")
# CHANGE ME: There are 372 SVs in this positive control. 
plot(unqs.PC, ylab="Number of seqs/SV, Neg Control", xlab="SVs", main="Number of sequences/SV in Positive Control") #CHANGE ME


library(phyloseq)
seqtab.nochim <- readRDS("C:\\Users\\Tessa.Goerish\\OneDrive - Bio-Techne Corporation\\Desktop\\AVS554\\Example_Data_Files\\seqtab.nochim.rds")
meta <- read.csv("C:\\Users\\Tessa.Goerish\\OneDrive - Bio-Techne Corporation\\Desktop\\AVS554\\Example_Data_Files\\Scallop_metadata.csv", header=T, row.names=1)
row.names(meta)
row.names(seqtab.nochim)
EX_ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
                  sample_data(meta))
EX_ps
# how many samples made it this far? 58 samples, 24 variables, 8734 taxa
library(ggplot2)
plot_richness(EX_ps, x="Sample_type", # CHANGE the x-axis to a factor of your choice
              measures=c("Observed","Chao1", "Shannon"), # these are some of the alpha diversity measures you can use 
              ) + # CHANGE the color to a factor of your choice
  theme_bw() # a ggplot theme to make the graph look nice

## Note: hopefully, diversity/richness is lower in the Negative controls than real samples
## SUMMARY: The negative controls for the most part are lower

plot(sort(taxa_sums(EX_ps), TRUE), 
     type="h", 
     ylim=c(0, 20))
EX.ord <- ordinate(EX_ps, #calculate similarities
                   method ="PCoA", #ordination type
                   "jaccard", binary = TRUE)
plot_ordination(EX_ps, EX.ord, type="samples", color="Sample_type")
sample_data(EX_ps)


# 1, Subset the controls that would pertain to all of the samples, like ones from reagents, PCR, or other whole-batch effects. prune to be only those taxa	
EX_ps_controls = subset_samples(EX_ps, Sample_type == "Negative_control")

EX_ps_controls <- prune_taxa(taxa_sums(EX_ps_controls) > 0, EX_ps_controls)
EX_ps_controls
#587 taxa in the 4 mock samples

# 2, Make the taxa names into a vector so you can remove them. 
control_vec <- as.vector(taxa_names(EX_ps_controls))
vec <- as.vector(taxa_names(EX_ps))
keep <- setdiff(vec, control_vec)

# 3. Use the keep vector for the prune taxa argument, because it wants the argument to be true (matching).
EX_ps_NCclean <- prune_taxa(keep, EX_ps)

# 4. Then remove the samples which are now empty, namely the NegCon_PCR and NegCon_swab.
EX_ps_NCclean <- prune_samples(sample_sums(EX_ps_NCclean) > 0, EX_ps_NCclean)

EX_ps_NCclean
# NOTE: If you have multiple batches of DNA extractions, you'll want to copy the code from steps 5,6,7 below and have a set for each batch. If you don't have negative controls for each DNA extraction, jump to step 10.
#not all extraction batches had negative controls, skip this part

## 9. clean out taxa/SV columns that are no longer present
EX_ps_NC_batch_clean <- prune_taxa(taxa_sums(EX_ps_NCclean) > 0, EX_ps_NCclean)
EX_ps_NC_batch_clean <- prune_samples(sample_sums(EX_ps_NCclean) > 0, EX_ps_NCclean)
EX_ps_NC_batch_clean
# 54 number of samples and 24 SVs remaining, 8147 taxa

# 10. Did it work? Check your ordination again
EX_cleaner.ord <- ordinate(EX_ps_NC_batch_clean, #calculate similarities
                           method ="PCoA", #ordination type
                           "jaccard", binary = TRUE) #similarity type. Jaccard is binary, Bray can be binary (unweighted) or not (weighted)

plot_ordination(EX_ps_NC_batch_clean, EX_cleaner.ord, type="samples", color="Sample_type", title="After cleaning out negative controls")

# save this graph for later

### Summary of decontamination ----
# Used all negative controls, both ethanol and DI water
# CHANGE ME: Did you notice clustering by batch before and after removing contaminants?

#for Taxonomy use SILVA 16S reference database

seqtab.nochim.decontam = as(otu_table(EX_ps_NC_batch_clean), "matrix")

seqtab.nochim.decontam = as.matrix(seqtab.nochim.decontam)

all.taxa <- assignTaxonomy(seqtab.nochim.decontam, 
                           "C:\\Users\\Tessa.Goerish\\OneDrive - Bio-Techne Corporation\\Desktop\\AVS554\\Example_Data_Files\\SILVA\\silva_nr99_v138.2_toSpecies_trainset.fa.gz",        # CHANGE file path
                           tryRC = TRUE, # didn't get any IDs the first time? Use this to try the reverse complement of your sequences instead
                           minBoot = 75, verbose = TRUE); beep() 

saveRDS(all.taxa, "C:\\Users\\Tessa.Goerish\\OneDrive - Bio-Techne Corporation\\Desktop\\AVS554\\Example_Data_Files\\Data\\taxa_silva.rds") #CHANGE ME
write.csv(all.taxa,"C:\\Users\\Tessa.Goerish\\OneDrive - Bio-Techne Corporation\\Desktop\\AVS554\\Example_Data_Files\\Data\\taxa_silva.csv") #CHANGE ME

all.taxa.sp <- addSpecies(all.taxa, "C:\\Users\\Tessa.Goerish\\OneDrive - Bio-Techne Corporation\\Desktop\\AVS554\\Example_Data_Files\\SILVA\\silva_v138.2_assignSpecies.fa.gz", allowMultiple = FALSE, verbose = FALSE); beep() #CHANGE file path

saveRDS(all.taxa.sp, "C:\\Users\\Tessa.Goerish\\OneDrive - Bio-Techne Corporation\\Desktop\\AVS554\\Example_Data_Files\\Data\\taxa_species_silva.rds") #CHANGE ME
write.csv(all.taxa.sp, "C:\\Users\\Tessa.Goerish\\OneDrive - Bio-Techne Corporation\\Desktop\\AVS554\\Example_Data_Files\\Data\\taxa_species_silva.csv") #CHANGE ME

otu_t <- otu_table(EX_ps_NC_batch_clean)

all.taxa <- readRDS("C:\\Users\\Tessa.Goerish\\OneDrive - Bio-Techne Corporation\\Desktop\\AVS554\\Example_Data_Files\\Data\\taxa_silva.rds")
all.taxa.sp <- readRDS("C:\\Users\\Tessa.Goerish\\OneDrive - Bio-Techne Corporation\\Desktop\\AVS554\\Example_Data_Files\\Data\\taxa_species_silva.rds")


EX_ps_NC_batch_clean <- phyloseq(otu_table(otu_t, taxa_are_rows=FALSE), 
                                 sample_data(meta),
                                 tax_table(all.taxa.sp))

saveRDS(EX_ps_NC_batch_clean, "Example_phyloseq_object_cleaned.rds")

require(dplyr)

EX_ps_clean <- EX_ps_NC_batch_clean %>% # CHANGE ME to your phyloseq object name 
  subset_taxa(Kingdom == "Bacteria" & Family != "Mitochondria" & Order != "Chloroplast") # CHANGE ME to taxa you want to remove


EX_ps_clean <- prune_taxa(taxa_sums(EX_ps_clean) > 0, EX_ps_clean)
EX_ps_clean <- prune_samples(sample_sums(EX_ps_clean) > 0, EX_ps_clean)

EX_ps_clean
# CHANGE ME: 54 number of samples and 5267 SVs left 

saveRDS(EX_ps_clean, 'EX_ps_clean_phyloseq_object.RDS') #CHANGE ME
EX_ps_clean <-readRDS('EX_ps_clean_phyloseq_object.RDS') 

library(vegan)
# Extract OTU table as a matrix
otu_mat <- as(otu_table(EX_ps_clean), "matrix")

# If taxa are rows (common in phyloseq), transpose it
if (taxa_are_rows(EX_ps_clean)) {
  otu_mat <- t(otu_mat)
}
EX_rarec <- rarecurve(otu_mat, step = 10, cex = 0.5, label = TRUE)
sort(rowSums(otu_table(EX_ps_clean)))
# CHANGE ME: smallest reads (sequences) in a sample 353
# CHANGE ME: largest in a sample 128595
# CHANGE ME: number of samples with <5000 23

EX_ps_clean.rar <- rarefy_even_depth(EX_ps_clean, 
                                     sample.size=3000, # CHANGE ME to the sequences/sample you want. 5-10k is a good amount, more is better
                                     replace=FALSE, #sampling with or without replacement
                                     trimOTUs=TRUE, #remove SVs left empty (called OTUs here but really they are SVs) 
                                     rngseed=711, 
                                     verbose=TRUE)

# CHANGE ME: how many 16 samples and 2780 SVs were removed
saveRDS(EX_ps_clean.rar, 'EX_ps_clean_rarefied_phyloseq_object.RDS')
write.csv(otu_table(EX_ps_clean.rar), 'EX_ps_clean_rarefied_phyloseq_object.csv')

EX_ps_clean.rar <- readRDS('EX_ps_clean_rarefied_phyloseq_object.RDS')

library(ggplot2)
library(RColorBrewer)

# optional, if you want factors to graph in a specific order, you can set that manually, and relabel them so they are more readable in the graph
sample_data(EX_ps_clean.rar)$Diet <- factor(sample_data(EX_ps_clean.rar)$Diet, levels=c("LooseAlfalfa", "PelletedAlfalfa"), labels=c("Loose Alfalfa", "Pelleted Alfalfa")) #CHANGE ME

# plot alpha diversity with phyloseq: https://www.rdocumentation.org/packages/phyloseq/versions/1.16.2/topics/plot_richness. 
# measures include c("Observed", "Chao1", "ACE", "Shannon", "Simpson", "InvSimpson", "Fisher")
sample_variables(EX_ps_clean.rar)

plot_richness(EX_ps_clean.rar, 
              x="Clean_or_dirty", #CHANGE ME, A is whatever factor you want on x-axis
              measures="Observed",#CHANGE ME, whatever alpha diversity measure you want. to have multiple, use: = c("Observed","Shannon")
              title = NULL) + 
  theme_set(theme_minimal(base_size = 14)) + #make it look pretty
  #theme(axis.text.x=element_blank()) + #example, get rid of x axis labels
  geom_violin(trim=TRUE, aes(fill=Tank_contents)) + #optional. CHANGE ME, A is whatever factor to color violins
  geom_boxplot(width = 0.1, aes(group=Clean_or_dirty)) + #optional. CHANGE ME, A is whatever factor to group boxplots
  facet_grid(.~Tank_contents, space="free") +
  # theme(legend.position = "none") + #use to get rid of your legend
  ylab("Alpha Diversity (Observed & Shannon") +
  xlab("Sample Type")

plot_richness(EX_ps_clean.rar, 
              x="Tank_setup", #CHANGE ME, A is whatever factor you want on x-axis
              measures=c("Observed","Shannon"),#CHANGE ME, whatever alpha diversity measure you want. to have multiple, use: = c("Observed","Shannon")
              title = NULL) + 
  theme_set(theme_minimal(base_size = 14)) + #make it look pretty
  #theme(axis.text.x=element_blank()) + #example, get rid of x axis labels
  geom_violin(trim=TRUE, aes(fill=Clean_or_dirty)) + #optional. CHANGE ME, A is whatever factor to color violins
  geom_boxplot(width = 0.1, aes(group=Tank_setup)) + #optional. CHANGE ME, A is whatever factor to group boxplots
  facet_grid(.~Clean_or_dirty, space="free") +
  # theme(legend.position = "none") + #use to get rid of your legend
  ylab("Alpha Diversity (Observed & Shannon)") +
  xlab("Tank Setup")

# Use geoms to add different shapes/styles of plots: https://ggplot2-book.org/individual-geoms.html#:~:text=3.1%20Basic%20plot%20types,plot%20has%20a%20special%20name.

# Use multiple geoms to create gaph overlays in a single plot, and rearrange their order to move to the foreground or background.

# want to map factors to geoms to get it to set colors by factor? use geom_point(aes(color=Factor)). otherwise, use geom_point(color=Black) to make them all the same.


# use phyloseq to measure alpha diversity
EX_ps_clean.rar.rich <- estimate_richness(EX_ps_clean.rar, measure=c("Observed", "Shannon")) #change to whatever measures you want
EX_ps_clean.rar.rich
# OPTIONAL: use phyloseq to calculate Faith's Diversity metric, https://rdrr.io/github/twbattaglia/btools/man/estimate_pd.html
install.packages("remotes")
library(remotes)
remotes::install_github("twbattaglia/btools")
library(btools)
seqs <- taxa_names(EX_ps_clean.rar)

# 2. Create a DNAStringSet
dna <- DNAStringSet(seqs)

# 3. Align sequences
alignment <- AlignSeqs(dna)

# 4. Make a distance matrix
phydist <- dist.ml(as.phyDat(alignment))

# 5. Build a tree (NJ or ML)
tree <- NJ(phydist)          # neighbor-joining
tree <- ladderize(tree)      # optional for plotting
taxa_names(EX_ps_clean.rar)[1:10]
tree$tip.label[1:10]
tree$tip.label <- taxa_names(EX_ps_clean.rar)
plot(tree)



# 6. Add tree to phyloseq object
EX_ps_clean.rar <- merge_phyloseq(EX_ps_clean.rar, phy_tree(tree))

EX_faith <- estimate_pd(EX_ps_clean.rar)
EX_faith

# measure evenness for each sample
EX_ps_clean.rar.even <- EX_ps_clean.rar.rich$Shannon/log(EX_ps_clean.rar.rich$Observed)

# Coerce to data.frame and add the metadata for these samples
EX_ps_clean.rar.sd = as(sample_data(EX_ps_clean.rar), "matrix")
EX_ps_clean.rar.sd = as.data.frame(EX_ps_clean.rar.sd)
EX_ps_clean.rar.rich.df <- cbind(EX_ps_clean.rar.rich, EX_ps_clean.rar.even, EX_ps_clean.rar.sd)

dim(EX_ps_clean.rar.rich.df)

saveRDS(EX_ps_clean.rar, 'EX_ps_clean_rarefied_phyloseq_object.RDS')

Ex_ps_clean.rar <- readRDS('EX_ps_clean_rarefied_phyloseq_object.RDS')


summary(aov(Shannon ~ Sample_type, data = EX_ps_clean.rar.rich.df))
ggplot(EX_ps_clean.rar.rich.df, aes(x = Sample_type, y = Observed)) +
  geom_boxplot() +
  theme_bw()


# make a graph using that dataframe. this is a generic example, you will want to personalize this.
ggplot(data=EX_ps_clean.rar.rich.df, aes(x=Sample_type, y=Observed)) + 
  theme_minimal() + 
  geom_point(aes(color=Location), size = 3) + # sets colors by group and a set point size 
  xlab("Sample Type") + 
  ylab("Observed Bacterial Richness (SVs)") + 
  theme(text = element_text(size = 15)) # increases font size


# make that same graph but drop any samples that lack data for that FactorA (replace FactorA in the code with your factor name).
ggplot(data=subset(EX_ps_clean.rar.rich.df, !is.na(Tank_setup)), aes(x=Sample_type, y=Observed)) + 
  theme_minimal() + 
  geom_point(aes(color=Location), size = 3) + # sets colors by group and a set point size 
  xlab("Sample Type") + 
  ylab("Observed Bacterial Richness (SVs)") + 
  theme(text = element_text(size = 15)) # increases font size



# Make an alpha diversity table: https://www.geeksforgeeks.org/create-table-from-dataframe-in-r/
EX_table = as.table(table(EX_ps_clean.rar.rich.df$Sample_type, EX_ps_clean.rar.rich.df$Location)) 
EX_table


# use phyloseq to measure alpha diversity
EX_ps_clean.rar.rich <- estimate_richness(EX_ps_clean.rar, measure=c("Observed", "Shannon")) #change to whatever measures you want

# use phyloseq to calculate Faith's Diversity metric (optional), https://rdrr.io/github/twbattaglia/btools/man/estimate_pd.html
EX_faith <- estimate_pd(EX_ps_clean.rar)
library(asbio)
library(microbiome)
EX_ps_clean.rar.even_SV <- evenness(otu_table(EX_ps_clean.rar))

# measure evenness for each sample
EX_ps_clean.rar.even <- EX_ps_clean.rar.rich$Shannon/log(EX_ps_clean.rar.rich$Observed)
EX_ps_clean.rar.sd = as(sample_data(EX_ps_clean.rar), "matrix")
EX_ps_clean.rar.sd = as.data.frame(EX_ps_clean.rar.sd)
EX_ps_clean.rar.rich.df <- cbind(EX_ps_clean.rar.rich, EX_ps_clean.rar.even, EX_ps_clean.rar.sd)

# make a histogram to look at the shape of the data (bell curve? skew?). You can save this graph for your own benefit if you want.
hist(EX_ps_clean.rar.rich.df$Observed)
hist(EX_ps_clean.rar.rich.df$Shannon)
# Want to know how much skew there is? Measure Kurtosis (a.k.a. tailedness). You can save this graph for your own benefit if you want.
install.packages("PerformanceAnalytics")
library(PerformanceAnalytics)
kurtosis(EX_ps_clean.rar.rich.df$Observed)
# output is 0.1545592, looking for whether +/- and if a large/small number
kurtosis(EX_ps_clean.rar.rich.df$Shannon)
# output is 0.02291181, looking for whether +/- and if a large/small number


head(EX_ps_clean.rar.rich.df)
## Use Shannon - normally distributed
library(vegan)
library(lme4)
library(lmerTest)
library(emmeans)
library(dplyr)
  
  # Subset tank samples
  tank_df <- subset(EX_ps_clean.rar.rich.df, Environment == "Tank")
  
  # Recode flow type for clarity
tank_df$Flow_type <- factor(tank_df$Tank_setup, levels = c("static", "flow_through"))
library(dplyr)
  
tank_df %>%
    group_by(Flow_type) %>%
    summarise(
      n = n(),
      Shannon_mean = mean(Shannon, na.rm = TRUE),
      Shannon_sd = sd(Shannon, na.rm = TRUE),
      Observed_mean = mean(Observed, na.rm = TRUE),
      Observed_sd = sd(Observed, na.rm = TRUE)
    )
saveRDS(tank_df, "tank_alpha_diversity_df.RDS")
tank_df <- readRDS("tank_alpha_diversity_df.RDS")

tank_df_clean <- subset(tank_df, Clean_or_dirty != "n/a")
library(ggplot2)
  
  # Shannon diversity
ggplot(tank_df_clean, aes(x = Flow_type, y = Shannon, fill = Flow_type)) +
    geom_boxplot() +
    geom_jitter(width = 0.1, alpha = 0.7) +
    theme_bw() +
    labs(title = "Shannon Diversity by Tank Flow Type",
         x = "Flow Type", y = "Shannon Index") +
    scale_fill_manual(values = c("static" = "#56B4E9", "flow_through" = "#D55E00"))
  
ggplot(tank_df_clean, aes(x = Clean_or_dirty, y = Shannon, fill = Clean_or_dirty)) +
  geom_boxplot() +
  geom_jitter(width = 0.1, alpha = 0.7) +
  theme_bw() +
  labs(title = "Shannon Diversity by Tank State",
       x = "State of Tank", y = "Shannon Index") +
  scale_fill_manual(values = c("after_cleaned_filled_with_water" = "Green" , "after_drain_down" = "Pink"))
  # Observed richness (optional)
ggplot(tank_df_clean, aes(x = Flow_type, y = Observed, fill = Flow_type)) +
    geom_boxplot() +
    geom_jitter(width = 0.1, alpha = 0.7) +
    theme_bw() +
    labs(title = "Observed Richness by Tank Flow Type",
         x = "Flow Type", y = "Observed Richness") +
    scale_fill_manual(values = c("static" = "#56B4E9", "flow_through" = "#D55E00"))
  
shapiro.test(tank_df$Shannon)
shapiro.test(tank_df$Observed)
shapiro.test(tank_df$EX_ps_clean.rar.even)

wilcox.test(Shannon ~ Flow_type, data = tank_df)
wilcox.test(Observed ~ Flow_type, data = tank_df)
wilcox.test(EX_ps_clean.rar.even ~ Flow_type, data = tank_df)



### ANOVA for multiple categorical factors -------
# ANOVA is used when your factor has 2+ variables in it. If you have 3+ variables, you should run a post-hoc test to correct the p-value for multiple group comparisons


fit_shannon <- aov(Shannon ~ Flow_type, data = tank_df_clean) #does flow type influence microbial diversity?
summary(fit_shannon)
#p=0.0903, marginally non-significant, trending but cant reject null hypothesis
wilcox.test(Shannon ~ Flow_type, data=tank_df_clean)
#p=0.1807 -> No evidence shannon diversity differs between the two flow types for the tank samples
ggplot(tank_df_clean, aes(x = Flow_type, y = Shannon, color = Flow_type)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.5) +     
  geom_jitter(width = 0.15, size = 3, alpha = 0.8) +  
  stat_summary(fun = mean, geom = "point", shape = 18, size = 4, color = "black")
  labs(
    title = "Shannon Diversity Across Tank Flow Types",
    x = "Flow Type",
    y = "Shannon Diversity (Alpha Diversity)"
  ) +
  theme_minimal(base_size = 15) +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, face = "bold")
  )


fit_shannon_cd <- aov(Shannon ~ Clean_or_dirty, data=tank_df_clean) #what is the effect of tank cleaning?
summary(fit_shannon_cd)
#p=0.94, shannon diversity not affected by the state of the tank when sample is taken
kruskal.test(Shannon ~ Clean_or_dirty, data = tank_df_clean)
#p=0.9624, significant results, run Post-hoc
pairwise.wilcox.test(tank_df_clean$Shannon, tank_df_clean$Clean_or_dirty, p.adjust.method = "BH")

fit_observed <- aov(Observed ~ Flow_type, data = tank_df_clean)
summary(fit_observed)
#P=0.255
wilcox.test(Observed ~ Flow_type, data=tank_df_clean)
#p=0.3553

ggplot(tank_df_clean, aes(x = Flow_type, y = Observed, color = Flow_type)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.5) +     
  geom_jitter(width = 0.15, size = 3, alpha = 0.8) +  
  stat_summary(fun = mean, geom = "point", shape = 18, size = 4, color = "black")+
  labs(
  title = "Observed Diversity Across Tank Flow Types",
  x = "Flow Type",
  y = "Observed Diversity (Alpha Diversity)"
  ) +
  theme_minimal(base_size = 15) +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, face = "bold")
  )


fit_observed_cd <- aov(Observed ~ Clean_or_dirty, data = tank_df_clean)
summary(fit_observed_cd)
#p=0.303
kruskal.test(Observed ~ Clean_or_dirty, data = tank_df_clean)
#p=0.5554
pairwise.wilcox.test(tank_df_clean$Observed, tank_df_clean$Clean_or_dirty, p.adjust.method = "BH")




# Two Way Factorial Design
fit_twoway <- aov(Shannon ~ Flow_type * Clean_or_dirty, data=tank_df_clean) 
summary(fit_twoway)
#p flow=0.103, clean or dirty=0.942, compared=0.746 -> not significant

# Analysis of Covariance
fit_co <- aov(Shannon ~ Flow_type + Days_since_full_moon, data=tank_df_clean) # x is a series of numerical values, like pH)
summary(fit_co)

ggplot(tank_df_clean, aes(x = Days_since_full_moon, y = Shannon, color = Flow_type)) +
  geom_point(size = 3, alpha = 0.8) +
  geom_smooth(method = "lm", se = TRUE) +  # regression line for each flow type
  labs(
    title = "Shannon Diversity vs. Flow Type with Lunar Phase as Covariate",
    x = "Days Since Full Moon",
    y = "Shannon Diversity (Alpha Diversity)"
  ) +
  theme_minimal(base_size = 15) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))


saveRDS(tank_df_clean, "tank_alpha_diversity_df_cleaned.RDS")


summary(lm(Shannon ~ Flow_type, data= tank_df_clean))
ggplot(tank_df_clean, aes(x = Flow_type, y = Shannon, color = Flow_type)) +
  geom_jitter(width = 0.1, size = 2, alpha = 0.7) +       # raw points
  stat_summary(fun = mean, geom = "point", size = 4) +     # mean
  stat_summary(fun.data = mean_cl_normal, geom = "errorbar", width = 0.2) + # 95% CI
  labs(
    title = "Shannon Diversity by Flow Type",
    subtitle = "Linear model trend: p = 0.090",
    x = "Flow Type",
    y = "Shannon Diversity"
  ) +
  theme_minimal(base_size = 15) +
  theme(legend.position = "none")


summary(lm(Shannon ~ Clean_or_dirty, data= tank_df_clean))
ggplot(tank_df_clean, aes(x = Clean_or_dirty, y = Shannon, color = Clean_or_dirty)) +
  geom_jitter(width = 0.1, size = 2, alpha = 0.7) +       # raw points
  stat_summary(fun = mean, geom = "point", size = 4) +     # mean
  stat_summary(fun.data = mean_cl_normal, geom = "errorbar", width = 0.2) + # 95% CI
  labs(
    title = "Shannon Diversity by State of Tank",
    subtitle = "Linear model trend: p = 0.090",
    x = "State of Tank",
    y = "Shannon Diversity"
  ) +
  theme_minimal(base_size = 15) +
  theme(legend.position = "none")

#Observed
linear_observed <- (glmer(Observed ~ Flow_type + Clean_or_dirty + (1|Tank_ID_match_to_the_isolate_Tank_ID), family=poisson, data=tank_df_clean))
summary(linear_observed)

linear_shannon <- lmer(Shannon ~ Flow_type + Clean_or_dirty + (1|Tank_ID_match_to_the_isolate_Tank_ID), data=tank_df_clean)
summary(linear_shannon)

plot_heatmap(EX_ps_clean.rar)
# follow the tutorial to make prettier versions in phyloseq: https://joey711.github.io/phyloseq/plot_heatmap-examples.html


topNOTUs = names(sort(taxa_sums(EX_ps_clean.rar), TRUE)[1:15])
top15 = prune_taxa(topNOTUs, EX_ps_clean.rar)

# combine with your metadata and create one dataframe object. you can include other info that you created for a previous dataframe, as long as those objects are still in your R environment. Reminder, you can only use numeric data in a correlation matrix, so you will have to drop certain columns or make them numeric instead.
# Coerce to data.frame and add the metadata for these samples
top15.sd = as(otu_table(top15), "matrix")
top15.sd = as.data.frame(top15.sd)

# add Genus names in place of the full sequence name that is in the SV columns
top15.tax <- as.data.frame(tax_table(top15))
## if the Genus is empty, replace with the Family
top15.tax$Genus[is.na(top15.tax$Genus)] <- top15.tax$Family[is.na(top15.tax$Genus)]
colnames(top15.sd) = top15.tax$Genus

# paste all the components together
EX_ps_clean.rar.corr.df <- cbind(top15.sd, EX_ps_clean.rar.rich, EX_ps_clean.rar.even, EX_ps_clean.rar.sd)

# check header to make sure it looks good
head(EX_ps_clean.rar.corr.df)

# change any column factor names to make them prettier
names(EX_ps_clean.rar.corr.df)[names(EX_ps_clean.rar.corr.df) == "EX_ps_clean.rar.even"] <- "Evenness"

# check header to make sure it looks good
head(EX_ps_clean.rar.corr.df)

# for example, drop a column with factorial data
# Keep only numeric columns
EX_ps_clean.rar.corr.df_numeric <- EX_ps_clean.rar.corr.df[, sapply(EX_ps_clean.rar.corr.df, is.numeric)]

# Check the structure
str(EX_ps_clean.rar.corr.df_numeric)

# run correlations
ex_corr_calc <- cor(EX_ps_clean.rar.corr.df_numeric, 
                    use="complete.obs", # use=Specifies the handling of missing data. Options are all.obs (assumes no missing data - missing data will produce an error), complete.obs (listwise deletion), and pairwise.complete.obs (pairwise deletion)
                    method="kendall") # correlation method=pearson, spearman, or kendall

### Note, if you have too few samples, you may receive an error about too few observations. You may ignore the error message, or remove than column



### test significance of correlations, copy and paste this whole chunk together
cor.mtest <- function(mat, conf.level = 0.95){
  mat <- as.matrix(mat)
  n <- ncol(mat)
  p.mat <- lowCI.mat <- uppCI.mat <- matrix(NA, n, n)
  diag(p.mat) <- 0
  diag(lowCI.mat) <- diag(uppCI.mat) <- 1
  for(i in 1:(n-1)){
    for(j in (i+1):n){
      tmp <- cor.test(mat[,i], mat[,j], conf.level = conf.level)
      p.mat[i,j] <- p.mat[j,i] <- tmp$p.value
      lowCI.mat[i,j] <- lowCI.mat[j,i] <- tmp$conf.int[1]
      uppCI.mat[i,j] <- uppCI.mat[j,i] <- tmp$conf.int[2]
    }
  }
  return(list(p.mat, lowCI.mat, uppCI.mat))
}

res1 <- cor.mtest(EX_ps_clean.rar.corr.df_numeric,0.95)
res2 <- cor.mtest(EX_ps_clean.rar.corr.df_numeric,0.99)
# Note: if you come up with a warning message, it means that one or more of your columns generated correlations of value 0, which makes it angry.  Visualize the plot, and "?" will come up in those columns.  To fix, remove the column or add a data transformation, to the dataframe you created to make this.  Set the corr.mtest function and run corr.mtest again.
# Use all samples in your OTU table

# 5. Compute Spearman correlations with Hmisc::rcorr
library(Hmisc)
rcorr_res <- rcorr(as.matrix(EX_ps_clean.rar.corr.df_numeric), type = "kendall")
ex_corr_calc <- rcorr_res$r
ex_corr_pval <- rcorr_res$P



library(corrplot)
library(RColorBrewer)
## plot only those correlations with a significant p-value <0.05, using hierarchical clustering
corrplot(ex_corr_calc, 
         type="lower", #shape of the plot itself: full, upper, lower
         method="circle", #shape of the data: circle, square, ellipse, number, shade, color, pie 
         order="hclust", #how to cluster samples: AOE, hclust, FPC, alphabet, or leave blank
         p.mat = ex_corr_pval, #which correlations to use
         sig.level=0.05, #sets significance cutoff
         insig="blank", #leaves p > 0.05 blank
         tl.col="black", # text color
         tl.cex=.9, #text size
         col=brewer.pal(n=10, name="RdYlBu")) #specify a color palette and number of shades needed

plot_bar(EX_ps_clean.rar, fill="Phylum")
EX_ps_clean.rar.glom = tax_glom(EX_ps_clean.rar, "Phylum")

# WARNING: if you group along the x-axis phyloseq just stacks seqence counts, it DOES NOT calculate show average for that group
plot_bar(EX_ps_clean.rar.glom, fill="Phylum")
EX_ps_clean.rar.stacked = transform_sample_counts(EX_ps_clean.rar, function(x) x / sum(x) )
EX_ps_clean.rar.glom.stacked = tax_glom(EX_ps_clean.rar.stacked, "Phylum")
plot_bar(EX_ps_clean.rar.glom.stacked, fill="Phylum") 

library(dplyr)
library(microbiome)
# grab the taxonomy to clean up first
tax.df <- data.frame(tax_table(EX_ps_clean.rar), stringsAsFactors=FALSE)

## if the Genus is empty, replace with the Family
tax.df$Genus[is.na(tax.df$Genus)] <- tax.df$Family[is.na(tax.df$Genus)]

# put the new taxa back in
tax_table(EX_ps_clean.rar) <- as.matrix(tax.df)

# check it worked
as.data.frame(tax_table(EX_ps_clean.rar))$Genus

#calculate compositional version of the data (relative abundances)
core_biomeW.rel <- microbiome::transform(EX_ps_clean.rar, "compositional") #CHANGE ME to your phyloseq object

#This returns the taxa that exceed the given prevalence and minimum abundance detection thresholds. Set your preferred thresholds.
core_taxa_standardW <- core_members(core_biomeW.rel, detection = 1/10000, prevalence = 60/100) #CHANGE ME

#A full phyloseq object of the core microbiota at those limits is obtained as follows
phylo.coreW <- core(core_biomeW.rel, detection = 1/100000, prevalence = .2)

###retrieving the associated taxa names from the phyloseq object and add it to what you just made.
core.taxaW <- taxa(phylo.coreW)
class(core.taxaW)

# get the taxonomy data and assign it to what you just made.
tax.matW <- tax_table(phylo.coreW)
tax.dfW <- as.data.frame(tax.matW)

# add the SVs to last column of what you just made.
tax.dfW$SV <- rownames(tax.dfW)
summary(tax.dfW)

## 30 taxa shared

# select taxonomy of only those OTUs that are core members based on the thresholds that were used.
core.taxa.classW <- dplyr::filter(tax.dfW, rownames(tax.dfW) %in% core.taxaW)
knitr::kable(head(core.taxa.classW))

# save the list so you can access later, can report just the list
write.csv(core.taxa.classW, "core.taxa.example.csv")


# graph the abundance of those shared taxa, here are some example: https://microbiome.github.io/tutorials/Core.html

#aggregate the genera so we don't get a lot of lines separating all the SVs
plot.gen <- aggregate_taxa(phylo.coreW, "Genus")
plot.gen <- microbiome::transform(plot.gen, "compositional")

# load libraries as needed
library(ggplot2)
library(RColorBrewer)
library(viridis)

prevalences <- seq(.05, 1, .05)
detections <- round(10^seq(log10(1e-4), log10(.2), length = 10), 3)

plot_core(plot.gen, 
          plot.type = "heatmap", 
          prevalences = prevalences, 
          detections = detections, min.prevalence = .1) + #CHANGE min prevalence
  xlab("Detection Threshold (Relative Abundance (%))") + ylab("Bacterial SVs") +
  theme_minimal() + scale_fill_viridis()

library(phyloseq) 
library(dplyr) 
library(tidyr)  
library(ggplot2) 

# transform to relative abundance
relabun.ps <- transform_sample_counts(EX_ps_clean.rar,function(x) x / sum(x)) 

# glom ASVs by genus
ps_genus <- tax_glom(relabun.ps, taxrank = "Genus", NArm = FALSE) 

#subset by the genus of choice. Taxon is not present if you get error: "length of 'dimnames' [1] not equal to array extent" 
ps_genusP <- subset_taxa(ps_genus, Genus %in% c("Streptococcus"))                      

# melt the data into a different configuration
genus.df <- psmelt(ps_genusP) 

# grab that abundance data. CHANGE ME where it says Factor1, Factor2
MySummary <- genus.df %>% group_by(Tank_setup, Clean_or_dirty, Genus) %>% summarize(mean_abund = mean(Abundance, na.rm=TRUE)) 

#check it
head(MySummary)

# graph it
ggplot(data=MySummary, aes(x = Clean_or_dirty, y = mean_abund)) +
  geom_point(aes(color=Tank_setup)) + ylab("Mean relative abundance of reads")


library(DESeq2)
packageVersion("DESeq2")

# OPTIONAL if you need to subset 
EX_ps_clean_subsetA = subset_samples(EX_ps_clean, FACTOR == "A") #CHANGE ME 
EX_ps_clean_subsetA <- prune_samples(sample_sums(EX_ps_clean_subsetA) > 0, EX_ps_clean_subsetA)
EX_ps_clean_subsetA <- prune_taxa(taxa_sums(EX_ps_clean_subsetA) > 0, EX_ps_clean_subsetA)

EX_ps_clean_noNA <- subset_samples(EX_ps_clean, Tank_setup != "n/a")
sample_data(EX_ps_clean_noNA)$Tank_setup <- as.factor(sample_data(EX_ps_clean_noNA)$Tank_setup)
levels(sample_data(EX_ps_clean_noNA)$Tank_setup)

saveRDS(EX_ps_clean_noNA, 'EX_ps_clean_noNA.RDS')
EX_ps_clean_noNA <- readRDS('EX_ps_clean_noNA.RDS')

# grab phyloseq data for use in deseq
diagdds = phyloseq_to_deseq2(EX_ps_clean_noNA, ~ Tank_setup) #CHANGE ME to an factor with only 2 levels

saveRDS(diagdds, 'DESEQ_obj.RDS')

# calculate differential abundance
gm_mean =  function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
geoMeans = apply(counts(diagdds), 1, gm_mean)
diagdds = estimateSizeFactors(diagdds, geoMeans = geoMeans)
diagdds = DESeq(diagdds, fitType="local")

# calculate significance for those abundance calculations
res = results(diagdds)
res = res[order(res$padj, na.last=NA), ]
alpha = 0.01
sigtab = res[(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(EX_ps_clean_noNA)[rownames(sigtab), ], "matrix"))

head(sigtab)

dim(sigtab)
#13 SVs significantly different between the two tank conditions


# calculate log changes and set
sigtab = sigtab[sigtab[, "log2FoldChange"] < 0, ] #CHANGE ME. use this to select only positive (or negative) log changes
sigtab = sigtab[, c("baseMean", "log2FoldChange", "lfcSE", "padj", "Phylum", "Class", "Family", "Genus")] #CHANGE ME add Order or Species - if you have it

# Phylum order
x = tapply(sigtab$log2FoldChange, sigtab$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtab$Phylum = factor(as.character(sigtab$Phylum), levels=names(x))

# Genus order
x = tapply(sigtab$log2FoldChange, sigtab$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtab$Genus = factor(as.character(sigtab$Genus), levels=names(x))

sum(is.na(sigtab$Genus))
sigtab[is.na(sigtab$Genus), c("Family", "Genus")]

## if the Genus is empty, replace with the Family NOTE: this used to work but the syntax is broken for unknown reasons
sigtab.df$Genus[is.na(sigtab.df$Genus)] <- sigtab.df$Family[is.na(sigtab.df$Genus)]

# if the Genus is empty, replace with the Family NOTE: works as of Feb 2020
sigtab$Genus = ifelse(is.na(sigtab$Genus), paste(sigtab$Family),paste(sigtab$Genus));sigtab
sigtab$Genus <- ifelse(is.na(sigtab$Genus), sigtab$Family, sigtab$Genus)


## graph differential abundance
ggplot(sigtab, aes(y=Genus, x=log2FoldChange, color=Phylum)) + #play with aesthetics to make graph informative
  geom_vline(xintercept = 0.0, color = "gray", size = 0.5) +
  geom_point(aes(size=baseMean)) + #scale size by mean relative abundance
  theme(axis.text.x = element_text(hjust = 0, vjust=0.5, size=10), axis.text.y = element_text(size=10)) + xlab("log2 Fold Change") + labs(size = "Mean Sequence Abundance") + theme_minimal()


saveRDS(sigtab, 'Sigtab.RDS')
sigtab <- readRDS('Sigtab.RDS')


install.packages("devtools")
library(devtools)
install_github("vqv/ggbiplot") #takes a few minutes
library(ggbiplot)

# calculate the components
EX_clean.rar.pca <- prcomp(otu_table(EX_ps_clean.rar), center = TRUE, scale = TRUE) #CHANGE ME

# take a look
summary(EX_clean.rar.pca) #CHANGE ME

# graph it. add ggplot2 text to make it pretty
ggbiplot(EX_clean.rar.pca)
p <- ggbiplot(EX_clean.rar.pca,
              obs.scale = 1,
              var.scale = 1,
              groups = sample_data(EX_ps_clean.rar)$Clean_or_dirty, # CHANGE ME if using a different factor
              ellipse = TRUE,        # Draw confidence ellipses
              circle = TRUE)

p <- p +
  scale_color_manual(values = c("#1b9e77", "#d95f02", "blue")) +  # Custom colors for groups
  scale_fill_manual(values = c("#1b9e77", "#d95f02", "blue")) +
  theme_minimal(base_size = 14) +  # Clean background and larger text
  theme(
    panel.grid.major = element_line(color = "gray90"),
    panel.grid.minor = element_blank(),
    axis.text = element_text(color = "black"),
    axis.title = element_text(size = 14, face = "bold"),
    legend.position = "right",
    legend.title = element_text(face = "bold"),
    legend.background = element_blank()
  ) +
  labs(title = "PCA of Clean vs Dirty Samples",
       x = "PC1",
       y = "PC2",
       color = "Condition",
       fill = "Condition")

# Show it
print(p)

## PCoA in phyloseq -----------------

# use phyloseq to calculate ordination for use in the plot
EX_uJ_pcoa <- ordinate(EX_ps_clean.rar, #calculate similarities
                       method ="PCoA", #ordination type
                       "jaccard", binary = TRUE) #similarity type. Jaccard is binary, Bray can be binary (unweighted) or not (weighted) but is usually run as binary=FALSE.

# simple ordination
plot_ordination(EX_ps_clean.rar, EX_uJ_pcoa, type="samples", color="Tank_setup") #CHANGE ME but leave type as "samples"

# fancy ordination
library(ggplot2) # can add components on top of the phyloseq graph
install.packages("viridis")
library(viridis) #add a cool color palette

plot_ordination(EX_ps_clean.rar, EX_uJ_pcoa, type="samples") + #CHANGE ME
  geom_point(aes(size = as.factor(Clean_or_dirty), color=as.factor(Tank_setup), shape=as.factor(Clean_or_dirty))) + # resize the points based on a numeric factor, and make light/dark based on another factor
  theme_minimal()+
  scale_color_viridis(discrete = TRUE, option="A", begin = 0, end = 0.75) + #begin/end indicates where on the color scale to use, A refers to 'magma' color palette
  stat_ellipse(aes(group=as.factor(Clean_or_dirty), linetype=as.factor(Clean_or_dirty))) + #add circles around a particular grouping, and make the circle lines different
  labs(color = "Tank_setup", size = "Clean_or_dirty", shape = "Clean_vs_dirty", alpha = "Weight (kg)", linetype="Clean_vs_dirty") #rename the headers in the legend

colnames(sample_data(EX_ps_clean.rar))

## NMDS in phyloseq -----------------

# use phyloseq to calculate ordination for use in the plot
EX_uJ_nmds <- ordinate(EX_ps_clean.rar, #calculate similarities
                       method ="NMDS", #ordination type
                       "jaccard", binary = TRUE) #similarity type. Jaccard is binary, Bray can be binary (unweighted) or not (weighted) but is usually run as binary=FALSE.

# report the amount of stress at the end of the NMDS calculation, as it is typical to report that in your manuscript. >0.2 is not great, and >0.3 plot is meaningless
# Paste output here: Run 20 stress 0.1748307 #CHANGE ME

# simple ordination
plot_ordination(EX_ps_clean.rar, EX_uJ_nmds, type="samples", color="Tank_setup") #CHANGE ME but leave type as "samples"


# fancy ordination
library(ggplot2)
install.packages("viridis")
library(viridis) #add a cool color palette

plot_ordination(EX_ps_clean.rar, EX_uJ_nmds, type="samples", shape="Diet") + #CHANGE ME
  geom_point(aes(size = as.numeric(Weight_kg), color=as.factor(Week), alpha = as.factor(Week))) + # resize the points based on a numeric factor, and make light/dark based on another factor
  scale_color_viridis(discrete = TRUE, option="A", begin = 0, end = 0.75) + #begin/end indicates where on the color scale to use, A refers to 'magma' color palette
  stat_ellipse(aes(group=Diet, linetype=Diet)) + #add circles around a particular grouping, and make the circle lines different
  labs(color = "Weeks 0 to 2", size = "Weight in kg", shape = "Diet", linetype="Diet", alpha= "Week") #rename the headers in the legend


### ----- permanova stats for ordinations --------

library(vegan)

# Compute Bray-Curtis distance
otu_mat <- as(otu_table(EX_ps_clean.rar), "matrix")
if(taxa_are_rows(EX_ps_clean.rar)) otu_mat <- t(otu_mat)

dist_mat <- vegdist(otu_mat, method="bray")

adonis2(dist_mat ~ Tank_setup * Clean_or_dirty, data = as(sample_data(EX_ps_clean.rar), "data.frame"),
        permutations = 9999)

#Permutation test for adonis under reduced model
#Permutation: free
#Number of permutations: 9999
#adonis2(formula = dist_mat ~ Tank_setup * Clean_or_dirty, data = as(sample_data(EX_ps_clean.rar), "data.frame"), permutations = 9999)
#Df SumOfSqs      R2     F Pr(>F)  
#Model     5   2.6086 0.14484 1.084  0.029 *
#  Residual 32  15.4012 0.85516               
#Total    37  18.0098 1.00000               
#---
#  Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

# interpretation: The F-model/F statistic is a measure of importance, R2 may not be relevant for factorial data


#### run betadispersion test to see how tight/loose the clusters are ----
betadisp_EX <- betadisper(EX_uJ_nmds, Tank_setup$Clean_or_dirty) # CHANGE me so variable and factor reflect your data

otu_mat <- as(otu_table(EX_ps_clean.rar), "matrix")
if(taxa_are_rows(EX_ps_clean.rar)) otu_mat <- t(otu_mat)

dist_mat <- vegdist(otu_mat, method = "bray")  # numeric 'dist' object
group_factor <- sample_data(EX_ps_clean.rar)$Clean_or_dirty
betadisp_EX <- betadisper(dist_mat, group_factor)



betadisp_EX
# paste output here

# Example Call: betadisper(d = EX_BC, group = sampledf$Diet)

#No. of Positive Eigenvalues: 37
#No. of Negative Eigenvalues: 0

#Average distance to median:
#  after_cleaned_filled_with_water                after_drain_down                             n/a 
#0.6543                          0.6754                          0.6621 

#Eigenvalues for PCoA axes:
#  (Showing 8 of 37 eigenvalues)
#PCoA1  PCoA2  PCoA3  PCoA4  PCoA5  PCoA6  PCoA7  PCoA8 
#1.1713 0.8215 0.7962 0.6818 0.6298 0.5874 0.5726 0.5606 


# run ANOVA to see if clusters overlap or not
anova(betadisp_EX)
#Response: Distances
#Df    Sum Sq   Mean Sq F value  Pr(>F)  
#Groups     2 0.0030677 0.0015338  4.3048 0.02131 *
#  Residuals 35 0.0124706 0.0003563                  
#---
#  Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1


permutest(betadisp_EX)
#Permutation test for homogeneity of multivariate dispersions
#Permutation: free
#Number of permutations: 999

#Response: Distances
#Df    Sum Sq   Mean Sq      F N.Perm Pr(>F)  
#Groups     2 0.0030677 0.0015338 4.3048    999  0.022 *
#  Residuals 35 0.0124706 0.0003563                       
#---
#  Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1


#  correct for multiple comparisons
TukeyHSD(betadisp_EX)
#Tukey multiple comparisons of means
#95% family-wise confidence level

#Fit: aov(formula = distances ~ group, data = df)

#$group
#diff         lwr         upr     p adj
#after_drain_down-after_cleaned_filled_with_water  0.021086236  0.00254427 0.039628202 0.0228821
#n/a-after_cleaned_filled_with_water               0.007847021 -0.01392941 0.029623456 0.6551521
#n/a-after_drain_down                             -0.013239215 -0.03178118 0.005302751 0.2022815

install.packages('vegan3d')
library(vegan3d)
install.packages("scatterplot3d", dependencies = TRUE)
library(scatterplot3d)
# library(vegan)


# basic version of the plot, using the object made in the Betadispersion section above
ordiplot3d(betadisp_EX)


## CCA in phyloseq -------------------------------------------
# correspondence analysis, or optionally, constrained correspondence analysis (a.k.a. canonical correspondence analysis), usually used with bell-curve or unimodal relationships hypothesis-based testing does not need normally distributed data, but if you have a lot of outliers you may want to consider adding a data transformation.

# first create a distance ordination. BE SURE THERE ARE NO NAs in your factorial data!! If you need to alter something, see the "beta diversity components" section.

# reload as needed
library(phyloseq)
library(vegan)
library(ggplot2)

bray_not_na <- phyloseq::distance(physeq = EX_ps_clean.rar, method = "bray") #CHANGE ME

cca_ord <- ordinate(
  physeq = EX_ps_clean.rar, #CHANGE ME
  method = "CCA",
  distance = bray_not_na,
  formula = ~ Tank_setup * Clean_or_dirty # CHANGE ME add factors
)
colnames(sample_data(EX_ps_clean.rar))

cca_ord
#Call: cca(formula = OTU ~ Tank_setup * Clean_or_dirty, data = data)

#Inertia Proportion Rank
#Total         30.2996     1.0000     
#Constrained    4.2133     0.1391    5
#Unconstrained 26.0863     0.8609   32

#Inertia is scaled Chi-square

#-- NOTE:
#  Some constraints or conditions were aliased because they were redundant. This can happen if terms are
#constant or linearly dependent (collinear): 'Tank_setupn/a:Clean_or_dirtyafter_drain_down',
#'Tank_setupn/a:Clean_or_dirtyn/a', 'Tank_setupstatic:Clean_or_dirtyn/a'

#Eigenvalues for constrained axes:
#  CCA1   CCA2   CCA3   CCA4   CCA5 
#0.9844 0.9038 0.8209 0.7942 0.7100 

#Eigenvalues for unconstrained axes:
#  CA1    CA2    CA3    CA4    CA5    CA6    CA7    CA8 
#0.9756 0.9717 0.9556 0.9451 0.9346 0.9267 0.9148 0.9090 
#(Showing 8 of 32 unconstrained eigenvalues)


# NOTE: do you have CCA eigenvalues for each of the factors you put in? CCA1 = x axis, and CCA2 = y axis, so if you don't have both of those it will use CA which is the unconstrained axis instead.



# anova of whole model
anova(cca_ord, permu=1000)
#Permutation test for cca under reduced model
#Permutation: free
#Number of permutations: 999

#Model: cca(formula = OTU ~ Tank_setup * Clean_or_dirty, data = data)
#Df ChiSquare      F Pr(>F)
#Model     5    4.2133 1.0337  0.108
#Residual 32   26.0863  


# anova of the factors (terms) you specified
anova(cca_ord, by="terms", permu=1000)
#Permutation test for cca under reduced model
#Terms added sequentially (first to last)
#Permutation: free
#Number of permutations: 999

#Model: cca(formula = OTU ~ Tank_setup * Clean_or_dirty, data = data)
#Df ChiSquare      F Pr(>F)    
#Tank_setup                 2    1.7040 1.0451  0.165    
#Clean_or_dirty             2    1.7791 1.0912  0.001 ***
#  Tank_setup:Clean_or_dirty  1    0.7302 0.8957  0.955    
#Residual                  32   26.0863                  
#---
#  Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

# cca plot with FACTORA in the model
cca_plot <- plot_ordination(
  physeq = EX_ps_clean.rar, #CHANGE ME
  ordination = cca_ord, 
  color = "Tank_setup", #CHANGE ME
  axes = c(1,2)) + 
  theme_minimal() +
  aes(shape = as.factor(Clean_or_dirty)) + #CHANGE ME
  geom_point(aes(colour = Tank_setup), size = 4) + #CHANGE ME
  geom_point(colour = "grey90", size = 1.5) 


# Now add the environmental variables as arrows
arrowmat <- vegan::scores(cca_ord, display = "bp")

# Add labels, make a data.frame
arrowdf <- data.frame(labels = rownames(arrowmat), arrowmat)

# Define the arrow aesthetic mapping
arrow_map <- aes(xend = CCA1, 
                 yend = CCA2, 
                 x = 0, 
                 y = 0, 
                 shape = NULL, 
                 color = NULL, 
                 label = labels)

label_map <- aes(x = 1.3 * CCA1, 
                 y = 1.3 * CCA2, 
                 shape = NULL, 
                 color = NULL, 
                 label = labels)

arrowhead = arrow(length = unit(0.02, "npc"))

# Make a new graphic
cca_plot + 
  geom_segment(
    mapping = arrow_map, 
    linewidth = .5, 
    data = arrowdf, 
    color = "gray", 
    arrow = arrowhead
  ) + 
  geom_text(
    mapping = label_map, 
    size = 4,  
    data = arrowdf, 
    show.legend = FALSE
  ) 





## RDA in phyloseq -------------------------------------------
# Redundancy analysis for linear relationships

# first create a distance ordination. BE SURE THERE ARE NO NAs in your factorial data!! If you need to alter something, see the "beta diversity components" section.
bray_not_na <- phyloseq::distance(physeq = EX_ps_clean.rar, method = "bray")

rda_ord <- ordinate(
  physeq = EX_ps_clean.rar, #CHANGEME
  method = "RDA",
  distance = bray_not_na,
  formula = ~ Tank_setup * Clean_or_dirty # CHANGE ME add factors
   #CHANGE ME use for repeated measures
)

rda_ord
# paste the output here
#Inertia Proportion Rank
#Total         1.290e+06  1.000e+00     
#Constrained   2.098e+05  1.627e-01    5
#Unconstrained 1.080e+06  8.373e-01   32

#Inertia is variance

#-- NOTE:
#  Some constraints or conditions were aliased because they were redundant. This can happen
#if terms are constant or linearly dependent (collinear):
#  'Tank_setupn/a:Clean_or_dirtyafter_drain_down', 'Tank_setupn/a:Clean_or_dirtyn/a',
#'Tank_setupstatic:Clean_or_dirtyn/a'

#Eigenvalues for constrained axes:
#  RDA1   RDA2   RDA3   RDA4   RDA5 
#112567  28350  26634  24677  17611 

#Eigenvalues for unconstrained axes:
#  PC1    PC2    PC3    PC4    PC5    PC6    PC7    PC8 
#179059 122895  61627  60296  50673  46533  45199  41781 
#(Showing 8 of 32 unconstrained eigenvalue
# NOTE: do you have RDA eigenvalues for each of the factors you put in? RDA1 = x axis, and RDA2 = y axis, so if you don't have both of those it will use PC which is the unconstrained axis instead.



# anova of whole model
anova(rda_ord, permu=1000)
#Permutation test for rda under reduced model
#Permutation: free
#Number of permutations: 999

#Model: rda(formula = OTU ~ Tank_setup * Clean_or_dirty, data = data)
#Df Variance      F Pr(>F)
#Model     5   209838 1.2432  0.119
#Residual 32  1080271   

# anova of the factors (terms) you specified
anova(rda_ord, by="terms", permu=1000)
#Permutation test for rda under reduced model
#Terms added sequentially (first to last)
#Permutation: free
#Number of permutations: 999

#Model: rda(formula = OTU ~ Tank_setup * Clean_or_dirty, data = data)
#Df Variance      F Pr(>F)   
#Tank_setup                 2   137516 2.0368  0.002 **
#  Clean_or_dirty             2    50836 0.7529  0.945   
#Tank_setup:Clean_or_dirty  1    21487 0.6365  0.882   
#Residual                  32  1080271                 
#---
#  Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1             


# rda plot with FACTORA in the model
rda_plot <- plot_ordination(
  physeq = EX_ps_clean.rar, #CHANGE ME
  ordination = rda_ord, 
  color = "Tank_setup", #CHANGE ME
  axes = c(1,2)) + 
  theme_minimal() +
  aes(shape = as.factor(Clean_or_dirty)) + #CHANGE ME
  geom_point(aes(colour = Tank_setup), size = 4) + #CHANGE ME
  geom_point(colour = "grey90", size = 1.5) 



# Now add the environmental variables as arrows
arrowmat <- vegan::scores(rda_ord, display = "bp")

# Add labels, make a data.frame
arrowdf <- data.frame(labels = rownames(arrowmat), arrowmat)

# Define the arrow aesthetic mapping
arrow_map <- aes(xend = RDA1, 
                 yend = RDA2, 
                 x = 0, 
                 y = 0, 
                 shape = NULL, 
                 color = NULL, 
                 label = labels)

label_map <- aes(x = 1.3 * RDA1, 
                 #y = 1.3 * RDA2, 
                 shape = NULL, 
                 color = NULL, 
                 label = labels)

arrowhead = arrow(length = unit(0.02, "npc"))

# Make a new graphic
rda_plot + 
  geom_segment(
    mapping = arrow_map, 
    linewidth = .5, 
    data = arrowdf, 
    color = "gray", 
    arrow = arrowhead
  ) #+ 
  #geom_text(
    #mapping = label_map, 
    #size = 4,  
    #data = arrowdf, 
    #show.legend = FALSE
  #) 



## dbRDA in phyloseq -------------------------------------------
# Constrained Analysis of Principal Coordinates or distance-based RDA, for non-linear relationships

# first create a distance ordination. BE SURE THERE ARE NO NAs in your factorial data!! If you need to alter something, see the "beta diversity components" section.
bray_not_na <- phyloseq::distance(physeq = EX_ps_clean.rar, method = "bray")

cap_ord <- ordinate(
  physeq = EX_ps_clean.rar, #CHANGEME
  method = "CAP", #borrows capscale from the vegan package
  distance = bray_not_na,
  formula = ~ Tank_setup * Clean_or_dirty # CHANGE ME add factors
)

cap_ord
#Inertia Proportion Rank
#Total         18.0098     1.0000     
#Constrained    2.6086     0.1448    5
#Unconstrained 15.4012     0.8552   32

#Inertia is squared Bray distance

#-- NOTE:
#  Some constraints or conditions were aliased because they were redundant. This can happen
#if terms are constant or linearly dependent (collinear):
#  'Tank_setupn/a:Clean_or_dirtyafter_drain_down', 'Tank_setupn/a:Clean_or_dirtyn/a',
#'Tank_setupstatic:Clean_or_dirtyn/a'

#Eigenvalues for constrained axes:
#  CAP1   CAP2   CAP3   CAP4   CAP5 
#0.6409 0.5889 0.5117 0.4723 0.3947 

#Eigenvalues for unconstrained axes:
#  MDS1   MDS2   MDS3   MDS4   MDS5   MDS6   MDS7   MDS8 
#1.0925 0.7775 0.7048 0.6010 0.5819 0.5583 0.5502 0.5268 
#(Showing 8 of 32 unconstrained eigenvalues)

# NOTE: do you have CAP eigenvalues for each of the factors you put in? CAP1 = x axis, and CAP2 = y axis, so if you don't have both of those it will use MDS which is the unconstrained axis instead.

# anova of whole model
anova(cap_ord, permu=1000)
#Permutation test for capscale under reduced model
#Permutation: free
#Number of permutations: 999

#Model: capscale(formula = distance ~ Tank_setup * Clean_or_dirty, data = data)
#Df SumOfSqs     F Pr(>F)  
#Model     5   2.6086 1.084   0.03 *
#  Residual 32  15.4012               
#---
#  Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# anova of the factors (terms) you specified
anova(cap_ord, by="terms", permu=1000)
#Permutation test for capscale under reduced model
#Terms added sequentially (first to last)
#Permutation: free
#Number of permutations: 999

#Model: capscale(formula = distance ~ Tank_setup * Clean_or_dirty, data = data)
#Df SumOfSqs      F Pr(>F)  
#Tank_setup                 2   1.1000 1.1427  0.020 *
#  Clean_or_dirty             2   1.0901 1.1324  0.034 *
#  Tank_setup:Clean_or_dirty  1   0.4186 0.8697  0.960  
#Residual                  32  15.4012                
#---
#  Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1            


# CAP plot with FACTORA in the model
cap_plot <- plot_ordination(
  physeq = EX_ps_clean.rar, #CHANGE ME
  ordination = cap_ord, 
  axes = c(1,2)) + 
  theme_minimal() +
  aes(shape = as.factor(Clean_or_dirty)) + #CHANGE ME
  geom_point(aes(colour = Tank_setup), size = 4) + #CHANGE ME
  geom_point(colour = "grey90", size = 1.5) 




# Now add the environmental variables as arrows
arrowmat <- vegan::scores(cap_ord, display = "bp")

# Add labels, make a data.frame
arrowdf <- data.frame(labels = rownames(arrowmat), arrowmat)

# Define the arrow aesthetic mapping
arrow_map <- aes(xend = CAP1, 
                 yend = CAP2, 
                 x = 0, 
                 y = 0, 
                 shape = NULL, 
                 color = NULL, 
                 label = labels)

label_map <- aes(x = 1.3 * CAP1, 
                 y = 1.3 * CAP2, 
                 shape = NULL, 
                 color = NULL, 
                 label = labels)

arrowhead = arrow(length = unit(0.02, "npc"))

# Make a new graphic
cap_plot + 
  geom_segment(
    mapping = arrow_map, 
    linewidth = .5, 
    data = arrowdf, 
    color = "gray", 
    arrow = arrowhead
  ) #+ 
  #geom_text(
    #mapping = label_map, 
    #size = 4,  
    #data = arrowdf, 
    #show.legend = FALSE
  #) 



## variance partitioning --------------
# How much did any one factor contribute to the sample clustering?
# http://www.hiercourse.com/docs/variationPartitioning.html 

library(vegan)

# create a dataframe from your SV table
EX.df <- as.data.frame(otu_table(EX_ps_clean.rar))
taxa_are_rows(EX_ps_clean.rar)
head(rownames(EX.df))
head(rownames(env.df))
str(env.df)


# extract your sample data from the phyloseq object
env.df <- as(sample_data(EX_ps_clean.rar), "data.frame")
# Remove columns with only 1 unique value
env.df.clean <- env.df[, sapply(env.df, function(x) length(unique(x)) > 1)]
env.df.clean
sapply(env.df.clean, function(x) length(unique(x)))


# calculate Principal Coordinates Of Neighborhood Matrix, diversity distance data transformed into rectangular format
EX.pcnm <- pcnm(dist(bray_not_na))

# environmental variables as predictors of community similarity
cap.env <- capscale(EX.df ~ ., data=env.df.clean, distance='bray')
cap.env <- capscale(
  EX.df ~ Tank_setup + Clean_or_dirty + Sample_type + Days_since_full_moon,
  data = env.df.clean,
  distance = "bray"
)

# calculate CCA 
cap.pcnm <- capscale(EX.df ~ ., data=as.data.frame(scores(EX.pcnm)), distance='bray')


# make a model using SV ordination and sample data
mod0.env <- capscale(EX.df ~ 1, data=env.df.clean, distance='bray') # add + Condition(SAMPLE_ID) for repeated measures

# make a model using SV ordination and the component scores
mod0.pcnm <- capscale(EX.df ~ 1, data=as.data.frame(scores(EX.pcnm)), distance='bray') # add + Condition(SAMPLE_ID) for repeated measures

# select variables in each predictor table
step.env <- ordistep(mod0.env, scope=formula(cap.env))
#Start: EX.df ~ 1 

#Df    AIC      F Pr(>F)   
#+ Days_since_full_moon 11 118.36 1.0964  0.005 **
#  + Sample_type           5 114.85 1.0922  0.015 * 
#  + Tank_setup            2 112.45 1.1383  0.020 * 
#  + Clean_or_dirty        2 112.50 1.1149  0.075 . 
#---
#  Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

#Step: EX.df ~ Days_since_full_moon 

#Df    AIC      F Pr(>F)   
#- Days_since_full_moon 11 110.84 1.0964  0.005 **
#  ---
#  Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

#Df    AIC      F Pr(>F)
#+ Sample_type     3 119.68 1.0045  0.565
#+ Tank_setup      2 119.37 0.9827  0.630
#+ Clean_or_dirty  2 119.54 0.9231  0.740


# check variance inflation factors, higher number = data are redundant to another factor, 1= data are unique.  If factors are redundant (conflated) to each other (for example they basically report the same thing, go back and remove one and mention they were redundant/conflated in your manuscript Methods)
vif.cca(step.env)
#Days_since_full_moon13 Days_since_full_moon14 Days_since_full_moon15 Days_since_full_moon16 
#1.168421               1.610526               1.168421               1.168421 
#Days_since_full_moon18 Days_since_full_moon24 Days_since_full_moon26  Days_since_full_moon4 
#1.168421               1.326316               1.736842               1.473684 
#Days_since_full_moon6  Days_since_full_moon7  Days_since_full_moon9 
#1.610526               1.736842               1.852632 

step.pcnm <- ordistep(mod0.pcnm, scope=formula(cap.pcnm))

step.env$anova
#Df    AIC      F Pr(>F)   
#+ Days_since_full_moon 11 118.36 1.0964  0.005 **
#  ---
#  Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

step.pcnm$anova
#Df    AIC      F Pr(>F)   
#+ PCNM1   1 110.68 2.1036  0.005 **
#  + PCNM2   1 110.61 1.9590  0.005 **
#  + PCNM3   1 110.92 1.5525  0.005 **
#  + PCNM4   1 111.35 1.3927  0.005 **
#  + PCNM9   1 111.75 1.3770  0.005 **
#  + PCNM5   1 112.13 1.3434  0.005 **
#  + PCNM6   1 112.50 1.3174  0.005 **
#  + PCNM10  1 112.83 1.3069  0.005 **
#  + PCNM8   1 113.10 1.3032  0.005 **
#  + PCNM7   1 113.35 1.2715  0.025 * 
#  ---
#  Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
EX.pcnm.sub <- scores(EX.pcnm, 
                      choices=c(1:10)) #CHANGE ME to include the significant pcnm numbers from previous command. If you want them all, write 1:16. to select only a few of them, write something like 1,3,6,7:9. 

# partition variation among four predictor tables:
EX.var <- varpart(EX.df, 
                  ~ Tank_setup, #CHANGE ME
                  ~ Clean_or_dirty, #CHANGE ME
                  EX.pcnm.sub, data=env.df)

#plot 
par(mfrow=c(1, 2))
showvarparts(4)
plot(EX.var)



EX.var
#Partition of variance in RDA 

#Call: varpart(Y = EX.df, X = ~Tank_setup, ~Clean_or_dirty, EX.pcnm.sub, data =
#                env.df)

#Explanatory tables:
#  X1:  ~Tank_setup
#X2:  ~Clean_or_dirty
#X3:  EX.pcnm.sub 

#No. of explanatory tables: 3 
#Total variation (SS): 47734029 
#Variance: 1290109 
#No. of observations: 38 

#Partition table:
#  Df R.square Adj.R.square Testable
#[a+d+f+g] = X1         2  0.10659      0.05554     TRUE
#[b+d+e+g] = X2         2  0.06548      0.01208     TRUE
#[c+e+f+g] = X3        10  0.34263      0.09916     TRUE
#[a+b+d+e+f+g] = X1+X2  4  0.14600      0.04248     TRUE
#[a+c+d+e+f+g] = X1+X3 12  0.40992      0.12669     TRUE
#[b+c+d+e+f+g] = X2+X3 12  0.40877      0.12498     TRUE
#[a+b+c+d+e+f+g] = All 14  0.44675      0.10999     TRUE
#Individual fractions                                   
#[a] = X1 | X2+X3       2              -0.01499     TRUE
#[b] = X2 | X1+X3       2              -0.01670     TRUE
#[c] = X3 | X1+X2      10               0.06751     TRUE
#[d]                    0               0.04252    FALSE
#[e]                    0               0.00364    FALSE
#[f]                    0               0.04540    FALSE
#[g]                    0              -0.01739    FALSE
#[h] = Residuals                        0.89001    FALSE
#Controlling 1 table X                                  
#[a+d] = X1 | X3        2               0.02753     TRUE
#[a+f] = X1 | X2        2               0.03041     TRUE
#[b+d] = X2 | X3        2               0.02582     TRUE
#[b+e] = X2 | X1        2              -0.01306     TRUE
#[c+e] = X3 | X1       10               0.07115     TRUE
#[c+f] = X3 | X2       10               0.11291     TRUE
#---
#  Use function 'rda' to test significance of fractions of interest

anova(rda(EX.df  ~ Tank_setup + Condition(EX.pcnm.sub), data=env.df)) # add + Condition(env.df$SAMPLE_ID) for repeated measures
#Permutation test for rda under reduced model
#Permutation: free
#Number of permutations: 999

#Model: rda(formula = EX.df ~ Tank_setup + Condition(EX.pcnm.sub), data = env.df)
#Df Variance      F Pr(>F)  
#Model     2    86819 1.4256  0.071 .
#Residual 25   761261                
#---
#  Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1











## Mantel test -------

bray_not_na <- phyloseq::distance(physeq = EX_ps_clean.rar, method = "bray")

# make a data frame that can be used in vegan from the sample_data in the phyloseq object
meta.df <- data.frame(sample_data(EX_ps_clean.rar))

# check the columns in your meta dataframe, need numeric data for mantel test
str(meta.df)

# make numeric any columns that are a number but as listed as "chr", for example 
meta.df$DNA_extraction_batch <- as.numeric(meta.df$DNA_extraction_batch) # CHANGE ME to any columns that should be numeric

# replace any binary factors with numbers, for example
meta.df["Diet"][meta.df["Diet"] == "LooseAlfalfa"] <- "1"
meta.df["Diet"][meta.df["Diet"] == "PelletedAlfalfa"] <- "2"
meta.df$Diet <- as.numeric(meta.df$Diet)

#remove any columns that you can't make numeric
meta.df <- subset(meta.df, select = -c(Sample_type, Sheep_ID, Treatment))

# Make a distance matrix for your meta data
bray.dist.meta <- vegdist(meta.df, method = "bray") 

# Mantel test to see if SV and meta distance matrixes correlate
bray.mantel <- mantel(bray_not_na, bray.dist.meta, permutations = 999)

# look at the output
bray.mantel

# Mantel statistic based on Pearson's product-moment correlation 
# 
# Call:
# mantel(xdis = bray_not_na, ydis = bray.dist.meta, permutations = 999) 
# 
# Mantel statistic r: 0.06023 
#       Significance: 0.224 
# 
# Upper quantiles of permutations (null model):
#   90%   95% 97.5%   99% 
# 0.108 0.143 0.166 0.189 
# Permutation: free
# Number of permutations: 999



