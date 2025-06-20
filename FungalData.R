#install.packages('phyloseq')
library(phyloseq)
library(dplyr)
#install.packages('stringr')
library(stringr)
#install.packages("writexl")
library(writexl)
library(ggplot2)

#install.packages("BiocManager")
#BiocManager::install("phyloseq")
#BiocManager::install("decontam")
library(decontam)
library(phyloseq)
#install.packages("vegan")
library(vegan)
#install.packages("pheatmap")
library(pheatmap)
library(tidyverse)
#install.packages("Matrix")
library(Matrix)


getwd()
setwd("/Users/caoyang/Desktop/Tetel Lab")

data <- readRDS("/Users/caoyang/Desktop/Tetel Lab/Walther-Antonio_Project_022_ITS2.rds") #the data is reading an email forwarded by Alice to Helena 

otu_table(data)
?otu_table
tax_table(data)

#View(otuF)

taxaF <- tax_table(data)
otuF <- otu_table(data)
otuF

otuF@.Data[1,]
dim(otuF@.Data)
dim(taxaF@.Data)

#sample_data(data)

#rows in both datasets are genetic sequences, the otu dataset says 1/0 whether it has that genetic sequence
#rownames in taxa table : specific sequences; whether a sample have this sequence: indicated by 0/1

length(table(taxaF@.Data[,"Species"]))
#so each row is not species, since the number of species is much shorter than #rows, maybe it's substring?
#there should be somewhere that associate the id numbers in otu with the specific participant 

# colnamesOtu <- colnames(otuF@.Data)
# colnamesOtu
# 
# rowOtu <- rownames(otuF@.Data)
# rowOtu

dim(otuF@.Data)
#dim(taxa_matrix)
str(taxaF)
taxa_matrix <- taxaF@.Data
taxa_matrix[1,]
rownames <- rownames(taxaF@.Data)
rownames
colnames <- colnames(taxaF@.Data)
colnames
taxa_matrix[, "SH"]
sequence <- taxaF@.Data[, "sequence"]
# sequence[1]
# sequence[2]
# sequence[5]

species <- taxaF@.Data[, "Species"]
dim(taxaF)
# species[1]
# species[2]
# species[3]
# species[4]
# species[5]

table(taxaF@.Data[, "Species"])
names(which.max(table(taxaF@.Data[, "Species"])))


species <- taxaF@.Data[, "Species"]
nonemptyspecies <- species[species != ""]# Remove empty strings
emptystrings <- species[species == ""]
length(emptystrings) #8560/15293 are empty
most_frequent_species <- names(which.max(table(nonemptyspecies)))
most_frequent_species


#spot check random rows
set.seed(123)
randomrows <- sample(1:nrow(taxaF@.Data), 25)
randomTaxa <- taxaF@.Data[randomrows, ]
randomTaxa[, c("sequence", "Species", "Species_exact")] #shows sequence and species

sequence_a0487137f959e09afcf2c8bbed1d3f0a <- taxaF@.Data["a0487137f959e09afcf2c8bbed1d3f0a", "sequence"] #species = "Metazoa_sp" 
sequence_a0487137f959e09afcf2c8bbed1d3f0a
sequence_2d9a26aa49b44db02454205f87f7229b <- taxaF@.Data["2d9a26aa49b44db02454205f87f7229b", "sequence"] #species = "Metazoa_sp"
sequence_2d9a26aa49b44db02454205f87f7229b
sequence_4e6a63a218b15e079e2be7182f463528 <- taxaF@.Data["4e6a63a218b15e079e2be7182f463528", "sequence"] #species = "Malassezia_globosa"
sequence_4e6a63a218b15e079e2be7182f463528
sequence_8fc4a855921af39c16de659ab850b298 <- taxaF@.Data["8fc4a855921af39c16de659ab850b298", "sequence"] #species = "Fusarium_mangiferae"
sequence_8fc4a855921af39c16de659ab850b298
sequence_3418b23c3ff6957730f17e4fd7e1f11a <- taxaF@.Data["3418b23c3ff6957730f17e4fd7e1f11a", "sequence"] #species = "Cercozoa_sp" 
sequence_3418b23c3ff6957730f17e4fd7e1f11a
sequence_ac6b493205438564c3895a162ba06d1a <- taxaF@.Data["ac6b493205438564c3895a162ba06d1a", "sequence"] #species = "Metazoa_sp"
sequence_ac6b493205438564c3895a162ba06d1a
sequence_e69ecb46b691be29b77534c32028225c <- taxaF@.Data["e69ecb46b691be29b77534c32028225c", "sequence"] #species = "subulatus"
sequence_e69ecb46b691be29b77534c32028225c
sequence_df61fec922e8ebfc48f748002c71926c <- taxaF@.Data["df61fec922e8ebfc48f748002c71926c", "sequence"] #species = "Metazoa_sp"
sequence_df61fec922e8ebfc48f748002c71926c
sequence_866e0176108bba26deba664c1a0ab8ff <- taxaF@.Data["866e0176108bba26deba664c1a0ab8ff", "sequence"] #species = "Metazoa_sp"
sequence_866e0176108bba26deba664c1a0ab8ff
sequence_f9d7bb29bdd41360e3a98318a78a9a9a <- taxaF@.Data["f9d7bb29bdd41360e3a98318a78a9a9a", "sequence"] #species = "Brassica_carinata"
sequence_f9d7bb29bdd41360e3a98318a78a9a9a

#Species, Confidence, Participant ID
#For each person, what is the most common specie/Shannon Diversity Index???
#what does each row represents? There is 15000+ rows

length(table(rownames)) #are they repeating rownames? no



#trying Alice's code #checking contamination
ntaxa(data)
nsamples(data)
sample_names(data)[1:5]
rank_names(data)
sample_variables(data, errorIfNULL=FALSE) # no colnames
otu_table(data)[1:5, 1:5]
tax_table(data)[1:5, 1:4]
phy_tree(data, errorIfNULL=FALSE) # no phylo tree
taxa_names(data)[1:10]

taxa_fungal.data <- as.data.frame(tax_table(data))
dim(taxa_fungal.data)
summary(as.numeric(taxa_fungal.data$confidence))
table(as.numeric(taxa_fungal.data$confidence))

otu_table_fungal.data <- as.data.frame(otu_table(data))

### Data Preprocessing

otu_table_fungal.data <- as.data.frame(t(otu_table_fungal.data))

metadata_fungal <- otu_table_fungal.data %>% 
  mutate(SampleID=rownames(otu_table_fungal.data),
         is_blank=as.logical(ifelse(str_detect(SampleID, "BLANK"), "TRUE", "FALSE")),
         SampleID= sub("\\..*", "", SampleID)) %>%
  select(SampleID, is_blank)

#The above code processes the otu_table_fungal.data object to:
#Extract row names (assumed to be sample IDs) into a new column called SampleID.
#Create a logical column (is_blank) that indicates whether the sample ID contains the word "BLANK".
#Clean up the SampleID column by removing everything after the first dot.
#Retain only the SampleID and is_blank columns in the final output, which is stored in metadata_fungal.

# Convert back to phyloseq obj
otu_table_obj <- otu_table(otu_table_fungal.data, taxa_are_rows = FALSE)
sample_data_obj <- sample_data(metadata_fungal)
fungal_physeq <- phyloseq(otu_table_obj, sample_data_obj, tax_table(data))

# Identify contaminants based on prevalence
contam_prev <- isContaminant(fungal_physeq, method = "prevalence", neg = "is_blank")
contaminants <- contam_prev$contaminant
contaminants

#filter the contaminants
fungal_physeq_no_contam <- prune_taxa(!contaminants, fungal_physeq)
fungal_physeq_subset <- subset_taxa(fungal_physeq_no_contam, Kingdom == "Fungi" & !is.na(Phylum) & Phylum != "")

#checking the dimension
#rows are the number of samples, columns are taxa
dim(otu_table_obj)
dim(otu_table(fungal_physeq_subset))
dim(otu_table(fungal_physeq_no_contam))
otu_table(fungal_physeq_subset)[1,]

#otutable spreadsheet

#readme

class(fungal_physeq_subset)

#3.26
dominant_spec <- apply(t(otu_table(fungal_physeq_subset)), 2, function(x) {
  spec <- tax_table(fungal_physeq_subset)[which.max(x), "Species"]
  # ifelse(is.na(spec), "Unknown", spec)
})

head(dominant_spec)


# cat("OTU table dimensions:", dim(otu_table(fungal_physeq_no_contam)), "\n")
# cat("Tax table dimensions:", dim(tax_table(fungal_physeq_no_contam)), "\n")
# 
# 
# 
# otu_corrected <- otu_table(t(otu_table(fungal_physeq_no_contam)), taxa_are_rows = TRUE)
# corrected_physeq <- phyloseq(
#   otu_corrected,
#   tax_table(fungal_physeq_no_contam),
#   sample_data(fungal_physeq_no_contam) # If you have sample data
# )
# 
# dominant_spec <- apply(otu_table(fungal_physeq_no_contam), 1, function(x) {
#   spec <- tax_table(fungal_physeq_no_contam)[which.max(x), "Species"]
#   ifelse(is.na(spec), "Unknown", spec)
# })

#creating dominant species column
sample_data(fungal_physeq_subset)$DominantSpecies <- dominant_spec
#check if the column is there
sample_data(fungal_physeq_subset)$DominantSpecies
#check the class of the object we created
summary(fungal_physeq_subset)
#the sample data part of our object
View(sample_data(fungal_physeq_subset))
#sort the most frequently appeared species in the sample data
sort(table(sample_data(fungal_physeq_subset)$DominantSpecies))  # Candida_albicans  the most frequently appeared

#4.2
#turning the otu table back to data frame
physeqOTU<- as.data.frame(otu_table(fungal_physeq_subset))
physeqOTU["F1376",1:10] #the first ten rows of this sample column

ncol(physeqOTU)
dim(physeqOTU)

# ones_indices <- which(physeqOTU["F1376", ] == 1) #looking for the sequences tagged as 1 
# ones_indices
# sum(physeqOTU["F1376", ] == 0)
# any(physeqOTU["F1376", ] == 1)
nonzero <- which(physeqOTU["F1376", ] != 0) 
nonzero #the column numbers of sequence ID
physeqOTU["F1376", nonzero] #the counts of the nonzero species
sequenceIDnonzero <- colnames(physeqOTU)[nonzero]
#> colnames(physeqOTU)[nonzero]
#"074f81db997e702d17c85be0c46b03ad" "9589a4186e70b56c852377d7562d4789" "f8f00151e8dd21d0f490343afb43d854"
#"f74f09973e588cd61df7e11c6620c2fd" "8d50d5e4d31b744eb9647746878c017e"
taxaF@.Data["074f81db997e702d17c85be0c46b03ad", "Species"] #it's Candida_albicans

#4.3
#install.packages('readxl')
library(readxl)
getwd()
sampleLabel<-read_excel("/Users/caoyang/Desktop/Tetel Lab/cleaned_samplesv2.xlsx")
dim(sampleLabel)
colnames(sampleLabel)
head(sampleLabel)
sampleLabel$sampleID <- sampleLabel$qr
sampleLabel <- sampleLabel[, c("sampleID", "biome_id", "sampleType")]
head(sampleLabel)

temp <- t(otu_table(fungal_physeq_subset))
dim(temp)
temp[1:10, 1:10]
temp2<-temp[, "F1376"]
head(temp2)
table(temp2)

which.max(temp2)
temp2[2024]
rownames(temp2[1:10])
rownames(temp2[temp2>0])
taxaF@.Data["f8f00151e8dd21d0f490343afb43d854", "Species"] 

# temp["Species", temp2>0] 


#merging sample label and sample data
?merge()



sum(duplicated(sample_names(fungal_physeq_subset)))  
sum(duplicated(sampleLabel$sampleID))  

head(sample_names(fungal_physeq_subset)) 

head(rownames(sampleLabel))  
head(sampleLabel$sampleID)

#removing the ones that are blank
fungal_physeq_subset <- subset_samples(
  fungal_physeq_subset,
  !str_detect(sample_names(fungal_physeq_subset), "^BLANK")
)
# Trim .ITS2
sample_names(fungal_physeq_subset) <- str_remove(sample_names(fungal_physeq_subset), "\\.ITS2$")

head(sample_data(fungal_physeq_subset))

temp <-sample_data(fungal_physeq_subset) 
sampleDataforMerge <- data.frame(temp$SampleID, temp$is_blank, temp$DominantSpecies)
head(sampleDataforMerge)
colnames(sampleDataforMerge) <- (c("SampleID", "is_blank", "DominantSpecies"))

labeled_sample <- merge(sampleDataforMerge, sampleLabel, by.x= "SampleID", by.y="sampleID", all=TRUE)
cl_labeled_sample <- merge(sampleDataforMerge, sampleLabel, by.x= "SampleID", by.y="sampleID", all=FALSE)
#View(cl_labeled_sample)

#the ones with my extra rows but not the extra rows in Alice's
dim(sampleDataforMerge)
sampleClean<-merge(sampleDataforMerge, sampleLabel, by.x= "SampleID", by.y="sampleID", all.x=TRUE, all.y=FALSE)
dim(sampleClean)
#View(sampleClean)

write.csv(cl_labeled_sample,"Cleaned Fungal Data")
write_xlsx(cl_labeled_sample, "Cleaned Fungal Data Sheet")

nrow(labeled_sample)

dim(cl_labeled_sample)
dim(labeled_sample)

#View(labeled_sample)
names(sample_data(fungal_physeq_subset))

#View(sample_data(fungal_physeq_subset))

head(sample_names(fungal_physeq_subset)) 
#head((sample_data(fungal_physeq_subset))
dim(sample_data(fungal_physeq_subset))
dim(sampleLabel)
#View(sampleLabel)

#table
table(is.na(labeled_sample$is_blank), is.na(labeled_sample$sampleType)) #1766 were matched, 105 were in fungal data but not Alice's data, 1249 were in Alice's data but not fungal data
table(is.na(labeled_sample$sampleType))


#April 17
#checking to see the sample names for the 105 and 1249 ones
subset1 <- labeled_sample[!is.na(labeled_sample$is_blank) & is.na(labeled_sample$sampleType), ]
dim(subset1)
View(subset1) #the ones not found in Alice's data
subset2 <- labeled_sample[is.na(labeled_sample$is_blank) & !is.na(labeled_sample$sampleType), ]
dim(subset2)
View(subset2) #the ones not found in my data

#check rownames
head(rownames(sampleLabel))  
head(sampleLabel$sampleID)  
rownames(sampleLabel) <- sampleLabel$sampleID 

clean_labeled_sample<-labeled_sample[, c("SampleID", "is_blank", "DominantSpecies", "sampleType")] #removing the collumns that I don't want
#View(sample_data(fungal_physeq_subset))
head(sample_data(fungal_physeq_subset))
head(clean_labeled_sample)
#View((clean_labeled_sample))

#rownames(clean_labeled_sample) <- sample_names(fungal_physeq_subset)
#nrow(clean_labeled_sample)

#sample_data(fungal_physeq_subset) <- clean_labeled_sample
head(sample_names(fungal_physeq_subset))
head(rownames(clean_labeled_sample))
nrow(clean_labeled_sample)
dim(sample_data(fungal_physeq_subset))


#randomly picking 5 fungus and BLAST
set.seed(345)
randomrowsFungus <- sample(1:nrow(tax_table(fungal_physeq_subset)), 5)
randomTaxaFungus <- tax_table(fungal_physeq_subset)[randomrowsFungus, ]
randomTaxaFungus[, c("sequence", "Species")]

#creating a new phyloseq object with updated sample table
fungal2.0 <- fungal_physeq_subset
rownames(cl_labeled_sample) <- cl_labeled_sample$SampleID
sample_data(fungal2.0) <- cl_labeled_sample

#get taxa table, otu table, sample data
fungalTaxa <- tax_table(fungal2.0)
fungalOTU <- otu_table(fungal2.0)
fungalSample <- as.data.frame(as.matrix(sample_data(fungal2.0)))

#filtering for just vaginal 
vaginal_sample_data_phyloseq <- fungalSample %>% 
  filter(sampleType=="vaginal")
vaginal_fungal_otu <- fungalOTU[rownames(vaginal_sample_data_phyloseq), , drop = FALSE]
vaginal_fungal_taxa <- fungalTaxa[colnames(vaginal_fungal_otu), , drop = FALSE]

#vaginal phyloseq obj
vaginal_phyloseq <- phyloseq(otu_table(vaginal_fungal_otu, taxa_are_rows=FALSE), 
                             sample_data(vaginal_sample_data_phyloseq), tax_table(vaginal_fungal_taxa))

#calculating for relative abundance
vaginal_phyloseq_rel <- transform_sample_counts(vaginal_phyloseq, function(x) x / sum(x))
#add C. albicans rel. abundance to sample data
ca_phy <- subset_taxa(vaginal_phyloseq_rel, Species == "Candida_albicans")
ca_abund <- rowSums( otu_table(ca_phy)[ , , drop = FALSE ] )
sample_data(vaginal_phyloseq_rel)$CA_abund <- ca_abund[ sample_names(vaginal_phyloseq_rel) ]



#Alpha Diversity
alpha_div <- estimate_richness(fungal_physeq_subset, measures = c("Shannon"))

summary(alpha_div[[1]]) #alpha_div is a list and we want the first element in the list


plot_richness(fungal_physeq_subset, x = "SampleID", measures = c("Shannon")) + theme_minimal()

plot_richness(fungal_physeq_subset, x = "SampleID", color = "DominantSpecies", measures = c("Shannon")) +
  theme_minimal() +
  ggtitle("Alpha Diversity Colored by Dominant Species")
#this looks chaotic

#creating the sampledata as a data frame
fungalSample <- sampleClean
fungalSample$Shannon <- alpha_div[[1]]
colnames(fungalSample)
#colnames(fungalSample)[4] <- "Shannon"
View(fungalSample)
summary(fungalSample)
hist(fungalSample$Shannon) 
summary(fungalSample$sampleType == "vaginal")

#strip rows with empty dominant species out so that when we're calculating shannon we don't take empty as a species
fungalSamplecl <- fungalSample %>%
  filter(DominantSpecies != "", !is.na(DominantSpecies))

#Shannon for vaginal
summary(fungalSample$Shannon[fungalSample$sampleType == "vaginal"])
boxplot(fungalSample$Shannon[fungalSample$sampleType == "vaginal"], fungalSample$Shannon[fungalSample$sampleType == "fecal"], names = c("vaginal", "fecal"))


VaginalSample <- fungalSamplecl %>% 
  filter(sampleType == "vaginal")

dim(VaginalSample)

#group vaginal sample by participants
groupedVaginal <- VaginalSample %>%
  group_by(biome_id) %>%
  summarise(
    avg_Shannon = mean(Shannon, na.rm = TRUE),
    n_samples = n(),
    dominantSpecies = names(sort(table(DominantSpecies), decreasing = TRUE))[1]
  )
#View(groupedVaginal)

#extracting fecal sample
FecalSample <- fungalSamplecl %>% 
  filter(sampleType == "fecal")

#grouping fecal sample by participants
groupedFecal <- FecalSample %>%
  group_by(biome_id) %>%
  summarise(
    avg_Shannon = mean(Shannon, na.rm = TRUE),
    n_samples = n(),
    dominantSpecies = names(sort(table(DominantSpecies), decreasing = TRUE))[1]
  )
#View(groupedFecal)

dim(FecalSample)

#boxplot by participants
?boxplot()
boxplot(groupedVaginal$avg_Shannon, groupedFecal$avg_Shannon, names = c("vaginal", "fecal"), main = "Shannon by Participants")

sort(table(VaginalSample$DominantSpecies)) 
#Candida_albicans 655 empty 193 globosa 94 restricta 80 arunalokei 18 Malassezia_globosa 17 Candida_parapsilosis 16

sort(table(FecalSample$DominantSpecies))
#Candida_albicans 267 empty 197 restricta 31 globosa 28 Malassezia_globosa 14 Candida_parapsilosis 11

#Diversity is different but not the dominant species 

#grouping sample by participants
groupedSampleP <- fungalSample %>%
  group_by(biome_id) %>%
  summarise(
    avg_Shannon = mean(Shannon, na.rm = TRUE),
    n_samples = n()
  )

#View(groupedSampleP) 

table(groupedVaginal$dominantSpecies)
table(groupedFecal$dominantSpecies)


#DASS importing
#file_path <- file.choose()
#file.info(file_path)

#Alice's codes
dass_data <- read.csv("/Users/caoyang/Desktop/Tetel Lab/datasets/Report 10-DASS-21.csv")
id_mapping <- read.csv("/Users/caoyang/Desktop/Tetel Lab/datasets/Original Study Mapping - Sheet3.csv")


### Map uid to study id

# Get unique study and biome health pairings
study_and_u_id <- unique(id_mapping %>% 
                           select(STUDY.ID, Biome.Health.App.ID))

# Match and join columns
study_and_u_id <- study_and_u_id %>% 
  rename("study_id" = "STUDY.ID") %>% 
  rename("biome_id" = "Biome.Health.App.ID")
dass_data <- dass_data %>% 
  rename("biome_id" = "Your.Biome.Health.App.ID")
study_and_u_id$study_id <- as.character(study_and_u_id$study_id)

# Map ids
dass_data <- dass_data %>%
  left_join(study_and_u_id, by = "biome_id") %>%
  mutate(biome_id = coalesce(study_id, biome_id)) %>%
  select(-study_id)

# Check missing ids
missing_list <- dass_data %>%
  filter(is.na(as.numeric(biome_id)))
print(unique(missing_list$biome_id))

dim(dass_data)
#View(dass_data)

### Save final data output
write.csv(dass_data,
          file = "/Users/caoyang/Desktop/Tetel Lab/datasets/cleaned_dass.csv",
          row.names = FALSE)

# ------------------------------------------------------------------------------
# calculate depression, anxiety, and stress scores 
# ------------------------------------------------------------------------------
dass <- dass_data
#Convert to Date with origin 1970-01-01 (Unix epoch, default in R)
dass$Timestamp <- as.Date(dass$Timestamp, format="%m/%d/%Y", tz="UTC") #the "%Y" must be capitalized
head(dass_data$Timestamp)
head(dass$Timestamp)
id_values <- unique(dass$biome_id) 
table(id_values)
# time_values <- sort(unique(dass$Timestamp))
# num_time_values <- as.numeric(time_values)
dass <- dass[order(dass$Timestamp),]

colnames(dass)

summary(dass$Timestamp)
#Remove any values that occurred after 12/16 (end of semester)
dass$Timestamp_numeric <- as.numeric(dass$Timestamp)
table(dass$Timestamp)
dass <- dass[dass$Timestamp_numeric <= 19342, ] 
View(dass)
?as.Date

dass$week <- rep(NA, nrow(dass))
dass$week[dass$Timestamp >= 19276 & dass$Timestamp <= 19280] <- 1
dass$week[dass$Timestamp >= 19281 & dass$Timestamp <= 19286] <- 2
dass$week[dass$Timestamp >= 19290 & dass$Timestamp <= 19293] <- 3
dass$week[dass$Timestamp >= 19296 & dass$Timestamp <= 19300] <- 4
dass$week[dass$Timestamp >= 19302 & dass$Timestamp <= 19308] <- 5
dass$week[dass$Timestamp >= 19312 & dass$Timestamp <= 19315] <- 6
dass$week[dass$Timestamp >= 19319 & dass$Timestamp <= 19321] <- 7
dass$week[dass$Timestamp >= 19323 & dass$Timestamp <= 19329] <- 8
dass$week[dass$Timestamp >= 19330 & dass$Timestamp <= 19336] <- 9
dass$week[dass$Timestamp >= 19339] <- 10

dass <- dass %>% 
filter(!is.na(week))

colnames(dass) 

#depression score used column name
dass$noPositiveFeeling <- dass$X3..I.couldn.t.seem.to.experience.any.positive.feeling.at.all
dass$initiative <- dass$X5..I.found.it.difficult.to.work.up.the.initiative.to.do.things
dass$noLookForward <- dass$X10..I.felt.that.I.had.nothing.to.look.forward.to
dass$downhearted <- dass$X13..I.felt.down.hearted.and.blue
dass$noEnthusiasm <- dass$X16..I.was.unable.to.become.enthusiastic.about.anything
dass$feelWorthless <- dass$X17..I.felt.I.wasn.t.worth.much.as.a.person
dass$lifeMeaningless <- dass$X21..I.felt.that.life.was.meaningless

#anxiety score used column name
dass$mouthDry <- dass$X2..I.was.aware.of.dryness.of.my.mouth
dass$difficultyBreathing <- dass$X4..I.experienced.breathing.difficulty..e.g..excessively.rapid.breathing..breathlessness.in.the.absence.of.physical.exertion..
dass$trembling <- dass$X7..I.experienced.trembling..e.g..in.the.hands.
dass$panicSituation <- dass$X9..I.was.worried.about.situations.in.which.I.might.panic.and.make.a.fool.of.myself
dass$closeToPanic <- dass$X15..I.felt.I.was.close.to.panic
dass$awareHeart <- dass$X19..I.was.aware.of.the.action.of.my.heart.in.the.absence.of.physical.exertion..e.g..sense.of.heart.rate.increase..heart.missing.a.beat.
dass$scared <- dass$X20..I.felt.scared.without.any.good.reason

#stress score used column name
dass$windDown <- dass$X1..I.found.it.hard.to.wind.down
dass$overreact <- dass$X6..I.tended.to.over.react.to.situations
dass$nervous <- dass$X8..I.felt.that.I.was.using.a.lot.of.nervous.energy
dass$agitated <- dass$X11..I.found.myself.getting.agitated
dass$difficultyRelax <- dass$X12..I.found.it.difficult.to.relax
dass$intolerant <- dass$X14..I.was.intolerant.of.anything.that.kept.me.from.getting.on.with.what.I.was.doing
dass$touchy <- dass$X18..I.felt.that.I.was.rather.touchy

colnames(dass)

dass$depression_score <- rep(NA, nrow(dass))
dass$anxiety_score <- rep(NA, nrow(dass))
dass$stress_score <- rep(NA, nrow(dass))

for(id in id_values){
  for(week in 1:10){
    dass$depression_score[dass$biome_id==id & dass$week==week] <- 2*sum(dass$noPositiveFeeling[dass$biome_id==id & dass$week==week] + dass$initiative[dass$biome_id==id & dass$week==week] + dass$noLookForward[dass$biome_id==id & dass$week==week] + 
                                                                          dass$downhearted[dass$biome_id==id & dass$week==week] + dass$noEnthusiasm[dass$biome_id==id & dass$week==week] + dass$feelWorthless[dass$biome_id==id & dass$week==week] + 
                                                                          dass$lifeMeaningless[dass$biome_id==id & dass$week==week], na.rm=TRUE)
    
    dass$anxiety_score[dass$biome_id==id & dass$week==week] <- 2*sum(dass$mouthDry[dass$biome_id==id & dass$week==week] + dass$difficultyBreathing[dass$biome_id==id & dass$week==week] + dass$trembling[dass$biome_id==id & dass$week==week] + 
                                                                       dass$panicSituation[dass$biome_id==id & dass$week==week] + dass$closeToPanic[dass$biome_id==id & dass$week==week] + 
                                                                       dass$awareHeart[dass$biome_id==id & dass$week==week] + dass$scared[dass$biome_id==id & dass$week==week], na.rm=TRUE)
    
    dass$stress_score[dass$biome_id==id & dass$week==week] <- 2*sum(dass$windDown[dass$biome_id==id  & dass$week==week] + dass$overreact[dass$biome_id==id & dass$week==week] + dass$nervous[dass$biome_id==id & dass$week==week] + 
                                                                      dass$agitated[dass$biome_id==id & dass$week==week] + dass$difficultyRelax[dass$biome_id==id & dass$week==week] + 
                                                                      dass$intolerant[dass$biome_id==id & dass$week==week] + dass$touchy[dass$biome_id==id & dass$week==week], na.rm=TRUE)
  }
}


#creation of variables to categorize scores into level of severity 
dass$depressionseverity <- rep(NA, nrow(dass))
dass$anxietyseverity <- rep(NA, nrow(dass))
dass$stressseverity<- rep(NA, nrow(dass))

dass$depressionseverity[dass$depression_score>=0 & dass$depression_score<=9] <- 0 
dass$depressionseverity[dass$depression_score>=10 & dass$depression_score<=13] <- 1 
dass$depressionseverity[dass$depression_score>=14 & dass$depression_score<=20] <- 2
dass$depressionseverity[dass$depression_score>=21 & dass$depression_score<=27] <- 3
dass$depressionseverity[dass$depression_score>=28] <- 4

dass$anxietyseverity[dass$anxiety_score>=0 & dass$anxiety_score<=7] <- 0
dass$anxietyseverity[dass$anxiety_score>=8 & dass$anxiety_score<=9] <- 1
dass$anxietyseverity[dass$anxiety_score>=10 & dass$anxiety_score<=14] <- 2
dass$anxietyseverity[dass$anxiety_score>=15 & dass$anxiety_score<=19] <- 3
dass$anxietyseverity[dass$anxiety_score>=20] <- 4

dass$stressseverity[dass$stress_score>=0 & dass$stress_score<=14] <- 0
dass$stressseverity[dass$stress_score>=15 & dass$stress_score<=18] <- 1
dass$stressseverity[dass$stress_score>=19 & dass$stress_score<=25] <- 2
dass$stressseverity[dass$stress_score>=26 & dass$stress_score<=33] <- 3
dass$stressseverity[dass$stress_score>=34] <- 4

View(dass)
plot(dass$Timestamp, dass$stress_score)
smoothingSpline = smooth.spline(dass$Timestamp, dass$stress_score, spar=0.35)
lines(smoothingSpline)


# average stress score
dass.avg <- dass %>% 
  group_by(biome_id) %>% 
  summarise(
    avg_depr=sum(depression_score)/n(),
    avg_anx=sum(anxiety_score)/n(),
    avg_stress=sum(stress_score)/n()
  )  

View(dass.avg)

shannonDassbyP <- merge(groupedSampleP, dass.avg, by.x="biome_id", by.y="biome_id", all=TRUE)

lmDepr <- lm(avg_Shannon ~ avg_depr, data = shannonDassbyP)
summary(lmDepr) #pval = 0.8, adj r sqr = -0.01

lmAnx <- lm(avg_Shannon ~ avg_anx, data = shannonDassbyP)
summary(lmAnx) #pval = 0.7, adj r sqr = -0.01

lmStr <- lm(avg_Shannon ~ avg_stress, data = shannonDassbyP)
summary(lmStr) #pval = 0.228, adj r sqr = 0.007



ggplot(shannonDassbyP) +
  geom_point(aes(x = avg_depr, y = avg_Shannon, color = "Depression"), size = 2) +
  geom_point(aes(x = avg_anx, y = avg_Shannon, color = "Anxiety"), size = 2) +
  geom_point(aes(x = avg_stress, y = avg_Shannon, color = "Stress"), size = 2) +
  geom_smooth(aes(x = avg_depr, y = avg_Shannon, color = "Depression"),
              method = "gam", formula = y ~ s(x, bs = "cs"), se = FALSE) +
  geom_smooth(aes(x = avg_anx, y = avg_Shannon, color = "Anxiety"),
              method = "gam", formula = y ~ s(x, bs = "cs"), se = FALSE) +
  geom_smooth(aes(x = avg_stress, y = avg_Shannon, color = "Stress"),
              method = "gam", formula = y ~ s(x, bs = "cs"), se = FALSE) +
  scale_color_manual(
    name = "DASS",
    values = c(
      "Depression" = "blue",  
      "Anxiety" = "red",     
      "Stress" = "green"       
    )
  ) +
  theme_classic() +
  labs(x = "DASS Score", y = "Shannon Diversity Index")

###########################################################################

#mapping c.albicans and dass scores
#convert sample_data to csv
vaginal_rel_metadata_df <- as(sample_data(vaginal_phyloseq_rel), "data.frame")
vaginal_rel_metadata_df$biome_id <- as.integer(vaginal_rel_metadata_df$biome_id)
View(vaginal_rel_metadata_df)

vaginal_rel_metadata_df$logDate <- as.Date(vaginal_rel_metadata_df$logDate)

dass$biome_id <- as.numeric(dass$biome_id)

####################
# in progress
# dass_clean <- dass %>%
#   select(biome_id, Timestamp, depression_score, anxiety_score, stress_score) %>%
#   distinct()
# #merge dass with csv
# vaginal_rel_metadata_dass_df <- vaginal_rel_metadata_df %>%
#   left_join(
#     dass %>% select(biome_id, depression_score, anxiety_score, stress_score),
#     by = "biome_id"
#   )





