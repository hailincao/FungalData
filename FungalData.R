library(phyloseq)
library(dplyr)
library(stringr)
library(writexl)
library(ggplot2)
library(decontam)
library(phyloseq)
library(vegan)
library(pheatmap)
library(tidyverse)
library(Matrix)
library(readxl)
library(magrittr)
library(lme4)
library(lmerTest)
library(performance)
library(ggforce)
library(scales)   
library(viridis)  
library(patchwork)
library(lubridate)
library(fuzzyjoin)
#remotes::install_github("david-barnett/microViz")
library(microViz)

data <- readRDS("/Users/caoyang/Desktop/Tetel Lab/Walther-Antonio_Project_022_ITS2.rds") #the data is reading an email forwarded by Alice to Helena 
############################################################################################################
#exploration on fungal data
taxaF <- tax_table(data)
otuF <- otu_table(data)

#rows in both datasets are genetic sequences, the otu dataset says 1/0 whether it has that genetic sequence
#rownames in taxa table : specific sequences; whether a sample have this sequence: indicated by 0/1

#so each row is not species, since the number of species is much shorter than #rows, maybe it's substring?
#there should be somewhere that associate the id numbers in otu with the specific participant 
taxa_matrix <- taxaF@.Data
rownames <- rownames(taxaF@.Data)
colnames <- colnames(taxaF@.Data)
taxa_matrix[, "SH"]
sequence <- taxaF@.Data[, "sequence"]

species <- taxaF@.Data[, "Species"]

table(taxaF@.Data[, "Species"])
names(which.max(table(taxaF@.Data[, "Species"])))

species <- taxaF@.Data[, "Species"]
nonemptyspecies <- species[species != ""]# Remove empty strings
emptystrings <- species[species == ""]
length(emptystrings) #8560/15293 are empty
most_frequent_species <- names(which.max(table(nonemptyspecies)))

############################################################################################################
#BLAST check
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


length(table(rownames)) #are there repeating rownames? no
############################################################################################################
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

#filter the contaminants
fungal_physeq_no_contam <- prune_taxa(!contaminants, fungal_physeq)
fungal_physeq_subset <- subset_taxa(fungal_physeq_no_contam, Kingdom == "Fungi" & !is.na(Phylum) & Phylum != "")

#checking the dimension
#rows are the number of samples, columns are taxa
dim(otu_table_obj)
dim(otu_table(fungal_physeq_subset))
dim(otu_table(fungal_physeq_no_contam))
otu_table(fungal_physeq_subset)[1,]


dominant_spec <- apply(t(otu_table(fungal_physeq_subset)), 2, function(x) {
  spec <- tax_table(fungal_physeq_subset)[which.max(x), "Species"]
  # ifelse(is.na(spec), "Unknown", spec)
})

#creating dominant species column
sample_data(fungal_physeq_subset)$DominantSpecies <- dominant_spec
#check if the column is there
sample_data(fungal_physeq_subset)$DominantSpecies
#check the class of the object we created
summary(fungal_physeq_subset)
#the sample data part of our object
sort(table(sample_data(fungal_physeq_subset)$DominantSpecies))  # Candida_albicans  the most frequently appeared

#turning the otu table back to data frame
physeqOTU<- as.data.frame(otu_table(fungal_physeq_subset))
physeqOTU["F1376",1:10] #the first ten rows of this sample column

nonzero <- which(physeqOTU["F1376", ] != 0) 
nonzero #the column numbers of sequence ID
physeqOTU["F1376", nonzero] #the counts of the nonzero species
sequenceIDnonzero <- colnames(physeqOTU)[nonzero]
#> colnames(physeqOTU)[nonzero]
#"074f81db997e702d17c85be0c46b03ad" "9589a4186e70b56c852377d7562d4789" "f8f00151e8dd21d0f490343afb43d854"
#"f74f09973e588cd61df7e11c6620c2fd" "8d50d5e4d31b744eb9647746878c017e"
taxaF@.Data["074f81db997e702d17c85be0c46b03ad", "Species"] #it's Candida_albicans

###################################################################
#creating sample data
sampleLabel<-read_excel("/Users/caoyang/Desktop/Tetel Lab/cleaned_samplesv2.xlsx")
sampleLabel$sampleID <- sampleLabel$qr
sampleLabel <- sampleLabel[, c("sampleID", "biome_id", "sampleType", "logDate")]

temp <- t(otu_table(fungal_physeq_subset))
#temp[1:10, 1:10]
temp2<-temp[, "F1376"]
#head(temp2)
#table(temp2)

which.max(temp2)
temp2[2024]
rownames(temp2[1:10])
rownames(temp2[temp2>0])
taxaF@.Data["f8f00151e8dd21d0f490343afb43d854", "Species"] 

###################################################################
#merging sample label and sample data
#removing the ones that are blank
fungal_physeq_subset <- subset_samples(
  fungal_physeq_subset,
  !str_detect(sample_names(fungal_physeq_subset), "^BLANK")
)
# Trim .ITS2
sample_names(fungal_physeq_subset) <- str_remove(sample_names(fungal_physeq_subset), "\\.ITS2$")
temp <-sample_data(fungal_physeq_subset) 
sampleDataforMerge <- data.frame(temp$SampleID, temp$is_blank, temp$DominantSpecies)
colnames(sampleDataforMerge) <- (c("SampleID", "is_blank", "DominantSpecies"))

#merging with or without all rows
labeled_sample <- merge(sampleDataforMerge, sampleLabel, by.x= "SampleID", by.y="sampleID", all=TRUE) #this is the one that did not omit anything from Alice or my data
cl_labeled_sample <- merge(sampleDataforMerge, sampleLabel, by.x= "SampleID", by.y="sampleID", all=FALSE) #this one only kept the ones that were matched in both Alice's and my data

#the ones with my extra rows but not the extra rows in Alice's
sampleClean<-merge(sampleDataforMerge, sampleLabel, by.x= "SampleID", by.y="sampleID", all.x=TRUE, all.y=FALSE)
dim(sampleClean)

#exporting out the sample data
write.csv(cl_labeled_sample,"Cleaned Fungal Data")
write_xlsx(cl_labeled_sample, "Cleaned Fungal Data Sheet")


#table of how many rows are missing from each merged dataset
table(is.na(labeled_sample$is_blank), is.na(labeled_sample$sampleType)) #1766 were matched, 105 were in fungal data but not Alice's data, 1249 were in Alice's data but not fungal data
table(is.na(labeled_sample$sampleType))


#checking to see the sample names for the 105 and 1249 ones
subset1 <- labeled_sample[!is.na(labeled_sample$is_blank) & is.na(labeled_sample$sampleType), ]
dim(subset1)
#View(subset1) #the ones not found in Alice's data
subset2 <- labeled_sample[is.na(labeled_sample$is_blank) & !is.na(labeled_sample$sampleType), ]
dim(subset2)
#View(subset2) #the ones not found in my data

#checking rownames
head(rownames(sampleLabel))  
head(sampleLabel$sampleID)  
rownames(sampleLabel) <- sampleLabel$sampleID 

#a sample data that I did not end up using -- cuz it only contained the ones in Alice's sheet?
clean_labeled_sample<-labeled_sample[, c("SampleID", "is_blank", "DominantSpecies", "sampleType")] #removing the collumns that I don't want
head(sample_data(fungal_physeq_subset))
head(clean_labeled_sample)

head(sample_names(fungal_physeq_subset))
head(rownames(clean_labeled_sample))
nrow(clean_labeled_sample)
dim(sample_data(fungal_physeq_subset))

###############################################################################################
#randomly picking 5 fungus and BLAST
set.seed(345)
randomrowsFungus <- sample(1:nrow(tax_table(fungal_physeq_subset)), 5)
randomTaxaFungus <- tax_table(fungal_physeq_subset)[randomrowsFungus, ]
randomTaxaFungus[, c("sequence", "Species")]

###############################################################################################
#creating a new phyloseq object with updated sample table
#I ended up with the sample data that were matched in both Alice's and my data
fungal2.0 <- fungal_physeq_subset
rownames(cl_labeled_sample) <- cl_labeled_sample$SampleID
sample_data(fungal2.0) <- cl_labeled_sample

#get taxa table, otu table, sample data
fungalTaxa <- tax_table(fungal2.0)
fungalOTU <- otu_table(fungal2.0)
fungalSample <- as.data.frame(as.matrix(sample_data(fungal2.0)))

###############################################################################################
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


###############################################################################################
#trying to plot shannon and colored by dominant species
#Alpha Diversity
alpha_div <- estimate_richness(fungal_physeq_subset, measures = c("Shannon"))

summary(alpha_div[[1]]) #alpha_div is a list and we want the first element in the list

###############################################################################################
#creating the sampledata as a data frame
fungalSample <- sampleClean #the sample data that has mine but not Alice's ID
fungalSample$Shannon <- alpha_div[[1]]
colnames(fungalSample)
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
###############################################################################################
#group vaginal sample by participants
groupedVaginal <- VaginalSample %>%
  group_by(biome_id) %>%
  summarise(
    avg_Shannon = mean(Shannon, na.rm = TRUE),
    n_samples = n(),
    dominantSpecies = names(sort(table(DominantSpecies), decreasing = TRUE))[1]
  )
#View(groupedVaginal)
###############################################################################################
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

dim(FecalSample)

#boxplot by participants
boxplot(groupedVaginal$avg_Shannon, groupedFecal$avg_Shannon, names = c("vaginal", "fecal"), main = "Shannon by Participants")
boxplot(groupedVaginal$avg_Shannon, groupedFecal$avg_Shannon,
        names = c("vaginal", "fecal"),
        main = "Shannon by Participants",
        ylab = "Shannon Index",
        col = c("#FF9999", "#9999FF"))

# Add jittered points
stripchart(list(groupedVaginal$avg_Shannon, groupedFecal$avg_Shannon),
           method = "jitter",
           pch = 16,       # filled circles
           col = rgb(0, 0, 0, 0.5), # semi-transparent black
           vertical = TRUE,
           add = TRUE)
###############################################################################################
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

table(groupedVaginal$dominantSpecies)
table(groupedFecal$dominantSpecies)

###############################################################################################
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

#how stress vary over time
plot(dass$Timestamp, dass$stress_score)
smoothingSpline = smooth.spline(dass$Timestamp, dass$stress_score, spar=0.35)
lines(smoothingSpline)

###############################################################################################
# average stress score
dass.avg <- dass %>% 
  group_by(biome_id) %>% 
  summarise(
    avg_depr=sum(depression_score)/n(),
    avg_anx=sum(anxiety_score)/n(),
    avg_stress=sum(stress_score)/n()
  )  
###############################################################################################
#mapping average shannon at both sites and dass together
shannonDassbyP <- merge(groupedSampleP, dass.avg, by.x="biome_id", by.y="biome_id", all=TRUE)

#relationship between Shannon and DASS -- no relationship seen
lmDepr <- lm(avg_Shannon ~ avg_depr, data = shannonDassbyP)
summary(lmDepr) #pval = 0.8, adj r sqr = -0.01

lmAnx <- lm(avg_Shannon ~ avg_anx, data = shannonDassbyP)
summary(lmAnx) #pval = 0.7, adj r sqr = -0.01

lmStr <- lm(avg_Shannon ~ avg_stress, data = shannonDassbyP)
summary(lmStr) #pval = 0.228, adj r sqr = 0.007

#relationship between DASS and Shannon, Plot -- no association seen
ggplot(shannonDassbyP) +
  geom_point(aes(x = avg_depr, y = avg_Shannon, color = "Depression"), size = 2) +
  geom_point(aes(x = avg_anx, y = avg_Shannon, color = "Anxiety"), size = 2) +
  geom_point(aes(x = avg_stress, y = avg_Shannon, color = "Stress"), size = 2) +
  geom_smooth(aes(x = avg_depr, y = avg_Shannon, color = "Depression"),
              method = "loess", se = FALSE, span = 0.75, size = 1.2) +
  geom_smooth(aes(x = avg_anx, y = avg_Shannon, color = "Anxiety"),
              method = "loess", se = FALSE, span = 0.75, size = 1.2) +
  geom_smooth(aes(x = avg_stress, y = avg_Shannon, color = "Stress"),
              method = "loess", se = FALSE, span = 0.75, size = 1.2) +
  scale_color_manual(
    name = "DASS",
    values = c(
      "Depression" = "blue",  
      "Anxiety" = "red",     
      "Stress" = "green"       
    )
  ) +
  theme_classic() +
  labs(
    title = "Relationship Between DASS Scores and Shannon Diversity",
    x = "DASS Score", 
    y = "Shannon Diversity Index"
  ) +
  xlim(0, 40)


###########################################################################
#mapping c.albicans and dass scores in vagina
#convert sample_data to csv
vaginal_rel_metadata_df <- as(sample_data(vaginal_phyloseq_rel), "data.frame")
vaginal_rel_metadata_df$biome_id <- as.integer(vaginal_rel_metadata_df$biome_id)

#reading in the lifestyle variables (merged with vaginal data) that Nicky created
hbcmerged_df <- read.csv("/Users/caoyang/Desktop/Tetel Lab/datasets/vaginal_rel_metadata_hbc_df_matched.csv")

#calculating average c.albicans rel.abundanxe
CAavgP <- hbcmerged_df %>%
  group_by(biome_id) %>%
  summarise(
    avg_CA = mean(CA_abund, na.rm = TRUE),
    n_samples = n()
  )

CAavgDass <- merge(CAavgP, dass.avg, by.x="biome_id", by.y="biome_id", all=TRUE)

#no relationship between avg CA and DASS seen
ggplot(CAavgDass) +
  geom_point(aes(x = avg_depr, y = avg_CA, color = "Depression"), size = 2) +
  geom_point(aes(x = avg_anx, y = avg_CA, color = "Anxiety"), size = 2) +
  geom_point(aes(x = avg_stress, y = avg_CA, color = "Stress"), size = 2) +
  geom_smooth(aes(x = avg_depr, y = avg_CA, color = "Depression"),
              method = "gam", formula = y ~ s(x, bs = "cs"), se = FALSE) +
  geom_smooth(aes(x = avg_anx, y = avg_CA, color = "Anxiety"),
              method = "gam", formula = y ~ s(x, bs = "cs"), se = FALSE) +
  geom_smooth(aes(x = avg_stress, y = avg_CA, color = "Stress"),
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
  labs(x = "DASS Score", y = "C.albicans Avg Relative Abundance by P")


#mapping c.albicans abundance by date 
dass$biome_id <- as.numeric(dass$biome_id)
hbcmerged_df$logDate <- as.Date(hbcmerged_df$logDate)

###########################################################################
#merging my dass data and the data that Nicky created
dass_CA <- left_join(dass, hbcmerged_df, by = c("Timestamp"="logDate", "biome_id"="biome_id"), relationship = "many-to-many")

####################
dass_CA2 <- dass_CA %>%
   select(biome_id, Timestamp, depression_score, anxiety_score, stress_score.x, CA_abund) %>%
   distinct()

dassCAplot <- dass_CA2 %>%
  pivot_longer(cols = c(depression_score, anxiety_score, stress_score.x),
               names_to = "ScoreType",
               values_to = "Score")

ggplot(dassCAplot, aes(x = Score, y = CA_abund, color = ScoreType)) +
  geom_point(alpha = 0.7, size = 2) +
  geom_smooth(method = "loess", se = FALSE, span = 0.75, size = 1.2) +
  scale_color_manual(values = c("depression_score" = "blue", 
                                "anxiety_score" = "red", 
                                "stress_score.x" = "green")) +
  labs(x = "DASS Score",
       y = "C. albicans Relative Abundance",
       color = NULL) +
  xlim(0, 50) +
  theme_minimal(base_size = 13) +
  theme(
    legend.position = "top",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )


ggplot(dassCAplot, aes(
  x = Score,
  y = CA_abund,
  color = recode(ScoreType,
                 "depression_score" = "Depression",
                 "anxiety_score" = "Anxiety",
                 "stress_score.x" = "Stress"))
) +
  geom_point(alpha = 0.7, size = 2) +
  geom_smooth(method = "loess", se = FALSE, span = 0.75, size = 1.2) +
  scale_color_manual(
    values = c(
      "Depression" = "blue", 
      "Anxiety" = "red", 
      "Stress" = "green"
    )
  ) +
  labs(
    title = "Relationship Between DASS Scores and C. albicans Abundance",
    x = "DASS Score",
    y = "C. albicans Relative Abundance",
    color = "Score Type"
  ) +
  xlim(0, 50) +
  theme_minimal(base_size = 13) +
  theme(
    legend.position = "right",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_text(hjust = 0.5, face = "bold")
  )


#checking p-val, non of them significant
dassCAplot %>%
  group_by(ScoreType) %>%
  summarise(correlation = cor(Score, CA_abund, use = "complete.obs"),
            p_value = cor.test(Score, CA_abund)$p.value,
            method = "Pearson")
#merge dass with csv
# vaginal_rel_metadata_dass_df <- vaginal_rel_metadata_df %>%
#   left_join(
#     dass %>% select(biome_id, depression_score, anxiety_score, stress_score),
#     by = "biome_id"
#   )

###############################################################################################
#adding a column that stands for c.albican relative abundance calculated by rate out of total reads per sample
otu_mat <- otu_table(fungal2.0)
# Transpose only if taxa are not rows
otu_mat <- if (!taxa_are_rows(fungal2.0)) {
  t(otu_mat)
} else {
  otu_mat
}
otu_table(fungal2.0) <- otu_mat

#adding total reads to sample data
total_reads<- colSums(otu_table(fungal2.0))
sample_data(fungal2.0)$total_reads <- total_reads

#how many ASVs are identified
asv_counts <- apply(otu_table(fungal2.0), 2, function(x) sum(x > 0))
sample_data(fungal2.0)$ASV_counts <- asv_counts

#adding C.albicans abudnance
sample_data(fungal2.0)$Candida_abundance <- fungal2.0 %>%
  {
    otu <- otu_table(.)
    tax <- tax_table(.)
    candida_asvs <- rownames(tax)[tax[, "Species"] == "Candida_albicans"]
    candida_counts <- otu[candida_asvs, , drop = FALSE] %>%
      apply(2, sum)
    candida_counts / sample_sums(.)
  }

#######################################################################
#adding a column for shannon diversity
shannon_df <- estimate_richness(fungal2.0, measures = "Shannon") %>%
  rownames_to_column(var = "SampleID")
meta_df <- sample_data(fungal2.0) %>%
  data.frame()
merged_df <- left_join(meta_df, shannon_df, by = "SampleID")
sample_data(fungal2.0) <- merged_df %>%
  column_to_rownames(var = "SampleID") %>%
  sample_data()

#######################################################################
#merging merging merging vaginal data

#reading in bacteria data
bacteria_abundance <- read.csv("/Users/caoyang/Desktop/Tetel Lab/datasets/microbiome_crosstalk_merged_abund.csv")

#separating gut and vaginal part for merge
fungal_sample_df <- data.frame(sample_data(fungal2.0))
# Subset to vaginal samples
fungal_vaginal <- fungal_sample_df %>%
  filter(sampleType == "vaginal")
# Subset to gut samples
fungal_gut <- fungal_sample_df %>%
  filter(sampleType == "fecal")

#mutating names
fungal_vaginal <- fungal_vaginal %>%
  rename(calbican_rel_abundance_vag = Candida_abundance) %>%
  rename(Shannon_vag = Shannon) %>%
  rename(dominant_species_fungal_vag = DominantSpecies)

fungal_gut <- fungal_gut %>%
  rename(calbican_rel_abundance_gut = Candida_abundance) %>%
  rename(Shannon_gut = Shannon) %>%
  rename(dominant_species_fungal_gut = DominantSpecies)


vag_abund <- fungal_vaginal %>%
  select(biome_id, logDate, calbican_rel_abundance_vag, Shannon_vag, dominant_species_fungal_vag)

gut_abund <- fungal_gut %>%
  select(biome_id, logDate, calbican_rel_abundance_gut, Shannon_gut, dominant_species_fungal_gut)

bacteria_abundance <- bacteria_abundance %>%
  mutate(logDate = as.Date(logDate))

vag_abund <- vag_abund %>%
  mutate(logDate = as.Date(logDate))

gut_abund <- gut_abund %>%
  mutate(logDate = as.Date(logDate))

#using the first sample for each participant on a certain day
vag_abund_unique <- vag_abund %>%
  arrange(biome_id, logDate) %>%
  group_by(biome_id, logDate) %>%
  slice(1) %>% ungroup()

gut_abund_unique <- gut_abund %>%
  arrange(biome_id, logDate) %>%
  group_by(biome_id, logDate) %>%
  slice(1) %>% ungroup()

bacteria_abundance_merged <- bacteria_abundance %>%
  left_join(vag_abund_unique, by = c("biome_id", "logDate")) %>%
  left_join(gut_abund_unique, by = c("biome_id", "logDate"))

#####################################################################
#correlation between c.albicans and lactobacillus
#548/1287 of the rows are NA for vaginal c.albicans abundance
cor_df <- bacteria_abundance_merged %>%
  filter(!is.na(calbican_rel_abundance_vag))
ggplot(cor_df, aes(x = lacto_rel_abundance_vag, y = calbican_rel_abundance_vag)) +
  geom_point(alpha = 0.7) +
  labs(x = "Vaginal Lactobacillus Abundance",
       y = "Vaginal C. albicans Relative Abundance",
       title = "Correlation between Lactobacillus and C. albicans in Vagina") +
  theme_minimal()

ggplot(cor_df, aes(x = lacto_rel_abundance_gut, y = calbican_rel_abundance_vag)) +
  geom_point(alpha = 0.7) +
  labs(x = "Gut Lactobacillus Abundance",
       y = "Vaginal C. albicans Relative Abundance",
       title = "Correlation between Lactobacillus in Gut and C. albicans in Vagina") +
  theme_minimal()

#correlation between Alistipes putredinis and C.albicans
ggplot(cor_df, aes(x = putredinis_rel_abundance_vag, y = calbican_rel_abundance_vag)) +
  geom_point(alpha = 0.7) +
  labs(x = "Vaginal A.putredinis Abundance",
       y = "Vaginal C. albicans Relative Abundance",
       title = "Correlation between A.putredinis and C. albicans in Vagina") +
  theme_minimal()

ggplot(cor_df, aes(x = putredinis_rel_abundance_gut, y = calbican_rel_abundance_gut)) +
  geom_point(alpha = 0.7) +
  labs(x = "Gut A.putredinis Abundance",
       y = "Gut C.albicans Relative Abundance",
       title = "Correlation between A.putredinis in Gut and C. albicans in Gut") +
  theme_minimal()


###############################################################################
#C.albicans abundance by CST
CST_df <- bacteria_abundance_merged %>%
  filter(!is.na(CST))

ggplot(CST_df, aes(x = CST, y = calbican_rel_abundance_vag, fill = CST)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.6) +
  geom_jitter(width = 0.2, alpha = 0.5, color = "black", size = 1) +
  labs(
    title = "Distribution of Vaginal C.albican Relative Abundance by CST per sample",
    x = "Vaginal CST",
    y = "Vaginal C.albicans Relative Abundance"
  ) +
  theme_minimal()


CST_ca <- lmer(calbican_rel_abundance_vag ~ CST + (1 | biome_id), data = CST_df)
summary(CST_ca)


#set CST V as the reference level
CST_df$CST <- factor(CST_df$CST)  # Convert to factor if not already
CST_df$CST <- relevel(CST_df$CST, ref = "V")  # Set CSTV as reference level
CST_ca_rev <- lmer(calbican_rel_abundance_vag ~ CST + (1 | biome_id), data = CST_df)
summary(CST_ca_rev)

#adding centered days 
all_days <- seq.Date(as.Date("2022-10-13"), as.Date("2022-12-16"), by = "day")
CST_df$study_day <- match(as.Date(CST_df$logDate), all_days) -1

CST_df <- CST_df %>%
  mutate(day_c = scale(study_day, center = TRUE, scale = FALSE))

CST_ca_rev_days <- lmer(calbican_rel_abundance_vag ~ CST + day_c + I(day_c^2) + (1 | biome_id), data = CST_df)
summary(CST_ca_rev_days)


ggplot(CST_df, aes(x = logDate, y = calbican_rel_abundance_vag, color = CST)) +
  geom_point(alpha = 0.6, size = 1.8) +
  geom_smooth(method = "loess", se = FALSE, size = 1) +
  scale_color_brewer(palette = "Set1") +
  labs(
    title = "Vaginal Candida albicans Abundance by CST",
    x = "Date",
    y = "Vaginal C.Albicans Relative Abundance"
  ) +
  theme_minimal()


########
#gut CA and CST
ggplot(CST_df, aes(x = CST, y = calbican_rel_abundance_gut, fill = CST)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.6) +
  geom_jitter(width = 0.2, alpha = 0.5, color = "black", size = 1) +
  labs(
    title = "Distribution of Gut C.albican Relative Abundance by Vaginal Microbiome CST per sample",
    x = "Vaginal CST",
    y = " Gut C.albicans Relative Abundance"
  ) +
  theme_minimal()

CST_ca_gut <- lmer(calbican_rel_abundance_gut ~ CST + (1 | biome_id), data = CST_df)
summary(CST_ca_gut)


ggplot(CST_df, aes(x = logDate, y = calbican_rel_abundance_gut, color = CST)) +
  geom_point(alpha = 0.6, size = 1.8) +
  geom_smooth(method = "loess", se = FALSE, size = 1) +
  scale_color_brewer(palette = "Set1") +
  labs(
    title = "Gut Candida albicans Abundance by CST",
    x = "Date",
    y = "Gut C.Albicans Relative Abundance"
  ) +
  theme_minimal()


################################################
#check who have CSTV
bacteria_abundance_merged %>%
  filter(CST == "V") %>%
  distinct(biome_id)
################################################
#case study of CA abundance overtime in ppl who has CST V
df_62 <- bacteria_abundance_merged%>%
  filter(biome_id == 62) 

df_62_long <- df_62 %>%
  select(logDate, CST, calbican_rel_abundance_vag, calbican_rel_abundance_gut) %>%
  pivot_longer(
    cols = starts_with("calbican_rel_abundance"),
    names_to = "Site",
    values_to = "Abundance"
  ) %>%
  filter(!is.na(logDate), !is.na(Abundance), is.finite(Abundance))

cst_rects62 <- df_62_long %>%
  distinct(logDate, CST)

# Plot
ggplot(df_62_long, aes(x = logDate, y = Abundance, color = Site, group = Site)) +
  geom_tile(
    data = cst_rects62,
    aes(x = logDate, y = 0, fill = CST),
    width = 0.9, height = Inf, alpha = 0.2, inherit.aes = FALSE
  ) +
  geom_point(size = 2) +
  geom_smooth(method = "loess", se = FALSE, linewidth = 1) +
  scale_color_manual(
    values = c(
      "calbican_rel_abundance_vag" = "#e41a1c",   # red-ish
      "calbican_rel_abundance_gut" = "#377eb8"    # blue-ish
    ),
    labels = c("Vagina", "Gut")
  ) +
  scale_fill_manual(
    values = c(
      "I" = "#1b9e77",
      "II" = "#d95f02",
      "III" = "#7570b3",
      "IV" = "#e7298a",
      "V" = "#66a61e"
    ),
    na.value = "white"
  ) +
  labs(
    title = "Candida albicans Abundance Over Time (Participant 62)",
    x = "Date",
    y = "Relative Abundance",
    color = "Site",
    fill = "CST"
  ) +
  theme_minimal()
###################################################################
#cross-site c.albicans abundance correlation -- no relationship seen
ggplot(bacteria_abundance_merged, aes(x = calbican_rel_abundance_gut, 
                     y = calbican_rel_abundance_vag)) +
  geom_point(alpha = 0.6, size = 2, color = "blue") +
  geom_smooth(method = "loess", se = FALSE, color = "black", linewidth = 1.2) +
  labs(
    title = "Overall Relationship Between Gut and Vaginal C.albicans Abundance",
    x = "Gut C.albicans Relative Abundance",
    y = "Vaginal C.albicans Relative Abundance"
  ) +
  theme_minimal()

###################################################################
#cross-site temporal trends
cross_long <- bacteria_abundance_merged %>%
  select(logDate, calbican_rel_abundance_gut, calbican_rel_abundance_vag, taken_antibiotics) %>%
  pivot_longer(cols = starts_with("calbican_rel_abundance"),
               names_to = "Site",
               values_to = "Abundance") %>%
  mutate(Site = recode(Site,
                       "calbican_rel_abundance_gut" = "Gut",
                       "calbican_rel_abundance_vag" = "Vagina"))


ggplot(cross_long, aes(x = logDate, y = Abundance, color = Site)) +
  geom_point(alpha = 0.5) +             # show points, maybe with some transparency
  geom_smooth(se = FALSE, method = "loess", span = 0.5) +  # smooth line
  labs(title = "C.albicans Abundance Over Time in Two Sites",
       x = "Time",
       y = "Relative Abundance",
       color = "Site") +
  theme_minimal()

###################################################################
#taking antibiotics into account
# Antibiotic group
cross_long_anti1 <- cross_long %>% filter(taken_antibiotics == 1)
# No antibiotic group
cross_long_anti0 <- cross_long %>% filter(taken_antibiotics == 0)

# Plot for antibiotic group
p1 <- ggplot(cross_long_anti1, aes(x = logDate, y = Abundance, color = Site)) +
  geom_point(alpha = 0.5) +
  geom_smooth(se = FALSE, method = "loess", span = 0.7) +
  labs(title = "Without Antibiotic Use",
       x = "Time", y = "Relative Abundance") +
  theme_minimal()

# Plot for no antibiotic group
p2 <- ggplot(cross_long_anti0, aes(x = logDate, y = Abundance, color = Site)) +
  geom_point(alpha = 0.5) +
  geom_smooth(se = FALSE, method = "loess", span = 0.7) +
  labs(title = "With Antibiotic Use",
       x = "Time", y = "Relative Abundance") +
  theme_minimal()

# Show them side-by-side
p1 + p2
###################################################################
#the ones who's taking SSRI case studies

df_60 <- bacteria_abundance_merged%>%
  filter(biome_id == 60) 

df_60_long <- df_60 %>%
  select(logDate, CST, calbican_rel_abundance_vag, calbican_rel_abundance_gut) %>%
  pivot_longer(
    cols = starts_with("calbican_rel_abundance"),
    names_to = "Site",
    values_to = "Abundance"
  ) %>%
  filter(!is.na(logDate), !is.na(Abundance), is.finite(Abundance))

cst_rects60 <- df_60_long %>%
  distinct(logDate, CST)

# Plot
ggplot(df_60_long, aes(x = logDate, y = Abundance, color = Site, group = Site)) +
  geom_tile(
    data = cst_rects60,
    aes(x = logDate, y = 0, fill = CST),
    width = 0.9, height = Inf, alpha = 0.2, inherit.aes = FALSE
  ) +
  geom_point(size = 2) +
  geom_smooth(method = "loess", se = FALSE, linewidth = 1) +
  scale_color_manual(
    values = c(
      "calbican_rel_abundance_vag" = "#e41a1c",   # red-ish
      "calbican_rel_abundance_gut" = "#377eb8"    # blue-ish
    ),
    labels = c("Gut", "Vagina")
  ) +
  scale_fill_manual(
    values = c(
      "I" = "#1b9e77",
      "II" = "#d95f02",
      "III" = "#7570b3",
      "IV" = "#e7298a",
      "V" = "#66a61e"
    ),
    na.value = "white"
  ) +
  labs(
    title = "Candida albicans Abundance Over Time (Participant 60 SSRI)",
    x = "Date",
    y = "Relative Abundance",
    color = "Site",
    fill = "CST"
  ) +
  theme_minimal()


df_44 <- bacteria_abundance_merged%>%
  filter(biome_id == 44) 
df_44_long <- df_44 %>%
  select(logDate, CST, calbican_rel_abundance_vag, calbican_rel_abundance_gut) %>%
  pivot_longer(
    cols = starts_with("calbican_rel_abundance"),
    names_to = "Site",
    values_to = "Abundance"
  ) %>%
  filter(!is.na(logDate), !is.na(Abundance), is.finite(Abundance))
cst_rects44 <- df_44_long %>%
  distinct(logDate, CST)
ggplot(df_44_long, aes(x = logDate, y = Abundance, color = Site, group = Site)) +
  geom_tile(
    data = cst_rects44,
    aes(x = logDate, y = 0, fill = CST),
    width = 0.9, height = Inf, alpha = 0.2, inherit.aes = FALSE
  ) +
  geom_point(size = 2) +
  geom_smooth(method = "loess", se = FALSE, linewidth = 1) +
  scale_color_manual(
    values = c(
      "calbican_rel_abundance_vag" = "#e41a1c",   # red-ish
      "calbican_rel_abundance_gut" = "#377eb8"    # blue-ish
    ),
    labels = c("Gut", "Vagina")
  ) +
  scale_fill_manual(
    values = c(
      "I" = "#1b9e77",
      "II" = "#d95f02",
      "III" = "#7570b3",
      "IV" = "#e7298a",
      "V" = "#66a61e"
    ),
    na.value = "white"
  ) +
  labs(
    title = "Candida albicans Abundance Over Time (Participant 44 non-SSRI)",
    x = "Date",
    y = "Relative Abundance",
    color = "Site",
    fill = "CST"
  ) +
  theme_minimal()





df_18 <- bacteria_abundance_merged%>%
  filter(biome_id == 18) 

df_18_long <- df_18 %>%
  select(logDate, CST, calbican_rel_abundance_vag, calbican_rel_abundance_gut) %>%
  pivot_longer(
    cols = starts_with("calbican_rel_abundance"),
    names_to = "Site",
    values_to = "Abundance"
  ) %>%
  filter(!is.na(logDate), !is.na(Abundance), is.finite(Abundance))

cst_rects18 <- df_18_long %>%
  distinct(logDate, CST)

# Plot
ggplot(df_18_long, aes(x = logDate, y = Abundance, color = Site, group = Site)) +
  geom_tile(
    data = cst_rects18,
    aes(x = logDate, y = 0, fill = CST),
    width = 0.9, height = Inf, alpha = 0.2, inherit.aes = FALSE
  ) +
  geom_point(size = 2) +
  geom_smooth(method = "loess", se = FALSE, linewidth = 1) +
  scale_color_manual(
    values = c(
      "calbican_rel_abundance_vag" = "#e41a1c",   # red-ish
      "calbican_rel_abundance_gut" = "#377eb8"    # blue-ish
    ),
    labels = c("Vagina", "Gut")
  ) +
  scale_fill_manual(
    values = c(
      "I" = "#1b9e77",
      "II" = "#d95f02",
      "III" = "#7570b3",
      "IV" = "#e7298a",
      "V" = "#66a61e"
    ),
    na.value = "white"
  ) +
  labs(
    title = "Candida albicans Abundance Over Time (Participant 18 adderall)",
    x = "Date",
    y = "Relative Abundance",
    color = "Site",
    fill = "CST"
  ) +
  theme_minimal()
################################################################################
#all the samples of the one's who's taking SSRI, attempt
ssri_users <- c(11,14,16,24,27,28,30,33,51,60,73)
ssridf <- bacteria_abundance_merged %>%
  mutate(
    SSRI_status = if_else(biome_id %in% ssri_users, "SSRI User", "Non-User")
  )

all_days <- seq.Date(as.Date("2022-10-13"), as.Date("2022-12-16"), by = "day")
ssridf$study_day <- match(as.Date(ssridf$logDate), all_days) -1

ssridf <- ssridf %>%
  mutate(day_c = scale(study_day, center = TRUE, scale = FALSE))

#SSRI and non-user, temporal trends of Vaginal CA
ggplot(ssridf, aes(x = logDate, y = calbican_rel_abundance_vag, color = SSRI_status)) +
  geom_jitter(width = 0.5, height = 0, alpha = 0.6, size = 2) +
  geom_smooth(se = FALSE, method = "loess") +
  labs(
    title = "Vaginal Candida albicans Abundance Over Time",
    x = "Date",
    y = "Relative Abundance",
    color = "SSRI Use"
  ) +
  theme_minimal()


# Make sure SSRI_status is a factor
ssridf$SSRI_status <- factor(ssridf$SSRI_status)

#linear mixed effect model of CA/Shannon and SSRI status
ssri_v <- lmer(calbican_rel_abundance_vag ~ SSRI_status + day_c + I(day_c^2) + (1| biome_id), data = ssridf)
summary(ssri_v)

ssri_shannonV <- lmer(Shannon_vag ~ SSRI_status + day_c + I(day_c^2) + (1| biome_id), data = ssridf)
summary(ssri_shannonV)

ssri_g <- lmer(calbican_rel_abundance_gut ~ SSRI_status + day_c + I(day_c^2) + (1| biome_id), data = ssridf)
summary(ssri_g)
r2(ssri_g)


#don't use rank sum test because each participant submit multiple samples
# wilcox.test(
#   calbican_rel_abundance_vag ~ SSRI_status,
#   data = ssridf
# )
#fluctuation of CA gut over time by SSRI
ggplot(ssridf, aes(x = logDate, y = calbican_rel_abundance_gut, color = SSRI_status)) +
  geom_jitter(width = 0.5, height = 0, alpha = 0.6, size = 2) +
  geom_smooth(se = FALSE, method = "loess") +
  labs(
    title = "Gut Candida albicans Abundance Over Time",
    x = "Date",
    y = "Relative Abundance",
    color = "SSRI Use"
  ) +
  theme_minimal()

# wilcox.test(
#   calbican_rel_abundance_gut ~ SSRI_status,
#   data = ssridf
# )


#boxplot of vaginal shannon diversity
ggplot(ssridf, aes(x = SSRI_status, y = Shannon_vag, fill = SSRI_status)) +
  geom_boxplot(alpha = 0.6, outlier.shape = NA) +  # hide default outliers
  geom_jitter(width = 0.15, size = 2, alpha = 0.7) +
  labs(
    title = "Vaginal Mycobiome Shannon Diversity by SSRI Use",
    x = "SSRI Status",
    y = "Shannon Diversity"
  ) +
  theme_minimal() +
  theme(legend.position = "none")

#sina plot
ggplot(ssridf, aes(x = SSRI_status, y = Shannon_vag, color = SSRI_status)) +
  geom_sina(alpha = 0.7, size = 2) +
  stat_summary(fun = median, geom = "crossbar", width = 0.3, color = "black") +  # show median line
  labs(
    title = "Vaginal Mycobiome Shannon Diversity by SSRI Use",
    x = "SSRI Status",
    y = "Shannon Diversity"
  ) +
  theme_minimal() +
  theme(legend.position = "none")

#fluctuation of Vaginal Shannon over time
ggplot(ssridf, aes(x = logDate, y = Shannon_vag, color = SSRI_status)) +
  geom_jitter(width = 0.5, height = 0, alpha = 0.6, size = 2) +
  geom_smooth(se = FALSE, method = "loess") +
  labs(
    title = "Vaginal Shannon Diversity Over Time",
    x = "Date",
    y = "Shannon Index",
    color = "SSRI Use"
  ) +
  theme_minimal()


#gut
#Shannon by SSRI boxplot
ggplot(ssridf, aes(x = SSRI_status, y = Shannon_gut, fill = SSRI_status)) +
  geom_boxplot(alpha = 0.6, outlier.shape = NA) +  # hide default outliers
  geom_jitter(width = 0.15, size = 2, alpha = 0.7) +
  labs(
    title = "Gut Mycobiome Shannon Diversity by SSRI Use",
    x = "SSRI Status",
    y = "Shannon Diversity"
  ) +
  theme_minimal() +
  theme(legend.position = "none")
#sina plot
ggplot(ssridf, aes(x = SSRI_status, y = Shannon_gut, color = SSRI_status)) +
  geom_sina(alpha = 0.7, size = 2) +
  stat_summary(fun = median, geom = "crossbar", width = 0.3, color = "black") +  # show median line
  labs(
    title = "Gut Mycobiome Shannon Diversity by SSRI Use",
    x = "SSRI Status",
    y = "Shannon Diversity"
  ) +
  theme_minimal() +
  theme(legend.position = "none")
#Gut Shannon fluctuation over time by SSRI
ggplot(ssridf, aes(x = logDate, y = Shannon_gut, color = SSRI_status)) +
  geom_jitter(width = 0.5, height = 0, alpha = 0.6, size = 2) +
  geom_smooth(se = FALSE, method = "loess") +
  labs(
    title = "Gut Shannon Diversity Over Time",
    x = "Date",
    y = "Shannon Index",
    color = "SSRI Use"
  ) +
  theme_minimal() 

#CST and SSRI
ggplot(ssridf, aes(x = SSRI_status, fill = CST)) +
  geom_bar(position = "fill") +  # proportion bars
  labs(title = "Distribution of CST by SSRI Use",
       x = "SSRI Status",
       y = "Proportion",
       fill = "CST") +
  theme_minimal()

ssridf_clean <- ssridf %>% 
  filter(!is.na(CST))

ssri_counts <- ssridf_clean %>%
  count(SSRI_status)

# Plot
ggplot(ssridf_clean, aes(x = SSRI_status, fill = CST)) +
  geom_bar(position = "fill") +  # proportion bars
  geom_text(data = ssri_counts,
            aes(x = SSRI_status, y = 1.05, label = paste0("n = ", n)),
            inherit.aes = FALSE,
            size = 4.5) +
  labs(title = "Distribution of CST by SSRI Use",
       x = "SSRI Status",
       y = "Proportion",
       fill = "CST") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "right"
  )


#CST and gut shannon
ggplot(ssridf_clean, aes(x = CST, y = Shannon_gut, fill = CST)) +
  geom_boxplot(alpha = 0.6, outlier.shape = NA) +  # Hide default outliers
  geom_jitter(width = 0.2, size = 2, alpha = 0.7, color = "black") +  # Overlay individual points
  scale_fill_brewer(palette = "Set3") +  # Or choose another palette
  labs(
    title = "Shannon Diversity in Gut Mycobiome by CST",
    x = "CST",
    y = "Shannon Diversity Index"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 30, hjust = 1)
  )

#CST and vaginal shannon
ggplot(ssridf_clean, aes(x = CST, y = Shannon_vag, fill = CST)) +
  geom_boxplot(alpha = 0.6, outlier.shape = NA) +  # Hide default outliers
  geom_jitter(width = 0.2, size = 2, alpha = 0.7, color = "black") +  # Overlay individual points
  scale_fill_brewer(palette = "Set3") +  # Or choose another palette
  labs(
    title = "Shannon Diversity in Vaginal Mycobiome by CST",
    x = "CST",
    y = "Shannon Diversity Index"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 30, hjust = 1)
  )

ggplot(ssridf_clean, aes(x = CST, y = Shannon_gut, fill = CST)) +
  geom_boxplot(alpha = 0.6, outlier.shape = NA) +  # Hide default outliers
  geom_jitter(width = 0.2, size = 2, alpha = 0.7, color = "black") +  # Overlay individual points
  scale_fill_brewer(palette = "Set3") +  # Or choose another palette
  labs(
    title = "Shannon Diversity in Gut Mycobiome by CST",
    x = "CST",
    y = "Shannon Diversity Index"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 30, hjust = 1)
  )




#############################################################
#case studies of all participants in a glance
ssridf <- ssridf %>%
  arrange(biome_id, logDate)

#tried to put everyone onto a panel but too many plots
# ggplot(ssridf, aes(x = logDate, y = calbican_rel_abundance_vag, group = biome_id, color = SSRI_status)) +
#   geom_smooth(se = FALSE, method = "loess") +
#   ylim(0,1) +
#   geom_point(size = 1.5, alpha = 0.7) +
#   facet_wrap(~ biome_id, scales = "free_x") +  # One panel per participant
#   scale_color_manual(values = c("Non-User" = "#1b9e77", "SSRI User" = "#d95f02")) +  # customize colors
#   labs(
#     title = "Candida albicans Vaginal Abundance Over Time",
#     x = "Date",
#     y = "Relative Abundance",
#     color = "SSRI Status"
#   ) +
#   theme_minimal(base_size = 12) +
#   theme(
#     strip.text = element_text(size = 7),
#     axis.text.x = element_text(angle = 45, hjust = 1),
#     legend.position = "bottom"
#   )

#panel of vag CA for each participant
example_ids <- unique(ssridf$biome_id)[51:76] #start from [1:24], then [25:50]

ssridf %>%
  filter(biome_id %in% example_ids) %>%
  ggplot(aes(x = logDate, y = calbican_rel_abundance_vag, group = biome_id, color = SSRI_status)) +
  geom_smooth(se = FALSE, method = "loess") +
  ylim(0,1) +
  geom_point(size = 1.5) +
  facet_wrap(~ biome_id, scales = "free_x", ncol = 4) +
  theme_minimal()

#panel of gut CA for each participant
ssridf %>%
  filter(biome_id %in% example_ids) %>%
  ggplot(aes(x = logDate, y = calbican_rel_abundance_gut, group = biome_id, color = SSRI_status)) +
  geom_smooth(se = FALSE, method = "loess") +
  ylim(0,1) +
  geom_point(size = 1.5) +
  facet_wrap(~ biome_id, scales = "free_x", ncol = 4) +
  theme_minimal()



#############################################################
#joining DASS into ssridf

cleanDass <- read.csv("/Users/caoyang/Desktop/Tetel Lab/datasets/cleaned_dass.csv")

ssridf <- ssridf %>%
  mutate(logDate = as.Date(logDate))

cleanDass <- cleanDass %>%
  mutate(Timestamp = as.Date(Timestamp))
# Step 1: Create all combinations of samples and mood data with same participant
combo <- ssridf %>%
  mutate(biome_id = as.character(biome_id)) %>%
  inner_join(
    cleanDass %>% mutate(biome_id = as.character(biome_id)),
    by = "biome_id"
  ) %>%
  mutate(
    date_diff = as.numeric(difftime(logDate, Timestamp, units = "days")),
    abs_diff = abs(date_diff)
  ) %>%
  #Step2: keep only those within ±7 days
  filter(abs_diff <= 7)

# Step 3: For each sample (logDate + biome_id), keep only the mood survey with smallest time difference
closest_match <- combo %>%
  group_by(biome_id, logDate) %>%
  slice_min(abs_diff, with_ties = FALSE) %>%
  ungroup()

# Step 4: Select relevant mood columns to join back
mood_matched <- closest_match %>%
  select(biome_id, logDate, depression_score, anxiety_score, stress_score, depressionseverity, anxietyseverity, stressseverity)

# Step 5: Left join to original sample data — unmatched samples get NA
ssridf_with_mood <- ssridf %>%
  mutate(biome_id = as.character(biome_id)) %>%
  left_join(
    mood_matched %>% mutate(biome_id = as.character(biome_id)),
    by = c("biome_id", "logDate")
  )
################################################################
#SSRI and DASS and CA analysis plots
#plotting dass and ssri vaginal CA
ssridf_with_mood <- ssridf_with_mood %>%
  mutate(logDate = as.Date(logDate)) %>%
  filter(!is.na(depression_score), !is.na(anxiety_score), !is.na(stress_score))

ssridf_with_mood <- ssridf_with_mood %>%
  mutate(
    depressionseverity = factor(depressionseverity,
                                 levels = 0:4,
                                 labels = c("Normal", "Mild", "Moderate", "Severe", "Extremely Severe")
    )
  )

ssridf_with_mood <- ssridf_with_mood %>%
  mutate(
    anxietyseverity = factor(anxietyseverity,
                                levels = 0:4,
                                labels = c("Normal", "Mild", "Moderate", "Severe", "Extremely Severe")
    )
  )

ssridf_with_mood <- ssridf_with_mood %>%
  mutate(
    stressseverity = factor(stressseverity,
                             levels = 0:4,
                             labels = c("Normal", "Mild", "Moderate", "Severe", "Extremely Severe")
    )
  )

#how dass score depression/anxiety/stress fluctuate over time by SSRI use
plot_mood <- function(data, mood_var) {
  # mood_var should be one of "depression_score", "anxiety_score", or "stress_score"
  data %>%
    select(biome_id, logDate, SSRI_status, all_of(mood_var)) %>%
    rename(Score = all_of(mood_var)) %>%
    ggplot(aes(x = logDate, y = Score, color = SSRI_status)) +
    geom_point(alpha = 0.6, size = 1.5) +
    geom_smooth(method = "loess", se = FALSE, span = 0.75) +
    labs(
      title = paste0(gsub("_", " ", tools::toTitleCase(mood_var)), " Over Time by SSRI Use"),
      x = "Date",
      y = "Score",
      color = "SSRI Status"
    ) +
    scale_color_manual(values = c("lightsalmon", "cornflowerblue")) +
    theme_minimal(base_size = 13) +
    theme(legend.position = "top")
}


plot_mood(ssridf_with_mood, "depression_score")
plot_mood(ssridf_with_mood, "anxiety_score")
plot_mood(ssridf_with_mood, "stress_score")

#boxplots comparing vaginal CA by each category of DASS severity
ggplot(ssridf_with_mood, aes(x = depressionseverity, y = calbican_rel_abundance_vag, fill = SSRI_status)) +
  geom_boxplot(
    position = position_dodge(width = 0.75),
    outlier.shape = NA,
    alpha = 0.6,
    width = 0.5
  ) +
  geom_jitter(
    aes(color = SSRI_status),
    position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.75),
    size = 1.5,
    alpha = 0.6
  ) +
  labs(
    title = "Vaginal Candida albicans Abundance by Depression Severity and SSRI Use",
    x = "Depression Severity",
    y = "Relative Abundance of C. albicans",
    fill = "SSRI Status",
    color = "SSRI Status"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1),
    plot.title = element_text(face = "bold")
  )

ggplot(ssridf_with_mood, aes(x = anxietyseverity, y = calbican_rel_abundance_vag, fill = SSRI_status)) +
  geom_boxplot(
    position = position_dodge(width = 0.75),
    outlier.shape = NA,
    alpha = 0.6,
    width = 0.5
  ) +
  geom_jitter(
    aes(color = SSRI_status),
    position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.75),
    size = 1.5,
    alpha = 0.6
  ) +
  labs(
    title = "Vaginal Candida albicans Abundance by Anxiety Severity and SSRI Use",
    x = "Anxiety Severity",
    y = "Relative Abundance of C. albicans",
    fill = "SSRI Status",
    color = "SSRI Status"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1),
    plot.title = element_text(face = "bold")
  )

ggplot(ssridf_with_mood, aes(x = stressseverity, y = calbican_rel_abundance_vag, fill = SSRI_status)) +
  geom_boxplot(
    position = position_dodge(width = 0.75),
    outlier.shape = NA,
    alpha = 0.6,
    width = 0.5
  ) +
  geom_jitter(
    aes(color = SSRI_status),
    position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.75),
    size = 1.5,
    alpha = 0.6
  ) +
  labs(
    title = "Vaginal Candida albicans Abundance by Stress Severity and SSRI Use",
    x = "Stress Severity",
    y = "Relative Abundance of C. albicans",
    fill = "SSRI Status",
    color = "SSRI Status"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1),
    plot.title = element_text(face = "bold")
  )

#vaginal CA and DASS
ggplot(ssridf_with_mood, aes(x = depression_score, y = calbican_rel_abundance_vag, color = SSRI_status)) +
  geom_jitter(width = 0.5, height = 0, alpha = 0.6, size = 2) +
  geom_smooth(se = FALSE, method = "loess") +
  labs(
    title = "Vaginal Candida albicans Abundance and Depression Score",
    x = "Depression Score",
    y = "Relative Abundance",
    color = "SSRI Use"
  ) +
  theme_minimal()

ggplot(ssridf_with_mood, aes(x = anxiety_score, y = calbican_rel_abundance_vag, color = SSRI_status)) +
  geom_jitter(width = 0.5, height = 0, alpha = 0.6, size = 2) +
  geom_smooth(se = FALSE, method = "loess") +
  labs(
    title = "Vaginal Candida albicans Abundance and Anxiety Score",
    x = "Anxiety Score",
    y = "Relative Abundance",
    color = "SSRI Use"
  ) +
  theme_minimal()

ggplot(ssridf_with_mood, aes(x = stress_score, y = calbican_rel_abundance_vag, color = SSRI_status)) +
  geom_jitter(width = 0.5, height = 0, alpha = 0.6, size = 2) +
  geom_smooth(se = FALSE, method = "loess") +
  labs(
    title = "Vaginal Candida albicans Abundance and Stress Score",
    x = "Stress Score",
    y = "Relative Abundance",
    color = "SSRI Use"
  ) +
  theme_minimal()



#########################################################################
#gut CA and DASS
ggplot(ssridf_with_mood, aes(x = depression_score, y = calbican_rel_abundance_gut, color = SSRI_status)) +
  geom_jitter(width = 0.5, height = 0, alpha = 0.6, size = 2) +
  geom_smooth(se = FALSE, method = "loess") +
  labs(
    title = "Gut Candida albicans Abundance and Depression Score",
    x = "Depression Score",
    y = "Relative Abundance",
    color = "SSRI Use"
  ) +
  theme_minimal()

ggplot(ssridf_with_mood, aes(x = anxiety_score, y = calbican_rel_abundance_gut, color = SSRI_status)) +
  geom_jitter(width = 0.5, height = 0, alpha = 0.6, size = 2) +
  geom_smooth(se = FALSE, method = "loess") +
  labs(
    title = "Gut Candida albicans Abundance and Anxiety Score",
    x = "Anxiety Score",
    y = "Relative Abundance",
    color = "SSRI Use"
  ) +
  theme_minimal()

ggplot(ssridf_with_mood, aes(x = stress_score, y = calbican_rel_abundance_gut, color = SSRI_status)) +
  geom_jitter(width = 0.5, height = 0, alpha = 0.6, size = 2) +
  geom_smooth(se = FALSE, method = "loess") +
  labs(
    title = "Gut Candida albicans Abundance and Stress Score",
    x = "Stress Score",
    y = "Relative Abundance",
    color = "SSRI Use"
  ) +
  theme_minimal()


ggplot(ssridf_with_mood, aes(x = depressionseverity, y = calbican_rel_abundance_gut, fill = SSRI_status)) +
  geom_boxplot(
    position = position_dodge(width = 0.75),
    outlier.shape = NA,
    alpha = 0.6,
    width = 0.5
  ) +
  geom_jitter(
    aes(color = SSRI_status),
    position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.75),
    size = 1.5,
    alpha = 0.6
  ) +
  labs(
    title = "Gut Candida albicans Abundance by Depression Severity and SSRI Use",
    x = "Depression Severity",
    y = "Relative Abundance of C. albicans",
    fill = "SSRI Status",
    color = "SSRI Status"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1),
    plot.title = element_text(face = "bold")
  )


ggplot(ssridf_with_mood, aes(x = anxietyseverity, y = calbican_rel_abundance_gut, fill = SSRI_status)) +
  geom_boxplot(
    position = position_dodge(width = 0.75),
    outlier.shape = NA,
    alpha = 0.6,
    width = 0.5
  ) +
  geom_jitter(
    aes(color = SSRI_status),
    position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.75),
    size = 1.5,
    alpha = 0.6
  ) +
  labs(
    title = "Gut Candida albicans Abundance by Anxiety Severity and SSRI Use",
    x = "Depression Severity",
    y = "Relative Abundance of C. albicans",
    fill = "SSRI Status",
    color = "SSRI Status"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1),
    plot.title = element_text(face = "bold")
  )

ggplot(ssridf_with_mood, aes(x = stressseverity, y = calbican_rel_abundance_gut, fill = SSRI_status)) +
  geom_boxplot(
    position = position_dodge(width = 0.75),
    outlier.shape = NA,
    alpha = 0.6,
    width = 0.5
  ) +
  geom_jitter(
    aes(color = SSRI_status),
    position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.75),
    size = 1.5,
    alpha = 0.6
  ) +
  labs(
    title = "Gut Candida albicans Abundance by Stress Severity and SSRI Use",
    x = "Stress Severity",
    y = "Relative Abundance of C. albicans",
    fill = "SSRI Status",
    color = "SSRI Status"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1),
    plot.title = element_text(face = "bold")
  )
################################################################################
#ssri and gut most abundant species
#stacked barplot of dominant species in gut

plot_df <- ssridf_with_mood %>%
  mutate(
    species_group = case_when(
      str_detect(dominant_species_fungal_gut, "Candida_albicans") ~ "Candida_albicans",
      str_detect(dominant_species_fungal_gut, "globosa") ~ "globosa",
      str_detect(dominant_species_fungal_gut, "restricta") ~ "restricta",
      str_detect(dominant_species_fungal_gut, "Malassezia_globosa") ~ "Malassezia_globosa",
      str_detect(dominant_species_fungal_gut, "Candida_parapsilosis") ~ "Candida_parapsilosis",
      str_detect(dominant_species_fungal_gut, "Rhodotorula_mucilaginosa") ~ " Rhodotorula_mucilaginosa",
      TRUE ~ "Other"
    )
  )

species_freq <- plot_df %>%
  group_by(SSRI_status, species_group) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(SSRI_status) %>%
  mutate(percent = count / sum(count) * 100)

ggplot(species_freq, aes(x = SSRI_status, y = percent, fill = species_group)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(
    title = "Distribution of Dominant Gut Microbiome Species by SSRI (samples)",
    x = "SSRI status",
    y = "Percentage of Samples",
    fill = "Gut Species"
  ) +
  theme_minimal() +
  scale_fill_manual(
    values = c(
      "Candida_albicans" = "#1f77b4",
      "globosa" = "#ff7f0e",
      "Malassezia_globosa" = "#ff7f0e",
      "restricta" = "#d62728",
      "Candida_parapsilosis" = "#9467bd",
      "Rhodotorula_mucilaginosa" = "#8c564b",
      "Other" = "gray80"
    )
  )


#########################################################################
#PCA
#PCA of Vag CA by SSRI status
sample_df <- data.frame(sample_data(fungal2.0))
# Create the new column
sample_df$SSRI_status <- ifelse(sample_df$biome_id %in% ssri_users, "SSRI", "non-SSRI")
# Assign it back to fungal2.0
sample_data(fungal2.0) <- sample_data(sample_df)



#all samples
ord_plot(fungal2.0 %>% 
           tax_transform("clr") %>%
           ord_calc("PCA"), 
         color = "SSRI_status", size = 2) +
  theme_minimal()
#vaginal samples
fungal_vaginal <- subset_samples(fungal2.0, sampleType == "vaginal")
fungal_vaginal <- fungal_vaginal %>% tax_transform("clr")
fungal_vaginal <- fungal_vaginal %>% ord_calc(method = "PCA")
ord_plot(fungal_vaginal, color = "SSRI_status", size = 2) +
  theme_minimal()
#gut samples
fungal_gut <- subset_samples(fungal2.0, sampleType == "fecal")
fungal_gut <- fungal_gut %>% tax_transform("clr")
fungal_gut <- fungal_gut %>% ord_calc(method = "PCA")
ord_plot(fungal_gut, color = "SSRI_status", size = 2) +
  stat_ellipse(aes(color = SSRI_status), type = "norm", size = 0.6, linetype = "dashed")+
  theme_minimal()


#######################################################################
#general exploration of difference between gut and vaginal fungus
ord_plot(fungal2.0 %>% 
           tax_transform("clr") %>%
           ord_calc("PCA"), 
         color = "sampleType", size = 2) +
  theme_minimal()
#############################################################
#everyone's fluctuation over time
ssri_ids <- unique(ssridf$biome_id[ssridf$SSRI_status == "SSRI User"])
nonssri_ids <- unique(ssridf$biome_id[ssridf$SSRI_status == "Non-User"])

# Generate distinct colors for each group using base R colorRampPalette
warm_palette <- colorRampPalette(c("#E69F00", "#F0E442")) # oranges/yellows
cold_palette <- colorRampPalette(c("#56B4E9", "#0072B2")) # blues

warm_colors <- warm_palette(length(nonssri_ids))
cold_colors <- cold_palette(length(ssri_ids))

# Map colors to biome_ids
biome_colors <- c(setNames(cold_colors, ssri_ids), setNames(warm_colors, nonssri_ids))

# Add color column to dataframe by matching biome_id
ssridf <- ssridf %>%
  mutate(color = biome_colors[as.character(biome_id)])

ggplot(ssridf, aes(x = logDate, y = calbican_rel_abundance_vag, group = biome_id)) +
  geom_smooth(aes(color = color), method = "loess", se = FALSE, alpha = 0.1, size = 0.5) +
  geom_point(aes(color = color), alpha = 0.9, size = 2) +
  scale_color_identity() +
  scale_y_continuous(limits = c(0, 1)) +  # set y-axis from 0 to 1
  theme_minimal() +
  labs(
    title = "Vaginal C.albicans Abundance Over Time by Biome ID and SSRI Status",
    x = "Log Date",
    y = "C. albicans Abundance (vaginal)"
  )

ggplot(ssridf, aes(x = logDate, y = calbican_rel_abundance_gut, group = biome_id)) +
  #geom_smooth(aes(color = color), method = "loess", se = FALSE, alpha = 0.1, size = 0.5) +
  geom_point(aes(color = color), alpha = 0.9, size = 2) +  # semi-transparent points
  scale_color_identity() +
  scale_y_continuous(limits = c(0, 1)) +
  theme_minimal() +
  labs(
    title = "Gut C. albicans Abundance Over Time by Biome ID and SSRI Status",
    x = "Log Date",
    y = "C. albicans Abundance (gut)"
  )

#################################################################################
#CA and L.iners
Liners_CA_vag <- ssridf %>%
  select(logDate, SSRI_status, calbican_rel_abundance_vag, Liners_rel_abundance_vag) %>%
  pivot_longer(
    cols = c(calbican_rel_abundance_vag, Liners_rel_abundance_vag),
    names_to = "Species",
    values_to = "Abundance"
  ) %>%
  mutate(
    Species = recode(Species,
                     calbican_rel_abundance_vag = "C. albicans",
                     Liners_rel_abundance_vag = "L. iners")
  )

ggplot(Liners_CA_vag, aes(x = logDate, y = Abundance, color = Species, linetype = SSRI_status)) +
  geom_smooth(method = "loess", se = FALSE, size = 1.2) +
  geom_point(alpha = 0.4, size = 1.0) + 
  scale_color_manual(values = c("C. albicans" = "purple", "L. iners" = "darkgreen")) +
  scale_linetype_manual(values = c("SSRI User" = "solid", "Non-User" = "dashed")) +
  labs(
    x = "Date",
    y = "Relative Abundance",
    color = "Species",
    linetype = "SSRI Status",
    title = "Vaginal C. albicans and L. iners Abundance Over Time by SSRI Status"
  ) +
  theme_minimal()

Liners_CA_gut <- ssridf %>%
  select(logDate, SSRI_status, calbican_rel_abundance_gut, Liners_rel_abundance_gut) %>%
  pivot_longer(
    cols = c(calbican_rel_abundance_gut, Liners_rel_abundance_gut),
    names_to = "Species",
    values_to = "Abundance"
  ) %>%
  mutate(
    Species = recode(Species,
                     calbican_rel_abundance_gut = "C. albicans",
                     Liners_rel_abundance_gut = "L. iners")
  )

ggplot(Liners_CA_gut, aes(x = logDate, y = Abundance, color = Species, linetype = SSRI_status)) +
  geom_smooth(method = "loess", se = FALSE, size = 1.2) +
  geom_point(alpha = 0.4, size = 1.0) + 
  scale_color_manual(values = c("C. albicans" = "purple", "L. iners" = "darkgreen")) +
  scale_linetype_manual(values = c("SSRI User" = "solid", "Non-User" = "dashed")) +
  labs(
    x = "Date",
    y = "Relative Abundance",
    color = "Species",
    linetype = "SSRI Status",
    title = "Gut C. albicans and L. iners Abundance Over Time by SSRI Status"
  ) +
  theme_minimal()

lgasseri_CA_vag <- ssridf %>%
  select(logDate, SSRI_status, calbican_rel_abundance_vag, lgasseri_rel_abundance_vag) %>%
  pivot_longer(
    cols = c(calbican_rel_abundance_vag, lgasseri_rel_abundance_vag),
    names_to = "Species",
    values_to = "Abundance"
  ) %>%
  mutate(
    Species = recode(Species,
                     calbican_rel_abundance_vag = "C. albicans",
                     lgasseri_rel_abundance_vag = "L. gasseri")
  )

ggplot(lgasseri_CA_vag, aes(x = logDate, y = Abundance, color = Species, linetype = SSRI_status)) +
  geom_smooth(method = "loess", se = FALSE, size = 1.2) +
  geom_point(alpha = 0.4, size = 1.0) + 
  scale_color_manual(values = c("C. albicans" = "purple", "L. gasseri" = "darkgreen")) +
  scale_linetype_manual(values = c("SSRI User" = "solid", "Non-User" = "dashed")) +
  labs(
    x = "Date",
    y = "Relative Abundance",
    color = "Species",
    linetype = "SSRI Status",
    title = "Vaginal C. albicans and L. gasseri Abundance Over Time by SSRI Status"
  ) +
  theme_minimal()

lgasseri_CA_gut <- ssridf %>%
  select(logDate, SSRI_status, calbican_rel_abundance_gut, lgasseri_rel_abundance_gut) %>%
  pivot_longer(
    cols = c(calbican_rel_abundance_gut, lgasseri_rel_abundance_gut),
    names_to = "Species",
    values_to = "Abundance"
  ) %>%
  mutate(
    Species = recode(Species,
                     calbican_rel_abundance_gut = "C. albicans",
                     lgasseri_rel_abundance_gut = "L. gasseri")
  )

ggplot(lgasseri_CA_gut, aes(x = logDate, y = Abundance, color = Species, linetype = SSRI_status)) +
  geom_smooth(method = "loess", se = FALSE, size = 1.2) +
  geom_point(alpha = 0.4, size = 1.0) + 
  scale_color_manual(values = c("C. albicans" = "purple", "L. gasseri" = "darkgreen")) +
  scale_linetype_manual(values = c("SSRI User" = "solid", "Non-User" = "dashed")) +
  labs(
    x = "Date",
    y = "Relative Abundance",
    color = "Species",
    linetype = "SSRI Status",
    title = "Gut C. albicans and L. gasseri Abundance Over Time by SSRI Status"
  ) +
  theme_minimal()


lcrispatus_CA_vag <- ssridf %>%
  select(logDate, SSRI_status, calbican_rel_abundance_vag, lcrispatus_rel_abundance_vag) %>%
  pivot_longer(
    cols = c(calbican_rel_abundance_vag, lcrispatus_rel_abundance_vag),
    names_to = "Species",
    values_to = "Abundance"
  ) %>%
  mutate(
    Species = recode(Species,
                     calbican_rel_abundance_vag = "C. albicans",
                     lcrispatus_rel_abundance_vag = "L. crispatus")
  )

ggplot(lcrispatus_CA_vag, aes(x = logDate, y = Abundance, color = Species, linetype = SSRI_status)) +
  geom_smooth(method = "loess", se = FALSE, size = 1.2) +
  geom_point(alpha = 0.4, size = 1.0) + 
  scale_color_manual(values = c("C. albicans" = "purple", "L. crispatus" = "darkgreen")) +
  scale_linetype_manual(values = c("SSRI User" = "solid", "Non-User" = "dashed")) +
  labs(
    x = "Date",
    y = "Relative Abundance",
    color = "Species",
    linetype = "SSRI Status",
    title = "Vaginal C. albicans and L. crispatus Abundance Over Time by SSRI Status"
  ) +
  theme_minimal()

ljensenii_CA_vag <- ssridf %>%
  select(logDate, SSRI_status, calbican_rel_abundance_vag, ljensenii_rel_abundance_vag) %>%
  pivot_longer(
    cols = c(calbican_rel_abundance_vag, ljensenii_rel_abundance_vag),
    names_to = "Species",
    values_to = "Abundance"
  ) %>%
  mutate(
    Species = recode(Species,
                     calbican_rel_abundance_vag = "C. albicans",
                     ljensenii_rel_abundance_vag = "L. jensenii")
  )


ggplot(ljensenii_CA_vag, aes(x = logDate, y = Abundance, color = Species, linetype = SSRI_status)) +
  geom_smooth(method = "loess", se = FALSE, size = 1.2) +
  geom_point(alpha = 0.4, size = 1.0) + 
  scale_color_manual(values = c("C. albicans" = "purple", "L. crispatus" = "darkgreen")) +
  scale_linetype_manual(values = c("SSRI User" = "solid", "Non-User" = "dashed")) +
  labs(
    x = "Date",
    y = "Relative Abundance",
    color = "Species",
    linetype = "SSRI Status",
    title = "Vaginal C. albicans and L.jensenii Abundance Over Time by SSRI Status"
  ) +
  theme_minimal()

#####################################################################################
#Menses and SSRI and CA and Depression
menses_df <- read.csv("/Users/caoyang/Desktop/Tetel Lab/datasets/imputed_menstruation_data_3_11.csv")
menses_df <- menses_df %>% 
  rename_with(~gsub("X2022.", "2022.", .), starts_with("X2022.")) %>% 
  rename_with(~gsub("\\.", "-", .))

menses_df_long <- menses_df %>% 
  pivot_longer(cols=starts_with("2022-"), names_to="logDate", values_to="menses_status")

menses_df_long$logDate <- as.Date(menses_df_long$logDate)
menses_df_long$biome_id <- as.character(menses_df_long$biome_id)

ssridf_mood_menses <- ssridf_with_mood %>% 
  left_join(menses_df_long, by=c("biome_id", "logDate"))

base_names <- gsub("\\.x$|\\.y$", "", colnames(ssridf_mood_menses))
duplicated_names <- base_names[duplicated(base_names)]
duplicated_names

colnames(ssridf_mood_menses)[endsWith(colnames(ssridf_mood_menses), ".y")]
ssridf_mood_menses <- ssridf_mood_menses %>%
  select(-ends_with(".y"))
names(ssridf_mood_menses) <- gsub("\\.x$", "", names(ssridf_mood_menses))

ssridf_mood_menses <- ssridf_mood_menses %>% 
  mutate(menses_day = ifelse(menses_status %in% c(1,2,3,7,9,78), "menses", 
                             ifelse(menses_status %in% c(4,5,6,10), "not_menses", NA)))


#boxplot
boxplot_menses_df <- ssridf_mood_menses %>%
  filter(!is.na(menses_day)) %>%  # Remove NA menses upfront
  mutate(
    menses_label = case_when(
      menses_day == "menses" ~ "On Menses",
      menses_day == "not_menses" ~ "Not on Menses"
    ),
    group = paste(SSRI_status, menses_label, sep = " - ")
  )


ggplot(boxplot_menses_df, aes(x = depressionseverity, y = calbican_rel_abundance_vag, fill = group)) +
  geom_boxplot(position = position_dodge(width = 0.8), outlier.shape = NA, width = 0.7, alpha = 0.5) +
  geom_jitter(aes(color = group), 
              position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8),
              size = 1.5, alpha = 0.6) +
  labs(
    x = "Depression Severity",
    y = "Vaginal C. albicans Relative Abundance",
    fill = "Group",
    color = "Group"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))
#violin
ggplot(boxplot_menses_df, aes(x = depressionseverity, y = calbican_rel_abundance_vag, fill = group)) +
  geom_violin(position = position_dodge(width = 0.8), trim = FALSE, alpha = 0.7) +
  geom_jitter(aes(color = group),
              position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8),
              size = 1.5, alpha = 0.6) +
  labs(
    x = "Depression Severity",
    y = "Vaginal C. albicans Relative Abundance",
    fill = "Group",
    color = "Group"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

#################gut

ggplot(boxplot_menses_df, aes(x = depressionseverity, y = calbican_rel_abundance_gut, fill = group)) +
  geom_boxplot(position = position_dodge(width = 0.8), outlier.shape = NA, width = 0.7, alpha = 0.5) +
  geom_jitter(aes(color = group), 
              position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8),
              size = 1.5, alpha = 0.6) +
  labs(
    x = "Depression Severity",
    y = "Gut C. albicans Relative Abundance",
    fill = "Group",
    color = "Group"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))



#with only SSRI and menses
ggplot(boxplot_menses_df, aes(x = group, y = calbican_rel_abundance_vag, fill = group)) +
  geom_boxplot(outlier.shape = NA, width = 0.7, alpha = 0.5) +
  geom_jitter(aes(color = group),
              position = position_jitter(width = 0.2),
              size = 1.5, alpha = 0.6) +
  labs(
    x = "Group",
    y = "Vaginal C. albicans Relative Abundance",
    fill = "Group",
    color = "Group"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 20, hjust = 1))


#gut
ggplot(boxplot_menses_df, aes(x = group, y = calbican_rel_abundance_gut, fill = group)) +
  geom_boxplot(outlier.shape = NA, width = 0.7, alpha = 0.5) +
  geom_jitter(aes(color = group),
              position = position_jitter(width = 0.2),
              size = 1.5, alpha = 0.6) +
  labs(
    x = "Group",
    y = "Gut C. albicans Relative Abundance",
    fill = "Group",
    color = "Group"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 20, hjust = 1))

################################
#vaginal and gut lactobacillus
ggplot(boxplot_menses_df, aes(x = group, y = lacto_rel_abundance_vag, fill = group)) +
  geom_boxplot(outlier.shape = NA, width = 0.7, alpha = 0.5) +
  geom_jitter(aes(color = group),
              position = position_jitter(width = 0.2),
              size = 1.5, alpha = 0.6) +
  labs(
    x = "Group",
    y = "Vaginal Lactobacillus Relative Abundance",
    fill = "Group",
    color = "Group"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 20, hjust = 1))


ggplot(boxplot_menses_df, aes(x = depressionseverity, y = lacto_rel_abundance_vag, fill = group)) +
  geom_boxplot(position = position_dodge(width = 0.8), outlier.shape = NA, width = 0.7, alpha = 0.5) +
  geom_jitter(aes(color = group), 
              position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8),
              size = 1.5, alpha = 0.6) +
  labs(
    x = "Depression Severity",
    y = "Vaginal Lactobaciilus Relative Abundance",
    fill = "Group",
    color = "Group"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

ggplot(boxplot_menses_df, aes(x = group, y = lacto_rel_abundance_vag, fill = group)) +
  geom_boxplot(outlier.shape = NA, width = 0.7, alpha = 0.5) +
  geom_jitter(aes(color = group),
              position = position_jitter(width = 0.2),
              size = 1.5, alpha = 0.6) +
  labs(
    x = "Group",
    y = "Vaginal Lactobacillus Relative Abundance",
    fill = "Group",
    color = "Group"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 20, hjust = 1))

ggplot(boxplot_menses_df, aes(x = group, y = Shannon_vag, fill = group)) +
  geom_boxplot(outlier.shape = NA, width = 0.7, alpha = 0.5) +
  geom_jitter(aes(color = group),
              position = position_jitter(width = 0.2),
              size = 1.5, alpha = 0.6) +
  labs(
    x = "Group",
    y = "Vaginal Fungal Shannon Diversity",
    fill = "Group",
    color = "Group"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 20, hjust = 1))

ggplot(boxplot_menses_df, aes(x = group, y = Shannon_gut, fill = group)) +
  geom_boxplot(outlier.shape = NA, width = 0.7, alpha = 0.5) +
  geom_jitter(aes(color = group),
              position = position_jitter(width = 0.2),
              size = 1.5, alpha = 0.6) +
  labs(
    x = "Group",
    y = "Gut Fungal Shannon Diversity",
    fill = "Group",
    color = "Group"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 20, hjust = 1))

#Linear Mixed Effect Model
vagCA_SSRI_Menses <- lmer(calbican_rel_abundance_vag ~ SSRI_status * menses_day + (1 | biome_id) , data = ssridf_mood_menses)
summary(vagCA_SSRI_Menses)

vagCA_SSRI <- lmer(calbican_rel_abundance_vag ~ SSRI_status + (1 | biome_id) , data = ssridf_mood_menses)
summary(vagCA_SSRI)

gutCA_SSRI_Menses <- lmer(calbican_rel_abundance_gut ~ SSRI_status * menses_day + (1 | biome_id) , data = ssridf_mood_menses)
summary(gutCA_SSRI_Menses)

gutCA_SSRI <- lmer(calbican_rel_abundance_gut ~ SSRI_status + (1 | biome_id) , data = ssridf_mood_menses)
summary(gutCA_SSRI)

vagLacto_SSRI_Menses <- lmer(lacto_rel_abundance_vag ~ SSRI_status * menses_day + (1 | biome_id) , data = ssridf_mood_menses)
summary(vagLacto_SSRI_Menses)

vagLacto_SSRI <- lmer(lacto_rel_abundance_vag ~ SSRI_status + (1 | biome_id) , data = ssridf_mood_menses)
summary(vagLacto_SSRI)

######################################################################
#with depression
#creating a new dataset where na are already dropped
ssridf_lacto_mm_clean <- ssridf_mood_menses %>%
  dplyr::select(lacto_rel_abundance_vag, SSRI_status, depression_score, menses_day, biome_id) %>%
  tidyr::drop_na()

vagCA_SSRI_Menses_d <- lmer(calbican_rel_abundance_vag ~ SSRI_status * menses_day * depression_score + (1 | biome_id) , data = ssridf_mm_clean)
summary(vagCA_SSRI_Menses_d)
r2(vagCA_SSRI_Menses_d)

gutCA_SSRI_Menses_d <- lmer(calbican_rel_abundance_gut ~ SSRI_status * menses_day * depression_score + (1 | biome_id) , data = ssridf_mood_menses)
summary(gutCA_SSRI_Menses_d)
r2(gutCA_SSRI_Menses_d)

vagLacto_SSRI_Menses_d <- lmer(lacto_rel_abundance_vag ~ SSRI_status + menses_day + depression_score +
                                 SSRI_status:menses_day + SSRI_status:depression_score + menses_day:depression_score +
                                 (1 | biome_id), data = ssridf_lacto_mm_clean) #took out the three-way interaction
summary(vagLacto_SSRI_Menses_d)
BIC(vagLacto_SSRI_Menses_d)
r2(vagLacto_SSRI_Menses_d)

withoutmenses <- lmer(lacto_rel_abundance_vag ~ SSRI_status * depression_score + (1 | biome_id) , data = ssridf_lacto_mm_clean)
BIC(withoutmenses)
anova(vagLacto_SSRI_Menses_d, withoutmenses)
##################################################################################
#full model actual comparisons for poster
model_lacto_mm_full <- lmer(lacto_rel_abundance_vag ~ SSRI_status + menses_day + depression_score +
                     SSRI_status:menses_day + SSRI_status:depression_score + menses_day:depression_score +
                     (1 | biome_id), data = ssridf_lacto_mm_clean, REML = FALSE)
summary(model_lacto_mm_full)

model_lacto_mm_noSSRI <- lmer(lacto_rel_abundance_vag ~ menses_day + depression_score +
                       menses_day:depression_score +
                       (1 | biome_id), data = ssridf_lacto_mm_clean, REML = FALSE)

model_lacto_mm_noMenses <- lmer(lacto_rel_abundance_vag ~ SSRI_status + depression_score +
                                                    SSRI_status:depression_score +
                                                    (1 | biome_id), data = ssridf_lacto_mm_clean, REML = FALSE)

model_lacto_mm_noDepression <- lmer(lacto_rel_abundance_vag ~ SSRI_status + menses_day +
                             SSRI_status:menses_day +
                             (1 | biome_id), data = ssridf_lacto_mm_clean, REML = FALSE)

anova(model_lacto_mm_noSSRI, model_lacto_mm_full)       # Test for SSRI main effect p = 0.1
anova(model_lacto_mm_noMenses, model_lacto_mm_full)     # Test for menses_day main effect
anova(model_lacto_mm_noDepression, model_lacto_mm_full) # Test for depression_score main effect
##################################################################################
#full model for poster - vaginal CA
#creating a new dataset where na are already dropped
ssridf_vagCA_mm_clean <- ssridf_mood_menses %>%
  dplyr::select(calbican_rel_abundance_vag, SSRI_status, depression_score, menses_day, biome_id) %>%
  tidyr::drop_na()

model_vagCA_mm_full <- lmer(calbican_rel_abundance_vag ~ SSRI_status + menses_day + depression_score +
                              SSRI_status:menses_day + SSRI_status:depression_score + menses_day:depression_score +
                              (1 | biome_id), data = ssridf_vagCA_mm_clean, REML = FALSE)

model_vagCA_mm_noSSRI <- lmer(calbican_rel_abundance_vag ~ menses_day + depression_score +
                                menses_day:depression_score +
                                (1 | biome_id), data = ssridf_vagCA_mm_clean, REML = FALSE)

model_vagCA_mm_noMenses <- lmer(calbican_rel_abundance_vag ~ SSRI_status + depression_score +
                                  SSRI_status:depression_score +
                                  (1 | biome_id), data = ssridf_vagCA_mm_clean, REML = FALSE)

model_vagCA_mm_noDepression <- lmer(calbican_rel_abundance_vag ~ SSRI_status + menses_day +
                                      SSRI_status:menses_day +
                                      (1 | biome_id), data = ssridf_vagCA_mm_clean, REML = FALSE)

anova(model_vagCA_mm_noSSRI, model_vagCA_mm_full)       # Test for SSRI main effect
anova(model_vagCA_mm_noMenses, model_vagCA_mm_full)     # Test for menses_day main effect
anova(model_vagCA_mm_noDepression, model_vagCA_mm_full) # Test for depression_score main effect
##################################################################################
ssridf_gutCA_mm_clean <- ssridf_mood_menses %>%
  dplyr::select(calbican_rel_abundance_gut, SSRI_status, depression_score, menses_day, biome_id) %>%
  tidyr::drop_na()

model_gutCA_mm_full <- lmer(calbican_rel_abundance_gut ~ SSRI_status + menses_day + depression_score +
                              SSRI_status:menses_day + SSRI_status:depression_score + menses_day:depression_score +
                              (1 | biome_id), data = ssridf_gutCA_mm_clean, REML = FALSE)

model_gutCA_mm_noSSRI <- lmer(calbican_rel_abundance_gut ~ menses_day + depression_score +
                                menses_day:depression_score +
                                (1 | biome_id), data = ssridf_gutCA_mm_clean, REML = FALSE)

model_gutCA_mm_noMenses <- lmer(calbican_rel_abundance_gut ~ SSRI_status + depression_score +
                                  SSRI_status:depression_score +
                                  (1 | biome_id), data = ssridf_gutCA_mm_clean, REML = FALSE)

model_gutCA_mm_noDepression <- lmer(calbican_rel_abundance_gut ~ SSRI_status + menses_day +
                                      SSRI_status:menses_day +
                                      (1 | biome_id), data = ssridf_gutCA_mm_clean, REML = FALSE)

anova(model_gutCA_mm_noSSRI, model_gutCA_mm_full)       # Test for SSRI main effect
anova(model_gutCA_mm_noMenses, model_gutCA_mm_full)     # Test for menses_day main effect
anova(model_gutCA_mm_noDepression, model_gutCA_mm_full) # Test for depression_score main effect
##################################################################################
#looking at depression, menses, and SSRI in a continuous way
lineplot_df <- ssridf_mood_menses %>%
  filter(!is.na(menses_day), !is.na(depression_score)) %>%  # remove NAs
  mutate(
    menses_label = case_when(
      menses_day == "menses" ~ "On Menses",
      menses_day == "not_menses" ~ "Not on Menses"
    ),
    group = paste(SSRI_status, menses_label, sep = " - ")
  )

ggplot(lineplot_df, aes(x = depression_score, y = calbican_rel_abundance_vag, color = group)) +
  geom_point(alpha = 0.5, size = 1.8) +  # Dots
  geom_smooth(method = "loess", se = FALSE, size = 0.8) +  # Smooth solid lines
  labs(
    x = "Depression Score",
    y = "C. albicans Relative Abundance (Vaginal)",
    color = "Group"
  ) +
  theme_minimal()

ggplot(lineplot_df, aes(x = depression_score, y = calbican_rel_abundance_gut, color = group)) +
  geom_point(alpha = 0.5, size = 1.8) +  # Dots
  geom_smooth(method = "loess", se = FALSE, size = 0.8) +  # Smooth solid lines
  labs(
    x = "Depression Score",
    y = "C. albicans Relative Abundance (Gut)",
    color = "Group"
  ) +
  theme_minimal()

ggplot(lineplot_df, aes(x = depression_score, y = lacto_rel_abundance_vag, color = group)) +
  geom_point(alpha = 0.5, size = 1.8) +  # Dots
  geom_smooth(method = "loess", se = FALSE, size = 0.8) +  # Smooth solid lines
  labs(
    x = "Depression Score",
    y = "Lactobacillus Relative Abundance (Vaginal)",
    color = "Group"
  ) +
  theme_minimal()




#################################################################################
#PCA by lots of categorical variable
forPCAmerge <- sample_data(fungal2.0) %>%
  data.frame() %>%
  rownames_to_column(var = "SampleID")

columns_to_join <- c("biome_id", "logDate", "menses_day", "CST", "ethnicity", "sexuality",
                     "activity_level", "birthControl", "field_hockey",
                     "depressionseverity", "anxietyseverity", "stressseverity")

ssri_subset <- ssridf_mood_menses %>%
  select(all_of(columns_to_join))

forPCAmerge$biome_id <- as.character(sample_df$biome_id)

joined_df <- forPCAmerge %>%
  left_join(ssri_subset, by = c("biome_id", "logDate"))

new_sample_data <- joined_df %>%
  column_to_rownames(var = "SampleID") %>%
  sample_data()



############
#creating study phase column
# Ensure logDate is in Date format
new_sample_data_df <- data.frame(new_sample_data)
new_sample_data_df$logDate <- as.Date(new_sample_data_df$logDate)

# Remove NAs in logDate before breaking into phases
valid_dates <- new_sample_data_df %>%
  filter(!is.na(logDate))

numeric_dates <- as.numeric(valid_dates$logDate)

quantile_cutoffs <- quantile(numeric_dates, probs = c(0, 1/3, 2/3, 1), na.rm = TRUE)

date_cutoffs <- as.Date(quantile_cutoffs, origin = "1970-01-01")

new_sample_data_df$study_phase <- cut(
  as.Date(new_sample_data_df$logDate),
  breaks = date_cutoffs,
  labels = c("Phase 1", "Phase 2", "Phase 3"),
  include.lowest = TRUE
)

updated_sample_data <- sample_data(new_sample_data_df)

# Replace the sample_data in the phyloseq object
phyleqforPCA <- phyloseq(
  otu_table(fungal2.0),
  tax_table(fungal2.0),
  updated_sample_data
)
##############
#separating vaginal phyloseq
vaginal_phy <- subset_samples(phyleqforPCA, sampleType == "vaginal")

vaginal_phy <- vaginal_phy %>% tax_transform("clr")
vaginal_phy <- vaginal_phy %>% ord_calc(method = "PCA")
ord_plot(vaginal_phy, color = "menses_day", size = 2) +
  stat_ellipse(aes(color = menses_day), type = "norm", size = 0.6, linetype = "dashed") +
  theme_minimal()

ord_plot(vaginal_phy, color = "CST", size = 2) +
  stat_ellipse(aes(color = CST), type = "norm", size = 0.6, linetype = "dashed") +
  theme_minimal()

ord_plot(vaginal_phy, color = "activity_level", size = 2) +
  #stat_ellipse(aes(color = activity_level), type = "norm", size = 0.6, linetype = "dashed") +
  theme_minimal()

ord_plot(vaginal_phy, color = "birthControl", size = 2) +
  stat_ellipse(aes(color = birthControl), type = "norm", size = 0.6, linetype = "dashed") +
  theme_minimal()

ord_plot(vaginal_phy, color = "field_hockey", size = 2) +
  stat_ellipse(aes(color = field_hockey), type = "norm", size = 0.6, linetype = "dashed") +
  theme_minimal()

ord_plot(vaginal_phy, color = "stressseverity", size = 2) +
  stat_ellipse(aes(color = stressseverity), type = "norm", size = 0.6, linetype = "dashed") +
  theme_minimal()

ord_plot(vaginal_phy, color = "study_phase", size = 2) +
  stat_ellipse(aes(color = study_phase), type = "norm", size = 0.6, linetype = "dashed") +
  theme_minimal()

ord_plot(vaginal_phy, color = "ethnicity", size = 2) +
  stat_ellipse(aes(color = ethnicity), type = "norm", size = 0.6, linetype = "dashed") +
  theme_minimal()

ord_plot(vaginal_phy, color = "menses_day", size = 2) +
  stat_ellipse(aes(color = menses_day), type = "norm", size = 0.6, linetype = "dashed") +
  theme_minimal()

ord_plot(vaginal_phy, color = "biome_id", size = 2) +
#   #stat_ellipse(aes(color = biome_id), type = "norm", size = 0.6, linetype = "dashed") +
  theme_minimal()
#none of the results are significant.....



############################################################################
biome1to10 <- updated_sample_data %>%
  data.frame() %>%
  mutate(biome_label = ifelse(biome_id %in% 1:10, as.character(biome_id), "Other"))

# Create phyloseq object with updated sample data
phyleqforPCA1to10 <- phyloseq(
  otu_table(fungal2.0),
  tax_table(fungal2.0),
  sample_data(biome1to10)
)

# Subset vaginal samples
vaginal_phy1to10 <- subset_samples(phyleqforPCA1to10, sampleType == "vaginal")

# Transform with CLR, calculate ordination, and plot with biome_label coloring
vaginal_phy1to10 <- vaginal_phy1to10 %>% tax_transform("clr") %>% ord_calc(method = "PCA")

ord_plot(vaginal_phy1to10, color = "biome_label", size = 2) +
  stat_ellipse(aes(color = biome_label), type = "norm", size = 0.6, linetype = "dashed") +
  scale_color_manual(values = c(
    "1" = "red", "2" = "blue", "3" = "green", "4" = "purple", "5" = "orange",
    "6" = "cyan", "7" = "pink", "8" = "brown", "9" = "yellow", "10" = "darkgreen",
    "Other" = "gray70"
  )) +
  theme_minimal()
############################################################################
biome11to20 <- updated_sample_data %>%
  data.frame() %>%
  mutate(biome_label = ifelse(biome_id %in% 11:20, as.character(biome_id), "Other"))

# Create phyloseq object with updated sample data
phyleqforPCA11to20 <- phyloseq(
  otu_table(fungal2.0),
  tax_table(fungal2.0),
  sample_data(biome11to20)
)

# Subset vaginal samples
vaginal_phy11to20 <- subset_samples(phyleqforPCA11to20, sampleType == "vaginal")

# Transform with CLR, calculate ordination, and plot with biome_label coloring
vaginal_phy11to20 <- vaginal_phy11to20 %>% tax_transform("clr") %>% ord_calc(method = "PCA")

ord_plot(vaginal_phy11to20, color = "biome_label", size = 2) +
  stat_ellipse(aes(color = biome_label), type = "norm", size = 0.6, linetype = "dashed") +
  scale_color_manual(values = c(
    "11" = "red", "12" = "blue", "13" = "green", "14" = "purple", "15" = "orange",
    "16" = "cyan", "17" = "pink", "18" = "brown", "19" = "yellow", "20" = "darkgreen",
    "Other" = "gray70"
  )) +
  theme_minimal()
############################################################################
biome21to30 <- updated_sample_data %>%
  data.frame() %>%
  mutate(biome_label = ifelse(biome_id %in% 21:30, as.character(biome_id), "Other"))

# Create phyloseq object with updated sample data
phyleqforPCA21to30 <- phyloseq(
  otu_table(fungal2.0),
  tax_table(fungal2.0),
  sample_data(biome21to30)
)

# Subset vaginal samples
vaginal_phy21to30 <- subset_samples(phyleqforPCA21to30, sampleType == "vaginal")

# Transform with CLR, calculate ordination, and plot with biome_label coloring
vaginal_phy21to30 <- vaginal_phy21to30 %>% tax_transform("clr") %>% ord_calc(method = "PCA")

ord_plot(vaginal_phy21to30, color = "biome_label", size = 2) +
  stat_ellipse(aes(color = biome_label), type = "norm", size = 0.6, linetype = "dashed") +
  scale_color_manual(values = c(
    "21" = "red", "22" = "blue", "23" = "green", "24" = "purple", "25" = "orange",
    "26" = "cyan", "27" = "pink", "28" = "brown", "29" = "yellow", "30" = "darkgreen",
    "Other" = "gray70"
  )) +
  theme_minimal()
############################################################################
biome31to40 <- updated_sample_data %>%
  data.frame() %>%
  mutate(biome_label = ifelse(biome_id %in% 31:40, as.character(biome_id), "Other"))

# Create phyloseq object with updated sample data
phyleqforPCA31to40 <- phyloseq(
  otu_table(fungal2.0),
  tax_table(fungal2.0),
  sample_data(biome31to40)
)

# Subset vaginal samples
vaginal_phy31to40 <- subset_samples(phyleqforPCA31to40, sampleType == "vaginal")

# Transform with CLR, calculate ordination, and plot with biome_label coloring
vaginal_phy31to40 <- vaginal_phy31to40 %>% tax_transform("clr") %>% ord_calc(method = "PCA")

ord_plot(vaginal_phy31to40, color = "biome_label", size = 2) +
  stat_ellipse(aes(color = biome_label), type = "norm", size = 0.6, linetype = "dashed") +
  scale_color_manual(values = c(
    "31" = "red", "32" = "blue", "33" = "green", "34" = "purple", "35" = "orange",
    "36" = "cyan", "37" = "pink", "38" = "brown", "39" = "yellow", "40" = "darkgreen",
    "Other" = "gray70"
  )) +
  theme_minimal()
############################################################################
biome41to50 <- updated_sample_data %>%
  data.frame() %>%
  mutate(biome_label = ifelse(biome_id %in% 41:50, as.character(biome_id), "Other"))

# Create phyloseq object with updated sample data
phyleqforPCA41to50 <- phyloseq(
  otu_table(fungal2.0),
  tax_table(fungal2.0),
  sample_data(biome41to50)
)

# Subset vaginal samples
vaginal_phy41to50 <- subset_samples(phyleqforPCA41to50, sampleType == "vaginal")

# Transform with CLR, calculate ordination, and plot with biome_label coloring
vaginal_phy41to50 <- vaginal_phy41to50 %>% tax_transform("clr") %>% ord_calc(method = "PCA")

ord_plot(vaginal_phy41to50, color = "biome_label", size = 2) +
  stat_ellipse(aes(color = biome_label), type = "norm", size = 0.6, linetype = "dashed") +
  scale_color_manual(values = c(
    "41" = "red", "42" = "blue", "43" = "green", "44" = "purple", "45" = "orange",
    "46" = "cyan", "47" = "pink", "48" = "brown", "49" = "yellow", "50" = "darkgreen",
    "Other" = "gray70"
  )) +
  theme_minimal()
############################################################################
biome51to60 <- updated_sample_data %>%
  data.frame() %>%
  mutate(biome_label = ifelse(biome_id %in% 51:60, as.character(biome_id), "Other"))

# Create phyloseq object with updated sample data
phyleqforPCA51to60 <- phyloseq(
  otu_table(fungal2.0),
  tax_table(fungal2.0),
  sample_data(biome51to60)
)

# Subset vaginal samples
vaginal_phy51to60 <- subset_samples(phyleqforPCA51to60, sampleType == "vaginal")

# Transform with CLR, calculate ordination, and plot with biome_label coloring
vaginal_phy51to60 <- vaginal_phy51to60 %>% tax_transform("clr") %>% ord_calc(method = "PCA")

ord_plot(vaginal_phy51to60, color = "biome_label", size = 2) +
  stat_ellipse(aes(color = biome_label), type = "norm", size = 0.6, linetype = "dashed") +
  scale_color_manual(values = c(
    "51" = "red", "52" = "blue", "53" = "green", "54" = "purple", "55" = "orange",
    "56" = "cyan", "57" = "pink", "58" = "brown", "59" = "yellow", "60" = "darkgreen",
    "Other" = "gray70"
  )) +
  theme_minimal()
############################################################################
biome61to70 <- updated_sample_data %>%
  data.frame() %>%
  mutate(biome_label = ifelse(biome_id %in% 61:70, as.character(biome_id), "Other"))

# Create phyloseq object with updated sample data
phyleqforPCA61to70 <- phyloseq(
  otu_table(fungal2.0),
  tax_table(fungal2.0),
  sample_data(biome61to70)
)

# Subset vaginal samples
vaginal_phy61to70 <- subset_samples(phyleqforPCA61to70, sampleType == "vaginal")

# Transform with CLR, calculate ordination, and plot with biome_label coloring
vaginal_phy61to70 <- vaginal_phy61to70 %>% tax_transform("clr") %>% ord_calc(method = "PCA")

ord_plot(vaginal_phy61to70, color = "biome_label", size = 2) +
  stat_ellipse(aes(color = biome_label), type = "norm", size = 0.6, linetype = "dashed") +
  scale_color_manual(values = c(
    "61" = "red", "62" = "blue", "63" = "green", "64" = "purple", "65" = "orange",
    "66" = "cyan", "67" = "pink", "68" = "brown", "69" = "yellow", "70" = "darkgreen",
    "Other" = "gray70"
  )) +
  theme_minimal()


############################################################################
#case study of SSRI user and menstruation and CA abundance
menses_11 <- ssridf_mood_menses %>%
  filter(biome_id == 11)

ggplot(menses_11, aes(x = logDate, y = calbican_rel_abundance_vag)) +
  geom_point(aes(shape = menses_day, color = menses_day), size = 3) +
  scale_shape_manual(values = c("menses" = 17, "not_menses" = 16)) +  # 17 = triangle, 16 = circle
  scale_color_manual(values = c("menses" = "red", "not_menses" = "black")) +
  labs(
    x = "Date",
    y = "C. albicans Abundance (Vaginal)",
    title = "C. albicans Abundance Over Time (biome_id = 11, SSRI user)"
  ) +
  theme_minimal()

ggplot(menses_11, aes(x = logDate, y = Shannon_vag)) +
  geom_point(aes(shape = menses_day, color = menses_day), size = 3) +
  scale_shape_manual(values = c("menses" = 17, "not_menses" = 16)) +  # 17 = triangle, 16 = circle
  scale_color_manual(values = c("menses" = "red", "not_menses" = "black")) +
  labs(
    x = "Date",
    y = "Shannon Diversity (Vaginal)",
    title = "Shannon Diversity Over Time (biome_id = 11, SSRI user)"
  ) +
  theme_minimal()

menses_33 <- ssridf_mood_menses %>%
  filter(biome_id == 33)
ggplot(menses_33, aes(x = logDate, y = calbican_rel_abundance_vag)) +
  geom_point(aes(shape = menses_day, color = menses_day), size = 3) +
  scale_shape_manual(values = c("menses" = 17, "not_menses" = 16)) +  # 17 = triangle, 16 = circle
  scale_color_manual(values = c("menses" = "red", "not_menses" = "black")) +
  labs(
    x = "Date",
    y = "C. albicans Abundance (Vaginal)",
    title = "C. albicans Abundance Over Time (biome_id = 33, SSRI user)"
  ) +
  theme_minimal()

ggplot(menses_33, aes(x = logDate, y = Shannon_vag)) +
  geom_point(aes(shape = menses_day, color = menses_day), size = 3) +
  scale_shape_manual(values = c("menses" = 17, "not_menses" = 16)) +  # 17 = triangle, 16 = circle
  scale_color_manual(values = c("menses" = "red", "not_menses" = "black")) +
  labs(
    x = "Date",
    y = "Shannon Diversity (Vaginal)",
    title = "Shannon Diversity Over Time (biome_id = 33, SSRI user)"
  ) +
  theme_minimal()

menses_60 <- ssridf_mood_menses %>%
  filter(biome_id == 60)
ggplot(menses_60, aes(x = logDate, y = calbican_rel_abundance_vag)) +
  geom_point(aes(shape = menses_day, color = menses_day), size = 3) +
  scale_shape_manual(values = c("menses" = 17, "not_menses" = 16)) +  # 17 = triangle, 16 = circle
  scale_color_manual(values = c("menses" = "red", "not_menses" = "black")) +
  labs(
    x = "Date",
    y = "C. albicans Abundance (Vaginal)",
    title = "C. albicans Abundance Over Time (biome_id = 60, SSRI user)"
  ) +
  theme_minimal()

ggplot(menses_60, aes(x = logDate, y = Shannon_vag)) +
  geom_point(aes(shape = menses_day, color = menses_day), size = 3) +
  scale_shape_manual(values = c("menses" = 17, "not_menses" = 16)) +  # 17 = triangle, 16 = circle
  scale_color_manual(values = c("menses" = "red", "not_menses" = "black")) +
  labs(
    x = "Date",
    y = "Shannon Diversity (Vaginal)",
    title = "Shannon Diversity Over Time (biome_id = 60, SSRI user)"
  ) +
  theme_minimal()
##############################
#non-user
menses_1 <- ssridf_mood_menses %>%
  filter(biome_id == 1)

ggplot(menses_1, aes(x = logDate, y = calbican_rel_abundance_vag)) +
  geom_point(aes(shape = menses_day, color = menses_day), size = 3) +
  scale_shape_manual(values = c("menses" = 17, "not_menses" = 16)) +  # 17 = triangle, 16 = circle
  scale_color_manual(values = c("menses" = "red", "not_menses" = "black")) +
  labs(
    x = "Date",
    y = "C. albicans Abundance (Vaginal)",
    title = "C. albicans Abundance Over Time (biome_id = 1, Non-User)"
  ) +
  theme_minimal()
ggplot(menses_1, aes(x = logDate, y = Shannon_vag)) +
  geom_point(aes(shape = menses_day, color = menses_day), size = 3) +
  scale_shape_manual(values = c("menses" = 17, "not_menses" = 16)) +  # 17 = triangle, 16 = circle
  scale_color_manual(values = c("menses" = "red", "not_menses" = "black")) +
  labs(
    x = "Date",
    y = "Shannon Diversity (Vaginal)",
    title = "Shannon Diversity Over Time (biome_id = 1, Non-User)"
  ) +
  theme_minimal()


menses_61 <- ssridf_mood_menses %>%
  filter(biome_id == 61)

ggplot(menses_61, aes(x = logDate, y = calbican_rel_abundance_vag)) +
  geom_point(aes(shape = menses_day, color = menses_day), size = 3) +
  scale_shape_manual(values = c("menses" = 17, "not_menses" = 16)) +  # 17 = triangle, 16 = circle
  scale_color_manual(values = c("menses" = "red", "not_menses" = "black")) +
  labs(
    x = "Date",
    y = "C. albicans Abundance (Vaginal)",
    title = "C. albicans Abundance Over Time (biome_id = 61, Non-User)"
  ) +
  theme_minimal()
ggplot(menses_61, aes(x = logDate, y = Shannon_vag)) +
  geom_point(aes(shape = menses_day, color = menses_day), size = 3) +
  scale_shape_manual(values = c("menses" = 17, "not_menses" = 16)) +  # 17 = triangle, 16 = circle
  scale_color_manual(values = c("menses" = "red", "not_menses" = "black")) +
  labs(
    x = "Date",
    y = "Shannon Diversity (Vaginal)",
    title = "Shannon Diversity Over Time (biome_id = 61, Non-User)"
  ) +
  theme_minimal()
#################################################################
#investigating the two clouds
#creating a column in the phyloseq to stand for the clouds
pca_obj <- ord_calc(vaginal_phy, method = "PCA", comp = 3)
#vegan_ord <- attr(pca_obj, "vegan_ord")
pca_scores <- as.data.frame(pca_obj@ord$CA$u)
set.seed(123)
kmeans_res <- kmeans(pca_scores[, 1:2], centers = 2)
pca_scores$PCA_cluster <- factor(kmeans_res$cluster)

#merging back to sample data
# Match sample IDs
pca_scores$SampleID <- rownames(pca_scores)

# Get sample_data as data.frame
samp_df_pca <- data.frame(sample_data(vaginal_phy))
samp_df_pca$SampleID <- rownames(samp_df_pca)

# Merge
merged_pca_df <- merge(samp_df_pca, pca_scores[, c("SampleID", "PCA_cluster")], by = "SampleID")
rownames(merged_pca_df) <- merged_pca_df$SampleID

# Replace sample_data in phyloseq object
sample_data(vaginal_phy) <- sample_data(merged_pca_df)

vaginal_phy <- vaginal_phy %>% ord_calc(method = "PCA")

ord_plot(vaginal_phy, color = "PCA_cluster", size = 2) +
  stat_ellipse(aes(color = PCA_cluster), type = "norm", size = 0.6, linetype = "dashed") +
  theme_minimal()
#################################################################
#investigate when the cluster change







#################################################################
#average stress, depression, and CA abundance for gut and vagina
avgDASS_CA_vag <- ssridf_mood_menses %>%
  group_by(biome_id, SSRI_status) %>%
  summarise(
    avg_stress = mean(stress_score, na.rm = TRUE),
    avg_calbicans = mean(calbican_rel_abundance_vag, na.rm = TRUE),
    .groups = "drop"
  )


ggplot(avgDASS_CA_vag, aes(x = avg_stress, y = avg_calbicans, color = SSRI_status)) +
  geom_point(size = 2, alpha = 0.7) +
  #geom_abline(intercept = 0, slope = 1)
  abline(coef = c(0,1)) +
  labs(
    x = "Average Stress Score for Each Participant Throughout Study",
    y = "Average Vaginal C. albicans Abundance Throughout Study",
    title = "Average Stress vs. Vaginal C. albicans Abundance by SSRI Use for Each Participant"
  ) +
  scale_color_manual(values = c("SSRI User" = "cyan3", "Non-User" = "coral1")) +
  theme_minimal()

avgDASS_CA_gut <- ssridf_mood_menses %>%
  group_by(biome_id, SSRI_status) %>%
  summarise(
    avg_depression = mean(depression_score, na.rm = TRUE),
    avg_calbicans = mean(calbican_rel_abundance_gut, na.rm = TRUE),
    .groups = "drop"
  )

ggplot(avgDASS_CA_gut, aes(x = avg_depression, y = avg_calbicans, color = SSRI_status)) +
  geom_point(size = 2, alpha = 0.7) +
  #geom_abline(intercept = 0, slope = 1)
  abline(coef = c(0,1)) +
  labs(
    x = "Average Depression Score for Each Participant Throughout Study",
    y = "Average Gut C. albicans Abundance Throughout Study",
    title = "Average Depression vs. C. albicans Abundance by SSRI Use for Each Participant"
  ) +
  scale_color_manual(values = c("SSRI User" = "cyan3", "Non-User" = "coral1")) +
  theme_minimal()


#if the average is normally distributed
ggplot(avgDASS_CA_vag, aes(x = avg_calbicans, fill = SSRI_status)) +
  geom_density(alpha = 0.5) +
  theme_minimal()

wilcox.test(avg_calbicans ~ SSRI_status, data = avgDASS_CA_vag) #not significant

t.test(avg_calbicans ~ SSRI_status, data = avgDASS_CA_vag)

boxplot(avg_calbicans ~ SSRI_status, data = avgDASS_CA_vag,
        main = "Average C. albicans Abundance by SSRI Status",
        ylab = "Average C. albicans Abundance",
        xlab = "SSRI Status",
        col = c("lightblue", "pink"))





