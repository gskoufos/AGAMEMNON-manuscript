#import libraries
library(dplyr)
library(Metrics)
library(ggplot2)
library(ggpubr)
library(stringr)
library(reshape2)
library(data.table)
library(RColorBrewer)
library(wavethresh)


###################################################################################################################################
# IMPORT RESULTS AND GROUND TRUTH INTO R  #
###################################################################################################################################

#directories to results files and ground truth
agamemnon.file <- ""
kraken2.file <- ""
bracken.species.file = ""
bracken.genus.file = ""
ground.truth.file = ""
kaiju.genus.file = ""
kaiju.species.file = ""
metaphlan.file <- ""


#function that reads results into dataframes
read_data <- function(filename, headers){
  temp.dataframe <- read.csv(file = filename,
                             sep = "\t",
                             header = headers,
                             stringsAsFactors = FALSE,
                             quote = "",
                             fill = FALSE)
  temp.dataframe <- temp.dataframe %>% mutate(across(where(is.character), str_trim))
  return(temp.dataframe)
}


ground.truth <- read_data(ground.truth.file, TRUE)
ground.truth$NumCounts <- as.integer(as.character(ground.truth$NumCounts))
ground.truth <- ground.truth[ground.truth$Genus != "---", ]
total.reads.in.dataset <- sum(ground.truth$NumCounts)

agamemnon <- read_data(agamemnon.file, TRUE)
agamemnon$NumCounts <- as.integer(as.character(agamemnon$NumCounts))

kraken2 <- read_data(kraken2.file, FALSE)
colnames(kraken2) <- c("percent", "reads_total", "reads_direct", "taxonomic_rank", "taxid", "scientific_name")
kraken2$reads_total <- as.numeric(as.character(kraken2$reads_total))
kraken2$reads_direct <- as.numeric(as.character(kraken2$reads_direct))

bracken.species <- read_data(bracken.species.file, TRUE)
bracken.species$new_est_reads <- as.numeric(as.character(bracken.species$new_est_reads))

bracken.genus <- read_data(bracken.genus.file, TRUE)
bracken.genus$new_est_reads <- as.numeric(as.character(bracken.genus$new_est_reads))

kaiju.species <- read_data(kaiju.species.file, TRUE)
kaiju.species$reads <- as.numeric(as.character(kaiju.species$reads))

kaiju.genus <- read_data(kaiju.genus.file, TRUE)
kaiju.genus$reads <- as.numeric(as.character(kaiju.genus$reads))

metaphlan3 <- read_data(metaphlan.file, TRUE)
metaphlan3$taxid <- as.character(lapply(strsplit(as.character(metaphlan3$clade_taxid), split = "\\|"), tail, n = 1))
metaphlan3$taxon <- as.character(lapply(strsplit(as.character(metaphlan3$clade_name), split = "\\|"), tail, n = 1))
metaphlan3 <- metaphlan3 %>% mutate(across(where(is.character), str_trim))
metaphlan3 <- metaphlan3 %>% dplyr::select(taxon, taxid, estimated_number_of_reads_from_the_clade)
#use kraken2 scientific names and taxids for consistency and bring output to the taxon "\t" counts format
metaphlan3 <- merge(metaphlan3, kraken2, by = "taxid", all.x = TRUE)
metaphlan3 <- metaphlan3[(metaphlan3$taxonomic_rank %in% c("G", "G1", "S", "S1") | 
                            is.na(metaphlan3$taxonomic_rank)), ]
metaphlan3$temp <- ifelse(as.character(lapply(strsplit(as.character(metaphlan3$taxon), split = "__"), head, n = 1)) %in% 
                            c("g", "s"), 1, 0)
metaphlan3 <- metaphlan3[metaphlan3$temp == 1, ]
metaphlan3$temp <- NULL
metaphlan3$scientific_name <- ifelse(is.na(metaphlan3$taxonomic_rank), metaphlan3$taxon, metaphlan3$scientific_name)
metaphlan3$taxonomic_rank <- ifelse(substring(metaphlan3$scientific_name, 1, 1) %in% c("s", "g"), 
                                    substring(metaphlan3$scientific_name, 1, 1), metaphlan3$taxonomic_rank)
metaphlan3$taxonomic_rank <- toupper(metaphlan3$taxonomic_rank)
metaphlan3$scientific_name <- as.character(lapply(strsplit(as.character(metaphlan3$scientific_name), split = "__"), tail, n = 1))
metaphlan3$temp <- gsub("_", " ", metaphlan3$scientific_name)
metaphlan3$scientific_name <- metaphlan3$temp
metaphlan3$temp <- NULL
metaphlan3 <- metaphlan3 %>% dplyr::select(scientific_name, taxonomic_rank, estimated_number_of_reads_from_the_clade)
metaphlan3_genus <- metaphlan3[metaphlan3$taxonomic_rank %in% c("G", "G1"), ]
metaphlan3_species <- metaphlan3[metaphlan3$taxonomic_rank %in% c("S", "S1"), ]
metaphlan3_genus$taxonomic_rank <- NULL; colnames(metaphlan3_genus) <- c("genus", "counts")
metaphlan3_species$taxonomic_rank <- NULL; colnames(metaphlan3_species) <- c("species", "counts")
#divide counts by 2 because of the way metaphlan3 is handling paired-end datasets
metaphlan3_species$counts <- metaphlan3_species$counts / 2
metaphlan3_genus$counts <- metaphlan3_genus$counts / 2



###################################################################################################################################
# EDIT DATAFRAMES TO GET GENUS, SPECIES AND STRAIN LEVEL COUNTS  #
###################################################################################################################################

#GROUND TRUTH
ground.truth.species <- aggregate(ground.truth$NumCounts, by = list(Species = ground.truth$Species), FUN = sum)
ground.truth.genus <- aggregate(ground.truth$NumCounts, by = list(Genus = ground.truth$Genus), FUN = sum)
ground.truth.strain <- aggregate(ground.truth$NumCounts, by = list(Strain = ground.truth$Scientific_Name), FUN = sum)
colnames(ground.truth.species) <- c("species", "counts"); colnames(ground.truth.genus) <- c("genus", "counts"); colnames(ground.truth.strain) <- c("strain", "counts")


#AGAMEMNON
agamemnon.species <- aggregate(agamemnon$NumCounts, by = list(Species = agamemnon$Species), FUN = sum)
agamemnon.genus <- aggregate(agamemnon$NumCounts, by = list(Genus = agamemnon$Genus), FUN = sum)
agamemnon.strain <- aggregate(agamemnon$NumCounts, by = list(Strain = agamemnon$Scientific_Name), FUN = sum)
colnames(agamemnon.species) <- c("species", "counts"); colnames(agamemnon.genus) <- c("genus", "counts"); colnames(agamemnon.strain) <- c("strain", "counts")


#KRAKEN 2
kraken2.species <- kraken2[(kraken2$taxonomic_rank == "S"), ]
kraken2.species <- kraken2.species %>% dplyr::select(scientific_name, reads_total)
colnames(kraken2.species) <- c("species", "counts")

kraken2.genus <- kraken2[(kraken2$taxonomic_rank == "G"), ]
kraken2.genus <- kraken2.genus %>% dplyr::select(scientific_name, reads_total)
colnames(kraken2.genus) <- c("genus", "counts")

kraken2.strain <- kraken2[kraken2$taxid %in% agamemnon$TaxID, ]
kraken2.strain <- kraken2.strain %>% dplyr::select(scientific_name, reads_total)
colnames(kraken2.strain) <- c("strain", "counts")


#BRACKEN
bracken.species <- bracken.species %>% dplyr::select(name, new_est_reads)
colnames(bracken.species) <- c("species", "counts")

bracken.genus <- bracken.genus %>% dplyr::select(name, new_est_reads)
colnames(bracken.genus) <- c("genus", "counts")


#KAIJU
kaiju.species <- kaiju.species %>% dplyr::select(taxon_name, reads)
colnames(kaiju.species) <- c("species", "counts")

kaiju.genus <- kaiju.genus %>% dplyr::select(taxon_name, reads)
colnames(kaiju.genus) <- c("genus", "counts")


#merge species-level dataframes and replace NULL values with zero
species.merged <- Reduce(function(x, y) merge(x, y, by = "species", all = TRUE), list(ground.truth.species, agamemnon.species, kraken2.species, bracken.species))
colnames(species.merged) <- c("species", "ground_truth", "agamemnon", "kraken2", "bracken")
species.merged[is.na(species.merged)] <- 0
#remove "---" entries from species-level results
species.merged <- species.merged[species.merged$species != "---", ]

#merge genus-level dataframes and replace NULL values with zero
genus.merged <- Reduce(function(x, y) merge(x, y, by = "genus", all = TRUE), list(ground.truth.genus, agamemnon.genus, kraken2.genus, bracken.genus))
colnames(genus.merged) <- c("genus", "ground_truth", "agamemnon", "kraken2", "bracken")
genus.merged[is.na(genus.merged)] <- 0
#remove "---" entries from genus-level results
genus.merged <- genus.merged[genus.merged$genus != "---", ]

#merge strain-level dataframes and replace NULL values with zero
strain.merged <- Reduce(function(x, y) merge(x, y, by = "strain", all = TRUE), list(ground.truth.strain, agamemnon.strain, kraken2.strain))
colnames(strain.merged) <- c("strain", "ground_truth", "agamemnon", "kraken2")
strain.merged[is.na(strain.merged)] <- 0
#remove "---" entries from genus-level results
strain.merged <- strain.merged[strain.merged$strain != "---", ]


#create species/genus dataframes with true counts and metaphlan3 counts
true.metaphlan.species <- merge(ground.truth.species, metaphlan3_species, by = "species", all = TRUE)
true.metaphlan.species[is.na(true.metaphlan.species)] <- 0
colnames(true.metaphlan.species) <- c("species", "ground_truth", "metaphlan3")

true.metaphlan.genus <- merge(ground.truth.genus, metaphlan3_genus, by = "genus", all = TRUE)
true.metaphlan.genus[is.na(true.metaphlan.genus)] <- 0
colnames(true.metaphlan.genus) <- c("genus", "ground_truth", "metaphlan3")


#create speceis/genus dataframes with true counts and kaiju counts
true.kaiju.species <- merge(ground.truth.species, kaiju.species, by = "species", all = TRUE)
true.kaiju.species[is.na(true.kaiju.species)] <- 0
colnames(true.kaiju.species) <- c("species", "ground_truth", "kaiju")

true.kaiju.genus <- merge(ground.truth.genus, kaiju.genus, by = "genus", all = TRUE)
true.kaiju.genus[is.na(true.kaiju.genus)] <- 0
colnames(true.kaiju.genus) <- c("genus", "ground_truth", "kaiju")


#create only true positive dataframes
species.tp <- species.merged[species.merged$ground_truth > 0, ]
genus.tp <- genus.merged[genus.merged$ground_truth > 0, ]
strains.tp <- strain.merged[strain.merged$ground_truth > 0, ]
metaphlan.species.tp <- true.metaphlan.species[true.metaphlan.species$ground_truth > 0, ]
metaphlan.genus.tp <- true.metaphlan.genus[true.metaphlan.genus$ground_truth > 0, ]
kaiju.species.tp <- true.kaiju.species[true.kaiju.species$ground_truth > 0, ]
kaiju.genus.tp <- true.kaiju.genus[true.kaiju.genus$ground_truth > 0, ]

species.tp <- Reduce(function(x, y) merge(x, y, by = "species", all = TRUE), list(species.tp, metaphlan.species.tp, kaiju.species.tp))
genus.tp <- Reduce(function(x, y) merge(x, y, by = "genus", all = TRUE), list(genus.tp, metaphlan.genus.tp, kaiju.genus.tp))
species.tp <- species.tp %>% dplyr::select(species, ground_truth.x, agamemnon, kraken2, bracken, metaphlan3, kaiju)
genus.tp <- genus.tp %>% dplyr::select(genus, ground_truth.x, agamemnon, kraken2, bracken, metaphlan3, kaiju)
colnames(species.tp)[2] <- "ground_truth"; colnames(genus.tp)[2] <- "ground_truth"


###################################################################################################################################
# COMPUTE METRICS  #
###################################################################################################################################

#function to calculate manhattan distance
manhattan_dist <- function(a, b){
  m_distance <- abs(a-b)
  m_distance <- sum(m_distance)
  return(m_distance)
}

#function to calculate canberra distance
canberra_dist <- function(x, y){
  a1 <-  abs(x - y)
  a2 <- abs(x) + abs(y)
  c_distance <- sum(a1 / a2, na.rm = T)
  return(c_distance)
}

#function to calculate mean squared log error
msle.continuous <- function(num, table) {
  table <- table[rowSums(table) >= num, ]
  return(msle(table$GROUND.TRUTH, table$APPLICATION))
}


#calculate species-level MSLE at different read thresholds
i <- seq(0, 1, 1)
comparison.table <- data.frame(GROUND.TRUTH = species.merged$ground_truth, APPLICATION = species.merged$agamemnon)
AGAMEMNON.msle.species <- unlist(lapply(1:length(i), function(x) msle.continuous(num = i[x], table = comparison.table)))
comparison.table <- data.frame(GROUND.TRUTH = species.merged$ground_truth, APPLICATION = species.merged$kraken2)
KRAKEN2.msle.species <- unlist(lapply(1:length(i), function(x) msle.continuous(num = i[x], table = comparison.table)))
comparison.table <- data.frame(GROUND.TRUTH = species.merged$ground_truth, APPLICATION = species.merged$bracken)
BRACKEN.msle.species <- unlist(lapply(1:length(i), function(x) msle.continuous(num = i[x], table = comparison.table)))
comparison.table <- data.frame(GROUND.TRUTH = species.merged$ground_truth, APPLICATION = species.merged$bracken)
BRACKEN.msle.species <- unlist(lapply(1:length(i), function(x) msle.continuous(num = i[x], table = comparison.table)))
comparison.table <- data.frame(GROUND.TRUTH = true.metaphlan.species$ground_truth, APPLICATION = true.metaphlan.species$metaphlan3)
METAPHLAN3.msle.species <- unlist(lapply(1:length(i), function(x) msle.continuous(num = i[x], table = comparison.table)))
comparison.table <- data.frame(GROUND.TRUTH = true.kaiju.species$ground_truth, APPLICATION = true.kaiju.species$kaiju)
KAIJU.msle.species <- unlist(lapply(1:length(i), function(x) msle.continuous(num = i[x], table = comparison.table)))

MSLE.species <- data.frame(c(i), c(AGAMEMNON.msle.species), c(KRAKEN2.msle.species), c(BRACKEN.msle.species), c(METAPHLAN3.msle.species), c(KAIJU.msle.species))
colnames(MSLE.species) <- c("Reads_threshold", "AGAMEMNON", "Kraken 2", "Bracken", "Metaphlan 3", "Kaiju")
melted.species.msle <- reshape2::melt(MSLE.species, id.vars = c('Reads_threshold'))


#calculate genus-level MSLE at different read thresholds
i <- seq(0, 1, 1)
comparison.table <- data.frame(GROUND.TRUTH = genus.merged$ground_truth, APPLICATION = genus.merged$agamemnon)
AGAMEMNON.msle.genus <- unlist(lapply(1:length(i), function(x) msle.continuous(num = i[x], table = comparison.table)))
comparison.table <- data.frame(GROUND.TRUTH = genus.merged$ground_truth, APPLICATION = genus.merged$kraken2)
KRAKEN2.msle.genus <- unlist(lapply(1:length(i), function(x) msle.continuous(num = i[x], table = comparison.table)))
comparison.table <- data.frame(GROUND.TRUTH = genus.merged$ground_truth, APPLICATION = genus.merged$bracken)
BRACKEN.msle.genus <- unlist(lapply(1:length(i), function(x) msle.continuous(num = i[x], table = comparison.table)))
comparison.table <- data.frame(GROUND.TRUTH = true.metaphlan.genus$ground_truth, APPLICATION = true.metaphlan.genus$metaphlan3)
METAPHLAN3.msle.genus <- unlist(lapply(1:length(i), function(x) msle.continuous(num = i[x], table = comparison.table)))
comparison.table <- data.frame(GROUND.TRUTH = true.kaiju.genus$ground_truth, APPLICATION = true.kaiju.genus$kaiju)
KAIJU.msle.genus <- unlist(lapply(1:length(i), function(x) msle.continuous(num = i[x], table = comparison.table)))

MSLE.genus <- data.frame(c(i), c(AGAMEMNON.msle.genus), c(KRAKEN2.msle.genus), c(BRACKEN.msle.genus), c(METAPHLAN3.msle.genus), KAIJU.msle.genus)
colnames(MSLE.genus) <- c("Reads_threshold", "AGAMEMNON", "Kraken 2", "Bracken", "Metaphlan 3", "Kaiju")
melted.genus.msle <- reshape2::melt(MSLE.genus, id.vars = c('Reads_threshold'))


#calculate strain-level MSLE at different read thresholds
i <- seq(0, 1, 1)
comparison.table <- data.frame(GROUND.TRUTH = strain.merged$ground_truth, APPLICATION = strain.merged$agamemnon)
AGAMEMNON.msle.strain <- unlist(lapply(1:length(i), function(x) msle.continuous(num = i[x], table = comparison.table)))
comparison.table <- data.frame(GROUND.TRUTH = strain.merged$ground_truth, APPLICATION = strain.merged$kraken2)
KRAKEN2.msle.strain <- unlist(lapply(1:length(i), function(x) msle.continuous(num = i[x], table = comparison.table)))

MSLE.strain <- data.frame(c(i), c(AGAMEMNON.msle.strain), c(KRAKEN2.msle.strain))
colnames(MSLE.strain) <- c("Reads_threshold", "AGAMEMNON", "Kraken 2")
melted.strain.msle <- reshape2::melt(MSLE.strain, id.vars = c('Reads_threshold'))


#calculate species-level number of false positives at different read thresholds
i <- seq(0, 10, 1)
AGAMEMNON.false.positives <- unlist(lapply(i, function(x) nrow(species.merged[(species.merged$ground_truth == 0 & species.merged$agamemnon - species.merged$ground_truth > x), ])))
KRAKEN2.false.positives <- unlist(lapply(i, function(x) nrow(species.merged[(species.merged$ground_truth == 0 & species.merged$kraken2 - species.merged$ground_truth > x), ])))
BRACKEN.false.positives <- unlist(lapply(i, function(x) nrow(species.merged[(species.merged$ground_truth == 0 & species.merged$bracken - species.merged$ground_truth > x), ])))
METAPHLAN3.false.positives <- unlist(lapply(i, function(x) nrow(true.metaphlan.species[(true.metaphlan.species$ground_truth == 0 & true.metaphlan.species$metaphlan3 - true.metaphlan.species$ground_truth > x), ])))
KAIJU.false.positives <- unlist(lapply(i, function(x) nrow(true.kaiju.species[(true.kaiju.species$ground_truth == 0 & true.kaiju.species$kaiju - true.kaiju.species$ground_truth > x), ])))
FP.species <- data.frame(c(i), c(AGAMEMNON.false.positives), c(KRAKEN2.false.positives), c(BRACKEN.false.positives), c(METAPHLAN3.false.positives), c(KAIJU.false.positives))
colnames(FP.species) <- c("Reads_threshold", "AGAMEMNON", "Kraken 2", "Bracken", "Metaphlan 3", "Kaiju")
melted.species.fp <- reshape2::melt(FP.species, id.vars = c('Reads_threshold'))


#calculate genus-level number of false positives at different read thresholds
i <- seq(0, 10, 1)
AGAMEMNON.false.positives <- unlist(lapply(i, function(x) nrow(genus.merged[(genus.merged$ground_truth == 0 & genus.merged$agamemnon - genus.merged$ground_truth > x), ])))
KRAKEN2.false.positives <- unlist(lapply(i, function(x) nrow(genus.merged[(genus.merged$ground_truth == 0 & genus.merged$kraken2 - genus.merged$ground_truth > x), ])))
BRACKEN.false.positives <- unlist(lapply(i, function(x) nrow(genus.merged[(genus.merged$ground_truth == 0 & genus.merged$bracken - genus.merged$ground_truth > x), ])))
METAPHLAN3.false.positives <- unlist(lapply(i, function(x) nrow(true.metaphlan.genus[(true.metaphlan.genus$ground_truth == 0 & true.metaphlan.genus$metaphlan3 - true.metaphlan.genus$ground_truth > x), ])))
KAIJU.false.positives <- unlist(lapply(i, function(x) nrow(true.kaiju.genus[(true.kaiju.genus$ground_truth == 0 & true.kaiju.genus$kaiju - true.kaiju.genus$ground_truth > x), ])))
FP.genus <- data.frame(c(i), c(AGAMEMNON.false.positives), c(KRAKEN2.false.positives), c(BRACKEN.false.positives), c(METAPHLAN3.false.positives), c(KAIJU.false.positives))
colnames(FP.genus) <- c("Reads_threshold", "AGAMEMNON", "Kraken 2", "Bracken", "Metaphlan 3", "Kaiju")
melted.genus.fp <- reshape2::melt(FP.genus, id.vars = c('Reads_threshold'))


#calculate strain-level number of false positives at different read thresholds
i <- seq(0, 10, 1)
AGAMEMNON.false.positives <- unlist(lapply(i, function(x) nrow(strain.merged[(strain.merged$ground_truth == 0 & strain.merged$agamemnon - strain.merged$ground_truth > x), ])))
KRAKEN2.false.positives <- unlist(lapply(i, function(x) nrow(strain.merged[(strain.merged$ground_truth == 0 & strain.merged$kraken2 - strain.merged$ground_truth > x), ])))
FP.strains <- data.frame(c(i), c(AGAMEMNON.false.positives), c(KRAKEN2.false.positives))
colnames(FP.strains) <- c("Reads_threshold", "AGAMEMNON", "Kraken 2")
melted.strains.fp <- reshape2::melt(FP.strains, id.vars = c('Reads_threshold'))



#calculate species-level total number  of classifications and number of correct classifications
agamemnon.total.species <- sum(species.merged$agamemnon)
kraken2.total.species <- sum(species.merged$kraken2)
bracken.total.species <- sum(species.merged$bracken)
metaphlan.total.species <- sum(true.metaphlan.species$metaphlan3)
kaiju.total.species <- sum(true.kaiju.species$kaiju)

species.merged$agamemnon_counts <- ifelse(species.merged$ground_truth > species.merged$agamemnon, 
                                          species.merged$agamemnon, ifelse(species.merged$ground_truth > 0,
                                                                           species.merged$ground_truth, 0))

species.merged$kraken2_counts <- ifelse(species.merged$ground_truth > species.merged$kraken2, 
                                        species.merged$kraken2, ifelse(species.merged$ground_truth > 0,
                                                                       species.merged$ground_truth, 0))

species.merged$bracken_counts <- ifelse(species.merged$ground_truth > species.merged$bracken, 
                                        species.merged$bracken, ifelse(species.merged$ground_truth > 0,
                                                                       species.merged$ground_truth, 0))

true.metaphlan.species$metaphlan_counts <- ifelse(true.metaphlan.species$ground_truth > true.metaphlan.species$metaphlan3, 
                                                  true.metaphlan.species$metaphlan3, ifelse(true.metaphlan.species$ground_truth > 0,
                                                                                            true.metaphlan.species$ground_truth, 0))

true.kaiju.species$kaiju_counts <- ifelse(true.kaiju.species$ground_truth > true.kaiju.species$kaiju,
                                          true.kaiju.species$kaiju, ifelse(true.kaiju.species$ground_truth > 0,
                                                                                true.kaiju.species$ground_truth, 0))


agamemnon.correct.species <- sum(species.merged$agamemnon_counts)
kraken2.correct.species <- sum(species.merged$kraken2_counts)
bracken.correct.species <- sum(species.merged$bracken_counts)
metaphlan.correct.species <- sum(true.metaphlan.species$metaphlan_counts)
kaiju.correct.species <- sum(true.kaiju.species$kaiju_counts)



#calculate genus-level total number  of classifications and number of correct classifications
agamemnon.total.genus <- sum(genus.merged$agamemnon)
kraken2.total.genus <- sum(genus.merged$kraken2)
bracken.total.genus <- sum(genus.merged$bracken)
metaphlan.total.genus <- sum(true.metaphlan.genus$metaphlan3)
kaiju.total.genus <- sum(true.kaiju.genus$kaiju)

genus.merged$agamemnon_counts <- ifelse(genus.merged$ground_truth > genus.merged$agamemnon, 
                                        genus.merged$agamemnon, ifelse(genus.merged$ground_truth > 0,
                                                                       genus.merged$ground_truth, 0))

genus.merged$kraken2_counts <- ifelse(genus.merged$ground_truth > genus.merged$kraken2, 
                                      genus.merged$kraken2, ifelse(genus.merged$ground_truth > 0,
                                                                   genus.merged$ground_truth, 0))

genus.merged$bracken_counts <- ifelse(genus.merged$ground_truth > genus.merged$bracken, 
                                      genus.merged$bracken, ifelse(genus.merged$ground_truth > 0,
                                                                   genus.merged$ground_truth, 0))

true.metaphlan.genus$metaphlan_counts <- ifelse(true.metaphlan.genus$ground_truth > true.metaphlan.genus$metaphlan3, 
                                                true.metaphlan.genus$metaphlan3, ifelse(true.metaphlan.genus$ground_truth > 0,
                                                                                        true.metaphlan.genus$ground_truth, 0))

true.kaiju.genus$kaiju_counts <- ifelse(true.kaiju.genus$ground_truth > true.kaiju.genus$kaiju,
                                        true.kaiju.genus$kaiju, ifelse(true.kaiju.genus$ground_truth > 0,
                                                                       true.kaiju.genus$ground_truth, 0))

agamemnon.correct.genus <- sum(genus.merged$agamemnon_counts)
kraken2.correct.genus <- sum(genus.merged$kraken2_counts)
bracken.correct.genus <- sum(genus.merged$bracken_counts)
metaphlan.correct.genus <- sum(true.metaphlan.genus$metaphlan_counts)
kaiju.correct.genus <- sum(true.kaiju.genus$kaiju_counts)



#calculate strain-level total number  of classifications and number of correct classifications
agamemnon.total.strain <- sum(strain.merged$agamemnon)
kraken2.total.strain <- sum(strain.merged$kraken2)
kallisto.total.strain <- sum(strain.merged$kallisto)

strain.merged$agamemnon_counts <- ifelse(strain.merged$ground_truth > strain.merged$agamemnon, 
                                         strain.merged$agamemnon, ifelse(strain.merged$ground_truth > 0,
                                                                         strain.merged$ground_truth, 0))

strain.merged$kraken2_counts <- ifelse(strain.merged$ground_truth > strain.merged$kraken2, 
                                       strain.merged$kraken2, ifelse(strain.merged$ground_truth > 0,
                                                                     strain.merged$ground_truth, 0))


agamemnon.correct.strain <- sum(strain.merged$agamemnon_counts)
kraken2.correct.strain <- sum(strain.merged$kraken2_counts)

total.correct <- data.frame(c("AGAMEMNON", "Kraken 2", "Bracken", "Metaphlan 3", "Kaiju",
                              "AGAMEMNON", "Kraken 2", "Bracken", "Metaphlan 3", "Kaiju",
                              "AGAMEMNON", "Kraken 2"),
                            c(rep("Genus", 5), rep("Species", 5), rep("Strain", 2)),
                            c(agamemnon.total.genus,  kraken2.total.genus, bracken.total.genus, metaphlan.total.genus, kaiju.total.genus,
                              agamemnon.total.species, kraken2.total.species, bracken.total.species, metaphlan.total.species, kaiju.total.species,
                              agamemnon.total.strain, kraken2.total.strain),
                            c(agamemnon.correct.genus, kraken2.correct.genus, bracken.correct.genus, metaphlan.correct.genus, kaiju.correct.genus,
                              agamemnon.correct.species, kraken2.correct.species, bracken.correct.species, metaphlan.correct.species, kaiju.correct.species,
                              agamemnon.correct.strain, kraken2.correct.strain))
colnames(total.correct) <- c("Method", "Taxonomic rank", "Total_reads", "Correct_reads")
total.correct$precision <- (total.correct$Correct_reads / total.correct$Total_reads) * 100
total.correct$sensitivity <- (total.correct$Correct_reads / total.reads.in.dataset) * 100
colnames(total.correct) <- c("Method", "Taxonomic rank", "Total_reads", "Correct_reads", "Precision", "Sensitivity")


#calculate manhattan and canberra distances for genus and species levels using only true positive results
s.aga.man <- manhattan_dist(species.tp$ground_truth, species.tp$agamemnon)
s.kra.man <- manhattan_dist(species.tp$ground_truth, species.tp$kraken2)
s.bra.man <- manhattan_dist(species.tp$ground_truth, species.tp$bracken)
s.kai.man <- manhattan_dist(species.tp$ground_truth, species.tp$kaiju)
s.met.man <- manhattan_dist(species.tp$ground_truth, species.tp$metaphlan3)

g.aga.man <- manhattan_dist(genus.tp$ground_truth, genus.tp$agamemnon)
g.kra.man <- manhattan_dist(genus.tp$ground_truth, genus.tp$kraken2)
g.bra.man <- manhattan_dist(genus.tp$ground_truth, genus.tp$bracken)
g.kai.man <- manhattan_dist(genus.tp$ground_truth, genus.tp$kaiju)
g.met.man <- manhattan_dist(genus.tp$ground_truth, genus.tp$metaphlan3)

str.aga.man <- manhattan_dist(strains.tp$ground_truth, strains.tp$agamemnon)
str.kra.man <- manhattan_dist(strains.tp$ground_truth, strains.tp$kraken2)

s.aga.can <- canberra_dist(species.tp$ground_truth, species.tp$agamemnon)
s.kra.can <- canberra_dist(species.tp$ground_truth, species.tp$kraken2)
s.bra.can <- canberra_dist(species.tp$ground_truth, species.tp$bracken)
s.kai.can <- canberra_dist(species.tp$ground_truth, species.tp$kaiju)
s.met.can <- canberra_dist(species.tp$ground_truth, species.tp$metaphlan3)

g.aga.can <- canberra_dist(genus.tp$ground_truth, genus.tp$agamemnon)
g.kra.can <- canberra_dist(genus.tp$ground_truth, genus.tp$kraken2)
g.bra.can <- canberra_dist(genus.tp$ground_truth, genus.tp$bracken)
g.kai.can <- canberra_dist(genus.tp$ground_truth, genus.tp$kaiju)
g.met.can <- canberra_dist(genus.tp$ground_truth, genus.tp$metaphlan3)

str.aga.can <- canberra_dist(strains.tp$ground_truth, strains.tp$agamemnon)
str.kra.can <- canberra_dist(strains.tp$ground_truth, strains.tp$kraken2)

manhattan.dataframe <- data.frame(rep(c("AGAMEMNON", "Kraken 2", "Bracken", "Metaphlan 3", "Kaiju"), 3),
                                  c(rep("Genus", 5), rep("Species", 5), rep("Strains", 5)),
                                  c(g.aga.man, g.kra.man, g.bra.man, g.met.man, g.kai.man,
                                    s.aga.man, s.kra.man, s.bra.man, s.met.man, s.kai.man,
                                    str.aga.man, str.kra.man, NA, NA, NA))
colnames(manhattan.dataframe) <- c("Method", "Taxonomic rank", "Manhattan distance")

canberra.dataframe <- data.frame(rep(c("AGAMEMNON", "Kraken 2", "Bracken", "Metaphlan 3", "Kaiju"), 3),
                                 c(rep("Genus", 5), rep("Species", 5), rep("Strains", 5)),
                                 c(g.aga.can, g.kra.can, g.bra.can, g.met.can, g.kai.can,
                                   s.aga.can, s.kra.can, s.bra.can, s.met.can, s.kai.can,
                                   str.aga.can, str.kra.can, NA, NA, NA))
colnames(canberra.dataframe) <- c("Method2", "Taxonomic rank2", "Canberra distance")

distances <- cbind(manhattan.dataframe, canberra.dataframe)
distances$Method2 <- NULL; distances$`Taxonomic rank2` <- NULL
distances <- distances[complete.cases(distances), ]


###################################################################################################################################
# CREATE FIGURES AND TABLES  #
###################################################################################################################################

p.species.msle <- ggplot(melted.species.msle, aes(Reads_threshold, value, colour = variable)) + geom_point() +
  geom_line() + ylim(c(0, max(melted.species.msle$value))) + scale_x_continuous(breaks = seq(0, 1, 1), lim = c(0, 1)) +
  ylab(NULL) +
  xlab(NULL) +
  scale_color_manual("Methods:", values = c("#619CFF", "#F8766D", "#00BA38", "#F564E3", "#B79F00", "#FF8F00")) +
  theme(text = element_text(size = 13), axis.text.y = element_text(angle = 90)) +
  theme(plot.title = element_text(size = 12, face = "bold"))


p.genus.msle <- ggplot(melted.genus.msle, aes(Reads_threshold, value, colour = variable)) + geom_point() +
  geom_line() + ylim(c(0, max(melted.genus.msle$value))) + scale_x_continuous(breaks = seq(0, 1, 1), lim = c(0, 1)) +
  ylab(NULL) +
  xlab(NULL) +
  scale_color_manual("Methods:", values = c("#619CFF", "#F8766D", "#00BA38", "#F564E3", "#B79F00", "#FF8F00")) +
  theme(text = element_text(size = 13), axis.text.y = element_text(angle = 90)) +
  ggtitle("Synthetic") +
  theme(plot.title = element_text(size = 12, face = "bold"))


p.strain.msle <- ggplot(melted.strain.msle, aes(Reads_threshold, value, colour = variable)) + geom_point() +
  geom_line() + ylim(c(0, max(melted.strain.msle$value))) + scale_x_continuous(breaks = seq(0, 1, 1), lim = c(0, 1)) +
  ylab(NULL) +
  xlab("Reads threshold") +
  scale_color_manual("Methods:", values = c("#619CFF", "#F8766D", "#00BA38", "#F564E3", "#B79F00", "#FF8F00")) +
  theme(text = element_text(size = 13), axis.text.y = element_text(angle = 90)) +
  theme(plot.title = element_text(size = 12, face = "bold"))


p.species <- ggplot(melted.species.fp, aes(Reads_threshold, value, colour = variable)) +
  geom_line() + ylim(0, max(melted.species.fp$value)) + scale_x_continuous(breaks = seq(0, 10, 1), lim = c(0, 10)) + geom_point() +
  ylab(NULL) +
  xlab(NULL) +
  scale_color_manual("Methods:", values = c("#619CFF", "#F8766D", "#00BA38", "#F564E3", "#B79F00", "#FF8F00")) +
  theme(text = element_text(size = 13)) +
  theme(plot.title = element_text(size = 12, face = "bold"))


p.genus <- ggplot(melted.genus.fp, aes(Reads_threshold, value, colour = variable)) +
  geom_line() + ylim(c(0, max(melted.genus.fp$value))) + scale_x_continuous(breaks = seq(0, 10, 1), lim = c(0, 10)) + geom_point() +
  ylab(NULL) +
  xlab(NULL) +
  scale_color_manual("Methods:", values = c("#619CFF", "#F8766D", "#00BA38", "#F564E3", "#B79F00", "#FF8F00")) +
  theme(text = element_text(size = 13)) +
  ggtitle("Synthetic") +
  theme(plot.title = element_text(size = 12, face = "bold"))


p.strain <- ggplot(melted.strains.fp, aes(Reads_threshold, value, colour = variable)) +
  geom_line() + ylim(c(0, max(melted.strains.fp$value))) + scale_x_continuous(breaks = seq(0, 10, 1), lim = c(0, 10)) + geom_point() +
  ylab(NULL) +
  xlab("Reads threshold") +
  scale_color_manual("Methods:", values = c("#619CFF", "#F8766D", "#FF8F00", "#00BA38", "#F564E3")) +
  theme(text = element_text(size = 13)) +
  theme(plot.title = element_text(size = 12, face = "bold"))


tp.distances <- ggplot(distances, aes(x = `Manhattan distance`, y = `Canberra distance`, color = Method, shape = `Taxonomic rank`)) + 
  geom_point(size = 4) + scale_color_manual("Methods", values = c("#619CFF", "#00BA38", "#B79F00", "#F8766D", "#F564E3", "#FF8F00", "#B79F00")) +
  theme(text = element_text(size = 15)) +
  ylab("Canberra distance (only TP)") + xlab("Manhattan distance (only TP)") +
  scale_x_continuous(labels = function(x) format(x, scientific = TRUE)) +
  ggtitle("Synthetic") +
  theme(plot.title = element_text(size = 14, face = "bold"))


sensitivity.precision <- ggplot(total.correct, aes(x = Precision, y = Sensitivity, color = Method, shape = `Taxonomic rank`)) + 
  geom_point(size = 4) + scale_color_manual("Methods", values = c("#619CFF", "#00BA38", "#B79F00", "#F8766D", "#F564E3", "#FF8F00", "#B79F00")) +
  xlim(c(0, 100)) + ylim(c(0, 100)) +
  theme(text = element_text(size = 15)) +
  ylab("Sensitivity") + xlab("Specificity")


tp.distances <- tp.distances + ggtitle("Title") +
  theme(plot.title = element_text(size = 14, face = "bold"))

