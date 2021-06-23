#import libraries
library(dplyr)
library(Metrics)
library(ggplot2)
library(ggpubr)
library(stringr)
library(reshape2)
library(data.table)
library(RColorBrewer)
library(VennDiagram)
library(corrplot)


###################################################################################################################################
# IMPORT RESULTS INTO R  #
###################################################################################################################################


#import results files
agamemnon.file.one <- ""
kraken2.file.one <- ""
bracken.species.file.one = ""
bracken.genus.file.one = ""
metaphlan.file.one <- ""
kaiju.species.file.one = ""
kaiju.genus.file.one = ""


#import AGAMEMNON results and store into dataframe
agamemnon <- read.csv(file = agamemnon.file.one,
                      sep = "\t",
                      header = TRUE,
                      stringsAsFactors = FALSE,
                      quote = "", fill = FALSE)
#clean from leading/trailing white characters
agamemnon <- agamemnon %>% mutate(across(where(is.character), str_trim))
#convert counts column into numeric
agamemnon$NumCounts <- as.numeric(as.character(agamemnon$NumCounts))


#import KRAKEN 2 results and store into dataframe
kraken2 <- read.csv(file = kraken2.file.one,
                    sep = "\t",
                    header = FALSE,
                    stringsAsFactors = FALSE,
                    quote = "", fill = FALSE)
#create header names
colnames(kraken2) <- c("percent", "reads_total", "reads_direct", "taxonomic_rank", "taxid", "scientific_name")
#clean from leading/trailing white characters
kraken2 <- kraken2 %>% mutate(across(where(is.character), str_trim))
#convert counts column into numeric
kraken2$reads_total <- as.numeric(as.character(kraken2$reads_total))
kraken2$reads_direct <- as.numeric(as.character(kraken2$reads_direct))


#import bracken species/genus results and store into dataframes
bracken.species <- read.csv(file = bracken.species.file.one,
                            sep = "\t",
                            header = TRUE,
                            stringsAsFactors = FALSE,
                            quote = "", fill = FALSE)
#clean from leading/trailing white characters
bracken.species <- bracken.species %>% mutate(across(where(is.character), str_trim))
#convert counts column into numeric
bracken.species$new_est_reads <- as.numeric(as.character(bracken.species$new_est_reads))


bracken.genus <- read.csv(file = bracken.genus.file.one,
                          sep = "\t",
                          header = TRUE,
                          stringsAsFactors = FALSE)
#clean from leading/trailing white characters
bracken.genus <- bracken.genus %>% mutate(across(where(is.character), str_trim))
#convert counts column into numeric
bracken.genus$new_est_reads <- as.numeric(as.character(bracken.genus$new_est_reads))


#import kaiju species/genus results and store into dataframes
kaiju.species <- read.csv(file = kaiju.species.file.one,
                            sep = "\t",
                            header = TRUE,
                            stringsAsFactors = FALSE)
#clean from leading/trailing white characters
kaiju.species <- kaiju.species %>% mutate(across(where(is.character), str_trim))
#convert counts column into numeric
kaiju.species$reads <- as.numeric(as.character(kaiju.species$reads))


kaiju.genus <- read.csv(file = kaiju.genus.file.one,
                          sep = "\t",
                          header = TRUE,
                          stringsAsFactors = FALSE)
#clean from leading/trailing white characters
kaiju.genus <- kaiju.genus %>% mutate(across(where(is.character), str_trim))
#convert counts column into numeric
kaiju.genus$reads <- as.numeric(as.character(kaiju.genus$reads))


#import metaphlan3 results and store into dataframe
metaphlan3 <- read.csv(file = metaphlan.file.one,
                       sep = "\t",
                       header = TRUE,
                       stringsAsFactors = FALSE,
                       quote = "", fill = FALSE)
#clean metaphlan3 output from special characters
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
# EDIT DATAFRAMES TO GET GENUS AND SPECIES LEVEL COUNTS  #
###################################################################################################################################

#AGAMEMNON
agamemnon.species <- aggregate(agamemnon$NumCounts, by = list(Species = agamemnon$Species), FUN = sum)
agamemnon.genus <- aggregate(agamemnon$NumCounts, by = list(Genus = agamemnon$Genus), FUN = sum)
colnames(agamemnon.species) <- c("species", "counts"); colnames(agamemnon.genus) <- c("genus", "counts")


#KRAKEN 2
kraken2.species <- kraken2[(kraken2$taxonomic_rank == "S"), ]
kraken2.species <- kraken2.species %>% dplyr::select(scientific_name, reads_total)
colnames(kraken2.species) <- c("species", "counts")

kraken2.genus <- kraken2[(kraken2$taxonomic_rank == "G"), ]
kraken2.genus <- kraken2.genus %>% dplyr::select(scientific_name, reads_total)
colnames(kraken2.genus) <- c("genus", "counts")


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
species.merged <- Reduce(function(x, y) merge(x, y, by = "species", all = TRUE), list(agamemnon.species, kraken2.species,  bracken.species, kaiju.species, metaphlan3_species))
colnames(species.merged) <- c("species", "AGAMEMNON", "Kraken 2", "Bracken", "Kaiju", "Metaphlan 3")
species.merged[is.na(species.merged)] <- 0
#remove "---" entries from species-level results
species.merged <- species.merged[species.merged$species != "---", ]


#merge genus-level dataframes and replace NULL values with zero
genus.merged <- Reduce(function(x, y) merge(x, y, by = "genus", all = TRUE), list(agamemnon.genus, kraken2.genus, bracken.genus, kaiju.genus, metaphlan3_genus))
colnames(genus.merged) <- c("genus", "AGAMEMNON", "Kraken 2", "Bracken", "Kaiju", "Metaphlan 3")
genus.merged[is.na(genus.merged)] <- 0
#remove "---" entries from genus-level results
genus.merged <- genus.merged[genus.merged$genus != "---", ]


rownames(species.merged) <- species.merged$species
species.merged$species <- NULL
species.merged <- species.merged[rowSums(species.merged) > 0, ]

rownames(genus.merged) <- genus.merged$genus
genus.merged$genus <- NULL
genus.merged <- genus.merged[rowSums(genus.merged) > 0, ]


species.corr <- round(cor(species.merged, method = "spearman"), 2)
melted.corr.species <- melt(species.corr)

genus.corr <- round(cor(genus.merged, method = "spearman"), 2)
melted.corr.genus <- melt(genus.corr)


###################################################################################################################################
# CREATE FIGURES  #
###################################################################################################################################


#create Venn diagrams for species- and genus-level results
species.merged$species <- rownames(species.merged)
agamemnon.species.names <- species.merged[species.merged$AGAMEMNON > 0, ]$species
kaiju.species.names <- species.merged[species.merged$Kaiju > 0, ]$species
kraken2.species.names <- species.merged[species.merged$`Kraken 2` > 0, ]$species
bracken.species.names <- species.merged[species.merged$Bracken > 0, ]$species
metaphlan3.species.names <- species.merged[species.merged$`Metaphlan 3` > 0, ]$species

genus.merged$genus <- rownames(genus.merged)
agamemnon.genera.names <- genus.merged[genus.merged$AGAMEMNON > 0, ]$genus
kaiju.genera.names <- genus.merged[genus.merged$Kaiju > 0, ]$genus
kraken2.genera.names <- genus.merged[genus.merged$Kaiju > 0, ]$genus
bracken.genera.names <- genus.merged[genus.merged$Bracken > 0, ]$genus
metaphlan3.genera.names <- genus.merged[genus.merged$`Metaphlan 3` > 0, ]$genus


venn.species <- venn.diagram(
  x = list(agamemnon.species.names, kaiju.species.names, kraken2.species.names, bracken.species.names, metaphlan3.species.names),
  category.names = c("AGAMEMNON", "Kaiju", "Kraken 2", "Bracken", "Metaphlan 3"),
  filename = "",
  output = TRUE,
  
  imagetype = "tiff",
  height = 960, 
  width = 1280, 
  resolution = 300,
  compression = "lzw",
  main = "Species - SRS011084",
  main.fontfamily = "serif",
  main.cex = 0.6,
  main.pos = c(0.1, 1.05),
  
  lwd = 1,
  lty = "blank",
  fill = c(alpha("#440154ff", 0.7), alpha("#21908dff", 0.7), alpha("#fde725ff", 0.7), alpha("#31688eff", 0.7), alpha("#111111", 0.7)),
  
  cex = .6,
  cat.cex = .6,
  fontface = "bold",
  fontfamily = "sans",
  cat.fontface = "bold",
  cat.fontfamily = "sans",
  cat.default.pos = "outer",
  cat.pos = c(1, 1, -125, 130, 2)
)


venn.genus <- venn.diagram(
  x = list(agamemnon.genera.names, kaiju.genera.names, kraken2.genera.names, bracken.genera.names, metaphlan3.genera.names),
  category.names = c("AGAMEMNON", "Kaiju", "Kraken 2", "Bracken", "Metaphlan 3"),
  filename = "",
  output = TRUE,
  
  imagetype = "tiff",
  height = 960, 
  width = 1280, 
  resolution = 300,
  compression = "lzw",
  main = "Genus - SRS011084",
  main.fontfamily = "serif",
  main.cex = 0.6,
  main.pos = c(0.1, 1.05),
  
  lwd = 1,
  lty = "blank",
  fill = c(alpha("#e64b35b2", 0.7), alpha("#4dbbd5b2", 0.7), alpha("#00a087b2", 0.7), alpha("#7e6148b2", 0.7), alpha("#111111", 0.7)),
  
  cex = .6,
  cat.cex = .6,
  fontface = "bold",
  fontfamily = "sans",
  cat.fontface = "bold",
  cat.fontfamily = "sans",
  cat.default.pos = "outer",
  cat.pos = c(1, 1, -125, 130, 2)
)



#create corerlation matrix figures
species.corr.heatmap21 <- ggplot(data = melted.corr.species, aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile()
species.corr.heatmap21 <- species.corr.heatmap21 + 
  geom_text(aes(Var2, Var1, label = value), color = "black", size = 5) +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.ticks = element_blank(),
        axis.text.x = element_text(angle = 20, vjust = 0.65),
        text = element_text(size = 15, face = "bold")) +
  guides(fill = guide_colorbar(barwidth = 7, barheight = 1,
                               title.position = "top", title.hjust = 0.5)) + theme(legend.position = "none") +
  scale_fill_gradient2(mid = "#FBFEF9",low = "#0C6291",high = "#A63446", limits = c(-1, 1)) +
  ggtitle("Tilte")


genus.corr.heatmap21 <- ggplot(data = melted.corr.genus, aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile()
genus.corr.heatmap21 <- genus.corr.heatmap21 + 
  geom_text(aes(Var2, Var1, label = value), color = "black", size = 5) +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.ticks = element_blank(),
        axis.text.x = element_text(angle = 20, vjust = 0.65),
        text = element_text(size = 15, face = "bold")) +
  guides(fill = guide_colorbar(barwidth = 7, barheight = 1,
                               title.position = "top", title.hjust = 0.5)) + theme(legend.position = "none") +
  scale_fill_gradient2(mid = "#FBFEF9",low = "#0C6291",high = "#A63446", limits = c(-1, 1)) +
  ggtitle("Title")



