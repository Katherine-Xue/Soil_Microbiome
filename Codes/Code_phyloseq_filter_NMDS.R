library(phyloseq)
library(dplyr)
library(readr)
library(tidyr)
library(vegan)

setwd("./")

## filter ps
otu <- read_csv("./otu.Bac.csv" )
tax <- read_csv("./tax.Bac.csv")
sam <- read_csv("./metadata.csv")

library(phyloseq)
taxa <- tax_table(tax)
colnames(taxa) <- c("Kingdom","Phylum","Class","Order","Family","Genus","Species")
ps <- phyloseq(otu_table(otu,taxa_are_rows=TRUE),taxa,sample_data(sam))
ps

### Filter the phyloseq object
ps1 <- prune_samples(sample_sums(ps)>10000,ps)
ps1<- filter_taxa(ps1,function(x)sum(x)>=10,TRUE)
ps1 <- filter_taxa(ps1,function(x)sum(x>0)>2,TRUE)
ps1

###

set.seed(1)
Bac.NMDS <- ordinate(Bac1, "NMDS", "bray", k=3, try=20,trymax=100)
set.seed(1)
Fg.NMDS <- ordinate(Fg1, "NMDS", "bray", k=3, try=20,trymax=100)

save.image(file="ps_NMDS.RData")

###########

