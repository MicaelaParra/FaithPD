## PHYLOGENETIC DIVERSITY ##

# Parra, Micaela, Laufer, Natalia, Manrique, Julieta M., Jones, Leandro R., Quarleri, Jorge, Phylogenetic Diversity in Core Region of Hepatitis C Virus Genotype 1a as a Factor Associated 
# with Fibrosis Severity in HIV-1-Coinfected Patients, BioMed Research International, 2017, 1728456, 12 pages, 2017. https://doi.org/10.1155/2017/1728456 

# How to run:
# R
# setwd("/path/to/archives/")  # corrected 'setwc' to 'setwd'
# source("pd.R")
# analyze_pd("sample_file", "tree_file", "output_name_file")
# Example:
# analyze_pd("patient3.sample.txt", "patient3.tree", "results_example_patient3.csv")

# install.packages("picante")
# install.packages("ape")

library(picante)
library(ape)

# Function
analyze_pd <- function(sample_file, tree_file, output_file) {
  # Read sample data
  sample <- read.table(sample_file, sep = ",", header = TRUE, row.names=1)
  
  # Read trees
  t <- read.tree(tree_file)
  
  # Number of trees in tree archive
  num_trees <- length(t)
  print(paste("Number of trees:", num_trees))
  
  results <- data.frame(sample_category = character(), index = integer(), PD = numeric(), len = numeric(), nPD = numeric(), stringsAsFactors = FALSE)
  
  # Calculate nPD for each tree
  for (index in 1:num_trees) {
    tree <- t[[index]]
    
    # Verify valid tree
    if (is(tree, "phylo")) {
      # Calculate PD
      pd_int <- pd(sample, tree, include.root=FALSE)
      
      # Calculate tree length
      len <- sum(tree$edge.length)
      
      # Calculate nPD as PD/len
      npd_value <- pd_int$PD / len
      
      # Create data frame with results
      temp_df <- data.frame(sample_category = rownames(sample),
                            index = rep(index, nrow(sample)),
                            PD = pd_int$PD,
                            len = rep(len, nrow(sample)),
                            nPD = pd_int$PD / len)
      
      results <- rbind(results, temp_df)
      
    } else {
      print(paste("The object at index", index, "is not a valid tree."))
    }
  }
  
  # Save results in CSV file
  write.csv(results, file = output_file, row.names = FALSE)
}

#----------------------------------------

# Verify some results such as:

# pd(sample,t[[1]],include.root=FALSE)  # where t[[number]] is the number of bootstrap tree
#        PD SR
#TA 0.11393  5
#TB 0.01945  9
#TC 0.16930 14
#TD 0.14851  2
#TE 0.09794 15

# sum(t[[1]]$edge.length)   # where t[[number]] is the number of bootstrap tree
#[1] 0.38872

#----------------------------------------

## BOXPLOT ##

# Load library
#library(ggplot2)

# Read data in CSV file
#x <- read.table("results_example_patient3.csv", sep = ",", header = T)

#head(x)
#  sample_category index      PD     len        nPD     # each index corresponds to each tree
#1              TA     1 0.11393 0.38872 0.29309014
#2              TB     1 0.01945 0.38872 0.05003602
#3              TC     1 0.16930 0.38872 0.43553200
#4              TD     1 0.14851 0.38872 0.38204878
#5              TE     1 0.09794 0.38872 0.25195513
#6              TA     2 0.20447 0.45803 0.44641181

# Define fill and line colors
#fill <- "#4271AE"   # Example color
#line <- "#1F3552"   # Example color

# Boxplot
#pdf("pd_example.pdf")

#ggplot(x, aes(x=sample_category, y=nPD, fill=sample_category)) + 
#    geom_boxplot(colour = line, alpha=0.7, outlier.colour = "#1F3552", outlier.shape = 20) +
#    xlab("sample_category") +
#    ylab("nPD") +
#    ggtitle("Phylogenetic diversity by sample_category") +
#    theme_bw()

#dev.off()

