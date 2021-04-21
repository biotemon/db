#Plot_community 
#a script to plot stacked bar plots representing community analysis
#Version April 21 2021

rm(list=ls()) 
graphics.off() 
#Stacked bar for community struture GZ project.

library(plyr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(RColorBrewer)
library(data.table)

#Threshold of lowest percentage to show as an individual taxonomy
TH = SETTHRESHOLDHERE

#Set working directory
setwd("SETWORKINGDIRHERE")

#load data
taxonomyXcounts <- read.delim("SETTAXCOUNTSFILEHERE")

#Remove unknown terms
taxonomyXcounts<-taxonomyXcounts[!(taxonomyXcounts$SUPERKINGDOM=="Unknown"),]

desired_order = c("SETDESIREDBARORDERHERE")
sample_names = c("SETSAMPLENAMESHERE")
samples <- c("SETSAMPLENUMBERSHERE")

#Definition of my_level
#my_level ==> Superkingdom level 1
#my_level ==> Kingdom level 2
#my_level ==> Phylum level 3
#my_level ==> Class level 4
#my_level ==> Order level 5
#my_level ==> Family level 6
#my_level ==> Genus level 7
#my_level ==> Species level 8

#Find unique headers
taxones = c()
for (i in samples)
{
  taxonomyXcounts_sample <- filter(taxonomyXcounts, ASSEMBLY_ID == i)
  sums_by_taxon <- ddply(taxonomyXcounts_sample, .(TAX_ID), summarise, READ_COUNTS=sum(READ_COUNTS))
  vec <- as.vector(sums_by_taxon$TAX_ID)
  taxones = c(taxones, vec)
}

uniq_taxons = unique(taxones)

# Initialize matrix

n = length(samples)
m = length(uniq_taxons)

tax_counts = matrix(rep(0,n*m),nrow=n,ncol=m)


#Fill the matrix
for (i in 1:length(samples))
{
  taxonomyXcounts_sample <- filter(taxonomyXcounts, ASSEMBLY_ID == samples[i])
  sums_by_taxon <- ddply(taxonomyXcounts_sample, .(TAX_ID), summarise, READ_COUNTS=sum(READ_COUNTS))
  vec <- as.vector(sums_by_taxon$TAX_ID)
  vec2 <- as.vector(sums_by_taxon$READ_COUNTS)
  uniq_taxons[which(is.na(uniq_taxons))] <- 0
  vec[which(is.na(vec))] <- 0
  for (j in 1:length(vec)){
    for(k in 1:length(uniq_taxons)){
      if (vec[j] == uniq_taxons[k]){
        tax_counts[i,k] <- as.numeric(vec2[j])
      }
    }
    
  }
}

#Relative Abundance
SumSamples = rowSums(tax_counts)
relative_abundance = (tax_counts *100) / SumSamples

colnames(tax_counts) <- uniq_taxons
colnames(relative_abundance) <- uniq_taxons

my_relative = relative_abundance
my_absolute = tax_counts
my_level = 8
my_taxonomycounts = taxonomyXcounts


#We want to coalesce or merge those taxonomies with lower than a certain threshold (e.g. 0.5%) into a higher rank of taxonomy

#Subset higher ranks
high_ranks = setNames(data.frame(my_taxonomycounts[,2], my_taxonomycounts[,7], my_taxonomycounts[,8], my_taxonomycounts[,9], my_taxonomycounts[,10], my_taxonomycounts[,11], my_taxonomycounts[,12], my_taxonomycounts[,13], my_taxonomycounts[,14]), c("ID", "SUPERKINGDOM","KINGDOM", "PHYLUM", "CLASS", "ORDER", "FAMILY", "GENUS", "SPECIES"))
#remove repeated rows
high_ranks = high_ranks[!duplicated(high_ranks), ]
high_ranks = high_ranks[order(high_ranks$SUPERKINGDOM, high_ranks$KINGDOM, high_ranks$PHYLUM, high_ranks$CLASS, high_ranks$ORDER, high_ranks$FAMILY, high_ranks$GENUS), ]
#Set 0 in NA column as well as for empty taxons
if(any(is.na(high_ranks$ID))){
high_ranks[which(is.na(high_ranks$ID)),1] <- 0
}

#Extend last word to replace NAs and '' cells
for(i in 1:dim(high_ranks)[1]){
  last_term = ''
  for(j in 2:dim(high_ranks)[2]){
    if(!is.na(high_ranks[i,j]))
      {
      if(high_ranks[i,j] != "_"){
        if( !grepl("[[:space:]]+", high_ranks[i,j])){
          last_term <- as.character(high_ranks[i,j])
        }
      }
    }
   if(last_term != ''){
      my_levels <- levels(factor(high_ranks[,j]))
      my_levels[length(my_levels) + 1 ] <- last_term
      my_levels <- unique(my_levels)
      high_ranks[,j] <- factor(high_ranks[,j], levels = my_levels)
      high_ranks[i,j] <- last_term
   }
    
  }
} 

zerocols <- high_ranks[which(high_ranks[,2] == '' & high_ranks[,3] == '' & high_ranks[,4] == '' & high_ranks[,5] == '' & high_ranks[,6] == '' & high_ranks[,7] == '' & high_ranks[,8] == '' & high_ranks[,1] != 0),1]
zerovec <- as.vector(which(high_ranks[,2] == '' & high_ranks[,3] == '' & high_ranks[,4] == '' & high_ranks[,5] == '' & high_ranks[,6] == '' & high_ranks[,7] == '' & high_ranks[,8] == '' & high_ranks[,1] != 0))
if(any(zerovec)){
high_ranks <- high_ranks[-zerovec,]
}

if(length(levels(zerocols)) > 0){  
for(i in zerocols){
  colnames(my_absolute)[colnames(my_absolute) == i] <- 0
  colnames(my_relative)[colnames(my_relative) == i] <- 0
}
}

#Loop to coalesce columns names
my_level = my_level - 1
for(k in my_level:2)
{
#k = 2
  matrix_col_names = as.vector(colnames(my_relative))
  
  for(i in 1:dim(my_relative)[2])
  {
    #i = 8
    low_counter = 0
    for(j in 1:dim(my_relative)[1])
    {
      if(my_relative[j,i] < TH)
      {
        low_counter = low_counter + 1
      } 
    }
    if(low_counter == length(sample_names))
    {
      #Get the row in taxonomy where is the ID in the ID column.
      if(any(grep("^[0-9]+$",matrix_col_names[i], perl = TRUE, value=FALSE))){
        row_val = which(high_ranks[,1] == matrix_col_names[i])
        #Look for the 
        higher_name = as.character(high_ranks[row_val[1],k])
      }else{
        #higher_name = as.character(filter(high_ranks, CLASS == matrix_col_names[i])[[2]])
        my_k = k+1
        row_val = which(high_ranks[,my_k] == matrix_col_names[i])
        #Look for the 
        higher_name = as.character(high_ranks[row_val[1],k])
      }
      #Change the column ID
      if(!is.na(higher_name)){
        if(higher_name != ''){
          colnames(my_relative)[colnames(my_relative) == matrix_col_names[i]] <- higher_name
          colnames(my_absolute)[colnames(my_absolute) == matrix_col_names[i]] <- higher_name
        }
      }
      
    }
 }
  
  #Merging columns by the same name
  
  my_absolute = t(rowsum(t(my_absolute), colnames(my_absolute)))
  my_relative = t(rowsum(t(my_relative), colnames(my_relative)))
  
}

#Change names of genus
matrix_col_names = as.vector(colnames(my_relative))
for(i in 1:dim(my_relative)[2])
{
  if(any(grep("^[0-9]+$",matrix_col_names[i], perl = TRUE, value=FALSE))){
    row_val = which(high_ranks[,1] == matrix_col_names[i])
    higher_name = as.character(high_ranks[row_val[1],8])
    part2_name = as.character(high_ranks[row_val[1],9])
    if(!is.na(part2_name)){
    if(part2_name != ""){
    if(part2_name != higher_name){
      higher_name = paste(higher_name, part2_name, sep=' ')
    }}}
    colnames(my_relative)[colnames(my_relative) == matrix_col_names[i]] <- higher_name
    colnames(my_absolute)[colnames(my_absolute) == matrix_col_names[i]] <- higher_name
    
  }
}


#Merging columns by the same name... one more time just in case

my_absolute = t(rowsum(t(my_absolute), colnames(my_absolute)))
my_relative = t(rowsum(t(my_relative), colnames(my_relative)))

#Some times a column with a _ name gets here.
my_absolute = my_absolute[, colnames(my_absolute) != "_"]
my_relative = my_relative[, colnames(my_relative) != "_"]


###### Preparing Matrices for plotting ########

simple_absolute_matrix_2 = my_absolute
simple_relative_matrix_2 = my_relative

#Prepare Matrices for ploting
#Next line binds + create dataframe + keep numeric as numeric and dont change to text
simple_absolute_matrix_2 <- cbind.data.frame(sample_names, simple_absolute_matrix_2)
simple_relative_matrix_2 <- cbind.data.frame(sample_names, simple_relative_matrix_2)

#Reorder Columns by taxonomy

guide_vector = c()
for(i in 1:dim(high_ranks)[1])
  {
  for(j in 2:8)
    {
    my_value = as.character(high_ranks[i,j])
      if(is.na(match(my_value, guide_vector)) && (my_value != ""))
        {
          guide_vector = c(guide_vector, my_value)
        }
    }
  }

#Ordering the tables
headers_in_simple = colnames(simple_absolute_matrix_2)

#Remove 'sample_names' and 'V1' from headers_in_simple

headers_in_simple = headers_in_simple[headers_in_simple != "sample_names"]
headers_in_simple = headers_in_simple[headers_in_simple != "V1"]

#First we order the names later we append columns 
numeric_index=c()
names_index=c()

for(i in 1:length(headers_in_simple)){
  x_name <- as.vector(strsplit(headers_in_simple[i], " "))[[1]][1]
  xx = paste("\\b", x_name, "\\b", sep = '')
  if(any(grep(xx, guide_vector, perl = TRUE, value=FALSE))){
    x_val <- as.numeric(grep(xx, guide_vector, perl = TRUE, value=FALSE))
    numeric_index = c(numeric_index, x_val[1])
    names_index = c(names_index, headers_in_simple[i])
  }else{
    x_val = 1000000
    numeric_index = c(numeric_index, x_val)
    names_index = c(names_index, headers_in_simple[i])
  }
}

the_index_df = cbind.data.frame(numeric_index, names_index)
the_index_df = the_index_df[order(numeric_index),]


#READ THE NEXT THREE LINES#

suited_guide_vector <- as.vector(the_index_df$names_index)
if(colnames(simple_absolute_matrix_2)[2] == "V1"){
suited_guide_vector <- as.vector(c(c("sample_names", "V1"), suited_guide_vector)) #Run this if there is a NA / 0 column in the high ranks
}else{
  suited_guide_vector <- as.vector(c(c("sample_names"), suited_guide_vector))
}

suited_guide_vector <-suited_guide_vector[!is.na(suited_guide_vector)]

simple_absolute_matrix_3 = simple_absolute_matrix_2[suited_guide_vector]
simple_relative_matrix_3 = simple_relative_matrix_2[suited_guide_vector]

colnames(simple_absolute_matrix_3)[colnames(simple_absolute_matrix_3) == "V1"] <- "Others"
colnames(simple_relative_matrix_3)[colnames(simple_relative_matrix_3) == "V1"] <- "Others"

#Now dataframes are simplified and ordered
#Reshaping data frame
#Making the dataframe from wide to log format

simple_absolute_melt <- data.table::melt(as.data.table(simple_absolute_matrix_3), id.vars='sample_names')
simple_relative_melt <- data.table::melt(as.data.table(simple_relative_matrix_3), id.vars='sample_names')

#Colors
xval =  dim(simple_relative_matrix_3)[2] - 1
xvalA = ceiling(xval/3)
xvalB = xval - xvalA - xvalA
xvalC = xval - xvalA - xvalB
colfuncA <- colorRampPalette(brewer.pal(11,"RdYlBu"))
colfuncB <- colorRampPalette(brewer.pal(11,"PiYG"))
colfuncC <- colorRampPalette(brewer.pal(11,"BrBG"))
simple_color_vecA <- colfuncA(xvalA)
simple_color_vecB <- colfuncB(xvalB)
simple_color_vecC <- colfuncC(xvalC)
simple_color_vec <- c(simple_color_vecA,simple_color_vecB,simple_color_vecC)
#If you need to change a color simply call
#simple_color_vec[i] <- "#HEXCOLOR"
simple_color_vec <- c(simple_color_vecA,simple_color_vecB,simple_color_vecC)
#SETCOLORSHERE

#Sorting bars (x-axis wise) by sorting the levels of the "Treatment" column
simple_absolute_melt$sample_names = factor(simple_absolute_melt$sample_names,levels = desired_order)
simple_relative_melt$sample_names = factor(simple_relative_melt$sample_names,levels = desired_order)

pdf("taxonomy_abs_cutoff_SETTHRESHOLDHERE.pdf", width=12, height=7)
ggplot(data = simple_absolute_melt, aes(x = sample_names, y = value, fill = variable)) + geom_bar(colour="black", stat = "identity", size = 0.25) + theme_classic() + theme(axis.text.x = element_text(color="black", angle = 90, hjust = 1),axis.text.y = element_text(color="black")) + scale_fill_manual(values = simple_color_vec) + scale_y_continuous(name="Read Counts", labels = scales::comma, expand = c(0, 0)) + guides(fill=guide_legend(ncol=1)) + xlab("Treatments")
dev.off()

pdf("taxonomy_rel_cutoff_SETTHRESHOLDHERE.pdf", width=10, height=8)
ggplot(data = simple_relative_melt, aes(x = sample_names, y = value, fill = variable)) + geom_bar(colour="black", stat = "identity", size = 0.25) + theme_classic() + theme(axis.text.x = element_text(color="black", angle = 90, hjust = 1),axis.text.y = element_text(color="black")) +  scale_colour_manual("black") + scale_fill_manual(values = simple_color_vec) + scale_y_continuous(name="Proportion of Read Counts", labels = scales::comma, expand = c(0, 0)) + guides(fill=guide_legend(ncol=1)) + xlab("Treatments")
dev.off()

svg("taxonomy_abs_cutoff_SETTHRESHOLDHERE.svg", width=12, height=7)
ggplot(data = simple_absolute_melt, aes(x = sample_names, y = value, fill = variable)) + geom_bar(colour="black", stat = "identity", size = 0.25) + theme_classic() + theme(axis.text.x = element_text(color="black", angle = 90, hjust = 1),axis.text.y = element_text(color="black")) + scale_fill_manual(values = simple_color_vec) + scale_y_continuous(name="Read Counts", labels = scales::comma, expand = c(0, 0)) + guides(fill=guide_legend(ncol=1)) + xlab("Treatments")
dev.off()

svg("taxonomy_rel_cutoff_SETTHRESHOLDHERE.svg", width=10, height=8)
ggplot(data = simple_relative_melt, aes(x = sample_names, y = value, fill = variable)) + geom_bar(colour="black", stat = "identity", size = 0.25) + theme_classic() + theme(axis.text.x = element_text(color="black", angle = 90, hjust = 1),axis.text.y = element_text(color="black")) +  scale_colour_manual("black") + scale_fill_manual(values = simple_color_vec) + scale_y_continuous(name="Proportion of Read Counts", labels = scales::comma, expand = c(0, 0)) + guides(fill=guide_legend(ncol=1)) + xlab("Treatments")
dev.off()

#Write a csv file with the actual numbers shown in the simplified plots.
write.csv(simple_absolute_melt, file = "simple_absolute_melt.csv")
write.csv(simple_relative_melt, file = "simple_relative_melt.csv")
