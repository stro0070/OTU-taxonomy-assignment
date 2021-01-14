##Compare taxonomy assignment from amptk with taxonomy assignment from consensus Blast, and keep the more precise / deepest level of taxonomy, e.g. species-level taxonomy of all classifications.

#Clear R's brain.
rm(list = ls())

##Set working directory and import table with taxonomy assigned to OTUs by two or more methods of taxonomy assignment.  Order of columns should be in order of "preferred" taxonomy assignment method.
setwd("~/OneDrive - ORASURE TECHNOLOGIES/Noah_personal/Excel workbooks/")

taxonomy <- read.table("mk5-S3_taxonomy_assignments.txt", comment="", header=TRUE, sep="\t", as.is=TRUE, check.names=F) 


#view imported taxonomy dataframe.
head(taxonomy)


##Replace all other taxonomy assignment methods with SynMock taxonomy assignment if Synthetic Mock Community sequences were assigned in AMPTK.  (I added these sequences to the database "on the fly.")
#taxonomy[258,'AMPTK'] #View an example of a Synmock assignment

for(i in 1:nrow(taxonomy)){
  if(grepl('|SYN', taxonomy$AMPTK[i], fixed = TRUE)){
    synmock <- strsplit(taxonomy$AMPTK[i], '|', fixed = TRUE)[[1]][3]
    synmock <- strsplit(synmock, ';')[[1]][1]
    taxonomy$Consensus.20[i] <- paste0('k__', synmock, ';p__', synmock, ';c__', synmock, ';o__', synmock, ';f__', synmock, ';g__', synmock, ';s__', synmock)
    taxonomy$Consensus.10[i] <- paste0('k__', synmock, ';p__', synmock, ';c__', synmock, ';o__', synmock, ';f__', synmock, ';g__', synmock, ';s__', synmock)
    taxonomy$Consensus.1[i] <- paste0('k__', synmock, ';p__', synmock, ';c__', synmock, ';o__', synmock, ';f__', synmock, ';g__', synmock, ';s__', synmock)
    taxonomy$AMPTK[i] <- paste0('k__', synmock, ';p__', synmock, ';c__', synmock, ';o__', synmock, ';f__', synmock, ';g__', synmock, ';s__', synmock)
  }

}
#taxonomy[258,] #View an example of a Synmock assignment


###Reformat table

##Get rid of Qiime information that I do not need.

#unassigned example.
taxonomy[59,]
taxonomy <- gsub("Unassigned", "", as.matrix(taxonomy))  
#head(taxonomy)
taxonomy[59,]

#unidentified
taxonomy[60,5]
taxonomy <- gsub("..__unidentified*", "", as.matrix(taxonomy))
#head(taxonomy)
taxonomy[60,5]

#Eliminate species information if we only know a higher level of taxonomy
#example
taxonomy[54,5]
taxonomy <- gsub(";s__.*._sp", "", as.matrix(taxonomy)) 
#head(taxonomy)
taxonomy[54,5]


##Reformat AMPTK taxonomy to match Qiime.

#put an underscore in place of a space to match Qiime genus_species format
head(taxonomy[,'AMPTK'])
taxonomy <- gsub(" ", "_", as.matrix(taxonomy))   
head(taxonomy[,'AMPTK'])

#reformat labels for kingdom
taxonomy <- gsub(".*;k:", "k__", as.matrix(taxonomy))
head(taxonomy[,'AMPTK'])

#phylum
taxonomy <- gsub(",p:", ";p__", as.matrix(taxonomy))
head(taxonomy[,'AMPTK'])

#class
taxonomy <- gsub(",c:", ";c__", as.matrix(taxonomy))
head(taxonomy[,'AMPTK'])

#order
taxonomy <- gsub(",o:", ";o__", as.matrix(taxonomy))
head(taxonomy[,'AMPTK'])

#family
taxonomy <- gsub(",f:", ";f__", as.matrix(taxonomy))
head(taxonomy[,'AMPTK'])

#genus
taxonomy <- gsub(",g:", ";g__", as.matrix(taxonomy))
head(taxonomy[,'AMPTK'])

#species
taxonomy <- gsub(",s:", ";s__", as.matrix(taxonomy))
head(taxonomy[,'AMPTK'])


#show taxonomy assignment methods
labels(taxonomy)[[2]][2:ncol(taxonomy)]  

#number of taxonomy assignment methods
number.of.methods <- ncol(taxonomy) - 1   

#Save OTU_IDs for later.
OTU_IDs <- taxonomy[,1]  


##First, solve conflicts between taxonomy assignment methods, trusting assignments in left-most columns over right-most columns.  (e.g. Consensus 20 > AMPTK > Consensus 10 > Best Blast Hit (Consensus 1))

#example
taxonomy[29,]

solveConflicts <- function(taxonomicLevel){
  for(i in 1:nrow(taxonomy)){
    split.string <- list()
    split.string <- strsplit(taxonomy[i,], ';')  
    for(k in 2:(number.of.methods)){  #start with k=2 to avoid the OTU_ID column.  End with second to last column.  (compare with column to the right of this one.)
      for(j in k:(number.of.methods)){ #k begins with k=j then moves up with j
        if(length(split.string[[k]]) >= taxonomicLevel & length(split.string[[j+1]]) >= taxonomicLevel){ #split strings must be at least length 7 to include species
          if(split.string[[k]][taxonomicLevel] != split.string[[j+1]][taxonomicLevel]){ #e.g. 6 refers to the taxonomicLevel of genus. 7 refers to species.
            split.string[[j+1]][taxonomicLevel:7] <- 'N' #replace species names of taxonomy assignment methods to the right with 'N' if there is disagreement with taxonomy assignment methods to the left (left-most assignment methods are assumed to be more accurate.)
            taxonomy[i,j+1] <- paste(split.string[[j+1]], collapse = ';')  #replace in taxonomy vector
          }
        }
      }
      
    }
  }
  return(taxonomy)
}


##Solve taxonomic assignment conflicts. More reliable methods override taxonomic assignments by less reliable methods. (Methods are ordered by reliability in input Excel file.) Overriden classifications are replaced by "N."

#Species level conflict resolution.
#example
taxonomy[29,]
taxonomy <- solveConflicts(7) #7 refers to species level.
taxonomy[29,]

#Genus level conflict resolution.
#example
taxonomy[29,]
taxonomy <- solveConflicts(6) #6 refers to genus level.
taxonomy[29,]

#Family level conflict resolution.
#example
taxonomy[29,]
taxonomy <- solveConflicts(5) #5 refers to family level.
taxonomy[29,]

#Order level conflict resolution.
#example
taxonomy[9,]
taxonomy <- solveConflicts(4) #4 refers to order level.
taxonomy[9,]

#Class level conflict resolution.
#example
taxonomy[9,]
taxonomy <- solveConflicts(3) #3 refers to class level.
taxonomy[9,]

#Phylum level conflict resolution.
#example
taxonomy[9,]
taxonomy <- solveConflicts(2) #2 refers to phylum level.
taxonomy[9,]

#Kingdom level conflict resolution.
#example
taxonomy[9,]
taxonomy <- solveConflicts(1) #1 refers to kingdom level.
taxonomy[9,]



##Create and fill a vector with most complete taxonomies, beginning with those that contain species-level taxonomy, biased in order of methods I trust the most  (e.g. Consensus 20 > AMPTK > Consensus 10 > Best Blast Hit (Consensus 1)).

#Create a vector to store the most precise taxonomy.
longest.taxonomy <- vector()

#function to retain deepest level of non-conflicting taxonomy, with preference given to more reliable methods.
retainDeepest <- function(taxonomicLevel){ #e.g. taxonomicLevel = "s__"
  
  for(j in 2:(number.of.methods+1)){ #number.of.methods+1 is the number of columns in taxonomy
    for(i in 1:nrow(taxonomy)){
      #Retain taxonomicLevel (e.g. species), if best method has one.
      if(grepl(taxonomicLevel, taxonomy[i,j])){
        longest.taxonomy[i] <- taxonomy[i,j]
        
        #Change taxonomy from less reliable methods to "N"
        if(j == 2){ #j == 2 is the most reliable (left-most) column
          taxonomy[i,(j+1):(number.of.methods+1)] <- "N"
        }
        else{
          
          #If a taxonomicLevel (e.g. species) assignment is given by any method, keep it, and replace all other methods with "N".
          if(j < (number.of.methods+1)){ #j is not the first column
            taxonomy[i,(j+1):(number.of.methods+1)] <- "N"   #Change taxonomy of methods to right of current method to "N"
            taxonomy[i,2:(j-1)] <- "N"    #Change taxonomy of methods to left of current method to "N"
          }
          
          #As a last resort, retain the taxonomicLevel (e.g. species) assignment of the least reliable method.
          else{
            taxonomy[i,2:(j-1)] <- "N"  #Change taxonomy of methods to left of rightmost method to "N"
          }
        }
      }
    }
  }
  return(list(longest.taxonomy, taxonomy))
}

#Retain the species level taxonomy, wherever possible.
longest.taxonomy <- retainDeepest("s__")[[1]]
taxonomy <- retainDeepest("s__")[[2]]

head(longest.taxonomy)
head(taxonomy)


#Retain the genus level taxonomy, wherever possible.
longest.taxonomy <- retainDeepest("g__")[[1]]
taxonomy <- retainDeepest("g__")[[2]]

head(longest.taxonomy)
head(taxonomy)


#Retain the family level taxonomy, wherever possible.
longest.taxonomy <- retainDeepest("f__")[[1]]
taxonomy <- retainDeepest("f__")[[2]]

tail(longest.taxonomy)
tail(taxonomy)


#Retain the order level taxonomy, wherever possible.
longest.taxonomy <- retainDeepest("o__")[[1]]
taxonomy <- retainDeepest("o__")[[2]]

tail(longest.taxonomy)
tail(taxonomy)


#Retain the class level taxonomy, wherever possible.
longest.taxonomy <- retainDeepest("c__")[[1]]
taxonomy <- retainDeepest("c__")[[2]]

tail(longest.taxonomy)
tail(taxonomy)


#Retain the phylum level taxonomy, wherever possible.
longest.taxonomy <- retainDeepest("p__")[[1]]
taxonomy <- retainDeepest("p__")[[2]]

tail(longest.taxonomy)
tail(taxonomy)


#Retain the kingdom level taxonomy, wherever possible.
longest.taxonomy <- retainDeepest("k__")[[1]]
taxonomy <- retainDeepest("k__")[[2]]

tail(longest.taxonomy)
tail(taxonomy)


##Add column to taxonomy dataframe.
taxonomy <- as.data.frame(taxonomy)
taxonomy$longest.taxonomy <- longest.taxonomy
taxonomy$`#OTU ID` <- OTU_IDs
#row.names(taxonomy) <- OTU_IDs
head(taxonomy)

head(taxonomy$longest.taxonomy) #View longest taxonomy vector




###Modify taxonomy so it works in QIIME! QIIME generally reports all levels of taxonomy.

#make new column to hold qiime-formatted taxonomy
taxonomy$qiime.formatted.longest.taxonomy <- taxonomy$longest.taxonomy

#view deficient taxon.
taxonomy$qiime.formatted.longest.taxonomy[59]

#Function to begin to make taxonomies QIIME-formatted. This adds text like ";p__unidentified" to fill in taxonomic levels.
qiimeFormat <- function(taxonomicLevel){ #e.g. "p__"
  for(i in 1: length(taxonomy$qiime.formatted.longest.taxonomy)){
    if(grepl(taxonomicLevel, taxonomy$qiime.formatted.longest.taxonomy[i]) == 'FALSE'){
      taxonomy$qiime.formatted.longest.taxonomy[i] <- paste(taxonomy$qiime.formatted.longest.taxonomy[i], ";", taxonomicLevel, "unidentified", sep = "")
    }
  }
  return(taxonomy)
}

#add ;p_unidentified
taxonomy$qiime.formatted.longest.taxonomy[59]
taxonomy <- qiimeFormat("p__")
taxonomy$qiime.formatted.longest.taxonomy[59]

#add ;c_unidentified
taxonomy <- qiimeFormat("c__")
taxonomy$qiime.formatted.longest.taxonomy[59]

#add ;o_unidentified
taxonomy <- qiimeFormat("o__")
taxonomy$qiime.formatted.longest.taxonomy[59]

#add ;f_unidentified
taxonomy <- qiimeFormat("f__")
taxonomy$qiime.formatted.longest.taxonomy[59]

#add ;g_unidentified
taxonomy <- qiimeFormat("g__")
taxonomy$qiime.formatted.longest.taxonomy[59]


#Function to add ";s__" and "_sp" to the end of the taxonomy if it is not present.  Note:  This function is unique!  It adds "_sp" to the lowest assigned level of taxonomy.
addSp <- function(numericTaxonomicLevel, characterTaxonomicLevel){ #e.g. 6, "g"
  for(i in 1: length(taxonomy$qiime.formatted.longest.taxonomy)){
    if(strsplit(taxonomy$qiime.formatted.longest.taxonomy[i], split = ";")[[1]][numericTaxonomicLevel] != paste0(characterTaxonomicLevel, "__unidentified") & length(strsplit(taxonomy$qiime.formatted.longest.taxonomy[i], split = ";")[[1]]) == 6){
      taxonLevel <- strsplit(taxonomy$qiime.formatted.longest.taxonomy[i], split = ";")[[1]][numericTaxonomicLevel]
      taxonLevel <- sub(".*__", "", taxonLevel)
      taxonomy$qiime.formatted.longest.taxonomy[i] <- paste(taxonomy$qiime.formatted.longest.taxonomy[i], ";s__", taxonLevel, "_sp", sep = "")
    }
  }
  return(taxonomy)
}


#Add "_sp" to genera.
taxonomy$qiime.formatted.longest.taxonomy[57]
taxonomy <- addSp(6, "g")
taxonomy$qiime.formatted.longest.taxonomy[57]

#Add "_sp" to family
taxonomy$qiime.formatted.longest.taxonomy[48]
taxonomy <- addSp(5, "f")
taxonomy$qiime.formatted.longest.taxonomy[48]

#Add "_sp" to order
taxonomy$qiime.formatted.longest.taxonomy[55]
taxonomy <- addSp(4, "o")
taxonomy$qiime.formatted.longest.taxonomy[55]

#Add "_sp" to class
taxonomy <- addSp(3, "c")

#Add "_sp" to phylum
taxonomy <- addSp(2, "p")

#Add "_sp" to kingdom
taxonomy$qiime.formatted.longest.taxonomy[59]
taxonomy <- addSp(1, "k")
taxonomy$qiime.formatted.longest.taxonomy[59]


#Check entire taxonomy column for all taxonomic levels.  (All should add to zero.)
sum(grepl("k__", taxonomy$qiime.formatted.longest.taxonomy) != 'TRUE')
sum(grepl("p__", taxonomy$qiime.formatted.longest.taxonomy) != 'TRUE')
sum(grepl("c__", taxonomy$qiime.formatted.longest.taxonomy) != 'TRUE')
sum(grepl("o__", taxonomy$qiime.formatted.longest.taxonomy) != 'TRUE')
sum(grepl("f__", taxonomy$qiime.formatted.longest.taxonomy) != 'TRUE')
sum(grepl("g__", taxonomy$qiime.formatted.longest.taxonomy) != 'TRUE')
sum(grepl("s__", taxonomy$qiime.formatted.longest.taxonomy) != 'TRUE')

#View final product.
View(taxonomy)

#Write to Excel
write.csv(taxonomy, file = '~/OneDrive - ORASURE TECHNOLOGIES/Noah_personal/Excel workbooks/taxonomy.with.longest.taxonomy.2020.1.12.csv', row.names = TRUE)
