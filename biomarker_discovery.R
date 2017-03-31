# load needed libraries
library(vegan)
library(ggplot2)
library(reshape2)
library(mvabund)
library(randomForest)

# Useful functions for writing lefse tables.
prep.lefse <- function(m.mod, gvar){
  m.mod.labled <- rbind(gvar,t(m.mod))
  # add the colnames
  m.mod.labled <- rbind(colnames(m.mod.labled),m.mod.labled)
  row.names(m.mod.labled)[1] <- "fields"
  row.names(m.mod.labled) <- gsub(" ","_", row.names(m.mod.labled))
  # add the rownames
  m.mod.labled <- cbind(row.names(m.mod.labled), m.mod.labled)
  return(m.mod.labled)
}

write.lefse <- function(in.table, fname){
  write.table(in.table, fname, sep="\t", row.names=F, col.names=F, quote=F)
}


########################
setwd("C:\\Users\\WLab\\Desktop\\work\\A_projects\\a_wei_lab_microbiome\\a_projects\\microbiom_R_plotting")
#setwd("path/to/data")

# load test data
species.table <- read.delim("microbiome_test_data_species.tsv", sep="\t", header=T, row.names = 1)
# get the groups from the first row, and covert it to a simple list (called a vector)
groups <- as.vector(unlist(species.table[1,]))

# remove the first row to get the table and convert to a matrix
# remove row 1
tmp.matrix <- species.table[-1,]
# Convert data to numberic types
species.counts <- apply(as.matrix.noquote(tmp.matrix),2,as.numeric)
# reset the correct column and row names
rownames(species.counts) <- rownames(tmp.matrix)
colnames(species.counts) <- colnames(tmp.matrix)

# optional filtering
species.filtered <- species.counts
species.filtered[species.filtered < 2000 ] <- 0
# remove zero rows
rsums <- rowSums(species.filtered)
species.filtered <- species.filtered[-(which(rsums == 0)),]
species.counts <- species.filtered

# calculate relative abundance (vegan's decostand). This transposes the data.
species.relab <- decostand(t(species.counts), "total")


#############################
# Formatting tables for LEfSe
lefse.table <- prep.lefse(species.relab, groups)
write.lefse(lefse.table, "species_table_lefse.tsv")

# once formatted for LEfSe, you can call lefse with the following commands
#$ format_input.py species_table_lefse.tsv lefse_4group.in -c 2 -u 1 -o 1000000
#$ run_lefse.py lefse_4group.in lefse_4group.res
#$ plot_res.py --format pdf lefse_4group.res lefse_4group.pdf --width 10
# or:$ plot_res.py --format png lefse_4group.res lefse_4group.png --width 10

######################################################
# Random Forest classification and varaible importance
# Make a data frame where species are columns, t()
frame.species <- data.frame(species.relab,group=groups)
# train a random forest model.
set.seed(42)
rf.species <- randomForest(group ~ ., data=frame.species, importance=TRUE, do.trace=100)
rf.species

# make a new data.frame with the variable importance, mean decrease in accuracy
rf.frame <- data.frame(names=rownames(rf.species$importance),
                       rf.species$importance,
                       sdev=rf.species$importanceSD[,"MeanDecreaseAccuracy"])

# do some filtering to remove useless variables
rf.frame <- subset(rf.frame, MeanDecreaseAccuracy > 0)
sort.order <- order(rf.frame$MeanDecreaseAccuracy)
rf.frame <- within(rf.frame, 
                   names <- factor(names, 
                                   levels=names[sort.order]))

# Plot varaible importance for top 20 varaibles in the random forest
#png("plots/random_forest_varaible_importance.png",width=5, height=7,unit="in",res=300)
ggplot(rf.frame[rev(sort.order)[1:20],], aes(y=names, x=MeanDecreaseAccuracy))+
  labs(list(title = "Random Forest Variable Importance",
            x = "Mean Decrease in Accuracy",
            y = "Species")) +
  geom_point() +
  geom_errorbarh(aes(xmax=MeanDecreaseAccuracy+sdev, xmin=MeanDecreaseAccuracy-sdev, height=0)) +
  # Center the plot title
  theme(plot.title = element_text(hjust = 0.5))
#dev.off()  






