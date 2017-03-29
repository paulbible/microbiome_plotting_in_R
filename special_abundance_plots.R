# load needed libraries
library(vegan)
library(ggplot2)
library(reshape2)

########################
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
species.filtered <- as.matrix(species.counts)
species.filtered[species.filtered < 2000 ] <- 0
# remove zero rows
rsums <- rowSums(species.filtered)
species.filtered <- species.filtered[-(which(rsums == 0)),]
species.counts <- species.filtered

# calculate relative abundance (vegan's decostand). This transposes the data.
species.relab <- decostand(t(species.counts), "total")


###################################
# dot abundance plot across samples
# Select some high abundance species
species.of.interest.index <- which(colSums(species.relab) > 0.33)
# select only the abundance for the species of interest
species.relab.subset <- species.relab[,species.of.interest.index]
# Create a data.frame with the abundances, names, and groups
species.frame <- data.frame(sample=rownames(species.relab.subset),
                            species.relab.subset,
                            group=groups)

# use melt to restructure the data.frame in a way that ggplot2 likes.
species.melt.table <- melt(species.frame, id=c("sample", "group"))
colnames(species.melt.table) <- c("sample", "group","species", "abundance")
## Set factor order, this orders the sample labels to display in regular order (not alphabetically).
species.melt.table <- within(species.melt.table, 
                             sample <- factor(sample,
                                               levels=rownames(species.frame)))

# use ggplot2 to plot each abundance as a cricle
ggplot(species.melt.table,
       aes(x=sample,
           y=species,
           size=abundance)) +
  labs(list(title = "Dot Abundance Plot Title",
            x = "Sample",
            y = "Species",
            color="Legend Title")) +
  geom_point(aes(x=sample, color=factor(group))) +
  scale_color_discrete(breaks = c("group_1", "group_2", "group_3", "group_4"),
                      labels = c("This is G1", "This is G2", "This is G3", "This is G4")) +
  scale_size_area() +
  theme(axis.text.x = element_text(angle = 90, hjust=1, vjust=0.40),
        plot.title = element_text(hjust = 0.5))


###############################################
# Higher clade level abundance with nice colors

# A nice list of colors
nice.colors <- c("palegoldenrod","palegreen4","mistyrose3", "chartreuse2","palevioletred","pink",
                 "peru","sienna1","turquoise4", "slateblue4","yellow","violet",
                 "sienna4","yellowgreen","yellow4", "steelblue","red","tan2",
                 "steelblue4","gold4","tan", "forestgreen","gray","gray87",
                 "gold4","forestgreen","blueviolet", "brown","darkgoldenrod1","darkgray",
                 "darkgreen","darkmagenta", "darkolivegreen1","darksalmon",
                 "greenyellow","hotpink4","lemonchiffon4", "green4","blue","hotpink",
                 "lightgreen","ivory4","lightsalmon","mediumorchid4","lightsalmon4","lightseagreen",
                 "black","mediumseagreen","mediumvioletred","midnightblue","mistyrose3","lightsteelblue4",
                 "magenta","olivedrab","orange","orange4","mediumorchid","red4",
                 "orchid4")

# import the genus level abundance table
genus.table <- read.delim("microbiome_test_data_genus.tsv", sep="\t", header=T, row.names = 1)
# get the groups from the first row, and covert it to a simple list (called a vector)
groups <- as.vector(unlist(genus.table[1,]))

# remove the first row to get the table and convert to a matrix
# remove row 1
tmp.matrix <- genus.table[-1,]
# Convert data to numberic types
genus.counts <- apply(as.matrix.noquote(tmp.matrix),2,as.numeric)
# reset the correct column and row names
rownames(genus.counts) <- rownames(tmp.matrix)
colnames(genus.counts) <- colnames(tmp.matrix)

# calculate the abundance
genus.relab <- decostand(t(genus.counts), "total")

# establish clade order by average abundnace
clade.order <- order(colMeans(genus.relab), decreasing = TRUE)
# establish sample order by top clade value 
sample.order <- order(genus.relab[,clade.order[1]], decreasing = TRUE)
clade.ordered.names <- colnames(genus.relab)[clade.order]

# make the rownames with new order
sample.names.reorder <- factor(rownames(genus.relab),
                               levels=rownames(genus.relab)[sample.order])
tmp.table.frame <- data.frame(sample=sample.names.reorder, genus.relab, group=groups)

clade.melt.table <- melt(tmp.table.frame, id=c("sample", "group"))
colnames(clade.melt.table) <- c("sample", "group", "clade", "abundance")

# update the clade order in the table
clade.melt.table <- within(clade.melt.table,
                           clade <- factor(clade,
                                           levels=clade.ordered.names))

# move clades with less than 1% into an "Other category"
clade.table.hi <- subset(clade.melt.table, abundance >= 0.01)
clade.table.low <- subset(clade.melt.table, abundance < 0.01)
other.clade <- aggregate(abundance ~ sample + group, data=clade.table.low, FUN=sum)

# reformat the "other" clade and attach to the table
clade.table.other <- data.frame(sample=other.clade$sample,
                                group=other.clade$group,
                                clade="Other < 1%",
                                abundance=other.clade$abundance)

# combine the "other" categories back into to plot table
plot.table <- rbind(clade.table.hi, clade.table.other)
# rename group labels
levels(plot.table$group) <- c("G1","G2","G3","G4")
# abundance plot
ggplot(plot.table,
       aes(x=sample,
           y=abundance,
           ymax=1,
           fill=clade)) + 
  labs(list(title = "Genus level Relative Abundance",
            x = "Sample",
            y = "Relative Abundance",
            fill = "Lengend Title")) +
  facet_grid(.~group, scales = "free_x", space="free", labeller= as_labeller(group.names)) +
  geom_bar(stat="identity",width=0.8) + 
  scale_fill_manual(values=nice.colors) +
  theme(axis.text.x = element_text(angle = 90, hjust=1, vjust=0.40),
        plot.title = element_text(hjust = 0.5))









