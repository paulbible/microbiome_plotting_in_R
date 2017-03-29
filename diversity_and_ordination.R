# load needed libraries
library(vegan)
library(ggplot2)
library(reshape2)

########################
# Move R to the data directory
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


#################
# Alpha Diversity
#################
alpha.diversity <- diversity(species.relab, "shannon")
alpha.div.table <- data.frame(sample=rownames(species.relab),
                              alpha=alpha.diversity,
                              group=groups)

# Normal R Boxplot alpha diversity by group
boxplot(alpha.diversity ~ groups,
        xlab="Groups",
        ylab="Shannon Index",
        main="Alpha Diversity (Shannon Index) by Group")

# ggplot2 boxplot for alpha diversity
ggplot(alpha.div.table,
       aes(x=group, y=alpha, fill=group)) +
  # set labels for the plot
  labs(list(title = "Alpha Diversity (Shannon Index) by Group",
            x = "Groups",
            y = "Shannon Index",
            fill="Fill Legend Title")) +
  # use a bar plot
  geom_boxplot() +
  # Manually set colors here
  scale_fill_discrete(breaks = c("group_1", "group_2", "group_3", "group_4"),
                      labels = c("This is G1", "This is G2", "This is G3", "This is G4")) +
  # Manaully edit X axis labels
  scale_x_discrete(breaks = c("group_1", "group_2", "group_3", "group_4"),
                   labels = c("This is G1", "This is G2", "This is G3", "This is G4")) +
  # Rotate the x-axis sample labels, center the plot title
  theme(plot.title = element_text(hjust = 0.5),
        legend.position="none") # hide legend, legend is not necessary for this plot

# change sample order
# ggplot2 will order samples base on their order in the 'levels' of a 'factor'
# This command sets the order of sample to be the same as their normal order
alpha.div.table <- within(alpha.div.table,
                          sample <- factor(sample, 
                                    levels=sample))

# Plot alpha diversity for each sample
ggplot(alpha.div.table,
       aes(x=sample, y=alpha)) +
  # set labels for the plot
  labs(list(title = "Alpha Diversity (Shannon Index) by Group",
            x = "Samples",
            y = "Shannon Index",
            fill="Fill Legend Title")) +
  # use a bar plot
  geom_bar(aes(x=sample, fill=factor(group)), stat="identity") +
  # set the group name labels, replace these with your group labels.
  scale_fill_discrete(breaks = c("group_1", "group_2", "group_3", "group_4"),
                     labels = c("This is G1", "This is G2", "This is G3", "This is G4")) +
  # Rotate the x-axis sample labels, center the plot title
  theme(axis.text.x = element_text(angle = 90, hjust=1, vjust=0.40),
        plot.title = element_text(hjust = 0.5))

# Plot alpha diversity for each sample, with manually set colors.
ggplot(alpha.div.table,
       aes(x=sample, y=alpha)) +
  # set labels for the plot
  labs(list(title = "Alpha Diversity (Shannon Index) by Group",
            x = "Samples",
            y = "Shannon Index",
            fill="Fill Legend Title")) +
  # use a bar plot
  geom_bar(aes(x=sample, fill=factor(group)), stat="identity") +
  # Manually set colors here, replace breaks with your group labels.
  scale_fill_manual(values = c("red", "green", "blue", "orange"),
                    breaks = c("group_1", "group_2", "group_3", "group_4"),
                    labels = c("This is G1", "This is G2", "This is G3", "This is G4")) +
  # Rotate the x-axis sample labels, center the plot title
  theme(axis.text.x = element_text(angle = 90, hjust=1, vjust=0.40),
        plot.title = element_text(hjust = 0.5))


###############################
# Beta diversity and ordination
###############################
# calculate the Bray-Curtis distance (beta-diversity), using vegandist
bray.distance <- vegdist(species.relab, "bray")
# calculate the mutlidimentional scaling (mds) associated with Bray-Curtis distance
mds.bray <- cmdscale(bray.distance, 3, eig=TRUE)

# Plot or PCoA for Bray-Curtis beta-diversity
# create a data.frame for the PCoA data that ggplot2 will like
df = data.frame(PCoA1=mds.bray$points[,1], PCoA2=mds.bray$points[,2], group=groups)
ggplot(data=df, 
       aes(x=PCoA1,
           y=PCoA2,
           color=group)) + 
  geom_point() +
  labs(list(title = "PCoA of Relative abundance (Bray-Curtis Distance)",
            x = "PCo1",
            y = "PCo2")) +
  scale_color_discrete(name="Fill Group Title",
                       breaks = c("group_1", "group_2", "group_3", "group_4"),
                      labels = c("This is G1", "This is G2", "This is G3", "This is G4")) +
  theme(legend.position = 'bottom', plot.title = element_text(hjust = 0.5))

# PERMANOVA, Test for location effects (mean) of the groups (adonis function in vegan)
permanova.mod <- adonis(bray.distance ~ groups, perm=20000)
print(permanova.mod)

# Test for differences in dispersion effects (variance)
dispersion.mod <- betadisper(bray.distance, groups)
anova(dispersion.mod)

















