# Microbiome Plotting in R
Some notes on plotting microbiome data. This repository is intended to help other in making nice looking microbiome figures and perform some simple analysis.

## Brief note on how data was simulated

These data were simulated from a mouse gut microbiome study. The distribution of abundance values for each species was modeled and then sampled randomly with noise added. The data were then clustered using k-means with k = 4. The data points were then adjusted to be 80% closer to their k-means centers. This adjustment made the ordination look better. This dataset should be "like" a real dataset but any results/conclusions are not meaningful.

The genus data were constructed in a similar manner and are not related to the species data.

## Genus abundance bar graph

![bar-ab](/plots/genus_abundance_bar_graph.png)

## Dot abundance graph

![dot-ab](/plots/species_dot_abundance_graph.png)

## Principle Coordinate Analysis

![pcoa](/plots/bray_curtis_PCoA_ordination.png)

## Alpha diversity by sample

![a-div-bargraph](/plots/alpha_diversity_by_sample_bar_graph.png)

## Alpha diversity

![a-div-boxplot](/plots/alpha_diversity_species_boxplot.png)

## Random Forest Variable Importance
![random-forest](/plots/random_forest_varaible_importance.png)

