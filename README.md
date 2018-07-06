# TransferBoot
Transfer Bootstrap in R

### References
Code computes the Transfer Bootstrap Value defined in: Renewing Felsenstein's Phylogenetic Bootstrap in the Era of Big Data F. Lemoine, J.-B. Domelevo-Entfellner, E. Wilkinson, D. Correia, M. Davila Felipe, T. De Oliveira, O. Gascuel. #Nature 556, 452-456 (2018)

If you use this function please reference: Improving Statistical Binning Techniques for Species Tree Reconstruction: Rusinko, Vandenbussche, Lu (under review)


### System Requirements

This is R code which requires the following libraries packages: caper, ape,phangorn, ips.  

### Computing the Transfer Bootstrap

The primary function Tbootsrtap has two inputs.  The first a single phylogenetic tree and the second a list of estimated bootstrap trees. The output of this function is the input tree with Transfer Bootstrap values as the node labels.

### Useage

Tbootstrap(InputTree,BootstrapTrees)

### Example (Assumes you have downloaded and run the File TransferBootstrap.R and downloaded the file sampleDNA)

dna<-read.phy("sampleDNA") 

EstTree = nj(dist.gene(dna)) 

BootstrapTrees <- boot.phylo(EstTree, dna, function(xx) nj(dist.gene(xx)), B = 50,trees=TRUE) 

TransferBootstrap<-Tbootstrap(EstTree,BootstrapTrees[[2]])

plot(TransferBootstrap, main = "My Tree with Transfer Bootstrap Values")

drawSupportOnEdges(TransferBootstrap$node.label)

### Recommendations and Limitations
This code is designed for easy integration with R.  This code is slow for large datasets.  If your data is taking too long to run in this script we recommend you write your bootstrap trees to a file and use another implementation of this algorithm.  See https://github.com/evolbioinfo/booster for alternatives.  








