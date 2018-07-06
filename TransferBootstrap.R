###  Authors Qingyi Lu and Joseph Rusinko 
## July 2018


### Code computes the Transfer Bootstrap Value defined in: Renewing Felsenstein's Phylogenetic Bootstrap in the Era of Big Data
##F. Lemoine, J.-B. Domelevo-Entfellner, E. Wilkinson, D. Correia, M. Davila Felipe, T. De Oliveira, O. Gascuel.
#Nature 556, 452-456 (2018)

### If you use this function please reference: Improving Statistical Binning Techniques for Species Tree Reconstruction: Rusinko, Vandenbussche, Lu (under review)

### Required Packages
library(caper)
library(ape)
library(phangorn)
library(ips)

#library(phytools)
#library(phylobase)





### Unless you are editing the code on your own, you should only need to call the Tboostrap function.

#######this is the final function, which can directly used to compute the transfer bootstrap value 
#a set of bootstrap value for an estgene tree;
#EstG is the estimated gene tree we are trying to compute the Transfer BS values for
#EstBoot is the  collection of bootsrtap trees
Tbootstrap <- function(EstG,EstBoot){
  boot <- c()
  for(i in 1:EstG$Nnode){
    boot[i] <- TBE(i,EstG,EstBoot)
  } 
  boot<-lapply(boot,round,2)
  EstG$node.label<-boot
  return(EstG)
}


###Internal Functions


######this is the first function used to compute the transfer dist function;
#y is the index of bootstrap tree; 
#z is the index of branch of the current estgene;
#w is the index of branch of the current bootstrap  tree;
#EstG is the estimated gene tree we are trying to compute the Transfer BS values for
#EstBoot is the  collection of bootsrtap trees
transfer.dist<-function(y,z,w,EstG,EstBoot){
  tempa<-EstG$Nnode
  taxa<-length(EstG$tip.label)
  des1 <- vector("list",tempa) #des1 is the descend descendants for the first estgene 
  for(i in 1:tempa){
    des1[[i]] <- clade.members(taxa+i, EstG, tip.labels=TRUE)
  }
  tempb<-EstBoot[[y]]$Nnode
  des2 <- vector("list",tempb) #des2 is the descend descendants for one of the bootstrap tree
  for(i in 1:tempb){
    des2[[i]] <- clade.members(taxa+i, EstBoot[[y]], tip.labels=TRUE)
  }
  
  node <- 0:(taxa-1)
  b <- des1[[z]]
  b.res <- node[is.na(pmatch(node,b))]
  
  b_ <- des2[[w]]
  b_.res <- node[is.na(pmatch(node,b_))]
  
  #intersection x with y
  inter1 <- intersect(b,b_)
  len1 <- length(b) - 2*length(inter1)+length(b_) 
  
  #intersection x with total\y
  inter2 <- intersect(b,b_.res)
  len2 <- length(b) - 2*length(inter2)+length(b_.res)
  
  #intersection total\x with y
  inter3 <- intersect(b.res,b_)
  len3 <- length(b.res) - 2*length(inter3)+length(b_)
  
  #intersection total\x with total\y
  inter4 <- intersect(b.res,b_.res)
  len4 <- length(b.res) - 2*length(inter4)+length(b_.res) 
  
  return(transfer.dist <- min(len1,len2,len3,len4))
}



######this is the second function used to compute the transfer index;
#y is the index of bootstrap tree in the  bootstrap tree list; 
#z is the index of branch of the current estgene tree
#EstG is the estimated gene tree we are trying to compute the Transfer BS values for
#EstBoot is the  collection of bootsrtap trees
transfer.index <- function(y,z,EstG,EstBoot){
  
  taxa<-length(EstG$tip.label)
  transfer.dist.list <- vector("list",EstG$Nnode)
  for(i in 1:EstG$Nnode){
    transfer.dist.list[[i]] <- transfer.dist(y,z,i,EstG,EstBoot)
  }
  
  des1 <- vector("list",EstG$Nnode) #des1 is the descend descendants for the first estgene 
  tempa<-EstG$Nnode
  for(i in 1:tempa){
    des1[[i]] <- clade.members(taxa+i, EstG, tip.labels=TRUE)
  }
  
  node <- 0:(taxa-1)
  b <- des1[[z]]
  b.res <- node[is.na(pmatch(node,b))]
  
  p <- min(length(b),length(b.res))
  p <- p-1
  transfer.dist.list[[taxa-1]] <- p
  
  return(min(as.numeric(transfer.dist.list)))
}



######this is the third function used to compute the TBE value
#y is the index of the  branch being analyized in the current estgene
#EstG is the estimated gene tree we are trying to compute the Transfer BS values for
#EstBoot is the  collection of bootsrtap trees
TBE <- function(y,EstG,EstBoot){
  taxa<-length(EstG$tip.label)
  numboot<-length(EstBoot)
  transfer.index.list <- vector("list",numboot)
  
  transfer.index.list<-lapply(1:numboot,function(i){transfer.index(i,y,EstG,EstBoot)}) 
  
  avg.transfer.index <- do.call(sum,transfer.index.list)/numboot
  
  des1 <- vector("list",EstG$Nnode) #des1 is the descend descendants for the first estgene 
  tempa<-EstG$Nnode
  for(i in 1:tempa){
    des1[[i]] <- clade.members(taxa+i, EstG, tip.labels=TRUE)
  }
  taxa<-length(EstG$tip.label)
  node <- 0:(taxa-1)
  b <- des1[[y]]
  b.res <- node[is.na(pmatch(node,b))]
  
  p <- min(length(b),length(b.res))
  TBE <- 1-(avg.transfer.index/(p-1))
  
  return(TBE) 
}




