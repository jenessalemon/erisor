#this code is from VCF-kit:
#http://vcf-kit.readthedocs.io/en/latest/phylo/
#generate neighbor joining trees from newick format.
#get newick format from running ...vk phylo tree nj my_vcf > my_vcf.newick one directory above where Python is stored. (Pipe to emtpy file). 
library(tidyverse)
library(ape)
library(ggmap)
library(phyloseq)

tree <- ape::read.tree(paste0("p1.newick")) 

# Optionally set an outgroup.
# tree <- root(tree,outgroup = "outgroup", resolve.root = T)

treeSegs <- phyloseq::tree_layout(
  phyloseq::phy_tree(tree),
  ladderize = T
)

treeSegs$edgeDT <- treeSegs$edgeDT  %>% 
  dplyr::mutate(edge.length = 
                  ifelse(edge.length < 0, 0, edge.length)
                , xright = xleft + edge.length
  )
edgeMap = aes(x = xleft, xend = xright, y = y, yend = y)
vertMap = aes(x = x, xend = x, y = vmin, yend = vmax)
labelMap <- aes(x = xright+0.001, y = y, label = OTU)

ggplot(data = treeSegs$edgeDT) + geom_segment(edgeMap) + 
  geom_segment(vertMap, data = treeSegs$vertDT) +
  geom_text(labelMap, data = dplyr::filter(treeSegs$edgeDT, !is.na(OTU)), na.rm = TRUE, hjust = -0.05) +
  ggmap::theme_nothing() + 
  scale_x_continuous(limits = c(
    min(treeSegs$edgeDT$xleft)-0.3,
    max(treeSegs$edgeDT$xright)+0.3
  ),
  expand = c(0,0))