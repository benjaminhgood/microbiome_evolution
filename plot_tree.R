args = commandArgs(TRUE)
species=args[1]

# set the working directory
setwd('~/ben_nandita_hmp_data/species')

library(ape)
distance=read.table(paste('~/ben_nandita_hmp_analysis/fst.dist_', species, sep=''))
distance=as.dist(distance)
tree=njs(distance)

out_file=paste("~/ben_nandita_hmp_analysis/",species,"_tree.pdf",sep='')
pdf(out_file,width=4.25,height=4.86, title = "Figure1",paper="special")   
plot(tree,"u",show.tip.label = FALSE, main=species)
add.scale.bar(cex = 0.7, font = 2)
dev.off()
