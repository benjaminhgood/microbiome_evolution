
# set the working directory
setwd('~/ben_nandita_hmp_data/species')

##################################
# alpha diversity per acccesion: #
##################################  

#read in the data 
data=read.table('relative_abundance.txt.bz2', header=T)

# iterate through each accession (columns) to compute alpha (shannon) diversity for each accession

H_vector=c() #vector to story every accession's shannon diveristy (H=SUM(-p*ln(p)), summing over all species

for (accession in c(2:length(data[1,]))){
    p=as.numeric(data[,accession])
    H=sum((-1)*p*log(p), na.rm=T)  # natural log
    H_vector=c(H_vector,H)
    }

## try again rescaling coverage
#coverage=read.table('coverage.txt.bz2', header=T)
#H_vector_2=c() #vector to story every accession's shannon diveristy (H=SUM(-p*ln(p)), summing over all species

#for (accession in c(2:length(data[1,]))){
#    cov=as.numeric(coverage[,accession])  
#    p=cov/sum(cov)
#    H=sum((-1)*p*log(p), na.rm=T)
#    H_vector_2=c(H_vector_2,H)        
#    }


####################################### 
# Alpha diversity per combined sample #
#######################################

# repeat for merged samples using the merged coverage file that ben made
sample_coverage=read.table('sample_coverage.txt.bz2', header=T)

# need to convert to relative abundances by scaling the sample coverages
H_vector_sample=c()
for (accession in c(2:length(sample_coverage[1,]))){ 
    cov=as.numeric(sample_coverage[,accession])
    p=cov/sum(cov)
    H=sum((-1)*p*log(p), na.rm=T)
    H_vector_sample=c(H_vector_sample,H)
    }


###########################
# compute gamma diversity #
###########################

# first, we want the coverage across all samples (total poplation)
total_coverage=rowSums(sample_coverage[,2:length(sample_coverage[1,])])

# now compute the gamma, which is alpha of the total pop. 
p=total_coverage/sum(total_coverage)
gamma=sum((-1)*p*log(p), na.rm=T)


##########################
# Compute beta diversity
# Beta=gamma/alpha
#########################

beta=gamma/H_vector_sample


###### 
#plot#
######
# Question: do the alpha diveristy computed per accession vs per sample comapare?
# Answer: YES
out_file="~/ben_nandita_hmp_analysis/alpha_diversity.pdf"
pdf(out_file,width=4.25,height=4.86, title = "Figure1",paper="special")   
boxplot(H_vector, H_vector_sample, col='blue', names=c('Per accession', 'Per sample'), main='Alpha values in HMP data', ylab='alpha (Shannon diversity)')
dev.off()

# Question: what does alpha vs beta diversity look like for the samples?
# Answer: beta diversity is a little higher than 1, meaning than within-patient diversity at the species' level is less than the global pop. species' diversity.

out_file="~/ben_nandita_hmp_analysis/alpha_beta_diversity.pdf"
pdf(out_file,width=4.25,height=4.86, title = "Figure1",paper="special")   
boxplot(H_vector_sample, beta, col='blue', names=c('alpha', 'beta'), main='Alpha and beta values in HMP data', ylab='value')
dev.off()
