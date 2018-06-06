from scipy.cluster.hierarchy import dendrogram, linkage
from scipy.cluster.hierarchy import cophenet
from scipy.cluster.hierarchy import fcluster
from numpy.random import shuffle, normal

import diversity_utils

import parse_midas_data
import numpy

# min_d = pick only a single sample per cluster with distance below this value
# max_d = cut tree at this distance
def cluster_samples(distance_matrix, min_d=0, max_ds=[1e09]):
 
    # calculate compressed distance matrix suitable for agglomerative clustering
    Y = []
    for i in xrange(0,distance_matrix.shape[0]):
        for j in xrange(i+1,distance_matrix.shape[1]):
            Y.append(distance_matrix[i,j]) 
    Y = numpy.array(Y) 
     
    Z = linkage(Y, method='average')        
    
    # First coarse-grain things less than min_d apart:
    #NRG: what does it mean to coarse-grain?
    subcluster_assignments = fcluster(Z, min_d, criterion='distance')
    
    coarse_grained_idxs = []
    subcluster_idx_map = {}
    for i in xrange(0,len(subcluster_assignments)):
        if subcluster_assignments[i] not in subcluster_idx_map:
            subcluster_idx_map[subcluster_assignments[i]] = i
            coarse_grained_idxs.append(True)
        else:
            coarse_grained_idxs.append(False)
            
    coarse_grained_idxs = numpy.array(coarse_grained_idxs)
    
        
    sorted_final_clusterss = []
    for max_d in max_ds:
        
        cluster_assignments = fcluster(Z, max_d, criterion='distance')
        
        cluster_idx_map = {}
    
        for i in xrange(0,len(cluster_assignments)):
        
            if not coarse_grained_idxs[i]:
                continue
                
            if cluster_assignments[i] not in cluster_idx_map:
                cluster_idx_map[cluster_assignments[i]] = []
            cluster_idx_map[cluster_assignments[i]].append(i)
                
        cluster_labels = set(cluster_idx_map.keys())
        cluster_idxss = [set(cluster_idx_map[cluster_label]) for cluster_label in cluster_labels]
        cluster_sizes = [len(cluster_idxs) for cluster_idxs in cluster_idxss]
     
        # only return ones with more than one individual
        final_clusters = []
        final_cluster_sizes = []
      
        for cluster_idx_set in cluster_idxss:
         
            if len(cluster_idx_set)>1:
         
                cluster_idxs = numpy.array([(i in cluster_idx_set) for i in xrange(0,len(cluster_assignments))])
            
                final_clusters.append(cluster_idxs)
                final_cluster_sizes.append((cluster_idxs*1.0).sum())
        
        if len(final_cluster_sizes) > 0:
             
            final_cluster_idxs = [i for i in xrange(0,len(final_cluster_sizes))]
         
            final_cluster_sizes, final_cluster_idxs = zip(*sorted(zip(final_cluster_sizes, final_cluster_idxs),reverse=True))
    
        
            sorted_final_clusters = [final_clusters[idx] for idx in final_cluster_idxs]
            sorted_final_clusterss.append(sorted_final_clusters)
        else:
            sorted_final_clusterss.append([])
        
    return coarse_grained_idxs, sorted_final_clusterss
         

# Perform hierarchical clustering respecting the clade boundaries in clade_idxss
def cluster_samples_within_clades(distance_matrix, clade_idxss=[], d=1e09):

    
    if len(clade_idxss)==0:
        clade_idxss = [numpy.array([True for i in xrange(0,distance_matrix.shape[0])])]
        
    subcluster_sets = []

    for clade_idxs in clade_idxss:
        
        numeric_clade_idxs = numpy.nonzero(clade_idxs)[0]
        
        # get subset distance matrix
        sub_distance_matrix = distance_matrix[numpy.ix_(numeric_clade_idxs, numeric_clade_idxs)]
 
        # calculate compressed distance matrix suitable for agglomerative clustering
        Y = []
        for i in xrange(0,sub_distance_matrix.shape[0]):
            for j in xrange(i+1,sub_distance_matrix.shape[1]):
                Y.append(sub_distance_matrix[i,j]) 
        Y = numpy.array(Y) 
     
        Z = linkage(Y, method='average')        
    
        # First coarse-grain things less than min_d apart:
        subcluster_assignments = fcluster(Z, d, criterion='distance')
    
        new_subcluster_sets = {}
    
        for i in xrange(0,len(subcluster_assignments)):
            
            if subcluster_assignments[i] not in new_subcluster_sets:
                new_subcluster_sets[subcluster_assignments[i]] = set([])
                
            new_subcluster_sets[subcluster_assignments[i]].add( numeric_clade_idxs[i] )
            
        for subcluster_set in new_subcluster_sets.values():
            if len(subcluster_set) > 1:
                subcluster_sets.append(subcluster_set)
            
    return subcluster_sets
 


 
def calculate_phylogenetic_consistency(allele_counts_map, passed_sites_map, proposed_clusters, allowed_variant_types=set(['1D','2D','3D','4D']), allowed_genes=set([])):
 
    clusters = []
    anticlusters = []
    for cluster_idxs in proposed_clusters:
        
        #print cluster_idxs.sum(), numpy.logical_not(cluster_idxs).sum()
        
        if cluster_idxs.sum() > 1.5: # Need at least two guys in a cluster to look for polymorphisms
            
            anticluster_idxs = numpy.logical_not(cluster_idxs)
            
            if anticluster_idxs.sum() > 1.5: # Likewise for the anticluster
                
                clusters.append(cluster_idxs)
                anticlusters.append( anticluster_idxs )
 
    total_genes = set(passed_sites_map.keys())
 
    if len(allowed_genes)==0:
        allowed_genes = set(passed_sites_map.keys())
     
    allowed_genes = (allowed_genes & total_genes)     
      
    singleton_freqs = [] # actual freq value is meaningless..                 
    polymorphic_freqs = [] # non-singleton freqs -- only ones that can be inconsistent!
    inconsistent_freqs = []
    null_inconsistent_freqs = []
    
    singleton_variant_types = {variant_type: 0 for variant_type in allowed_variant_types}
    polymorphic_variant_types = {variant_type: 0 for variant_type in allowed_variant_types}
    inconsistent_variant_types = {variant_type: 0 for variant_type in allowed_variant_types}
    null_inconsistent_variant_types = {variant_type: 0 for variant_type in allowed_variant_types}
    
    if len(clusters)>0: # Can only do stuff if there are clusters!
        
     
        for gene_name in allowed_genes:
         
            for variant_type in passed_sites_map[gene_name].keys():
              
                if variant_type not in allowed_variant_types:
                    continue
         
                allele_counts = allele_counts_map[gene_name][variant_type]['alleles']                        
                if len(allele_counts)==0:
                    continue
                
                # good to go, let's get calculating
                
                # take consensus approximation
                genotype_matrix, passed_sites_matrix = diversity_utils.calculate_consensus_genotypes(allele_counts)
             
                population_prevalence = (genotype_matrix*passed_sites_matrix).sum(axis=1)
                population_max_prevalence = (passed_sites_matrix).sum(axis=1)
            
                population_minor_prevalence = numpy.fmin(population_prevalence, population_max_prevalence - population_prevalence)
            
                population_freqs = population_prevalence*1.0/(population_max_prevalence+10*(population_max_prevalence<0.5))
                population_freqs = numpy.fmin(population_freqs, 1-population_freqs)
     
                is_polymorphic = numpy.zeros(genotype_matrix.shape[0])
                is_inconsistent = numpy.zeros(genotype_matrix.shape[0])
     
                for cluster_idxs,anticluster_idxs in zip(clusters,anticlusters):
             
                
                    cluster_prevalence = (genotype_matrix[:,cluster_idxs]*passed_sites_matrix[:,cluster_idxs]).sum(axis=1)
                    cluster_min_prevalence = 1-1e-09
                    cluster_max_prevalence = (passed_sites_matrix[:,cluster_idxs]).sum(axis=1)-1+1e-09
                
                    cluster_freqs = cluster_prevalence*1.0/(cluster_max_prevalence+10*(cluster_max_prevalence<0.5))
                    cluster_freqs = numpy.fmin(cluster_freqs, 1-cluster_freqs)
             
                    anticluster_prevalence = (genotype_matrix[:,anticluster_idxs]*passed_sites_matrix[:,anticluster_idxs]).sum(axis=1)
                    anticluster_min_prevalence = 1-1e-09
                    anticluster_max_prevalence = (passed_sites_matrix[:,anticluster_idxs]).sum(axis=1) -1+1e-09
             
                    # Those that are polymorphic in the clade!
                    polymorphic_sites = (cluster_prevalence>=cluster_min_prevalence)*(cluster_prevalence<=cluster_max_prevalence)
                 
                    # Those that are also polymorphic in the remaining population!
                    inconsistent_sites = polymorphic_sites*(anticluster_prevalence>=anticluster_min_prevalence)*(anticluster_prevalence<=anticluster_max_prevalence)
             
                    is_polymorphic = numpy.logical_or(is_polymorphic, polymorphic_sites)
                    is_inconsistent = numpy.logical_or(is_inconsistent, inconsistent_sites)
            
                if is_polymorphic.sum() > 0:
            
                    is_singleton = (numpy.fabs(population_minor_prevalence-1)<1e-08)*is_polymorphic
                
                    is_polymorphic = (population_minor_prevalence>1.5)*is_polymorphic
                
                    singleton_freqs.extend( population_freqs[is_singleton] )
                    singleton_variant_types[variant_type] += is_singleton.sum()
                
                    polymorphic_freqs.extend( population_freqs[is_polymorphic] )
                    polymorphic_variant_types[variant_type] += is_polymorphic.sum()
                
                    if is_inconsistent.sum() > 0:
                        #inconsistent_freqs.extend( cluster_freqs[is_inconsistent] )
                        inconsistent_freqs.extend( population_freqs[is_inconsistent] )
                        inconsistent_variant_types[variant_type] += is_inconsistent.sum()
                
                    # now try to compute a null expectation for a completely unlinked genome
                    polymorphic_idxs = numpy.arange(0,genotype_matrix.shape[0])[is_polymorphic]
                    # Loop over sites that were polymorphic, generate a "null" draw for them
                    for site_idx in polymorphic_idxs:
                    
                        genotypes = genotype_matrix[site_idx,:]
                        passed_sites = passed_sites_matrix[site_idx,:]
                        population_freq = population_freqs[site_idx]
                    
                        permuted_idxs = numpy.arange(0,len(genotypes))
                    
                        is_polymorphic = False
                        is_inconsistent = False
                        # loop until we find a polymorphic site
                        while not is_polymorphic:
                    
                            # permute indexes 
                            shuffle(permuted_idxs)
                        
                            permuted_genotypes = genotypes[permuted_idxs]
                            permuted_passed_sites = passed_sites[permuted_idxs]
                        
                            # loop through clusters
                            is_inconsistent = False
                            for cluster_idxs,anticluster_idxs in zip(clusters,anticlusters):
             
                            
                                cluster_prevalence = (permuted_genotypes[cluster_idxs]*permuted_passed_sites[cluster_idxs]).sum()
                                cluster_min_prevalence = 0.5
                                cluster_max_prevalence = (permuted_passed_sites[cluster_idxs]).sum()-0.5
                
                    
                                anticluster_prevalence = (permuted_genotypes[anticluster_idxs]*permuted_passed_sites[anticluster_idxs]).sum()
                                anticluster_min_prevalence = 0.5
                                anticluster_max_prevalence = (permuted_passed_sites[anticluster_idxs]).sum() - 0.5
             
                                polymorphic_in_cluster = ((cluster_prevalence>cluster_min_prevalence)*(cluster_prevalence<cluster_max_prevalence))
                                inconsistent_in_cluster = (polymorphic_in_cluster*(anticluster_prevalence>anticluster_min_prevalence)*(anticluster_prevalence<anticluster_max_prevalence))
                            
                                is_polymorphic = is_polymorphic or polymorphic_in_cluster
                                is_inconsistent = is_inconsistent or inconsistent_in_cluster
                    
                        if is_inconsistent:
                            null_inconsistent_freqs.append(population_freq)
                            null_inconsistent_variant_types[variant_type] += 1
                        
    singleton_freqs = numpy.array(singleton_freqs)            
    polymorphic_freqs = numpy.array(polymorphic_freqs)
    inconsistent_freqs = numpy.array(inconsistent_freqs)
    null_inconsistent_freqs = numpy.array(null_inconsistent_freqs)
         
    return singleton_freqs, polymorphic_freqs, inconsistent_freqs, null_inconsistent_freqs, singleton_variant_types, polymorphic_variant_types, inconsistent_variant_types, null_inconsistent_variant_types

def calculate_clade_allele_freqs(allele_counts_map, passed_sites_map, clusters, allowed_variant_types=set(['1D','2D','3D','4D']), allowed_genes=set([])):
 
    anticlusters = []
    for cluster_idxs in clusters:
        anticlusters.append( numpy.logical_not(cluster_idxs) )
 
    total_genes = set(passed_sites_map.keys())
 
    if len(allowed_genes)==0:
        allowed_genes = set(passed_sites_map.keys())
     
    allowed_genes = (allowed_genes & total_genes)     
      
    ktotals = []
    ntotals = []
    
    clade_ks = [[] for cluster in clusters]
    clade_ns = [[] for cluster in clusters]  
     
    for gene_name in allowed_genes:
        for variant_type in passed_sites_map[gene_name].keys():
            if variant_type not in allowed_variant_types:
                continue
         
            allele_counts = allele_counts_map[gene_name][variant_type]['alleles']                        
            if len(allele_counts)==0:
                continue
                
            # good to go, let's get calculating
                
            # take consensus approximation
            genotype_matrix, passed_sites_matrix = diversity_utils.calculate_consensus_genotypes(allele_counts)
            
            ktots = (genotype_matrix*passed_sites_matrix).sum(axis=1)
            ntots = (passed_sites_matrix).sum(axis=1)
            kminortots =  numpy.fmin(ktot, ntot - ktot)
            
            polymorphic_sites = (kminortots>0.5)
            
            genotype_matrix = genotype_matrix[polymorphic_sites,:]
            passed_sites_matrix = passed_sites_matrix[polymorhpic_sites,:]
            
            ktotals.extend(ktots)
            ntotals.extend(ntots)
            
            for cluster_idx in xrange(0,len(clusters)):
                
                cluster_idxs = clusters[cluster_idx]
                anticluster_idxs = anticlusters[cluster_idx]
                
                
                k0s = (genotype_matrix[:,cluster_idxs]*passed_sites_matrix[:,cluster_idxs]).sum(axis=1)
                n0s = cluster_max_prevalence = (passed_sites_matrix[:,cluster_idxs]).sum(axis=1)
                
                clade_ks.extend(k0s)
                clade_ns.extend(n0s) 
            
    return ktotals, ntotals, clade_ks, clade_ns

def calculate_phylogenetic_inconsistency_from_sfs(ktotals, ntotals, clade_ks, clade_ns):

    # first coarse-grain allele frequencies
    n = long(numpy.median(ntotals))
    frequency_bins = numpy.arange(1,n+1)*1.0/n
    frequency_bins = frequency_bins - (frequency_bins[1]-frequency_bins[0])/2.0 
    frequency_bins[0] = -1
    frequency_bins[-1] = 2
    fs = numpy.arange(1,n)*1.0/n
    pfs = numpy.histogram(ktotals*1.0/ntotals, bins=frequency_bins)[0]
    pfs = (pfs+pfs[::-1])/2
    
def load_manual_clade_divergence_threshold(species_name):
    file = open(parse_midas_data.scripts_directory+"manual_clade_thresholds.txt","r")
    file.readline()
    divergence_threshold = 1e-02
    for line in file:
        items=line.split("\t")
        print items[0].strip(), species_name
        if items[0].strip()==species_name:
            divergence_threshold = float(items[1])
            print "Setting divergence threshold!", divergence_threshold
            break
    file.close()
    return divergence_threshold
    
def load_manual_clades(species_name):
    file = open(parse_midas_data.scripts_directory+"manual_clade_definitions.txt","r")
    
    line = file.readline().strip()
    
    # Just put some default values there in case of issue
    current_species = ""
    clades = {"":set()}
    
    while line!="":
    
        items = line.split()
        
        if items[0].isdigit():
            # is a sample entry
            sample_name = items[1].strip()
            clades[current_species][-1].add(sample_name)
        elif items[0].startswith('-'):
            # delimits a clade
            clades[current_species].append(set())
        else:
            # new species
            current_species = items[0].strip()
            clades[current_species] = [set()]
    
        line = file.readline().strip()
    
    file.close()
    
    if species_name in clades:
        return clades[species_name]
    else:
        return []

def calculate_clade_idxs_from_clade_sets(samples, clade_sets):

    clade_idxss = []
    for clade_set in clade_sets:
        
        clade_idxs = numpy.array([sample in clade_set for sample in samples])
        clade_idxss.append(clade_idxs)
        
    return clade_idxss

def permute_idxs_within_clades(cluster_idxss):

    permuted_idxs = numpy.arange(0,cluster_idxss[0].shape[0])
    
    for cluster_idxs in cluster_idxss:
        idxs = list(permuted_idxs[cluster_idxs])
        shuffle(idxs)
        permuted_idxs[cluster_idxs] = numpy.array(idxs)
        
    return permuted_idxs  
    
    
###########
#
# Here is another set of algorithms that does it without tree structure
#
###########

###
#
# Finds SNPs. 
# Calculates min distance between the two alleles (i.e., approximate "date" if asexual)
# Then calculates average and max distance within each allele. 
#
###
def calculate_snp_distances(allele_counts_map, passed_sites_map, distance_matrix, allowed_variant_types=set(['1D','2D','3D','4D']), allowed_genes=set([])):
 
    total_genes = set(passed_sites_map.keys())
 
    if len(allowed_genes)==0:
        allowed_genes = set(passed_sites_map.keys())
     
    allowed_genes = (allowed_genes & total_genes)     
 
    snp_data = []
     
    for gene_name in allowed_genes:
         
        for variant_type in passed_sites_map[gene_name].keys():
              
            if variant_type not in allowed_variant_types:
                continue
         
            allele_counts = allele_counts_map[gene_name][variant_type]['alleles']                        
            if len(allele_counts)==0:
                continue
            locations = allele_counts_map[gene_name][variant_type]['locations']
            
            # good to go, let's get calculating
                
            # take consensus approximation
            genotype_matrix, passed_sites_matrix = diversity_utils.calculate_consensus_genotypes(allele_counts)
             
            derived_samples = numpy.logical_and(genotype_matrix, passed_sites_matrix)
            ancestral_samples = numpy.logical_and(numpy.logical_not(genotype_matrix), passed_sites_matrix)
                   
            derived_allele_counts = derived_samples.sum(axis=1)
            ancestral_allele_counts = ancestral_samples.sum(axis=1)
            
            # Get SNVs to look at (no singletons)
            polymorphic_idxs = numpy.nonzero( (derived_allele_counts>1.5)*(ancestral_allele_counts>1.5) )[0]
            
            for snp_idx in polymorphic_idxs:
                 
                # Get cross snp distance
                between_distances = distance_matrix[derived_samples[snp_idx,:],:][:,ancestral_samples[snp_idx,:]]
                
                min_between_d = between_distances.min()    

                # Get within snp distance
                within_derived_distances =  distance_matrix[derived_samples[snp_idx,:]][:,derived_samples[snp_idx,:]]
                within_derived_distances = within_derived_distances[numpy.triu_indices(within_derived_distances.shape[0], k = 1)]
                
                max_derived_d = within_derived_distances.max()
                avg_derived_d = within_derived_distances.mean()
                
        
                within_ancestral_distances =  distance_matrix[ancestral_samples[snp_idx,:],:][:,ancestral_samples[snp_idx,:]]
                within_ancestral_distances = within_ancestral_distances[numpy.triu_indices(within_ancestral_distances.shape[0], k = 1)]
                
                max_ancestral_d = within_ancestral_distances.max()
                avg_ancestral_d = within_ancestral_distances.mean()
                
                snp_data.append((locations[snp_idx],variant_type,derived_allele_counts[snp_idx],ancestral_allele_counts[snp_idx], min_between_d, avg_derived_d, max_derived_d, avg_ancestral_d, max_ancestral_d))
                
    return snp_data
    
