import parse_midas_data
import pylab
import sys
import numpy
import diversity_utils

species_names = ["Bacteroides_uniformis_57318"]

for species_name in species_names:

    # Load subject and sample metadata
    sys.stderr.write("Loading HMP metadata...\n")
    subject_sample_map = parse_midas_data.parse_subject_sample_map()
    sys.stderr.write("Done!\n")
    
    
    # Load SNP information for species_name
    sys.stderr.write("Loading %s...\n" % species_name)
    samples, allele_counts_map, passed_sites_map = parse_midas_data.parse_snps(species_name, combination_type="sample", debug=False)
    sys.stderr.write("Done!\n")

    locations = []
    pi_syns = []
    pi_nons = []
    for gene_name in passed_sites_map.keys():
        
        passed_sites = passed_sites_map[gene_name]['4D']['sites']
        passed_pairs = (passed_sites>0.5)
        passed_pairs[numpy.diag_indices_from(passed_pairs)] = False
        pi_matrix, avg_pi_matrix = diversity_utils.calculate_pi_matrix(allele_counts_map, passed_sites_map, variant_type='4D',allowed_genes=set([gene_name]))
        total_pairs = passed_pairs.sum()
        
        if total_pairs<0.5:
            continue
        
        avg_pi = (pi_matrix*passed_pairs).sum()/(total_pairs+(total_pairs<0.5))
        
        pi_syns.append(avg_pi)
        locations.append(passed_sites_map[gene_name]['4D']['location'])
        
        passed_sites = passed_sites_map[gene_name]['1D']['sites']
        passed_pairs = (passed_sites>0.5)
        passed_pairs[numpy.diag_indices_from(passed_pairs)] = False
        pi_matrix, avg_pi_matrix = diversity_utils.calculate_pi_matrix(allele_counts_map, passed_sites_map, variant_type='1D',allowed_genes=set([gene_name]))
        total_pairs = passed_pairs.sum()
        avg_pi = (pi_matrix*passed_pairs).sum()/(total_pairs+(total_pairs<0.5))
        
        pi_nons.append(avg_pi)
        
    locations, pi_syns, pi_nons = zip(*sorted(zip(locations, pi_syns, pi_nons)))
    
    xvalues = numpy.arange(0,len(pi_syns))
    pi_syns = numpy.array(pi_syns)
    pi_nons = numpy.array(pi_nons)
    
    # Done calculating... now plot figure!
    pylab.figure(figsize=(24,3))
    pylab.xlabel('Genes (ordered)')
    pylab.ylabel('Pairwise diversity between hosts')
    #pylab.ylim([0,1.1])
    pylab.title(species_name)
  
    pylab.plot(xvalues, pi_syns, 'b.-', label='Synonymous (4D)')
    pylab.plot(xvalues, pi_nons, 'r.-', label='Nonsynonymous (1D)')
    pylab.semilogy(xvalues, numpy.ones_like(xvalues)*1e-02,'k:')
    pylab.legend(loc='upper right',frameon=False)
    pylab.savefig('figures/%s_gene_pi.pdf' % species_name,bbox_inches='tight')
    #pylab.savefig('figures/%s_gene_pi.png' % species_name,bbox_inches='tight')
    
    pylab.figure()
    pylab.loglog(pi_syns, pi_nons,'b.',alpha=0.5)
    pylab.xlabel('Gene $\\pi_s$')
    pylab.ylabel('Gene $\\pi_n$')
    pylab.savefig('figures/%s_gene_pN_vs_pS.pdf' % species_name,bbox_inches='tight')
    pylab.savefig('figures/%s_gene_pN_vs_pS.png' % species_name,bbox_inches='tight')
    
    