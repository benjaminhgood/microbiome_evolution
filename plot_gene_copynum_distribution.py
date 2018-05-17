import numpy
import sys
import parse_midas_data
import gzip
import stats_utils
import pylab

species_name = sys.argv[1]
desired_gene_names = sys.argv[2:]

pylab.xlim([-0.1,2])
pylab.title('%s' % (species_name))
pylab.xlabel('Copynum, c')
pylab.ylabel('# samples >= c')


sys.stderr.write("Loading pangenome data for %s...\n" % species_name)
gene_samples, gene_names, gene_presence_matrix, gene_depth_matrix, marker_coverages, gene_reads_matrix = parse_midas_data.parse_pangenome_data(species_name)
sys.stderr.write("Done!\n")  

# calculating copynumbers
gene_names = numpy.array(gene_names)
gene_copynum_matrix = gene_depth_matrix*1.0/(marker_coverages+(marker_coverages==0))

desired_gene_idxs = numpy.array([gene_name in desired_gene_names for gene_name in gene_names])

min_copynum_distribution = gene_copynum_matrix[desired_gene_idxs,:].min(axis=0)

for gene_name in desired_gene_names:

    gene_idx = numpy.nonzero(gene_names==gene_name)[0][0]

    gene_copynum_distribution = gene_copynum_matrix[gene_idx,:]

    print gene_copynum_matrix.shape, gene_copynum_distribution.shape



    #print gene_copynum_distribution

    xvalues, ns = stats_utils.calculate_unnormalized_survival_from_vector(gene_copynum_distribution, min_x=0, max_x=gene_copynum_distribution.max())

    pylab.step(xvalues,ns,label=gene_name)
    #pylab.semilogy([4],[4])
    
xvalues, ns = stats_utils.calculate_unnormalized_survival_from_vector(min_copynum_distribution, min_x=0, max_x=min_copynum_distribution.max())

pylab.step(xvalues,ns,label='Both')

    
pylab.legend(loc='upper right',frameon=False)

pylab.savefig('../morteza_collaboration/ben_figures/Alistipes_onderdonkii_gene_gain_HMP_prevalence.pdf',bbox_inches='tight')

#pylab.show()
    