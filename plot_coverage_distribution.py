import pylab
import numpy
import parse_midas_data
import stats_utils

combination_type="sample"
#species = "Bacteroides_uniformis_57318"
#species = "Roseburia_intestinalis_56239"
#species = "Eubacterium_eligens_61678"
species = "Prevotella_copri_61740"
# Load marker gene coverages
species_coverage_matrix, sample_list, species_list = parse_midas_data.parse_marker_gene_coverages(species, combination_type=combination_type)
marker_coverages = species_coverage_matrix[species_list.index(species),:]

# Load genomic coverage distributions
sample_coverage_histograms, samples = parse_midas_data.parse_coverage_distribution(species, combination_type=combination_type)
median_coverages = numpy.array([stats_utils.calculate_median_from_histogram(sample_coverage_histogram) for sample_coverage_histogram in sample_coverage_histograms])
iqrs = numpy.array([stats_utils.calculate_IQR_from_histogram(sample_coverage_histogram) for sample_coverage_histogram in sample_coverage_histograms])

# Gene gene specific coverages
gene_coverages, samples = parse_midas_data.parse_gene_coverages(species, combination_type=combination_type)


peak_coverages = []
for i in xrange(0,len(sample_coverage_histograms)):
    
    sample_coverage_histogram = sample_coverage_histograms[i]
    x0 = median_coverages[i]
    
    if x0==0:
        peak_coverages.append(0)
        continue
        
    #x0 = marker_coverages[i]
    xs = sorted(sample_coverage_histogram.keys())
    pmfs = numpy.array([sample_coverage_histogram[x] for x in xs])
    pmfs /= pmfs.sum()
    
    xs = numpy.array(xs)
    
    xs = xs/x0
    
    dx = 0.1
    coarse_xs = numpy.linspace(0,4,400)
    coarse_pmfs = []
    
    #print "--"
    for x in coarse_xs:
        numerator = (pmfs*numpy.exp(-numpy.square(xs-x)/2/dx/dx)).sum()
        denominator = (numpy.exp(-numpy.square(xs-x)/2/dx/dx).sum())
        #print numerator, denominator
        coarse_pmfs.append( numerator / denominator )
    
    mode_idxs = []
    for j in xrange(1,len(coarse_xs)-1):
        if (coarse_pmfs[j]>coarse_pmfs[j-1]) and (coarse_pmfs[j]>coarse_pmfs[j+1]):
            mode_idxs.append(j)
            
    modes = numpy.array([coarse_xs[j] for j in mode_idxs])
    peak_pmfs = numpy.array([coarse_pmfs[j] for j in mode_idxs])
    
    x0_peak = modes[peak_pmfs.argmax()]
    x0_1 = modes[numpy.square(modes-1).argmin()]
    x0 = x0_peak
    x0 *= median_coverages[i]
    
    peak_coverages.append(x0)
peak_coverages = numpy.array(peak_coverages)

#pylab.figure(1)
#pylab.plot(marker_coverages, peak_coverages/marker_coverages,'bo')
#pylab.plot(marker_coverages, median_coverages/marker_coverages,'ro')

#pylab.semilogx([1,1000],[1,1],'k:')
#pylab.ylim([0,2])
#pylab.xlabel('Marker gene coverage')
#pylab.ylabel('Median genomic coverage')

pylab.figure(1)
sorted_median_coverages = numpy.array([x for x in sorted(median_coverages)])
#print sorted_median_coverages
#print len(sorted_median_coverages)
pylab.step(sorted_median_coverages, len(sorted_median_coverages)-numpy.arange(0,len(sorted_median_coverages)),where='post')
pylab.semilogx([1,1000],len(sorted_median_coverages)*numpy.ones(2),'k:')
pylab.xlim([1,1000])
pylab.ylim([0,250])
pylab.xlabel('Median coverage across genome')
pylab.ylabel('Samples >= coverage')
pylab.savefig('figures/%s_depth_distribution.pdf' % species, bbox_inches='tight')
pylab.savefig('figures/%s_depth_distribution.png' % species, bbox_inches='tight',dpi=300)


pylab.figure(2)
pylab.semilogx(median_coverages, iqrs/median_coverages, 'ko')
#pylab.loglog([1,1000],[1,1000],'k:')

pylab.xlabel('Median genomic coverage')
pylab.ylabel('Inter-quartile range / median')

pylab.savefig('figures/%s_iqr_vs_median_coverage.pdf' % species, bbox_inches='tight')
pylab.savefig('figures/%s_iqr_vs_median_coverage.png' % species, bbox_inches='tight',dpi=300)


pylab.figure(3)

for i in xrange(0,len(sample_coverage_histograms)):
    sample_coverage_histogram = sample_coverage_histograms[i]
    x0 = median_coverages[i]
    
    if x0<20:
        continue
    
    #x0 = marker_coverages[i]
    xs, CDFs = stats_utils.calculate_CDF_from_histogram(sample_coverage_histogram)
    pylab.plot(xs/x0, 1-CDFs, '-')

pylab.xlim([-0.1,4])
pylab.plot([0.25,0.25],[0,1],'k:')
pylab.plot([3,3],[0,1],'k:')
pylab.xlabel('Coverage relative to median')
pylab.ylabel('Fraction sites with >= coverage')


pylab.savefig('figures/%s_genome_depth_distribution.pdf' % species, bbox_inches='tight')
pylab.savefig('figures/%s_genome_depth_distribution.png' % species, bbox_inches='tight',dpi=300)

pylab.figure(4)

for i in xrange(0,len(sample_coverage_histograms)):
    sample_coverage_histogram = sample_coverage_histograms[i]
    x0 = peak_coverages[i]
    
    if x0<20:
        continue
    
    #x0 = marker_coverages[i]
    xs, CDFs = stats_utils.calculate_CDF_from_histogram(sample_coverage_histogram)
    pylab.plot(xs/x0, 1-CDFs, '-')

pylab.xlim([-0.1,4])
pylab.plot([0.25,0.25],[0,1],'k:')
pylab.plot([3,3],[0,1],'k:')
pylab.xlabel('Coverage relative to median')
pylab.ylabel('Fraction sites with >= coverage')


pylab.savefig('figures/%s_scaled_genome_depth_distribution.pdf' % species, bbox_inches='tight')
pylab.savefig('figures/%s_scaled_genome_depth_distribution.png' % species, bbox_inches='tight',dpi=300)


pylab.figure(5)

num_plotted=0
for gene_name in gene_coverages.keys():

    scale = peak_coverages
    #scale = marker_coverages
    scaled_coverages = numpy.array(sorted(gene_coverages[gene_name][scale>=20]/(scale[scale>=20])))
    
    m = numpy.median(scaled_coverages)
    
    #scaled_coverages /= numpy.median(scaled_coverages)
    
    p = (scaled_coverages>0.5).sum()*1.0/len(scaled_coverages)
    
    if False:
        pylab.step(scaled_coverages, 1-numpy.arange(0,len(scaled_coverages))*1.0/len(scaled_coverages),where='post',color='0.7',alpha=0.5,linewidth=0.5,zorder=0)
    else:
        num_plotted+=1
        pylab.step(scaled_coverages, 1-numpy.arange(0,len(scaled_coverages))*1.0/len(scaled_coverages),where='post',zorder=1)

print num_plotted        
pylab.xlim([-0.1,3])

pylab.xlabel('Avg coverage of gene / scale')
pylab.ylabel('Fraction of samples with >= scaled coverage')


pylab.savefig('figures/%s_scaled_gene_depth_distribution.pdf' % species, bbox_inches='tight')
pylab.savefig('figures/%s_scaled_gene_depth_distribution.png' % species, bbox_inches='tight',dpi=300)


pylab.figure(6)

num_plotted=0
for gene_name in gene_coverages.keys():

    #scale = peak_coverages
    scale = median_coverages
    #scale = marker_coverages
    raw_coverages = numpy.array(sorted(gene_coverages[gene_name][scale>=20]))
    
    p = (raw_coverages>=5).sum()*1.0/len(raw_coverages)
    
    if p>=0.95:
        pylab.step(raw_coverages, 1-numpy.arange(0,len(raw_coverages))*1.0/len(raw_coverages),where='post',color='0.7',alpha=0.5,linewidth=0.5,zorder=0)
    else:
        num_plotted+=1
        pylab.step(raw_coverages, 1-numpy.arange(0,len(raw_coverages))*1.0/len(raw_coverages),where='post',zorder=1)

pylab.semilogx([1],[-1])
pylab.ylim([0,1.1])        
pylab.xlim([1,1000])
pylab.xlabel('Avg coverage of gene')
pylab.ylabel('Fraction of samples with >= coverage')

print num_plotted
pylab.savefig('figures/%s_gene_depth_distribution.pdf' % species, bbox_inches='tight')
pylab.savefig('figures/%s_gene_depth_distribution.png' % species, bbox_inches='tight',dpi=300)


pylab.figure(7)
for i in xrange(0,len(sample_coverage_histograms)):
    sample_coverage_histogram = sample_coverage_histograms[i]
    x0 = median_coverages[i]
    
    if x0==0:
        continue
    #x0 = marker_coverages[i]
    xs = sorted(sample_coverage_histogram.keys())
    pmfs = numpy.array([sample_coverage_histogram[x] for x in xs])
    pmfs /= pmfs.sum()
    
    xs = numpy.array(xs)
    
    xs = xs/x0
    
    dx = 0.1
    coarse_xs = numpy.linspace(0,4,400)
    coarse_pmfs = []
    
    #print "--"
    for x in coarse_xs:
        numerator = (pmfs*numpy.exp(-numpy.square(xs-x)/2/dx/dx)).sum()
        denominator = (numpy.exp(-numpy.square(xs-x)/2/dx/dx).sum())
        #print numerator, denominator
        coarse_pmfs.append( numerator / denominator )
    
    mode_idxs = []
    for i in xrange(1,len(coarse_xs)-1):
        if (coarse_pmfs[i]>coarse_pmfs[i-1]) and (coarse_pmfs[i]>coarse_pmfs[i+1]):
            mode_idxs.append(i)
            
    modes = numpy.array([coarse_xs[i] for i in mode_idxs])
    peak_pmfs = numpy.array([coarse_pmfs[i] for i in mode_idxs])
    
    x0_peak = modes[peak_pmfs.argmax()]
    x0_1 = modes[numpy.square(modes-1).argmin()]
    x0 = x0_peak
    
    #if x0_peak > x0_1:
    pylab.plot(coarse_xs/x0, coarse_pmfs, '-')
    
    pylab.xlim([-0.1,4])


#pylab.show()