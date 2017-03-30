import matplotlib  
matplotlib.use('Agg') 
import parse_midas_data
import pylab
import sys
import numpy
import diversity_utils
import gene_diversity_utils
import stats_utils
import os

# Load time metadata
subject_sample_time_map = parse_midas_data.parse_subject_sample_time_map()

# store data in these variables
num_visnos={1:0,2:0,3:0}
num_samples_per_visno_aggregate={1:0}
num_samples_per_visno={1:{},2:{},3:{}}
days=[] # get a list of all days
days_by_visno={1:[],2:[],3:[]}
distance_between_days={'1-2':[],'1-3':[],'2-3':[]}


# iterate through data and store in variables above
for subject in subject_sample_time_map.keys():
    num_visnos[len(subject_sample_time_map[subject].keys())] +=1
    visnos=subject_sample_time_map[subject].keys()
    for vis in visnos:
        # compute num samples per vis
        num_samples_per_vis= len(subject_sample_time_map[subject][vis])
        if num_samples_per_vis not in num_samples_per_visno_aggregate.keys():
            num_samples_per_visno_aggregate[num_samples_per_vis]=1
        else:
            num_samples_per_visno_aggregate[num_samples_per_vis]+=1
            
        # store num samples per vis by visit number
        if num_samples_per_vis not in num_samples_per_visno[vis].keys():
            num_samples_per_visno[vis][num_samples_per_vis]=1
        else:
            num_samples_per_visno[vis][num_samples_per_vis]+=1    

        # store the days
        day=subject_sample_time_map[subject][vis][0][1]
        days.append(day)
        days_by_visno[vis].append(day)
        
    # compute differences between days
    if 1 in visnos and 2 in visnos:
        distance_between_days['1-2'].append(subject_sample_time_map[subject][2][0][1] - subject_sample_time_map[subject][1][0][1])
    if 1 in visnos and 3 in visnos:
        distance_between_days['1-3'].append(subject_sample_time_map[subject][3][0][1] - subject_sample_time_map[subject][1][0][1])
    if 2 in visnos and 3 in visnos:
        distance_between_days['2-3'].append(subject_sample_time_map[subject][3][0][1] - subject_sample_time_map[subject][2][0][1])



# plot the metadata

# plot distribution of days

pylab.figure()
pylab.xlabel('days')
pylab.ylabel('number of samples')
pylab.title('Distribution of days')
pylab.hist(days)
pylab.savefig('%s/metadata_days.png' % (parse_midas_data.analysis_directory),bbox_inches='tight', dpi=300)

# plot distribution of days for visnos 1-2, 1-3, 2-3:

fig=pylab.figure()

ax1=fig.add_subplot(311)
pylab.ylim(0,20)
pylab.xlim(0,500)
pylab.title('Days between visnos 1-2')
ax1.hist(distance_between_days['1-2'])

ax2=fig.add_subplot(312)
pylab.ylabel('number of samples')
pylab.ylim(0,20)
pylab.xlim(0,500)
pylab.title('Days between visnos 1-3')
ax2.hist(distance_between_days['1-3'])

ax3=fig.add_subplot(313)
pylab.ylim(0,20)
pylab.xlim(0,500)
pylab.xlabel('days')
pylab.title('Days between visnos 2-3')
ax3.hist(distance_between_days['2-3'])


fig.savefig('%s/metadata_days_by_visno.png' % (parse_midas_data.analysis_directory),bbox_inches='tight', dpi=300)


# Plot number of visnos per person
pylab.figure()       
pylab.xlabel('Number of visits/subject')     
pylab.ylabel('Number of subjects') 
pylab.bar([1,2,3],[num_visnos[1],num_visnos[2],num_visnos[3]])
pylab.plot([1,2,3],[num_visnos[1],num_visnos[2],num_visnos[3]], color='y')
pylab.xticks([1,2,3], ['1', '2', '3'])
pylab.savefig('%s/metadata_num_visnos_per_subject.png' % (parse_midas_data.analysis_directory),bbox_inches='tight', dpi=300)

# Plot number of samples per visno

pylab.figure()       
pylab.xlabel('Number of samples/visno/subject')     
pylab.ylabel('Number of visnos') 
pylab.bar([1,2,3,4],[num_samples_per_visno_aggregate[1],num_samples_per_visno_aggregate[2],num_samples_per_visno_aggregate[3],num_samples_per_visno_aggregate[4]])
pylab.plot([1,2,3,4],[num_samples_per_visno_aggregate[1],num_samples_per_visno_aggregate[2],num_samples_per_visno_aggregate[3],num_samples_per_visno_aggregate[4]], color='y')
pylab.xticks([1,2,3,4], ['1', '2', '3','4'])
pylab.savefig('%s/metadata_num_samples_per_visno.png' % (parse_midas_data.analysis_directory),bbox_inches='tight', dpi=300)


# plot distribution of num samples per visno
pylab.figure()       
pylab.xlabel('Number of replicate samples/visit/subject')     
pylab.ylabel('Number of subjects') 
pylab.title('Number of replicate samples/visit/subject')
width=0.1
pylab.bar([1,2,3,4],[num_samples_per_visno[1][1],num_samples_per_visno[1][2],num_samples_per_visno[1][3],num_samples_per_visno[1][4]],width)
pylab.bar([1+width,2+width,3+width,4+width],[num_samples_per_visno[2][1],num_samples_per_visno[2][2],0,0 ],width, color='r')
pylab.bar([1+2*width,2+2*width,3+2*width,4+2*width],[num_samples_per_visno[3][1],num_samples_per_visno[3][2],0,0 ],width, color='g')

pylab.xticks([1,2,3,4], ['1', '2', '3','4'])

pylab.legend(['Visit1','Visit2','Visit3'],'upper right',prop={'size':6})

pylab.savefig('%s/metadata_distribution_of_num_samples_visno_all.png' % (parse_midas_data.analysis_directory),bbox_inches='tight', dpi=300)
