import os
import sys

outFile=open(os.path.expanduser('~/ben_nandita_hmp_analysis/supplemental_table_sample_ids.txt'),'w')

inFN_hmp=open(os.path.expanduser('~/ben_nandita_hmp_scripts/HMP_ids_order.txt'),'r')
inFN_qin=open(os.path.expanduser('~/ben_nandita_hmp_scripts/qin_ids.txt'),'r')

inFN_hmp.readline() # header
inFN_qin.readline() # header
outFile.write('subject_id\tsample_id\trun_accession\tcountry\tcontinent\tVISNO\tstudy\n')

for line in inFN_hmp:
    outFile.write(line.strip() + '\tHMP\n')

for line in inFN_qin:
    outFile.write(('\t').join(line.strip().split('\t')[:-1]) + '\tNA\tQin et al. 2012\n')
