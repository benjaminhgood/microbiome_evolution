import os
import pylab

# Read in file matching sample id with visit number
visNo_FN= '/Users/nanditagarud/Documents/microbiome/BenNanditaProject/HMP_metadata/phs000228.v3.pht001193.v3.p1.c1.EMMES_HMP_GTV.HMP.txt'
visNo   = open(visNo_FN, "r")

visNo_dict={}

for i in range(0,10):
    visNo.readline()

header_visno=visNo.readline()

for line in visNo:
    day = line.strip().split('\t')[5] 
    VISNO=line.strip().split('\t')[7][1]
    sample_id=line.strip().split('\t')[6]
    RANDSID= line.strip().split('\t')[13]
    visNo_dict[sample_id]=[VISNO, day]

# outFile:
# open the HMP file to read in each line to match up the subject ID and sample ID to the metadata. 
HMP_ids_FN= '/Users/nanditagarud/Documents/microbiome/BenNanditaProject/HMP_metadata/HMP_ids.txt'
HMP_ids = open(HMP_ids_FN,"r")
HMP_header=HMP_ids.readline() #header

outFN='/Users/nanditagarud/Documents/microbiome/BenNanditaProject/HMP_metadata/HMP_ids_time.txt'
outFile=open(outFN,'w')
outFile.write(HMP_header.strip() +'\tVISNO\n')

for line in HMP_ids:  
    sample_id=line.split('\t')[1]
    VISNO=visNo_dict[sample_id][0]
    day=visNo_dict[sample_id][1]
    outFile.write(line.strip() +'\t' + VISNO + '\n' ) 


# other files to cross-check the visnos:


'''
# read in the accessions from the HMP1-2 data file from the web:
# for now I am not using this data, but can be used to get the visnos as well 
inFN_HMP2='/Users/nanditagarud/Documents/microbiome/BenNanditaProject/HMP_metadata/HMP_samples_website.txt'
HMP2=open(inFN_HMP2)
header_HMP2=HMP2.readline()

visNo_web_dict={} #sample_id -> visno, site
for line in HMP2:
    sample_id=line.strip().split('\t')[0]
    subject_id=line.strip().split('\t')[1]
    visno=line.strip().split('\t')[2]
    site=line.strip().split('\t')[3]
    accession=line.strip().split('\t')[6]
    visNo_web_dict[sample_id] = [visno, site]


# read in data from Owen White
SRA_mapping_FN = '/Users/nanditagarud/Documents/microbiome/BenNanditaProject/HMP_metadata/sra_mapping.txt'
SRA_mapping = open(SRA_mapping_FN, "r")

SRA_dict={}
SRA_mapping.readline()
for line in SRA_mapping:
    sample_id=line.strip().split()[2]
    SRS=line.strip().split()[1]
    SRA_dict[SRS]=sample_id


inFN_OW='/Users/nanditagarud/Documents/microbiome/BenNanditaProject/HMP_metadata/wgs_sample_info.csv'
inFile_OW=open(inFN_OW,'r')
header_OW=inFile_OW.readline()
visno_OW_dict = {} #sample_id ->visno

for line in inFile_OW:
    items=line.strip().split(',')
    SRS=items[3]
    visno=items[4]
    if SRS in SRA_dict.keys():
        sample_id=SRA_dict[SRS]
    else:
        print SRS +'\tnot in SRA dict'
    visno_OW_dict[sample_id]=visno

'''
