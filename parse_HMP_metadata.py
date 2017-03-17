import os
import pylab



# relevant input files
days_FN    ='/Users/nanditagarud/Documents/microbiome/BenNanditaProject/HMP_metadata/phs000228.v3.pht002157.v1.p1.c1.EMMES_HMP_DTP_DHX_DVD.HMP.txt'
HMP_ids_FN ='/Users/nanditagarud/Dropbox (Gladstone)/ben_nandita_collaboration/midas_python_code/HMP_ids.txt'
IDs_FN     ='/Users/nanditagarud/Documents/microbiome/BenNanditaProject/HMP_metadata/SRS2RAND.txt'
visNo_FN   ='/Users/nanditagarud/Documents/microbiome/BenNanditaProject/HMP_metadata/project_catalog.csv'
SRA_mapping_FN = '/Users/nanditagarud/Documents/microbiome/BenNanditaProject/HMP_metadata/sra_mapping_HMP.txt'
visNo2_FN= '/Users/nanditagarud/Documents/microbiome/BenNanditaProject/HMP_metadata/phs000228.v3.pht001193.v3.p1.c1.EMMES_HMP_GTV.HMP.txt'

days    = open(days_FN,"r")
HMP_ids = open(HMP_ids_FN,"r")
IDs     = open(IDs_FN,"r")
visNo   = open(visNo_FN,"r")
SRA_mapping = open(SRA_mapping_FN, "r")
visNo2     = open(visNo2_FN, "r")

# read HMP_ids into a dictionary with run_accession as the key
HMP_ids_dict={}
HMP_header=HMP_ids.readline() #header
for line in HMP_ids:
    run_accession=line.split('\t')[2]
    HMP_ids_dict[run_accession]=line.strip().split('\t')


# read days_FN into a dictionary with RANDSID and VISNO as the key and [dbGaP.SubjID, study_day] as the values


days_dict={}
days.readline() #header
for line in days:
    fields=line.strip().split('\t')
    RANDSID=fields[162]
    dbGaP_SubjID=fields[0]
    VISNO=fields[5]
    study_day=fields[6]
    if RANDSID not in days_dict.keys():
        days_dict[RANDSID]={}
    days_dict[RANDSID][VISNO]=[dbGaP_SubjID, study_day]

# read IDs into a dictionary with run_accession as the key
IDs_dict={}
IDs.readline() #header

for line in IDs:
    fields=line.strip().split('\t')
    run_accession=fields[0]
    RANDSID=fields[1]
    IDs_dict[run_accession]=RANDSID

# read the visit numbers into a dicationary with run_accession as the key
visNo_dict={}
visNo.readline() #header
for line in visNo:
    fields=line.strip().split(',')
    run_accession=fields[0]
    visno_string=fields[2].split()
    visit_number=0
    subject=0
    for i in range(0, len(visno_string)):
        if visno_string[i]=='number':
            visit_number=visno_string[i+1]
        elif visno_string[i]=='subject':
            subject=visno_string[i+1]  
    visNo_dict[run_accession]=[visit_number, subject]




# read the SRA_mapping to a dictionary with sample_accession as key
SRA_dict={}
SRA_mapping.readline()

for line in SRA_mapping:
    run_accession=line.strip().split()[0]
    sample_accession=line.strip().split()[1]
    SRA_dict[run_accession]=sample_accession

visNo2_dict = {}
visNo2_dict_alt={}

for i in range(0,10):
    visNo2.readline()

header_visno=visNo2.readline()

for line in visNo2:
    VISNO=line.split('\t')[7][1]
    sample_id=line.split('\t')[10]
    sample_id_alt=line.split('\t')[6]
    RANDSID= line.split('\t')[13]
    visNo2_dict[sample_id]=VISNO
    visNo2_dict_alt[sample_id_alt]=[VISNO, RANDSID]

# outFile:
outFN='/Users/nanditagarud/Documents/microbiome/BenNanditaProject/HMP_metadata/HMP_ids_time.txt'
outFile=open(outFN,'w')
outFile.write(HMP_header.strip() +'\tVISNO\tstudy_day\n')

# match all the dictionaries:
for key in HMP_ids_dict.keys():
    RANDSID=HMP_ids_dict[key][0]
    sample_id=HMP_ids_dict[key][1]
    if key in SRA_dict.keys():
        SRS=SRA_dict[key]
        if SRS in visNo_dict.keys():
            VISNO=str(visNo_dict[SRS][0])
            RANDSID=IDs_dict[SRS]
            study_day=days_dict[RANDSID][VISNO][1]
            outFile.write('\t'.join(HMP_ids_dict[key]) +'\t' + VISNO + '\t' + study_day +'\n')
        else:
#            outFile.write('\t'.join(HMP_ids_dict[key]) +'\tNA\tNA\n')
            if sample_id in visNo2_dict.keys():
                VISNO=visNo2_dict[sample_id]
                study_day=days_dict[RANDSID][VISNO][1]
                outFile.write('\t'.join(HMP_ids_dict[key]) +'\t' + VISNO + '\t' + study_day +'\n')
            elif sample_id in visNo2_dict_alt.keys():
                if visNo2_dict_alt[sample_id][1] == RANDSID:
                    VISNO=visNo2_dict_alt[sample_id][0]
                    study_day=days_dict[RANDSID][VISNO][1] 
                    outFile.write('\t'.join(HMP_ids_dict[key]) +'\t' + VISNO + '\t' + study_day +'\n')
                print sample_id
            else:
                print '\t'.join(HMP_ids_dict[key]) +'\tNA\tNA\n'





