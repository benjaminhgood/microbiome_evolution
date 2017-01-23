import parse_midas_data

subject_sample_map = parse_midas_data.parse_subject_sample_map()

for subject in subject_sample_map.keys(): # loop over subjects (hosts)
    for sample in subject_sample_map[subject].keys(): # loop over samples
        print "Accessions for", sample
        for accession in subject_sample_map[subject][sample]:
            print accession