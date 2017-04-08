while read PATRIC; do 
  file=`ls /pollard/shattuck0/snayfach/databases/PATRIC/genomes/${PATRIC}/${PATRIC}.PATRIC.features.tab*`

  if file --mime-type "$file" | grep -q gzip$; then
      zcat /pollard/shattuck0/snayfach/databases/PATRIC/genomes/${PATRIC}/${PATRIC}.PATRIC.features.tab | cut -f6,21 > ~/ben_nandita_hmp_data/kegg/${PATRIC}.kegg.txt
  else
      cat /pollard/shattuck0/snayfach/databases/PATRIC/genomes/${PATRIC}/${PATRIC}.PATRIC.features.tab | cut -f6,21 > ~/ben_nandita_hmp_data/kegg/${PATRIC}.kegg.txt  
	  
  fi
done < ~/ben_nandita_hmp_data/patric_genomes.txt


