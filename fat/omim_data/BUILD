#! /bin/zsh

wget 'https://data.omim.org/downloads/Le_sy8F_SBKL0hYy6Bek9g/genemap.txt'
wget 'https://data.omim.org/downloads/Le_sy8F_SBKL0hYy6Bek9g/mimTitles.txt'

catUncomment mimTitles.txt | awk -F '\t' '$1 == "Number Sign"' | cut -f 2 | sort -u > omim_diseases_with_molecular_basis.txt

cat omim_diseases_with_molecular_basis.txt | while read i
do
  wget "http://api.omim.org/api/entry?mimNumber=${i}&apiKey=Le_sy8F_SBKL0hYy6Bek9g&include=text&include=referenceList" -O molecularGenetics/${i}.xml
done

for i in molecularGenetics/*
do
  ./parseMolecularGenetics.py $i
done > diseases_to_molecular_genetics.tsv
sort -u diseases_to_molecular_genetics.tsv -o diseases_to_molecular_genetics.tsv
sort -t $'\t' -k2,2 diseases_to_molecular_genetics.tsv -o diseases_to_molecular_genetics.tsv
join -t $'\t' -1 2 -2 1 diseases_to_molecular_genetics.tsv \
  <(cat ../../../pmid_to_year_month.tsv) > all_papers_to_omim_to_date.tsv

for i in molecularGenetics/*
do
  ./parseMolecularGeneticsFirstPMID.py $i
done > diseases_to_molecular_genetics_first_pmid.tsv
sort -u diseases_to_molecular_genetics_first_pmid.tsv -o diseases_to_molecular_genetics_first_pmid.tsv
sort -t $'\t' -k2,2 diseases_to_molecular_genetics_first_pmid.tsv -o diseases_to_molecular_genetics_first_pmid.tsv
join -t $'\t' -1 2 -2 1 diseases_to_molecular_genetics_first_pmid.tsv \
  <(cat ../../../pmid_to_year_month.tsv) > paper_to_omim_to_date.tsv
