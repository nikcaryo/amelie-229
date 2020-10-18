#! /bin/zsh -x

download_data --num_threads 8 --omim_key "jDJzYZ0jQAuhiFF-c0s32g" --hgmd_location /cluster/data/labResources/medgen/hgmd/HGMD_PRO_2018.1/allmut-table-2018.1.tsv amelie_out_dir/

download_all_papers --elsevier_key "9375813031e16ecaa64295d9b45976e6" --sfx_server "http://sul-sfx.stanford.edu/sfxlcl41" --conda_path /cluster/u/jbirgmei/anaconda3/bin --negative_pmids_list negatives.json amelie_out_dir amelie_papers

process_all_papers --conda_path /cluster/u/jbirgmei/anaconda3/bin --pubmunch_data_dir /cluster/u/jbirgmei/shallowfloat/amelie_ethan/pubmunchData --delete_parasol_dir True amelie_out_dir amelie_papers amelie_process_dir

train_classifiers --forbidden_genesymbols_file ega215_forbidden_genesymbols.tsv amelie_out_dir amelie_process_dir amelie_ta_genephenos

extract_all_papers --pubmunch_data_dir /cluster/u/jbirgmei/shallowfloat/amelie_ethan/pubmunchData --conda_path /cluster/u/jbirgmei/miniconda3/bin --max_jobs 50 --elsevier_key 9375813031e16ecaa64295d9b45976e6 ${PWD}/amelie_out_dir ${PWD}/amelie_papers ${PWD}/amelie_process_dir

# for stupid stats that I forgot to collect:
extract_all_papers --pubmunch_data_dir /cluster/u/jbirgmei/shallowfloat/amelie_ethan/pubmunchData --conda_path /cluster/u/jbirgmei/miniconda3/bin --max_jobs 50 --elsevier_key 9375813031e16ecaa64295d9b45976e6 ${PWD}/amelie_out_dir ${PWD}/amelie_papers_DUMMY_SHOULD_BE_EMPTY ${PWD}/amelie_process_dir_DUMMY_SHOULD_BE_EMPTY --title_abstract_relevance_classify_only --no_pubtator_download --max_update_number 1280 --extracted_subdir extracted_DUMMY_SHOULD_BE_EMPTY

# cross-validation performance:
train_classifiers --forbidden_genesymbols_file ega215_forbidden_genesymbols.tsv amelie_out_dir amelie_process_dir amelie_ta_genephenos --do_not_save --cross_validate 2>&1 | tee train_classifiers.out.txt
train_title_abstract --forbidden_genesymbols_file ega215_forbidden_genesymbols.tsv amelie_out_dir amelie_process_dir amelie_ta_genephenos --do_not_save --cross_validate --l1_regularization 2>&1 | tee train_title_abstract_l1.out.txt

# to find median gene and pheno mentions, see gene_pheno_stats/


# to find stats about phenotype mentions:
sqlite3 amelie_knowledgebase_HPO_2018-07.sqlite3 'select pmid, count(distinct hpo_id) from amelie_paper p left outer join amelie_paper_phenotypes pph on (p.id = pph.paper_id) left outer join amelie_phenotype ph on (pph.phenotype_id = ph.id) group by pmid;' | tr '|' '\t' | datamash median 2 mean 2

# to find stats about variant mentions:
sqlite3 /dev/shm/amelie_knowledgebase_HPO_2018-07.sqlite3 'select pmid, count(pmid) from (select distinct pmid, chrom, vcf_pos, vcf_ref, vcf_alt from amelie_paper p join amelie_extractedvariant v on p.id = v.paper_id) group by pmid;' | tr '|' '\t' | datamash count 1 median 2

# train amelie:
train_amelie --patient_features_dir amelie_patient_features --test_patients_dir amelie/amelie/patient_data/amelie_real_EXAC_ALL --test_solutions_list amelie/amelie/patient_data/amelie_real.solutions.json /dev/shm/amelie_knowledgebase_HPO_2018-07.sqlite3 amelie_out_dir amelie_out_dir/1000g_EXAC_ALL amelie_out_dir/fake_patients_1000g_EXAC_ALL_fake_phenotypes_post_2011.tsv.json --count_filter --hmct_cutoff 1 --dominant_alct_cutoff 3 --rankings_file amelie_real_EXAC_ALL.rankings.txt --ignore_if_no_causative_gene
# test if AMELIE performs better if a gene has more papers:
test_amelie amelie_out_dir amelie/amelie/patient_data/amelie_real_EXAC_ALL amelie/amelie/patient_data/amelie_real.solutions.json --count_filter --hmct_cutoff 1 --dominant_alct_cutoff 3 /dev/shm/amelie_knowledgebase_HPO_2018-07.sqlite3 --rankings_file amelie_real_EXAC_ALL.rankings_with_num_papers.txt --output_num_papers_for_causative_gene 
amelie_rank_vs_num_articles amelie_real_EXAC_ALL.rankings_with_num_papers.txt amelie_real_EXAC_ALL.rank_vs_num_papers.png
join -t $'\t' -1 1 -2 1 <(sort -t $'\t' -k1,1 amelie_real_EXAC_ALL.rankings_with_num_papers.txt) <(sort -t $'\t' -k1,1 amelie_real_EXAC_ALL.exomiser.txt) | awk '$3 <= 10' | cut -f 1,5 > amelie_real_EXAC_ALL.exomiser_with_less_than_10_papers.txt
join -t $'\t' -1 1 -2 1 <(sort -t $'\t' -k1,1 amelie_real_EXAC_ALL.rankings_with_num_papers.txt) <(sort -t $'\t' -k1,1 amelie_real_EXAC_ALL.exomiser.txt) | awk '$3 <= 10' | cut -f 1,2 > amelie_real_EXAC_ALL.rankings_with_less_than_10_papers.txt
amelie_p_values amelie/amelie/patient_data/amelie_real_EXAC_ALL amelie/amelie/patient_data/amelie_real.solutions.json --count_filter --hmct_cutoff 1 --dominant_alct_cutoff 3 --rankings_files amelie_real_EXAC_ALL.rankings_with_less_than_10_papers.txt,amelie_real_EXAC_ALL.exomiser_with_less_than_10_papers.txt --rankings_files_names 'AMELIE,Exomiser' --nonscientific_format --subset_to_list1_patients

for i in {1..4}
do
  num_ranked=`cat amelie_real_EXAC_ALL.rankings_with_num_papers.txt | awk -v i=$i '$2 == i' | wc -l`
  echo -n "$i\t$num_ranked\t"
  cat amelie_real_EXAC_ALL.rankings_with_num_papers.txt | awk -v i=$i '$2 == i' | datamash median 3
  echo
done | awk '$0 != ""'
cat amelie_real_EXAC_ALL.rankings_with_num_papers.txt | awk '$2 > 4' | datamash median 3

# feature ablation:
# amelie_features_no_domrec.txt
# amelie_features_no_evs.txt
# amelie_features_no_phrank.txt
# amelie_features_no_ptv_ntv.txt
# amelie_features_no_relevance.txt
# amelie_features_no_scores.txt
train_amelie --test_patients_dir amelie/amelie/patient_data/amelie_real_EXAC_ALL --test_solutions_list amelie/amelie/patient_data/amelie_real.solutions.json /dev/shm/amelie_knowledgebase_HPO_2018-07.sqlite3 amelie_out_dir amelie_out_dir/1000g_EXAC_ALL amelie_out_dir/fake_patients_1000g_EXAC_ALL_fake_phenotypes_post_2011.tsv.json --count_filter --hmct_cutoff 1 --dominant_alct_cutoff 3 --rankings_file amelie_real_EXAC_ALL.no_domrec.rankings.txt --feature_names_list amelie/amelie/amelie/classifier_data/amelie_features_no_domrec.txt --classifier_filename amelie_out_dir/amelie_classifier_no_domrec.pkl --ignore_if_no_causative_gene --features_dict_dir features_dict_dir
train_amelie --test_patients_dir amelie/amelie/patient_data/amelie_real_EXAC_ALL --test_solutions_list amelie/amelie/patient_data/amelie_real.solutions.json /dev/shm/amelie_knowledgebase_HPO_2018-07.sqlite3 amelie_out_dir amelie_out_dir/1000g_EXAC_ALL amelie_out_dir/fake_patients_1000g_EXAC_ALL_fake_phenotypes_post_2011.tsv.json --count_filter --hmct_cutoff 1 --dominant_alct_cutoff 3 --rankings_file amelie_real_EXAC_ALL.no_evs.rankings.txt --feature_names_list amelie/amelie/amelie/classifier_data/amelie_features_no_evs.txt --classifier_filename amelie_out_dir/amelie_classifier_no_evs.pkl --ignore_if_no_causative_gene --features_dict_dir features_dict_dir
train_amelie --test_patients_dir amelie/amelie/patient_data/amelie_real_EXAC_ALL --test_solutions_list amelie/amelie/patient_data/amelie_real.solutions.json /dev/shm/amelie_knowledgebase_HPO_2018-07.sqlite3 amelie_out_dir amelie_out_dir/1000g_EXAC_ALL amelie_out_dir/fake_patients_1000g_EXAC_ALL_fake_phenotypes_post_2011.tsv.json --count_filter --hmct_cutoff 1 --dominant_alct_cutoff 3 --rankings_file amelie_real_EXAC_ALL.no_phrank.rankings.txt --feature_names_list amelie/amelie/amelie/classifier_data/amelie_features_no_phrank.txt --classifier_filename amelie_out_dir/amelie_classifier_no_phrank.pkl --ignore_if_no_causative_gene --no_phrank_interpolation --features_dict_dir features_dict_dir
train_amelie --test_patients_dir amelie/amelie/patient_data/amelie_real_EXAC_ALL --test_solutions_list amelie/amelie/patient_data/amelie_real.solutions.json /dev/shm/amelie_knowledgebase_HPO_2018-07.sqlite3 amelie_out_dir amelie_out_dir/1000g_EXAC_ALL amelie_out_dir/fake_patients_1000g_EXAC_ALL_fake_phenotypes_post_2011.tsv.json --count_filter --hmct_cutoff 1 --dominant_alct_cutoff 3 --rankings_file amelie_real_EXAC_ALL.no_ptv_ntv.rankings.txt --feature_names_list amelie/amelie/amelie/classifier_data/amelie_features_no_ptv_ntv.txt --classifier_filename amelie_out_dir/amelie_classifier_no_ptv_ntv.pkl --ignore_if_no_causative_gene --features_dict_dir features_dict_dir
train_amelie --test_patients_dir amelie/amelie/patient_data/amelie_real_EXAC_ALL --test_solutions_list amelie/amelie/patient_data/amelie_real.solutions.json /dev/shm/amelie_knowledgebase_HPO_2018-07.sqlite3 amelie_out_dir amelie_out_dir/1000g_EXAC_ALL amelie_out_dir/fake_patients_1000g_EXAC_ALL_fake_phenotypes_post_2011.tsv.json --count_filter --hmct_cutoff 1 --dominant_alct_cutoff 3 --rankings_file amelie_real_EXAC_ALL.no_relevance.rankings.txt --feature_names_list amelie/amelie/amelie/classifier_data/amelie_features_no_relevance.txt --classifier_filename amelie_out_dir/amelie_classifier_no_relevance.pkl --ignore_if_no_causative_gene --features_dict_dir features_dict_dir
train_amelie --test_patients_dir amelie/amelie/patient_data/amelie_real_EXAC_ALL --test_solutions_list amelie/amelie/patient_data/amelie_real.solutions.json /dev/shm/amelie_knowledgebase_HPO_2018-07.sqlite3 amelie_out_dir amelie_out_dir/1000g_EXAC_ALL amelie_out_dir/fake_patients_1000g_EXAC_ALL_fake_phenotypes_post_2011.tsv.json --count_filter --hmct_cutoff 1 --dominant_alct_cutoff 3 --rankings_file amelie_real_EXAC_ALL.no_scores.rankings.txt --feature_names_list amelie/amelie/amelie/classifier_data/amelie_features_no_scores.txt --classifier_filename amelie_out_dir/amelie_classifier_no_scores.pkl --ignore_if_no_causative_gene --features_dict_dir features_dict_dir

amelie_feature_ablation_table amelie_real_EXAC_ALL.rankings.txt,amelie_real_EXAC_ALL.no_domrec.rankings.txt,amelie_real_EXAC_ALL.no_evs.rankings.txt,amelie_real_EXAC_ALL.no_phrank.rankings.txt,amelie_real_EXAC_ALL.no_ptv_ntv.rankings.txt,amelie_real_EXAC_ALL.no_relevance.rankings.txt,amelie_real_EXAC_ALL.no_scores.rankings.txt 'Baseline (AMELIE),No inheritance mode,No AVADA variants,No Phrank,No variant types,No gene/article relevance,No pathogenicity scores' DDD,Stanford,MAN,UDN DDD,Stanford,Manton,UDN > amelie_feature_ablation_table.tsv
amelie_p_values amelie/amelie/patient_data/amelie_real_EXAC_ALL amelie/amelie/patient_data/amelie_real.solutions.json --count_filter --hmct_cutoff 1 --dominant_alct_cutoff 3 --rankings_files amelie_real_EXAC_ALL.rankings.txt,amelie_real_EXAC_ALL.no_domrec.rankings.txt,amelie_real_EXAC_ALL.no_evs.rankings.txt,amelie_real_EXAC_ALL.no_phrank.rankings.txt,amelie_real_EXAC_ALL.no_ptv_ntv.rankings.txt,amelie_real_EXAC_ALL.no_relevance.rankings.txt,amelie_real_EXAC_ALL.no_scores.rankings.txt --rankings_files_names 'Baseline (AMELIE),No inheritance mode,No AVADA variants,No Phrank,No variant types,No gene/article relevance,No pathogenicity scores' > feature_ablation.pvalues.txt

amelie_comparisons_info amelie/amelie/patient_data/amelie_real_EXAC_ALL amelie/amelie/patient_data/amelie_real.solutions.json --count_filter --hmct_cutoff 1 --dominant_alct_cutoff 3 --rankings_files amelie_real_EXAC_ALL.rankings.txt,amelie_real_EXAC_ALL.no_domrec.rankings.txt,amelie_real_EXAC_ALL.no_evs.rankings.txt,amelie_real_EXAC_ALL.no_phrank.rankings.txt,amelie_real_EXAC_ALL.no_ptv_ntv.rankings.txt,amelie_real_EXAC_ALL.no_relevance.rankings.txt,amelie_real_EXAC_ALL.no_scores.rankings.txt --rankings_files_names 'Baseline (AMELIE),Inheritance mode,AVADA variants,Phrank,Variant types,Gene/article relevance,Pathogenicity scores' "real patients" --output_image1 feature_ablation_output_image1.png --output_image2 feature_ablation_output_image2.png --output_image3 feature_ablation_output_image3.png --auto_color

amelie_feature_ablation_barchart amelie/amelie/patient_data/amelie_real.solutions.json amelie_real_EXAC_ALL.rankings.txt AMELIE amelie_real_EXAC_ALL.no_domrec.rankings.txt,amelie_real_EXAC_ALL.no_evs.rankings.txt,amelie_real_EXAC_ALL.no_phrank.rankings.txt,amelie_real_EXAC_ALL.no_ptv_ntv.rankings.txt,amelie_real_EXAC_ALL.no_relevance.rankings.txt,amelie_real_EXAC_ALL.no_scores.rankings.txt 'Inheritance mode,AVADA variants,Phrank,Variant types,Gene/article relevance,Pathogenicity scores' 'real patients' feature_ablation_barchart.png

# histogram of number of phenotypes per gene:
sqlite3 amelie_knowledgebase_HPO_2018-07.sqlite3 'select g.gene_symbol, ph.hpo_id from amelie_gene g join amelie_papergene pg on pg.gene_id = g.id join amelie_paper p on pg.paper_id = p.id join amelie_paper_phenotypes pph on pph.paper_id = p.id join amelie_phenotype ph on ph.id = pph.phenotype_id where right_gene_score >= 0.5' | tr '|' '\t' | sort -u > amelie_all_genephenos.tsv
sort -t $'\t' -k1,1 -k2,2 amelie_all_genephenos.tsv -o  amelie_all_genephenos.tsv
cat amelie_all_genephenos.tsv | datamash -g 1 count 2 > gene_to_num_phenos.tsv
sort -t $'\t' -k2,2nr gene_to_num_phenos.tsv -o gene_to_num_phenos.tsv
amelie_phenotypes_histogram gene_to_num_phenos.tsv gene_to_num_phenos.png

# cross-validate on fake patients:
train_amelie /dev/shm/amelie_knowledgebase_HPO_2018-07.sqlite3 amelie_out_dir amelie_out_dir/1000g_EXAC_ALL amelie_out_dir/fake_patients_1000g_EXAC_ALL_fake_phenotypes_post_2011.tsv.json --count_filter --hmct_cutoff 1 --dominant_alct_cutoff 3 --rankings_file fake_patients_EXAC_ALL.rankings.txt --do_not_save --cross_validate 5 --ignore_if_no_causative_gene

train_amelie /dev/shm/amelie_knowledgebase_HPO_2018-07.sqlite3 amelie_out_dir amelie_out_dir/1000g_EXAC_ALL amelie_out_dir/fake_patients_1000g_EXAC_ALL_fake_phenotypes_post_2011.tsv.json --count_filter --hmct_cutoff 1 --dominant_alct_cutoff 3 --rankings_file fake_patients_EXAC_ALL.rankings_with_is_causative_gene_ranked.txt --do_not_save --cross_validate 5 --ignore_if_no_causative_gene --output_if_causative_gene_ranked

# title/abstract only pipeline:
process_all_papers --conda_path /cluster/u/jbirgmei/miniconda3/bin --pubmunch_data_dir /cluster/u/jbirgmei/shallowfloat/amelie_ethan/pubmunchData --title_abstract_only ${PWD}/amelie_out_dir_ta_only ${PWD}/title_abstract_papers ${PWD}/amelie_process_dir_ta_only
train_classifiers --title_abstract_only --forbidden_genesymbols_file ega215_forbidden_genesymbols.tsv amelie_out_dir_ta_only amelie_process_dir_ta_only  # might have to train a couple individually, variant classifier won't train
# before, have 367809 papers in amelie_process_dir_ta_only
extract_all_papers --title_abstract_only --pubmunch_data_dir /cluster/u/jbirgmei/shallowfloat/amelie_ethan/pubmunchData --conda_path /cluster/u/jbirgmei/miniconda3/bin --max_jobs 50 --elsevier_key 9375813031e16ecaa64295d9b45976e6 ${PWD}/amelie_out_dir_ta_only ${PWD}/amelie_papers_ta_only ${PWD}/amelie_process_dir_ta_only
load_db_all_papers --create_db amelie_out_dir_ta_only amelie_knowledgebase_ta_only.sqlite3

train_amelie --test_patients_dir amelie/amelie/patient_data/amelie_real_EXAC_ALL --test_solutions_list amelie/amelie/patient_data/amelie_real.solutions.json /dev/shm/amelie_knowledgebase_ta_only.sqlite3 amelie_out_dir_ta_only amelie_out_dir_ta_only/1000g_EXAC_ALL amelie_out_dir_ta_only/fake_patients_1000g_EXAC_ALL_fake_phenotypes_post_2011.tsv.json --count_filter --hmct_cutoff 1 --dominant_alct_cutoff 3 --rankings_file amelie_real_EXAC_ALL.ta_only_rankings.txt --ignore_if_no_causative_gene

for i in amelie_real_EXAC_ALL.{ta_only_rankings,rankings}.txt
do
  for j in {DDD,Stanford,MAN,UDN,"(Stanford|MAN)"}
  do
  echo -n "$i $j "
  cat $i | egrep $j | awk '$2 == 1' | wc -l
  done
done

# Supplementary table about fake patients:
amelie_simulated_patients_table amelie_out_dir amelie_out_dir/1000g_EXAC_ALL amelie_out_dir/fake_patients_1000g_EXAC_ALL_fake_phenotypes_post_2011.tsv.json --count_filter --hmct_cutoff 1 --dominant_alct_cutoff 3 --ignore_if_no_causative_gene amelie/stm_submission_analyses/genepheno_recall_analysis/ALL_SOURCES_ALL_FREQUENCIES_diseases_to_genes_to_phenotypes.txt --output_file simulated_patients_table.tsv


############ DISGENET
# disgenet BEFREE rankings:
train_amelie --test_patients_dir amelie/amelie/patient_data/amelie_real_EXAC_ALL --test_solutions_list amelie/amelie/patient_data/amelie_real.solutions.json amelie/stm_submission_analyses/disgenet/befree_knowledgebase.sqlite3 amelie_out_dir amelie_out_dir/1000g_EXAC_ALL amelie_out_dir/fake_patients_1000g_EXAC_ALL_fake_phenotypes_post_2011.tsv.json --count_filter --hmct_cutoff 1 --dominant_alct_cutoff 3 --rankings_file amelie_real_EXAC_ALL.disgenet_befree.txt --ignore_if_no_causative_gene --do_not_save

# topical genes:
sqlite3 amelie/stm_submission_analyses/disgenet/befree_knowledgebase.sqlite3 'select pmid, count(distinct ensembl_id) from amelie_paper p left outer join amelie_papergene pg on (p.id = pg.paper_id) left outer join amelie_gene g on (pg.gene_id = g.id) group by pmid;' | tr '|' '\t' | datamash median 2 mean 2

# to find stats about phenotype mentions:
sqlite3 amelie/stm_submission_analyses/disgenet/befree_knowledgebase.sqlite3 'select pmid, count(distinct hpo_id) from amelie_paper p left outer join amelie_paper_phenotypes pph on (p.id = pph.paper_id) left outer join amelie_phenotype ph on (pph.phenotype_id = ph.id) group by pmid;' | tr '|' '\t' | datamash median 2 mean 2

# to find stats about variant mentions:
sqlite3 amelie/stm_submission_analyses/disgenet/befree_knowledgebase.sqlite3 'select pmid, count(pmid) from (select distinct pmid, chrom, vcf_pos, vcf_ref, vcf_alt from amelie_paper p join amelie_extractedvariant v on p.id = v.paper_id) group by pmid;' | tr '|' '\t' | datamash count 1 median 2

# disgenet curated rankings
train_amelie --test_patients_dir amelie/amelie/patient_data/amelie_real_EXAC_ALL --test_solutions_list amelie/amelie/patient_data/amelie_real.solutions.json amelie/stm_submission_analyses/disgenet/disgenet_curated_knowledgebase.sqlite3 amelie_out_dir amelie_out_dir/1000g_EXAC_ALL amelie_out_dir/fake_patients_1000g_EXAC_ALL_fake_phenotypes_post_2011.tsv.json --count_filter --hmct_cutoff 1 --dominant_alct_cutoff 3 --rankings_file amelie_real_EXAC_ALL.disgenet_curated.txt --ignore_if_no_causative_gene --do_not_save

# topical genes:
sqlite3 amelie/stm_submission_analyses/disgenet/disgenet_curated_knowledgebase.sqlite3 'select pmid, count(distinct ensembl_id) from amelie_paper p left outer join amelie_papergene pg on (p.id = pg.paper_id) left outer join amelie_gene g on (pg.gene_id = g.id) group by pmid;' | tr '|' '\t' | datamash count 1 median 2 mean 2

# to find stats about phenotype mentions:
sqlite3 amelie/stm_submission_analyses/disgenet/disgenet_curated_knowledgebase.sqlite3 'select pmid, count(distinct hpo_id) from amelie_paper p left outer join amelie_paper_phenotypes pph on (p.id = pph.paper_id) left outer join amelie_phenotype ph on (pph.phenotype_id = ph.id) group by pmid;' | tr '|' '\t' | datamash count 1 median 2 mean 2

# to find stats about variant mentions:
sqlite3 amelie/stm_submission_analyses/disgenet/disgenet_curated_knowledgebase.sqlite3 'select pmid, count(pmid) from (select distinct pmid, chrom, vcf_pos, vcf_ref, vcf_alt from amelie_paper p join amelie_extractedvariant v on p.id = v.paper_id) group by pmid;' | tr '|' '\t' | datamash count 1 median 2

# disgenet total rankings
train_amelie --test_patients_dir amelie/amelie/patient_data/amelie_real_EXAC_ALL --test_solutions_list amelie/amelie/patient_data/amelie_real.solutions.json amelie/stm_submission_analyses/disgenet/disgenet_total_knowledgebase.sqlite3 amelie_out_dir amelie_out_dir/1000g_EXAC_ALL amelie_out_dir/fake_patients_1000g_EXAC_ALL_fake_phenotypes_post_2011.tsv.json --count_filter --hmct_cutoff 1 --dominant_alct_cutoff 3 --rankings_file amelie_real_EXAC_ALL.disgenet_total.txt --ignore_if_no_causative_gene --do_not_save

# topical genes:
sqlite3 amelie/stm_submission_analyses/disgenet/disgenet_total_knowledgebase.sqlite3 'select pmid, count(distinct ensembl_id) from amelie_paper p left outer join amelie_papergene pg on (p.id = pg.paper_id) left outer join amelie_gene g on (pg.gene_id = g.id) group by pmid;' | tr '|' '\t' | datamash count 1 median 2 mean 2

# to find stats about phenotype mentions:
sqlite3 amelie/stm_submission_analyses/disgenet/disgenet_total_knowledgebase.sqlite3 'select pmid, count(distinct hpo_id) from amelie_paper p left outer join amelie_paper_phenotypes pph on (p.id = pph.paper_id) left outer join amelie_phenotype ph on (pph.phenotype_id = ph.id) group by pmid;' | tr '|' '\t' | datamash count 1 median 2 mean 2

# to find stats about variant mentions:
sqlite3 amelie/stm_submission_analyses/disgenet/disgenet_total_knowledgebase.sqlite3 'select pmid, count(pmid) from (select distinct pmid, chrom, vcf_pos, vcf_ref, vcf_alt from amelie_paper p join amelie_extractedvariant v on p.id = v.paper_id) group by pmid;' | tr '|' '\t' | datamash count 1 median 2

## p value comparison for disgenet:
amelie_p_values amelie/amelie/patient_data/amelie_real_EXAC_ALL amelie/amelie/patient_data/amelie_real.solutions.json --count_filter --hmct_cutoff 1 --dominant_alct_cutoff 3 --rankings_files amelie_real_EXAC_ALL.rankings.txt,amelie_real_EXAC_ALL.disgenet_befree.txt,amelie_real_EXAC_ALL.disgenet_curated.txt,amelie_real_EXAC_ALL.disgenet_total.txt --rankings_files_names 'AMELIE,BeFree,Curated,Total' > amelie_real_EXAC_ALL.disgenet_pvalues.txt

amelie_comparisons_info amelie/amelie/patient_data/amelie_real_EXAC_ALL amelie/amelie/patient_data/amelie_real.solutions.json --count_filter --hmct_cutoff 1 --dominant_alct_cutoff 3 --rankings_files amelie_real_EXAC_ALL.rankings.txt,amelie_real_EXAC_ALL.disgenet_befree.txt,amelie_real_EXAC_ALL.disgenet_curated.txt,amelie_real_EXAC_ALL.disgenet_total.txt --rankings_files_names 'AMELIE,BeFree,Curated,Total' 'real' --output_table disgenet_table.tsv --output_image1 disgenet_comparison1.png --output_image2 disgenet_comparison2.png --output_image3 disgenet_comparison3.png --auto_color

############# comparisons to other methods:
amelie_comparisons_info amelie/amelie/patient_data/ega_215_EXAC_ALL amelie/amelie/patient_data/ega_215.solutions.json "DDD" --count_filter --hmct_cutoff 1 --dominant_alct_cutoff 3 --rankings_files amelie_real_EXAC_ALL.rankings.txt,amelie_real_EXAC_ALL.exomiser.txt,amelie_real_EXAC_ALL.phenolyzer.txt,amelie_real_EXAC_ALL.extasy.txt,amelie_real_EXAC_ALL.phengen.txt,amelie_real_EXAC_ALL.pubcasefinder.txt --rankings_files_names AMELIE,Exomiser,Phenolyzer,eXtasy,Phen-Gen,PubCaseFinder --output_table ega_215.comparison_output_table.txt --output_image1 ega_215.comparison_output_image1.png --output_image2 ega_215.comparison_output_image2.png --output_image3 ega_215.comparison_output_image3.png --add_phenotypes

amelie_p_values amelie/amelie/patient_data/amelie_real_EXAC_ALL amelie/amelie/patient_data/ega_215.solutions.json --count_filter --hmct_cutoff 1 --dominant_alct_cutoff 3 --rankings_files amelie_real_EXAC_ALL.rankings.txt,amelie_real_EXAC_ALL.exomiser.txt,amelie_real_EXAC_ALL.extasy.txt,amelie_real_EXAC_ALL.phengen.txt,amelie_real_EXAC_ALL.phenolyzer.txt,amelie_real_EXAC_ALL.pubcasefinder.txt --rankings_files_names 'AMELIE,Exomiser,eXtasy,Phen-Gen,Phenolyzer,PubCaseFinder' > ega_215_EXAC_ALL.pvalues.txt

amelie_comparisons_info amelie/amelie/patient_data/amelie_2nd_test_set_EXAC_ALL amelie/amelie/patient_data/stanford_manton.solutions.json "Stanford/Manton" --count_filter --hmct_cutoff 1 --dominant_alct_cutoff 3 --rankings_files amelie_real_EXAC_ALL.rankings.txt,amelie_real_EXAC_ALL.exomiser.txt,amelie_real_EXAC_ALL.phenolyzer.txt,amelie_real_EXAC_ALL.extasy.txt,amelie_real_EXAC_ALL.phengen.txt,amelie_real_EXAC_ALL.pubcasefinder.txt --rankings_files_names AMELIE,Exomiser,Phenolyzer,eXtasy,Phen-Gen,PubCaseFinder --output_table stanford_manton.comparison_output_table.txt --output_image1 stanford_manton.comparison_output_image1.png --output_image2 stanford_manton.comparison_output_image2.png --output_image3 stanford_manton.comparison_output_image3.png

amelie_p_values amelie/amelie/patient_data/amelie_2nd_test_set_EXAC_ALL amelie/amelie/patient_data/stanford_manton.solutions.json --count_filter --hmct_cutoff 1 --dominant_alct_cutoff 3 --rankings_files amelie_real_EXAC_ALL.rankings.txt,amelie_real_EXAC_ALL.exomiser.txt,amelie_real_EXAC_ALL.phenolyzer.txt,amelie_real_EXAC_ALL.extasy.txt,amelie_real_EXAC_ALL.phengen.txt,amelie_real_EXAC_ALL.pubcasefinder.txt --rankings_files_names AMELIE,Exomiser,Phenolyzer,eXtasy,Phen-Gen,PubCaseFinder > stanford_manton_EXAC_ALL.pvalues.txt

amelie_comparisons_info amelie/amelie/patient_data/amelie_2nd_test_set_EXAC_ALL amelie/amelie/patient_data/stanford.solutions.json "Stanford" --count_filter --hmct_cutoff 1 --dominant_alct_cutoff 3 --rankings_files amelie_real_EXAC_ALL.rankings.txt,amelie_real_EXAC_ALL.exomiser.txt,amelie_real_EXAC_ALL.phenolyzer.txt,amelie_real_EXAC_ALL.extasy.txt,amelie_real_EXAC_ALL.phengen.txt,amelie_real_EXAC_ALL.pubcasefinder.txt --rankings_files_names AMELIE,Exomiser,Phenolyzer,eXtasy,Phen-Gen,PubCaseFinder --output_table stanford.comparison_output_table.txt --output_image1 stanford.comparison_output_image1.png --output_image2 stanford.comparison_output_image2.png --output_image3 stanford.comparison_output_image3.png

amelie_comparisons_info amelie/amelie/patient_data/amelie_2nd_test_set_EXAC_ALL amelie/amelie/patient_data/manton.solutions.json "Manton" --count_filter --hmct_cutoff 1 --dominant_alct_cutoff 3 --rankings_files amelie_real_EXAC_ALL.rankings.txt,amelie_real_EXAC_ALL.exomiser.txt,amelie_real_EXAC_ALL.phenolyzer.txt,amelie_real_EXAC_ALL.extasy.txt,amelie_real_EXAC_ALL.phengen.txt,amelie_real_EXAC_ALL.pubcasefinder.txt --rankings_files_names AMELIE,Exomiser,Phenolyzer,eXtasy,Phen-Gen,PubCaseFinder --output_table manton.comparison_output_table.txt --output_image1 manton.comparison_output_image1.png --output_image2 manton.comparison_output_image2.png --output_image3 manton.comparison_output_image3.png

amelie_comparisons_info amelie/amelie/patient_data/amelie_2nd_test_set_EXAC_ALL amelie/amelie/patient_data/udn.solutions.json "UDN" --count_filter --hmct_cutoff 1 --dominant_alct_cutoff 3 --rankings_files amelie_real_EXAC_ALL.rankings.txt,amelie_real_EXAC_ALL.exomiser.txt,amelie_real_EXAC_ALL.phenolyzer.txt,amelie_real_EXAC_ALL.extasy.txt,amelie_real_EXAC_ALL.phengen.txt,amelie_real_EXAC_ALL.pubcasefinder.txt --rankings_files_names AMELIE,Exomiser,Phenolyzer,eXtasy,Phen-Gen,PubCaseFinder --output_table udn.comparison_output_table.txt --output_image1 udn.comparison_output_image1.png --output_image2 udn.comparison_output_image2.png --output_image3 udn.comparison_output_image3.png

amelie_p_values amelie/amelie/patient_data/amelie_2nd_test_set_EXAC_ALL amelie/amelie/patient_data/udn.solutions.json --count_filter --hmct_cutoff 1 --dominant_alct_cutoff 3 --rankings_files amelie_real_EXAC_ALL.rankings.txt,amelie_real_EXAC_ALL.exomiser.txt,amelie_real_EXAC_ALL.phenolyzer.txt,amelie_real_EXAC_ALL.extasy.txt,amelie_real_EXAC_ALL.phengen.txt,amelie_real_EXAC_ALL.pubcasefinder.txt --rankings_files_names AMELIE,Exomiser,Phenolyzer,eXtasy,Phen-Gen,PubCaseFinder > udn_EXAC_ALL.pvalues.txt

amelie_comparisons_info amelie_out_dir/1000g_EXAC_ALL amelie/amelie/patient_data/fake_patients.solutions.json "simulated" --count_filter --hmct_cutoff 1 --dominant_alct_cutoff 3 --rankings_files fake_patients_EXAC_ALL.rankings.txt,fake_patients_EXAC_ALL.exomiser.txt,fake_patients_EXAC_ALL.phenolyzer.txt,fake_patients_EXAC_ALL.extasy.txt,fake_patients_EXAC_ALL.phengen.txt,fake_patients_EXAC_ALL.pubcasefinder.txt --rankings_files_names AMELIE,Exomiser,Phenolyzer,eXtasy,Phen-Gen,PubCaseFinder --output_table fake_patients.comparison_output_table.txt --output_image1 fake_patients.comparison_output_image1.png --output_image2 fake_patients.comparison_output_image2.png --output_image3 fake_patients.comparison_output_image3.png

amelie_p_values amelie_out_dir/1000g_EXAC_ALL amelie/amelie/patient_data/fake_patients.solutions.json --count_filter --hmct_cutoff 1 --dominant_alct_cutoff 3 --rankings_files fake_patients_EXAC_ALL.rankings.txt,fake_patients_EXAC_ALL.exomiser.txt,fake_patients_EXAC_ALL.phenolyzer.txt,fake_patients_EXAC_ALL.extasy.txt,fake_patients_EXAC_ALL.phengen.txt,fake_patients_EXAC_ALL.pubcasefinder.txt --rankings_files_names AMELIE,Exomiser,Phenolyzer,eXtasy,Phen-Gen,PubCaseFinder > fake_patients_EXAC_ALL.pvalues.txt

amelie_convert_rankings amelie/amelie/patient_data/amelie_real_EXAC_ALL amelie/amelie/patient_data/amelie_real.solutions.json --count_filter --hmct_cutoff 1 --dominant_alct_cutoff 3 --exomiser_output_dir amelie/stm_submission_analyses/comparisons/exomiser/results --phenolyzer_output_dir amelie/stm_submission_analyses/comparisons/phenolyzer/results --extasy_output_dir amelie/stm_submission_analyses/comparisons/extasy/results --phengen_output_dir amelie/stm_submission_analyses/comparisons/phen-gen/results --pubcasefinder_output_dir amelie/stm_submission_analyses/comparisons/pubcasefinder/results amelie_real_EXAC_ALL

amelie_convert_rankings amelie_out_dir/1000g_EXAC_ALL amelie/amelie/patient_data/fake_patients.solutions.json --count_filter --hmct_cutoff 1 --dominant_alct_cutoff 3 --exomiser_output_dir amelie/stm_submission_analyses/comparisons/exomiser/results --phenolyzer_output_dir amelie/stm_submission_analyses/comparisons/phenolyzer/results --extasy_output_dir amelie/stm_submission_analyses/comparisons/extasy/results --phengen_output_dir amelie/stm_submission_analyses/comparisons/phen-gen/results --pubcasefinder_output_dir amelie/stm_submission_analyses/comparisons/pubcasefinder/results fake_patients_EXAC_ALL

amelie_p_values amelie/amelie/patient_data/amelie_2nd_test_set_EXAC_ALL amelie/amelie/patient_data/stanford_manton.solutions.json --count_filter --hmct_cutoff 1 --dominant_alct_cutoff 3 --amelie_rankings_file stanford_manton_EXAC_ALL.rankings.txt --exomiser_output_dir amelie/stm_submission_analyses/comparisons/exomiser/results --phenolyzer_output_dir amelie/stm_submission_analyses/comparisons/phenolyzer/results --extasy_output_dir amelie/stm_submission_analyses/comparisons/extasy/results --phengen_output_dir amelie/stm_submission_analyses/comparisons/phen-gen/results --pubcasefinder_output_dir amelie/stm_submission_analyses/comparisons/pubcasefinder/results > stanford_manton_EXAC_ALL.pvalues.txt

amelie_p_values amelie/amelie/patient_data/amelie_2nd_test_set_EXAC_ALL amelie/amelie/patient_data/udn.solutions.json --count_filter --hmct_cutoff 1 --dominant_alct_cutoff 3 --filter_only_variants PASS,. --amelie_rankings_file udn_EXAC_ALL.PASS_ONLY.rankings.txt --exomiser_output_dir amelie/stm_submission_analyses/comparisons/exomiser/results --phenolyzer_output_dir amelie/stm_submission_analyses/comparisons/phenolyzer/results --extasy_output_dir amelie/stm_submission_analyses/comparisons/extasy/results --phengen_output_dir amelie/stm_submission_analyses/comparisons/phen-gen/results --pubcasefinder_output_dir amelie/stm_submission_analyses/comparisons/pubcasefinder/results > udn_EXAC_ALL.PASS_ONLY.pvalues.txt

amelie_p_values amelie_out_dir/1000g_EXAC_ALL amelie/amelie/patient_data/fake_patients.solutions.json --count_filter --hmct_cutoff 1 --dominant_alct_cutoff 3 --amelie_rankings_file fake_patients_EXAC_ALL.rankings.txt --exomiser_output_dir amelie/stm_submission_analyses/comparisons/exomiser/results --phenolyzer_output_dir amelie/stm_submission_analyses/comparisons/phenolyzer/results --extasy_output_dir amelie/stm_submission_analyses/comparisons/extasy/results --phengen_output_dir amelie/stm_submission_analyses/comparisons/phen-gen/results --pubcasefinder_output_dir amelie/stm_submission_analyses/comparisons/pubcasefinder/results > fake_patients_EXAC_ALL.pvalues.txt

# create feature importances:
test_amelie amelie_out_dir amelie/amelie/patient_data/amelie_real_EXAC_ALL amelie/amelie/patient_data/amelie_real.solutions.json --count_filter --hmct_cutoff 1 --dominant_alct_cutoff 3 /dev/shm/amelie_knowledgebase_HPO_2018-07.sqlite3 --rankings_file amelie_real_EXAC_ALL.rankings.txt --patient_features_dir stm_detailed_info --extended_feature_info

# most important features table:
amelie_most_important_features_table amelie/amelie/patient_data/amelie_real.solutions.json stm_detailed_info amelie_real_EXAC_ALL.rankings.txt --output_table most_important_featuresrank_leq1.tsv --output_image most_important_featuresrank_leq1.png --rank_leq 1

cat most_important_featuresrank_leq1.tsv | tail -n+2 | cut -f 2 | sort | uniq -c | sort -k1,1n
cat most_important_featuresrank_leq1.tsv | cut -f 2,3,4 | tail -n+2 | tr '\t' '\n' | sort | uniq -c | sort -k1,1n | wc -l

# xref phenos extraction:
extract_all_xref_phenos amelie_knowledgebase_HPO_2018-07.sqlite3 --conda_path /cluster/u/jbirgmei/miniconda3/bin --max_jobs 150 ${PWD}/amelie_out_dir ${PWD}/amelie_papers ${PWD}/amelie_process_dir_xref_phenos
amelie_add_xref_phenos_to_db amelie_knowledgebase_HPO_2018-07.sqlite3 amelie_out_dir amelie_knowledgebase_xref_phenos.sqlite3
sqlite3 amelie_knowledgebase_HPO_2018-07.sqlite3 "select avg(count) from (select p.pmid, count(distinct phenotype_id) count from amelie_paper p join amelie_paper_phenotypes pph on pph.paper_id = p.id group by p.pmid) a"
sqlite3 amelie_knowledgebase_xref_phenos.sqlite3 "select avg(count) from (select p.pmid, count(distinct phenotype_id) count from amelie_paper p join amelie_paper_phenotypes pph on pph.paper_id = p.id group by p.pmid) a"
train_amelie --test_patients_dir amelie/amelie/patient_data/amelie_real_EXAC_ALL --test_solutions_list amelie/amelie/patient_data/amelie_real.solutions.json amelie_knowledgebase_xref_phenos.sqlite3 amelie_out_dir amelie_out_dir/1000g_EXAC_ALL amelie_out_dir/fake_patients_1000g_EXAC_ALL_fake_phenotypes_post_2011.tsv.json --count_filter --hmct_cutoff 1 --dominant_alct_cutoff 3 --rankings_file amelie_real_EXAC_ALL.xref_phenos_rankings.txt --ignore_if_no_causative_gene --do_not_save
amelie_p_values amelie/amelie/patient_data/amelie_real_EXAC_ALL amelie/amelie/patient_data/amelie_real.solutions.json --count_filter --hmct_cutoff 1 --dominant_alct_cutoff 3 --rankings_files amelie_real_EXAC_ALL.rankings.txt,amelie_real_EXAC_ALL.xref_phenos_rankings.txt --rankings_files_names 'AMELIE,XREF' --nonscientific_format
amelie_pheno_precision_table amelie_knowledgebase_xref_phenos.sqlite3 amelie_process_dir_xref_phenos > xref_pheno_precision_table.tsv

# precision tables:
amelie_gene_precision_table amelie_knowledgebase_HPO_2018-07.sqlite3 amelie_process_dir > gene_precision_table.tsv
amelie_pheno_precision_table amelie_knowledgebase_HPO_2018-07.sqlite3 amelie_process_dir > pheno_precision_table.tsv
amelie_domrec_precision_table amelie_knowledgebase_HPO_2018-07.sqlite3 > domrec_precision_table.tsv

# recall of UniProt and HGNC gene recognition vs. everything (including PubTator):
amelie_get_num_pdf_pages amelie_papers extracted_uniprot_hgnc_genes.pmids.tsv > pmid_to_num_pages.tsv
amelie_create_all_extracted_genes_file amelie_out_dir/extracted_uniprot_hgnc_genes.tsv amelie_process_dir > all_extracted_genes.tsv
amelie_uniprot_hgnc_genes_hgmd_recall amelie_out_dir/extracted_uniprot_hgnc_genes.tsv all_extracted_genes.tsv --pmid_to_num_pdf_pages pmid_to_num_pages.tsv > amelie_genes_hgmd_recall.txt
amelie_uniprot_hgnc_genes_hgmd_recall amelie_out_dir/extracted_uniprot_hgnc_genes.tsv all_extracted_genes.tsv --pmid_to_num_pdf_pages pmid_to_num_pages.tsv --min_year 2013 --pmid_to_year pmid_to_year.tsv > amelie_genes_hgmd_recall_2013.txt

# get number of articles AMELIE analyzed per patient (ATTENTION: THIS IS FOR ALL, SUBSET IF NECESSARY):
amelie_num_articles_per_patient amelie/amelie/patient_data/amelie_real.solutions.json stm_detailed_info

# get number of candidate genes per patient:
amelie_num_candidate_genes_per_patient amelie/amelie/patient_data/amelie_real_EXAC_ALL amelie/amelie/patient_data/amelie_real.solutions.json --count_filter --hmct_cutoff 1 --dominant_alct_cutoff 3 --no_mcap
amelie_num_candidate_variants_per_patient amelie/amelie/patient_data/amelie_real_EXAC_ALL amelie/amelie/patient_data/amelie_real.solutions.json --count_filter --hmct_cutoff 1 --dominant_alct_cutoff 3 --no_mcap

amelie_num_candidate_genes_per_patient amelie/amelie/patient_data/amelie_2nd_test_set_EXAC_ALL amelie/amelie/patient_data/stanford_manton.solutions.json --count_filter --hmct_cutoff 1 --dominant_alct_cutoff 3 --no_mcap | datamash median 2

amelie_num_candidate_genes_per_patient amelie/amelie/patient_data/amelie_2nd_test_set_EXAC_ALL amelie/amelie/patient_data/udn.solutions.json --count_filter --hmct_cutoff 1 --dominant_alct_cutoff 3 --filter_only_variants PASS,. --no_mcap | datamash median 2

# number classified as relevant etc:
cat out_* | grep 'for relevance, ' > /cluster/u/jbirgmei/shallowfloat/amelie_ethan/tested_for_relevance.txt
cat tested_for_relevance.txt | cut -d ' ' -f 5,9 | tr ' ' '\t' | sort -u > tested_for_relevance.tsv
grep -w True tested_for_relevance.tsv | sort -u > relevant_pmids.txt
wc -l relevant_pmids.txt
cat relevant_pmids.txt | cut -f 1 | while read i
do
  ll amelie_papers/${i}.pdf | awk '$5 > 0'
done > downloaded_pmids.txt
cut -d '/' -f 2 downloaded_pmids.txt | cut -d '.' -f 1 | sort -u > downloaded_pmids.tsv
# download success rate:
join -t $'\t' -1 1 -2 1 relevant_pmids.txt downloaded_pmids.tsv | cut -f 1 | sort -u | wc -l

amelie_p_values amelie_out_dir/1000g_EXAC_ALL amelie/amelie/patient_data/fake_patients.solutions.json --count_filter --hmct_cutoff 1 --dominant_alct_cutoff 3 --rankings_files fake_patients_EXAC_ALL.no_phrank_interpolation.rankings.txt,fake_patients_EXAC_ALL.no_phrank_interpolation.no_allele_counts.rankings.txt,fake_patients_EXAC_ALL.no_phrank_interpolation.no_domrec.rankings.txt,fake_patients_EXAC_ALL.no_phrank_interpolation.no_evs.rankings.txt,fake_patients_EXAC_ALL.no_phrank_interpolation.no_mcap.rankings.txt,fake_patients_EXAC_ALL.no_phrank_interpolation.no_phrank.rankings.txt,fake_patients_EXAC_ALL.no_phrank_interpolation.no_right_gene.rankings.txt,fake_patients_EXAC_ALL.no_phrank_interpolation.no_rvis_pli.rankings.txt,fake_patients_EXAC_ALL.no_phrank_interpolation.no_strict_score.rankings.txt,fake_patients_EXAC_ALL.no_phrank_interpolation.no_vdestructive_het.rankings.txt --rankings_files_names "AMELIE,ExAC,Dom/rec,AVADA,M-CAP,Phenos,Topical,RVIS/pLI,Relevance,VInfo" --nonscientific_format
amelie_p_values amelie/amelie/patient_data/amelie_real_EXAC_ALL amelie/amelie/patient_data/amelie_real.solutions.json --count_filter --hmct_cutoff 1 --dominant_alct_cutoff 3 --rankings_files amelie_real_EXAC_ALL.no_phrank_interpolation.rankings.txt,amelie_real_EXAC_ALL.no_phrank_interpolation.no_allele_counts.rankings.txt,amelie_real_EXAC_ALL.no_phrank_interpolation.no_domrec.rankings.txt,amelie_real_EXAC_ALL.no_phrank_interpolation.no_evs.rankings.txt,amelie_real_EXAC_ALL.no_phrank_interpolation.no_mcap.rankings.txt,amelie_real_EXAC_ALL.no_phrank_interpolation.no_phrank.rankings.txt,amelie_real_EXAC_ALL.no_phrank_interpolation.no_right_gene.rankings.txt,amelie_real_EXAC_ALL.no_phrank_interpolation.no_rvis_pli.rankings.txt,amelie_real_EXAC_ALL.no_phrank_interpolation.no_strict_score.rankings.txt,amelie_real_EXAC_ALL.no_phrank_interpolation.no_vdestructive_het.rankings.txt --rankings_files_names "AMELIE,ExAC,Dom/rec,AVADA,M-CAP,Phenos,Topical,RVIS/pLI,Relevance,VInfo" --nonscientific_format

# ta_only comparisons:
amelie_comparisons_info --rankings_files amelie_real_EXAC_ALL.rankings.txt,amelie_real_EXAC_ALL.ta_only_rankings.txt,amelie_real_EXAC_ALL.exomiser.txt,amelie_real_EXAC_ALL.extasy.txt,amelie_real_EXAC_ALL.phengen.txt,amelie_real_EXAC_ALL.phenolyzer.txt,amelie_real_EXAC_ALL.pubcasefinder.txt --rankings_files_names "AMELIE,AMELIE from title/abstract only,Exomiser,eXtasy,Phen-Gen,Phenolyzer,PubCaseFinder" amelie/amelie/patient_data/amelie_real_EXAC_ALL amelie/amelie/patient_data/amelie_real.solutions.json --count_filter --hmct_cutoff 1 --dominant_alct_cutoff 3 --output_table amelie_ta_only_comparison_output_table.txt --output_image1 amelie_ta_only_comparison_output_image1.png --output_image2 amelie_ta_only_comparison_output_image2.png --output_image3 amelie_ta_only_comparison_output_image3.png "real"
# amelie_p_values --rankings_files amelie_real_EXAC_ALL.ta_only_rankings.txt,amelie_real_EXAC_ALL.exomiser.txt,amelie_real_EXAC_ALL.extasy.txt,amelie_real_EXAC_ALL.phengen.txt,amelie_real_EXAC_ALL.phenolyzer.txt,amelie_real_EXAC_ALL.pubcasefinder.txt --rankings_files_names "AMELIE from title/abstract only,Exomiser,eXtasy,Phen-Gen,Phenolyzer,PubCaseFinder" amelie/amelie/patient_data/amelie_real_EXAC_ALL amelie/amelie/patient_data/amelie_real.solutions.json --count_filter --hmct_cutoff 1 --dominant_alct_cutoff 3 > amelie_ta_only.pvalues.txt
amelie_p_values --rankings_files amelie_real_EXAC_ALL.rankings.txt,amelie_real_EXAC_ALL.ta_only_rankings.txt --rankings_files_names "AMELIE,AMELIE from title/abstract only" amelie/amelie/patient_data/amelie_real_EXAC_ALL amelie/amelie/patient_data/amelie_real.solutions.json --count_filter --hmct_cutoff 1 --dominant_alct_cutoff 3 > amelie_ta_only.pvalues.txt

#### NO PHENOTYPES ANALYSES
cp amelie_knowledgebase_HPO_2018-07.sqlite3 amelie_knowledgebase_WITHOUT_TOP_PAPER_PHENOS.sqlite3
for i in stm_detailed_info/*top_ranked*
do
  cat $i | head -n 1 | cut -f 1
done | sort -u | while read j
do
  sqlite3 amelie_knowledgebase_WITHOUT_TOP_PAPER_PHENOS.sqlite3 "delete from amelie_paper_phenotypes where paper_id in (select id from amelie_paper where pmid = '$j')"
done

test_amelie amelie_out_dir amelie/amelie/patient_data/amelie_real_EXAC_ALL amelie/amelie/patient_data/amelie_real.solutions.json --count_filter --hmct_cutoff 1 --dominant_alct_cutoff 3 amelie_knowledgebase_WITHOUT_TOP_PAPER_PHENOS.sqlite3 --rankings_file amelie_real_EXAC_ALL.drop_out_1paper_phenos_rankings.txt --drop_out_top_paper_phenos 1 --drop_out_ranking_info_dir stm_detailed_info

amelie_p_values amelie/amelie/patient_data/amelie_real_EXAC_ALL amelie/amelie/patient_data/amelie_real.solutions.json --count_filter --hmct_cutoff 1 --dominant_alct_cutoff 3 --rankings_files amelie_real_EXAC_ALL.drop_out_1paper_phenos_rankings.txt,amelie_real_EXAC_ALL.exomiser.txt,amelie_real_EXAC_ALL.extasy.txt,amelie_real_EXAC_ALL.phengen.txt,amelie_real_EXAC_ALL.phenolyzer.txt,amelie_real_EXAC_ALL.pubcasefinder.txt --rankings_files_names 'AMELIE without top paper,Exomiser,eXtasy,Phen-Gen,Phenolyzer,PubCaseFinder' > amelie_real_EXAC_ALL.drop_out_1paper.pvalues.txt

# replacing phenos with nothing:
train_classifiers --forbidden_genesymbols_file ega215_forbidden_genesymbols.tsv amelie_out_dir_replace_phenos_with_nothing amelie_process_dir --replace_phenos_with_nothing

extract_all_papers --pubmunch_data_dir /cluster/u/jbirgmei/shallowfloat/amelie_ethan/pubmunchData --conda_path /cluster/u/jbirgmei/miniconda3/bin --max_jobs 50 --elsevier_key 9375813031e16ecaa64295d9b45976e6 ${PWD}/amelie_out_dir_replace_phenos_with_nothing ${PWD}/amelie_papers ${PWD}/amelie_process_dir --replace_phenos_with_nothing --max_update_number 1280 --no_pubtator_download

cat /cluster/tmp/tmp5h94sk1q/out* | grep 'Tested title/abstract of' | cut -d ' ' -f 5,9 | tr ' ' '\t' > tested_for_relevance.NO_PHENOS.tsv

load_db_all_papers --create_db amelie_out_dir_replace_phenos_with_nothing /dev/shm/amelie_knowledgebase_no_phenos.sqlite3
cp /dev/shm/amelie_knowledgebase_no_phenos.sqlite3 amelie_knowledgebase_no_phenos.sqlite3

train_amelie --test_patients_dir amelie/amelie/patient_data/amelie_real_EXAC_ALL --test_solutions_list amelie/amelie/patient_data/amelie_real.solutions.json /dev/shm/amelie_knowledgebase_no_phenos.sqlite3 amelie_out_dir amelie_out_dir/1000g_EXAC_ALL amelie_out_dir/fake_patients_1000g_EXAC_ALL_fake_phenotypes_post_2011.tsv.json --count_filter --hmct_cutoff 1 --dominant_alct_cutoff 3 --no_phrank_interpolation --rankings_file amelie_real_EXAC_ALL.no_phrank.rankings.txt --feature_names_list amelie/amelie/amelie/classifier_data/amelie_features_no_phrank.txt --classifier_filename amelie_out_dir/amelie_classifier_no_phrank.pkl --ignore_if_no_causative_gene --no_phrank_interpolation

# TODO insert no phenotypes amelie training here

# no phenotypes figure:
amelie_comparisons_info amelie/amelie/patient_data/amelie_real_EXAC_ALL amelie/amelie/patient_data/amelie_real.solutions.json "real" --count_filter --hmct_cutoff 1 --dominant_alct_cutoff 3 --rankings_files amelie_real_EXAC_ALL.no_phrank.rankings.txt,amelie_real_EXAC_ALL.exomiser.txt,amelie_real_EXAC_ALL.phenolyzer.txt,amelie_real_EXAC_ALL.extasy.txt,amelie_real_EXAC_ALL.phengen.txt,amelie_real_EXAC_ALL.pubcasefinder.txt --rankings_files_names "AMELIE-all phenos,Exomiser,Phenolyzer,eXtasy,Phen-Gen,PubCaseFinder" --output_table all_real_no_phrank.comparison_output_table.txt --output_image1 all_real_no_phrank.comparison_output_image1.png --output_image2 all_real_no_phrank.comparison_output_image2.png --output_image3 all_real_no_phrank.comparison_output_image3.png

amelie_p_values amelie/amelie/patient_data/amelie_real_EXAC_ALL amelie/amelie/patient_data/amelie_real.solutions.json --count_filter --hmct_cutoff 1 --dominant_alct_cutoff 3 --rankings_files amelie_real_EXAC_ALL.no_phrank.rankings.txt,amelie_real_EXAC_ALL.exomiser.txt,amelie_real_EXAC_ALL.phenolyzer.txt,amelie_real_EXAC_ALL.extasy.txt,amelie_real_EXAC_ALL.phengen.txt,amelie_real_EXAC_ALL.pubcasefinder.txt --rankings_files_names "AMELIE-all phenos,Exomiser,Phenolyzer,eXtasy,Phen-Gen,PubCaseFinder" > all_real_no_phrank.pvalues.txt

# phrank-only ranking:
amelie_phrank_only amelie_out_dir amelie/amelie/patient_data/amelie_real_EXAC_ALL amelie/amelie/patient_data/amelie_real.solutions.json --count_filter --hmct_cutoff 1 --dominant_alct_cutoff 3 /dev/shm/amelie_knowledgebase_HPO_2018-07.sqlite3 --rankings_file amelie_real_EXAC_ALL.phrank_only_rankings.txt --features_dict_dir features_dict_dir

amelie_comparisons_info amelie/amelie/patient_data/amelie_real_EXAC_ALL amelie/amelie/patient_data/amelie_real.solutions.json "real" --count_filter --hmct_cutoff 1 --dominant_alct_cutoff 3 --rankings_files amelie_real_EXAC_ALL.rankings.txt,amelie_real_EXAC_ALL.phrank_only_rankings.txt --rankings_files_names "AMELIE,Phrank" --output_table only_phrank.comparison_output_table.txt --output_image1 only_phrank.comparison_output_image1.png --output_image2 only_phrank.comparison_output_image2.png --output_image3 only_phrank.comparison_output_image3.png --auto_color

amelie_p_values amelie/amelie/patient_data/amelie_real_EXAC_ALL amelie/amelie/patient_data/amelie_real.solutions.json --count_filter --hmct_cutoff 1 --dominant_alct_cutoff 3 --rankings_files amelie_real_EXAC_ALL.rankings.txt,amelie_real_EXAC_ALL.phrank_only_rankings.txt --rankings_files_names "AMELIE,Phrank"

# DDD papers cited in OMIM:
cat amelie_real_EXAC_ALL.rankings.txt | grep DDD | awk '$2 == 1' | cut -f 1 | while read i
do
head -n 1 stm_detailed_info/${i}.top_ranked_papers.txt | cut -f 1
done | sort -u > ddd_top_ranked_papers.tsv

cat ddd_top_ranked_papers.tsv | while read i
do
sqlite3 amelie_knowledgebase_HPO_2018-07.sqlite3 "select pmid, title from amelie_paper where pmid = '$i'"
done | tr '|' '\t' > omim_titles_table.tsv

# growth of literature picture:
amelie_growth_of_literature amelie_knowledgebase_HPO_2018-07.sqlite3 growth_of_literature.png --min_year 1980 --max_year 2017

# num genes per article histogram:
sqlite3 /dev/shm/amelie_knowledgebase_HPO_2018-07.sqlite3 'select pmid, count(*) from (select distinct pmid, gene_symbol from amelie_paper p join amelie_papergene pg on pg.paper_id = p.id join amelie_gene g on g.id = pg.gene_id) group by pmid' | tr '|' '\t' | sort -t $'\t' -k2,2nr > num_genes_histogram.tsv

amelie_relevant_genes_histogram num_genes_histogram.tsv num_genes_histogram.png
cat num_genes_histogram.tsv | datamash min 2 max 2 mean 2 median 2 pstdev 2 iqr 2


##### for website build, no forbidden genes, GNOMAD_MAX:
create_fake_patients --cohort 1000g_GNOMAD_MAX --real_solutions amelie/amelie/patient_data/empty.solutions.json --clinvar_exac_or_gnomad gnomad_max amelie/amelie/patient_data/fake_phenotypes_post_2011.tsv amelie_out_dir_GNOMAD_MAX
train_amelie /dev/shm/amelie_knowledgebase_HPO_2018-07.sqlite3 amelie_out_dir_GNOMAD_MAX amelie_out_dir_GNOMAD_MAX/1000g_GNOMAD_MAX amelie_out_dir_GNOMAD_MAX/fake_patients_1000g_GNOMAD_MAX_fake_phenotypes_post_2011.tsv.json --ignore_if_no_causative_gene --test_patients_dir amelie/amelie/patient_data/amelie_real_GNOMAD_MAX --test_solutions_list amelie/amelie/patient_data/amelie_real.solutions.json --rankings_file amelie_real_GNOMAD_MAX.no_count_filter.rankings.txt --output_if_causative_gene_ranked

# use --do_not_save for don't save
train_amelie /dev/shm/amelie_knowledgebase_2011_tfidf0195.sqlite3 amelie_out_dir_2011 amelie_out_dir_2011/1000g_GNOMAD_MAX amelie_out_dir_2011/fake_patients_1000g_GNOMAD_MAX_fake_phenotypes_pre_2011_post_2006.tsv.json --ignore_if_no_causative_gene --test_patients_dir amelie/amelie/patient_data/amelie_real_GNOMAD_MAX --test_solutions_list amelie/amelie/patient_data/amelie_real.solutions.json --rankings_file amelie_real_GNOMAD_MAX.rankings.2011.txt --output_if_causative_gene_ranked > train_amelie_classifier_2011.out.txt

test_amelie amelie_out_dir_GNOMAD_MAX amelie/amelie/patient_data/amelie_real_GNOMAD_MAX amelie/amelie/patient_data/amelie_real.solutions.json --count_filter --hmct_cutoff 1 --dominant_alct_cutoff 3 /dev/shm/amelie_knowledgebase_HPO_2018-07.sqlite3 --rankings_file amelie_real_GNOMAD_MAX.no_count_filter.rankings.txt --output_if_causative_gene_ranked --ignore_if_no_causative_gene

# gnomad max, retrain the whole BS without forbidden anything for website, including TAs in place of missing full text, tfidf 0.01-0.95 
train_classifiers amelie_out_dir_GNOMAD_MAX amelie_process_dir amelie_ta_genephenos

# journals with most gene-pheno relationships:
sqlite3 amelie_knowledgebase_HPO_2018-07.sqlite3 'select distinct journal, gene_symbol, hpo_id from amelie_gene g join amelie_papergene pg on pg.gene_id = g.id join amelie_paper p on p.id = pg.paper_id join amelie_paper_phenotypes pph on pph.paper_id = p.id join amelie_phenotype ph on pph.phenotype_id = ph.id' | tr '|' '\t' > journal_to_gene_to_pheno.tsv
sort -t $'\t' -k1,1 journal_to_gene_to_pheno.tsv -o journal_to_gene_to_pheno.tsv
cat journal_to_gene_to_pheno.tsv | datamash -g 1 count 2,3 > highest_gp_journals.tsv
sort -t $'\t' -k2,2nr highest_gp_journals.tsv -o highest_gp_journals.tsv
head highest_gp_journals.tsv
./plot_highest_gp_journals.py AMELIE_v77/highest_gp_journals.tsv AMELIE_v77/highest_gp_journals.png

# reanalysis all pubmed info db:
mkdir all_pubmed_info_dir
all_pubmed_info_db ${PWD}/all_pubmed_info.sqlite3 ${PWD}/all_pubmed_info_dir --conda_path /cluster/u/jbirgmei/miniconda3/bin
cat all_pubmed_info_dir/*.tsv | sort -t $'\t' -k1,1 --parallel 20 -S 64G > all_pubmed_info.tsv
sqlite3 /dev/shm/all_pubmed_info.sqlite3 'create table data (pmid integer primary key, title text, date text, journal text, authors text)'
# in sqlite3 /dev/shm/all_pubmed_info.sqlite3:
cat all_pubmed_info.tsv | sed 's/"//g' > /dev/shm/all_pubmed_info.tsv
.separator "\t"
.import /dev/shm/all_pubmed_info.tsv data
mv /dev/shm/all_pubmed_info.sqlite3 .

./create_patient_characteristics_table.py simulated_patients_table.tsv amelie_out_dir/fake_patients_1000g_EXAC_ALL_fake_phenotypes_post_2011.tsv.json > simulated_patient_high_level_characteristics.tsv
