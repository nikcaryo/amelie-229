
from setuptools import setup, find_packages
# To use a consistent encoding
from codecs import open
from os import path
import setuptools.command.build_ext
import subprocess
import shutil
import sys

here = path.abspath(path.dirname(__file__))

# class BuildExtCommand(setuptools.command.build_ext.build_ext):
#   """Custom build command."""

#   def run(self):
#     result = subprocess.run(['mvn', '-f', 'PDFTextToPos', 'compile', 'assembly:single'])
#     if result.returncode != 0:
#         print('Unable to compile PDFTextToPos', result)
#         sys.exit(1)
#     shutil.copyfile('PDFTextToPos/target/PDFTextToPos-1.0-jar-with-dependencies.jar', 'amelie/avada_data/PDFTextToPos.jar')
#     setuptools.command.build_ext.build_ext.run(self)

setup(
    name='amelie',
    version='1.0.0',
    description='Automatic Mendelian Literature Evaluation',

    # cmdclass={
    #     'build_ext': BuildExtCommand,
    # },

    packages=find_packages(exclude=['contrib', 'docs', 'tests', 'stm_submission_analyses']),

    install_requires=[
        'pytest==3.5.0',
        'pubMunch3==1.0.3',
        'nltk==3.2.5',
        'scikit-learn==0.20.0',
        'scipy==1.1.0',
        'numpy==1.15.2',
        'constdb==3.3.0',
        'pyfaidx==0.5.3.1',
        'sqlalchemy==1.2.0',
        'pytabix==0.1',
        'pandas==0.23.4',
        'matplotlib==3.0.0',
        'PyPDF2==1.26.0'
    ],

    zip_safe=False,

    package_data={
        'amelie': ['avada_data/PDFTextToPos.jar',
                   'pheno_data/hpo_ontology.obo',
                   'avada_data/clinvar20160705_pathogenic_variants_realigned_with_pmids.tsv',
                   'avada_data/PDFTextToPos.jar',
                   'variant_extractor/variant_data/regex.txt',
                   'gwas_data/gwas_catalog_v1.0-studies_r2018-04-10.tsv',
                   'gene_data/ensembl_genes_with_names.tsv',
                   'gene_data/gene.entrez.map',
                   'sql_schemas/knowledgebase.sql',
                   'gene_data/ensg_to_rvis_to_mcap_to_pli.tsv',
                   'gene_data/amelie_gene.csv',
                   'gene_data/nm_to_ensembl_to_np.tsv',
                   ]
    },

    entry_points={
        'console_scripts': [
            'download_data=amelie.download_main:download_data_program',
            'download_all_papers=amelie.download_main:download_all_papers_program',
            'process_all_papers=amelie.process_main:process_all_papers_program',

            'train_classifiers=amelie.train_main:train_classifiers_program',
            'train_full_text_relevance = amelie.train_main:train_full_text_relevance',
            'train_inheritance_mode = amelie.train_main:train_inheritance_mode',
            'train_title_abstract = amelie.train_main:train_title_abstract',
            'train_topical_gene = amelie.train_main:train_topical_gene',
            'train_variants = amelie.train_main:train_variants',
            'train_variant_type = amelie.train_main:train_variant_type',

            'extract_all_papers=amelie.extract_main:extract_all_papers_program',
            'extract_all_papers_single_thread=amelie.extract_main_single_thread:extract_all_papers_program',
            'load_db_all_papers=amelie.load_db_main:load_all_papers_program',
            'load_disgenet_db=amelie.stm_analyses.load_disgenet_db:load_disgenet_program',
            'amelie_extraction_report=amelie.extraction_report:extraction_report',

            'download_papers=amelie.download_main:download_papers_program',
            'process_papers=amelie.process_main:process_papers_program',
            'extract_papers=amelie.extract_main:extract_papers_program',
            'extract_papers_single_thread=amelie.extract_main_single_thread:extract_papers_program',

            'annotate_vcfs=amelie.amelie.annotate_vcf_main:annotate_vcfs_program',
            'annotate_vcf=amelie.amelie.annotate_vcf_main:annotate_vcf_program',
            'create_fake_patients=amelie.amelie.fake_patients_main:create_fake_patients_program',
            'train_amelie=amelie.amelie.amelie_classifier_main:train_amelie_classifier_program',
            'test_amelie=amelie.amelie.amelie_classifier_main:test_amelie_classifier_program',
            'amelie_reanalysis=amelie.amelie.reanalysis_main:reanalysis_program',

            'pmid_to_year=amelie.amelie.pmid_to_date_main:pmid_to_date_program',
            'amelie_filter_variants=amelie.amelie.variant_filtering:filter_variants_program',

            'canonicalize_phenotypes=amelie.stm_analyses.canonicalize_phenotypes_main:canonicalize_phenotypes_program',
            'amelie_simulated_patients_table=amelie.stm_analyses.create_simulated_patients_table_main:create_simulated_patients_table',
            'amelie_comparisons_info=amelie.stm_analyses.comparisons_info_main:comparisons_info',
            'amelie_most_important_features_table=amelie.stm_analyses.most_important_features_table_main:most_important_features_table',
            'amelie_gene_precision_table=amelie.stm_analyses.manual_eval_tables:gene_precision_table',
            'amelie_pheno_precision_table=amelie.stm_analyses.manual_eval_tables:pheno_precision_table',
            'amelie_domrec_precision_table=amelie.stm_analyses.manual_eval_tables:domrec_precision_table',

            'extract_all_xref_phenos=amelie.extract_xref_phenos_main:extract_all_xref_phenos_program',
            'extract_xref_phenos=amelie.extract_xref_phenos_main:extract_xref_phenos_program',
            'extract_all_uniprot_hgnc_genes=amelie.extract_uniprot_hgnc_genes_main:extract_all_uniprot_hgnc_genes_program',
            'extract_uniprot_hgnc_genes=amelie.extract_uniprot_hgnc_genes_main:extract_uniprot_hgnc_genes_program',

            'amelie_uniprot_hgnc_genes_hgmd_recall=amelie.stm_analyses.uniprot_hgnc_genes_hgmd_recall:compare_hgmd_genes_with_uniprot_hgnc_extracted_program',
            'amelie_get_num_pdf_pages=amelie.stm_analyses.get_num_pdf_pages:get_num_pdf_pages',
            'amelie_add_xref_phenos_to_db=amelie.stm_analyses.add_xref_phenos_to_db:add_xref_phenos_to_db',
            'amelie_create_all_extracted_genes_file=amelie.stm_analyses.create_all_extracted_genes_file:create_all_extracted_genes_file',
            'amelie_p_values=amelie.stm_analyses.p_values_main:p_values_program',
            'amelie_num_articles_per_patient=amelie.stm_analyses.analyzed_articles_analysis:num_articles_per_patient',
            'amelie_num_articles_per_causative_gene=amelie.stm_analyses.analyzed_articles_analysis:num_articles_per_causative_gene',
            'amelie_num_candidate_genes_per_patient=amelie.stm_analyses.num_candidate_gene_stats:num_candidate_genes_per_patient',
            'amelie_num_candidate_variants_per_patient=amelie.stm_analyses.num_candidate_gene_stats:num_candidate_variants_per_patient',
            'amelie_convert_rankings=amelie.stm_analyses.comparisons_info_main:convert_rankings',
            'amelie_rank_vs_num_articles=amelie.stm_analyses.rank_vs_num_articles:rank_vs_num_articles_program',
            'amelie_num_pheno_terms=amelie.pheno_extractor:num_pheno_terms',
            'amelie_feature_ablation_table=amelie.stm_analyses.feature_ablation_table:feature_ablation_table',
            'amelie_phrank_only=amelie.amelie.amelie_classifier_main:test_phrank_only',
            'amelie_phenotypes_histogram=amelie.stm_analyses.phenotypes_histogram:draw_phenotypes_histogram',
            'amelie_feature_ablation_barchart=amelie.stm_analyses.feature_ablation_barchart:feature_ablation_barchart',
            'amelie_growth_of_literature=amelie.stm_analyses.growth_of_literature:growth_of_literature',
            'amelie_anonymize_column=amelie.stm_analyses.anonymize_column:anonymize_column',
            'amelie_relevant_genes_histogram=amelie.stm_analyses.relevant_genes_histogram:relevant_genes_histogram',
            'amelie_allowed_genes_list=amelie.amelie.variant_filtering:allowed_genes_list',

            # reanalysis:
            'all_pubmed_info_db=amelie.all_pubmed_info_db_main:all_pubmed_info_db_program',
            'pubmed_info_db=amelie.all_pubmed_info_db_main:pubmed_info_db_program',
        ],
    },
)
