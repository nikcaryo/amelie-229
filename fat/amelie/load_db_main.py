import argparse
import sys, os
import json
import glob
import re
import datetime
import sqlite3
from collections import defaultdict

import pkg_resources

from .dataloaders import hpo
from amelie.sqlite_util import dotdict_factory


def get_db_loaded_updates(connection):
    cur = connection.cursor()
    results = cur.execute("SELECT pubmed_update, pubmed_update_index FROM amelie_pubmed_update")
    rv = []
    for x in results:
        rv.append(json.dumps({"update": x["pubmed_update"],
                              "update_index": x["pubmed_update_index"]}))
    return rv


def parse_update_from_filename(filename):
    filename = os.path.basename(filename)
    pattern = re.compile('(\w+)-(\w+)-(\d+).json')
    groups = pattern.match(filename)
    return json.dumps({"update": groups[1] + '-' + groups[2], "update_index": int(groups[3])})


def get_extract_papers_updates(args):
    extracted_dir = os.path.join(args.out_dir, args.out_extracteddir)
    extracted_files = glob.glob(extracted_dir + "/*-*.json")
    rv = []
    for f in extracted_files:
        rv.append(parse_update_from_filename(f))
    return rv


def get_simple_inheritance_mode(inheritance_mode_value):
    if inheritance_mode_value['dominant'] > 0.5:
        return 'DOMINANT'
    if inheritance_mode_value['recessive'] > 0.5:
        return 'RECESSIVE'
    return 'NONE'


def parse_publish_date(publish_date, pmid):
    year, month = publish_date_year_month(publish_date, pmid)
    return "%d %02d" % (year, month)


def publish_date_year_month(publish_date, pmid):
    try:
        year, month, day = [int(x) for x in publish_date.split('-')]
    except (ValueError, AttributeError):
        print("Publish date %s invalid for pmid %s, setting to today" % (publish_date, pmid))
        now = datetime.datetime.now()
        current_year, current_month = now.year, now.month
        return current_year, current_month
    return year, month


def get_gene_ids(eids, cur):
    rv = []
    for eid in eids:
        if eid not in get_gene_ids.cache:
            results = cur.execute("select id from amelie_gene where ensembl_id = ?", [eid])
            gene_db_ids = []
            for result in results:
                gene_db_ids.append(result['id'])
            get_gene_ids.cache[eid] = gene_db_ids
        rv.extend(get_gene_ids.cache[eid])
    return rv
get_gene_ids.cache = {}


def get_phenotype_id(hpo_id, cur):
    if hpo_id in get_phenotype_id.cache:
        return get_phenotype_id.cache[hpo_id]
    results = cur.execute("select id from amelie_phenotype where hpo_id = ?", [hpo_id])
    for result in results:
        db_id = result['id']
        get_phenotype_id.cache[hpo_id] = db_id
        return db_id
    get_phenotype_id.cache[hpo_id] = None
    return None
get_phenotype_id.cache = {}


def get_db_loaded_pmids(connection):
    rv = set()
    cur = connection.cursor()
    results = cur.execute("SELECT DISTINCT pmid FROM amelie_paper")
    for result in results:
        rv.add(result['pmid'])
    return rv


def get_next_id(table_name, cur):
    results = [x for x in cur.execute('SELECT MAX(id)+1 as max FROM %s' % table_name)]
    result = results[0]["max"]
    if result is None:
        return 0
    return result


def process_update(update, args, connection, previously_loaded_pmids,
                   ensembl_to_gene_ids, overwrite=False):
    update = json.loads(update)
    update_filename = os.path.join(args.out_dir, args.out_extracteddir,
                                   '%s-%d.json' % (update['update'], update['update_index']))
    print("Processing update %s" % json.dumps(update))
    cur = connection.cursor()
    cur.execute("PRAGMA synchronous = OFF")
    cur.execute("PRAGMA locking_mode = EXCLUSIVE")
    cur.execute("PRAGMA journal_mode = MEMORY")
    cur.execute("BEGIN TRANSACTION")
    with open(update_filename) as update_file:
        update_dict = json.load(update_file)
        next_paper_id = get_next_id('amelie_paper', cur=cur)
        try:
            for pmid, value in update_dict.items():
                if not value["is_valid"]:
                    continue
                year, month = publish_date_year_month(value['publish_date'], pmid=pmid)
                if args.max_year is not None and args.max_month is not None:
                    if (year, month) > (args.max_year, args.max_month):
                        print("Skipping PubMed ID; publish date %d-%02d is later than given max %d-%02d" %
                              (year, month, args.max_year, args.max_month))
                if pmid in previously_loaded_pmids:
                    print("Pubmed ID %s already loaded, skipping" % pmid)
                    continue
                if overwrite:
                    amelie_paper_id = [x for x in cur.execute('SELECT id FROM amelie_paper WHERE pmid = ?',
                                                               [pmid])][0]["id"]
                    if amelie_paper_id is None:
                        amelie_paper_id = next_paper_id
                        next_paper_id += 1
                    cur.execute('DELETE FROM amelie_paper WHERE pmid = ? CASCADE', [pmid])
                else:
                    amelie_paper_id = next_paper_id
                    next_paper_id += 1
                max_topical_gene_score = max(value['topical_genes'].values(), default=0)
                if max_topical_gene_score < 0.1:
                    # print("PubMed ID %s has a max topical gene score of less than .1, skipping" % pmid)
                    continue

                # print("Adding PubMed ID %s" % pmid)
                save_amelie_paper(amelie_paper_id, cur, pmid, value)
                next_paper_id += 1
                for classified_variant in value['variants']:
                    save_extracted_variant(amelie_paper_id, classified_variant,
                                           ensembl_to_gene_ids, cur,
                                           vscore_cutoff=args.vscore_cutoff)
                for ensg, score in value['topical_genes'].items():
                    save_topical_gene(amelie_paper_id, cur, ensg, score)
                for hpo_id in value['phenotypes']:
                    save_paper_phenotype(amelie_paper_id, cur, hpo_id)
                previously_loaded_pmids.add(pmid)
            cur.execute("INSERT INTO amelie_pubmed_update ("
                        "pubmed_update, "
                        "pubmed_update_index) "
                        "VALUES (?,?)",
                        (update['update'], update['update_index']))
            cur.execute("COMMIT")
        except:
            cur.execute("ROLLBACK")
            raise
    cur.close()


def save_paper_phenotype(amelie_paper_id, cur, hpo_id):
    phenotype_id = get_phenotype_id(hpo_id, cur)
    if phenotype_id is not None:
        cur.execute("INSERT INTO amelie_paper_phenotypes (paper_id, phenotype_id) "
                    "VALUES (?,?)", (amelie_paper_id, phenotype_id))


def save_topical_gene(amelie_paper_id, cur, ensg, score):
    if score >= 0.1:
        gene_ids = get_gene_ids([ensg], cur)
        for gene_id in gene_ids:
            cur.execute("INSERT INTO amelie_papergene (right_gene_score, gene_id, paper_id) "
                        "VALUES (?,?,?)", (score, gene_id, amelie_paper_id))


def save_extracted_variant(amelie_paper_id, classified_variant, ensembl_to_gene_ids, cur, vscore_cutoff):
    vscore = classified_variant['variant_score']
    v = classified_variant['variant']
    gene_ids = ensembl_to_gene_ids[v['gene']]
    modifications = v['modifications']
    assert len(gene_ids) > 0
    # print(str(v))
    if vscore >= vscore_cutoff:
        # print("Variant has score %.2f, adding" % vscore)
        for gene_id in gene_ids:
            for m in modifications:
                if m['vcf_pos'] is None or m['vcf_pos'] == "None" \
                        or m['vcf_ref'] is None or m['vcf_ref'] == "None" \
                        or m['vcf_alt'] is None or m['vcf_alt'] == "None" \
                        or m['variant_type'] == "synonymous" \
                        or len(m['vcf_ref']) > 20 \
                        or len(m['vcf_alt']) > 20:
                    continue
                cur.execute("INSERT INTO amelie_extractedvariant (chrom, bed_start, bed_end,"
                            "checked, variant_type, refseqprot_id, strand, orig_str, vcf_pos,"
                            "vcf_ref, vcf_alt, refseq_id, start_codon, classification_score,"
                            "gene_id, paper_id, regions) "
                            "VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)",
                            (m['chrom'], m['bed_start'], m['bed_end'], m['checked'], m['variant_type'],
                             m['grounded_refseq'], m['strand'], v['variant_text'], m['vcf_pos'], m['vcf_ref'],
                             m['vcf_alt'], m['nm_refseq'], m['codon_num'], vscore,
                             gene_id, amelie_paper_id, m['region']))
            # else:
            #     print("Variant has score %.2f, skipping" % vscore)


def save_amelie_paper(amelie_paper_id, cur, pmid, value):
    max_topical_gene_score = max(value['topical_genes'].values(), default=0)
    inheritance_mode = get_simple_inheritance_mode(value["inheritance_modes"])
    title = value["title"]
    medline_creation_date = parse_publish_date(value['publish_date'], pmid=pmid)
    cur.execute("INSERT INTO amelie_paper (id, pmid, inheritance_mode,"
                "title, medline_creation_date, authors, journal, destructive_score,"
                "not_destructive_score, dominant_score, recessive_score,"
                "strict_score, highest_right_gene_score) VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?)",
                (amelie_paper_id, pmid, inheritance_mode, title, medline_creation_date,
                 value['authors'], value['journal'], value['variant_types']['destructive'],
                 value['variant_types']['not_destructive'], value['inheritance_modes']['dominant'],
                 value['inheritance_modes']['recessive'], value['relevance'], max_topical_gene_score))


def load_hpo(cur, hpo_owl=None):
    hpo_info = hpo.HPO(hpo_owl=hpo_owl)
    term_id = 0
    term_ids = {}
    for term in hpo_info.all_descendants_mapping:
        assert term in hpo_info.terms
        assert 'name' in hpo_info.terms[term]
        name = hpo_info.terms[term]['name'][0]
        if 'def' in hpo_info.terms[term]:
            definition = hpo_info.terms[term]['def'][0]
        else:
            definition = None
        if 'alt_id' in hpo_info.terms[term]:
            alt_ids = ','.join(hpo_info.terms[term]['alt_id'])
        else:
            alt_ids = None
        cur.execute("INSERT INTO amelie_phenotype (id, hpo_id, name, definition, alt_ids) "
                    "VALUES (?,?,?,?,?)", (term_id, term, name, definition, alt_ids))
        term_ids[term] = term_id
        term_id += 1
    hpoedge_id = 0
    for parent, children in hpo_info.children_mapping.items():
        for child in children:
            cur.execute("INSERT INTO amelie_hpoedge (id, bottom_term_id, top_term_id) VALUES (?,?,?)",
                        (hpoedge_id, term_ids[child], term_ids[parent]))
            hpoedge_id += 1
    descendant_id = 0
    for parent, descendants in hpo_info.all_descendants_mapping.items():
        for descendant in descendants:
            cur.execute("INSERT INTO amelie_hpodescendants (id, bottom_term_id, top_term_id) VALUES (?,?,?)",
                        (descendant_id, term_ids[descendant], term_ids[parent]))
            descendant_id += 1


def load_genes(cur, genes_table=None):
    if genes_table is None:
        resource = 'gene_data/amelie_gene.csv'
    else:
        resource = genes_table
    with pkg_resources.resource_stream(__name__, resource) as gene_table:
        data = gene_table.read().decode('utf-8').split('\n')
        for line in data:
            if line.startswith('#'):
                continue
            if len(line) == 0:
                continue
            line = [x.strip('"') for x in line.strip().split(',')]
            # print("Row: %s" % str(line))
            id, ensembl_id, gene_symbol, rvis_score, mcap_gene, pli_score, created_on = line
            cur.execute("INSERT INTO amelie_gene (id, ensembl_id, gene_symbol,"
                        "rvis_score, mcap_gene, pli_score, created_on) "
                        "VALUES (?,?,?,?,?,?,?)", (id, ensembl_id,
                                                   gene_symbol, rvis_score,
                                                   mcap_gene, pli_score,
                                                   created_on))


def load_ensembl_id_to_gene_ids(connection):
    cur = connection.cursor()
    results = cur.execute("SELECT ensembl_id, id FROM amelie_gene")
    rv = defaultdict(lambda: [])
    for x in results:
        rv[x['ensembl_id']].append(x['id'])
    return rv


def create_knowledgebase_db(db_filename, hpo_owl=None, genes_table=None):
    connection = sqlite3.connect(db_filename)
    connection.row_factory = dotdict_factory
    cur = connection.cursor()
    cur.execute("PRAGMA synchronous = OFF")
    cur.execute("BEGIN TRANSACTION")
    with pkg_resources.resource_stream(__name__, 'sql_schemas/'
                                                 'knowledgebase.sql') as kb_schema:
        cur.executescript(kb_schema.read().decode('utf-8'))
    load_hpo(cur, hpo_owl=hpo_owl)
    load_genes(cur, genes_table=genes_table)
    cur.execute("COMMIT")
    connection.commit()
    cur.close()
    connection.close()


def load_all_papers_program():
    parser = argparse.ArgumentParser(description='The amelie tool loads to db')
    parser.add_argument('out_dir',  type=str, default='amelie_out_dir',
                        help='The name of the database to output')
    parser.add_argument('db_filename', type=str)
    parser.add_argument('--overwrite', action='store_true', default=False)
    parser.add_argument('--update_file', type=str, default=None)
    parser.add_argument('--create_db', action='store_true', default=False)
    parser.add_argument('--out_extracteddir', type=str, default='extracted')
    parser.add_argument('--vscore_cutoff', type=float, default=0.9)
    parser.add_argument('--max_year', type=int, default=None)
    parser.add_argument('--max_month', type=int, default=None)
    args = parser.parse_args()
    print(args)

    if not os.path.isdir(args.out_dir):
        print("Out dir %s doesn't exist, exiting" % args.out_dir)
        sys.exit(1)

    if args.create_db and not os.path.isfile(args.db_filename):
        create_knowledgebase_db(args.db_filename)

    if not os.path.isfile(args.db_filename):
        print("Knowledgebase database %s doesn't exist; create with --create_db" % args.db_filename)
        sys.exit(1)

    connection = sqlite3.connect(args.db_filename)
    connection.row_factory = dotdict_factory

    if args.overwrite:
        do_not_process_updates = []
    else:
        do_not_process_updates = get_db_loaded_updates(connection=connection)
    print("Have %d loaded updates (%s, ...)" % (len(do_not_process_updates), ', '.join(do_not_process_updates[:3])))

    if args.update_file is None:
        extract_papers_updates = get_extract_papers_updates(args=args)
    else:
        extract_papers_updates = [parse_update_from_filename(args.update_file)]
    print("Got %d updates ready (%s, ...)" % (len(extract_papers_updates), ', '.join(extract_papers_updates[:3])))

    updates_to_process = list(set(extract_papers_updates) - set(do_not_process_updates))
    print("Processing %d updates (%s, ...)" % (len(updates_to_process), ', '.join(updates_to_process[:3])))

    if args.overwrite:
        previously_loaded_pmids = set()
    else:
        previously_loaded_pmids = get_db_loaded_pmids(connection=connection)

    ensembl_to_gene_ids = load_ensembl_id_to_gene_ids(connection)

    updates_to_process.sort()
    for update in updates_to_process:
        process_update(update, args=args, connection=connection, previously_loaded_pmids=previously_loaded_pmids,
                       ensembl_to_gene_ids=ensembl_to_gene_ids)

    connection.close()
