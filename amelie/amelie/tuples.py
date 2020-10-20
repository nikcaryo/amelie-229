import collections

PhenoMentionCore = [("wordidxs", "int[]"), 
                    ("words", "text[]"), 
                    ("hpo_id", "text"), 
                    ("canon_name", "text")]
GeneMentionCore = [("gene_name", "text"), 
                   ("wordidxs_start", "int[]"), 
                   ("wordidxs_end", "int[]"), 
                   ("num_mentions", "int"), 
                   ("mapping_type", "text"),
                   ("surrounding_words_left", "text[]"),
                   ("surrounding_words_right", "text[]"),
                   ("orig_words", "text")]
GeneNameCore = [("eid", "text"),
                ("canonical_name", "text"),
                ("mapping_type", "text")]


def pull_fields(core):
    return [n[0] for n in core]

PhenoMention = collections.namedtuple("PhenoMention", pull_fields(PhenoMentionCore))
GeneMention = collections.namedtuple("GeneMention", pull_fields(GeneMentionCore))
GeneName = collections.namedtuple("GeneName", pull_fields(GeneNameCore))

DNAModification = collections.namedtuple("DNAModification", ["chrom",
                                                             "bed_start",
                                                             "bed_end",
                                                             "gene_symbol",
                                                             "entrez_id",
                                                             "checked",
                                                             "variant_type",
                                                             "grounded_refseq",
                                                             "strand",
                                                             "vcf_pos",
                                                             "vcf_ref",
                                                             "vcf_alt",
                                                             "nm_refseq",
                                                             "codon_num",
                                                             "region"])
VariantMention = collections.namedtuple("VariantMention", ["variant_text", "gene", "modifications"])
VCFCoords = collections.namedtuple("VCFCoords", ["chrom", "pos", "ref", "alt"])
