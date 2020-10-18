
from amelie.dataloaders.hpo import HPO
import argparse


def canonicalize_phenotypes_program():
    hpo = HPO()
    parser = argparse.ArgumentParser()
    parser.add_argument('input_file',
                        type=str)
    parser.add_argument('column',
                        help='1-indexed column containing HPO phenotypes to expand',
                        type=int)
    parser.add_argument('--add_info', action="store_true", default=False)

    args = parser.parse_args()

    with open(args.input_file) as f:
        for line in f:
            columns = line.rstrip('\n').split('\t')
            idx = args.column - 1
            hpo_id = columns[idx]
            if hpo_id not in hpo.terms:
                continue
            if args.add_info:
                info = "orig,%s" % [hpo.terms[hpo_id]['name'][0]]
            else:
                info = ""
            print('\t'.join(columns[:idx] + [hpo_id] + columns[idx+1:] + [info]))
            for ancestor in hpo.all_ancestors_mapping[hpo_id]:
                if ancestor == hpo_id:
                    continue  # adding that separately above
                if args.add_info:
                    info = "ancestor,%s" % [hpo.terms[ancestor]['name'][0]]
                else:
                    info = ""
                print('\t'.join(columns[:idx] + [ancestor] + columns[idx+1:] + [info]))
