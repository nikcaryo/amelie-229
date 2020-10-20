import json
from shutil import copyfile

def copy_training_data(omim_json, dataset_meta_json, og_process_dir, new_process_dir, new_out_dir):
    omim = json.load(open(omim_json))
    gwas_pmids = json.load(open(dataset_meta_json))['gwas_pmids']
    all_pmids = list(omim.keys()) + gwas_pmids

    copyfile(dataset_meta_json, f"{new_out_dir}/dataset_meta.json")
    copyfile(omim_json, f"{new_out_dir}/omim.json")
    
    succeded = 0
    failed = 0
    count = 0
    print(f"Attempting to copy {len(all_pmids)} pickled papers from {og_process_dir} to {new_process_dir}...")

    for pmid in all_pmids:
        count += 1
        try:
            if count % 1000 == 0:
                print(round(count/len(all_pmids), 2))
            succeded += 1
            process_path = f"{og_process_dir}/{pmid}.pkl"
            new_path = f"{new_process_dir}/{pmid}.pkl"
            copyfile(process_path, new_path)
        except Exception as e:
            failed += 1
    
    
    
    print("Done copying training data (OMIM and GWAS)!")
    print("PMIDs Succeded", succeded)
    print("PMIDs Failed:", failed)
      
if __name__ == "__main__":
    omim_json = "/data/amelie/amelie_out_dir/omim.json"
    dataset_meta_json = "/data/amelie/amelie_out_dir/dataset_meta.json"
    og_process_dir = "/data/amelie/amelie_process_dir"
    new_process_dir = "../amelie/amelie_process_dir"
    new_out_dir = "../amelie/amelie_out_dir"
      
    copy_training_data(omim_json, dataset_meta_json, og_process_dir, new_process_dir, new_out_dir)
