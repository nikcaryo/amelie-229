# amelie-229
Final project for Stanford CS229

## Installation
- Run `conda env create -f environment.yml`
- Activate the conda environment with `conda activate amelie`
- `cd` to `amelie/` and run `python -m pip install -e .`
- `cd` to `amelie/pheno_data` and run `./BUILD`
- Find the path to your conda environment (`should be something like /cluster/u/nikcaryo/software/miniconda/envs/amelie/`)
- Open the python terminal with `python -i`
- Do `import nltk`
- Then `nltk.download('wordnet', download_dir="{path to your conda}/lib/nltk_data}"`
- Then `nltk.download('stopwords', download_dir="{path to your conda}/lib/nltk_data}"`
- Exit by typing `exit()`
- Download the training data and unzip all the `pkl` files into amelie_process_dir (James can just `cp` it from my directory)
- Run `jupyter notebook` in the folder containing `train.ipynb`
- Open your browser to whatever it address it gives you (and do ssh tunneling if you're on lab machine)
- Run all the blocks and hopefully it works lol 
