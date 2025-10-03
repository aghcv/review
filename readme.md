
# PubMed Noise Index Script

This repo contains a reproducible environment and Python script to evaluate
the noisiness of PubMed search terms for systematic review protocols.

## Repository Structure
```
/pubmed-noise-index
├── environment.yml # Conda environment definition
├── requirements.txt # pip fallback
├── pubmed_noise_index.py # main script
├── keywords.txt # optional keyword list
├── filters.txt # optional filters
└── README.md # this file
```

## Setup with Conda
```bash
git clone https://github.com/aghcv/review.git
conda env create -f environment.yml
conda activate pubmed-search-env
```

