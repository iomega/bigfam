# Compute LDA of MiBIG (documents) and bigfam subdomains (words)

```
pipenv install
pipenv shell
sudo apt install hmmer

mkdir MIBiG
wget https://mibig.secondarymetabolites.org/mibig_gbk_1.4.tar.gz
tar -zxf mibig_gbk_1.4.tar.gz
cd ..

mkdir json
cd temp_scripts
python build_json.py ../MIBiG ../json
python apply_hmmscan.py ../json ../bigfam_models/BiGFAM-v1-domains.hmm
python apply_subdomainscan.py ../json ../bigfam_models/BiGFAM-v1-coredomains.txt ../bigfam_models/subdomains
```

```
jupyter notebook
```

Open and run lda.opynb notebook.
