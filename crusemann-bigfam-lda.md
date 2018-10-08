The http://doi.org/10.1021/acs.jnatprod.6b00722 article has genome and metabolomics data.

The genomes where downloaded from https://genome.jgi.doe.gov/portal/ and Genbank.
Run Anti-Smash on all genomes and store output in bgc_crusemann/ dir.

```
mkdir json_crusemann
cd temp_scripts
python build_json.py ../bgc_crusemann ../json_crusemann
python apply_hmmscan.py ../json_crusemann ../bigfam_models/BiGFAM-v1-domains.hmm
python apply_subdomainscan.py ../json_crusemann ../bigfam_models/BiGFAM-v1-coredomains.txt ../bigfam_models/subdomains
```
