# bwamapping
A snakemake workflows for aligning reads

Setting the correct environment (see below for the snakemake python3 virtual env) :
```bash
source bwaenv/bin/activate
module load bioinfo/fastp-0.19.4
module load bioinfo/bwa-0.7.17
module load bioinfo/samtools-1.9
module load bioinfo/qualimap-11-12-16
```

- **Warning** : The genome file must be bwa indexed

- Make appropriate changes to the config.yaml

- Provide a sample file

Executing the pipeline
```bash
snakemake --jobs 30 --cluster-config cluster.yaml --drmaa " --mem-per-cpu={cluster.mem}000 --mincpus={threads} --time={cluster.time} -J {cluster.name} -N 1=1" -p -n
```
In order to use snakemake 5.8.2 I use a python3 venv
``bash 
module load system/Python-3.6.3
python3 -m venv bwaenv
source bwaenv/bin/activate
pip install -r requirements.txt 
``



In order to use snakemake 5.8.2 I use a conda env (see environment.yaml) easily created by
```bash 
conda env create --prefix ./env  -f environment.yaml 
```



