# bwamapping
A snakemake workflows for aligning reads


Setting the correct environment :
 (I use a conda env with snakemake-minimal=5.8.2)
```bash
conda activate snakemake
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
