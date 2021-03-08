# bwamapping
A snakemake workflows for aligning reads

Setting the correct conda environment :
```bash
module load bioinfo/fastp-0.19.4
module load bioinfo/bwa-0.7.17
module load bioinfo/samtools-1.9
module load bioinfo/qualimap-11-12-16
module load system/Python-3.6.3
conda create -p ./bwaenv python=3.6.3
conda install fastqplitter
pip install -r requirements.txt
```

- **Warning** : The genome file must be bwa indexed

- Make appropriate changes to the config.yaml

- Provide a sample file

Executing the pipeline
```bash
snakemake --jobs 99 --cluster-config cluster.yaml --drmaa " --mem-per-cpu={cluster.mem-per-cpu}000 --mincpus={threads} --time={cluster.time} -J {cluster.name} -N 1=1" -p -n
```

- **MutiQC example**

Copy/paste multiqc_report.html path to https://htmlpreview.github.io/?



- **Generating the test data set from genome.fa.gz**

```python
import random
with open("random.fa","w") as fout:
  N=800000
  seq=''.join(random.choices(['A','T','C','G'], k=N))
  l_width= 60
  lines = [seq[i:i+l_width] for i in range(0, len(seq), l_width)]
  fastaseq = ">%s\n" %("random")
  fastaseq += "\n".join(lines) + "\n"
  fout.write(fastaseq)
```

```
wgsim -1 100 -2 100 -N 10000 genome.fa.gz tmp_reads_R1.fastq tmp_reads_R2.fastq
wgsim -1 100 -2 100 -N 500 random.fa tmp_randreads_R1.fastq tmp_randreads_R2.fastq
cat tmp_reads_R1.fastq tmp_randreads_R1.fastq | gzip -c > readsB_R1.fastq.gz
cat tmp_reads_R2.fastq tmp_randreads_R2.fastq | gzip -c > readsB_R2.fastq.gz
rm -f tmp_*.fastq
```
