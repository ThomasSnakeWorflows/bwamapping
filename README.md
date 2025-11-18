# bwamapping
A snakemake workflows for aligning reads

Setting the correct conda environment :
```bash
mamba create -n bwammaping -c bioconda snakemake=9.10.1
mamba activate bwammaping
mamba install bioconda::fastp
mamba install bioconda::bwa
mamba install bioconda::qualimap
mamba install bioconda::fastqsplitter
mamba install bioconda::samtools
pip install snakemake-executor-plugin-slurm
pip install termcolor==1.1.0
pip install multiQC==1.9
```
# Make sure to activate your conda environment
```
mamba activate bwammaping
```

- **Warning** : The genome file must be bwa indexed, if not run

```
bwa index REF.fa
```

- Make appropriate changes to the config.yaml

- Provide a sample file

Executing the pipeline
```bash
snakemake --configfile config.yaml --profile genotoul -j 40 -p -n
```

- **MutiQC example**

[MultiQC for the testdata](https://htmlpreview.github.io/?https://github.com/ThomasSnakeWorflows/bwamapping/blob/master/multiqc_report.html)

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
