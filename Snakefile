import sys
import pandas as pd
from snakemake.utils import validate

include: "rules/common.smk"

if not bwa_index_exists(config["ref"]["genome"]):
    print_error_exit("%s is not bwa-indexed" % config["ref"]["genome"])

# Modules
include: "rules/mapping.smk"
include: "rules/stats.smk"

workdir: config['workdir']

# Target rules
rule all:
    input:
        expand("mapping/{sample}/{sample}.bam.bai", sample=samples.index),
        expand("mapping/{sample}/fastp.html", sample=samples.index),
        expand("mapping/{sample}/qualimap/qualimapReport.html", sample=samples.index),
        expand("stats/{sample}_samtools.stats", sample=samples.index),
        "stats/multiqc/multiqc_report.html"


rule fastpall:
	input:
		expand("mapping/{sample}/fastp.html", sample=samples.index)
	output:
		touch("fastp.done")
