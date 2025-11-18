import sys
import pandas as pd
from snakemake.utils import validate

include: "rules/common.smk"

samples = pd.read_table(config["samples"],
                        dtype={"sample": str}).set_index("sample", drop=False)
#validate(samples, schema="../schemas/samples.schema.yaml")


if not bwa_index_exists(config["ref"]["genome"]):
    print_error_exit("%s is not bwa-indexed" % config["ref"]["genome"])

# Modules
include: "rules/mapping.smk"
include: "rules/stats.smk"

workdir: config['workdir']

# Wildcard constraints
wildcard_constraints:
    sample="|".join(samples.index),
    batch=r"batch\d+"


# Target rules
rule all:
    input:
        expand("mapping/{sample}/{sample}.bam.bai", sample=samples.index),
        expand("stats/{sample}/fastp.html", sample=samples.index),
        expand("mapping/{sample}/qualimap/qualimapReport.html", sample=samples.index),
        expand("stats/{sample}_samtools.stats", sample=samples.index),
        "stats/multiqc/multiqc_report.html"


rule fastpall:
	input:
		expand("stats/{sample}/fastp.html", sample=samples.index)
	output:
		touch("fastp.done")
