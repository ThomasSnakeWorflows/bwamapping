import sys
import pandas as pd
from snakemake.utils import validate

include: "rules/common.smk"

if not bwa_index_exists(config["ref"]["genome"]):
    print_error_exit("%s is not bwa-indexed" % config["ref"]["genome"])
# Modules
include: "rules/mapping.smk"

# Target rules
rule all:
    input:
        expand("mapping/{sample}.bam", sample=samples.index)