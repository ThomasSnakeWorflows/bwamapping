
import os
import sys
from termcolor import cprint
from collections import defaultdict
import math

# report: "../report/workflow.rst"

# Config file and sample sheets
configfile: "config.yaml"
validate(config, schema="../schemas/config.schema.yaml")


samples = pd.read_table(config["samples"],
                        dtype={"sample": str}).set_index("sample", drop=False)
validate(samples, schema="../schemas/samples.schema.yaml")


# Wildcard constraints
wildcard_constraints:
    sample="|".join(samples.index),
    batch="batch\d+"


def is_interleaved(sample):
    fastqs = samples.loc[sample, "reads"].split(",")
    return False if len(fastqs) > 1 else True


# Helper functions #####
def print_error_exit(message):
    """ Print soft error message in mangenta and exit """
    cprint("WARNING: " + message + ", exiting softly!", 'magenta',
           attrs=['bold'], file=sys.stderr)
    sys.exit(1)


def get_genome(wildcards):
    """ the genome file """
    return config["ref"]["genome"]


def get_fastq(wildcards):
    """Get fastq files of given sample-unit."""
    fastqs = samples.loc[wildcards.sample, "reads"].split(",")
    return fastqs


def get_fastq4fastp(wildcards):
    """Get fastq files of given sample-unit."""
    fastqs = samples.loc[wildcards.sample, "reads"].split(",")
    if len(fastqs) > 1:
        return "-i %s -I %s " % (fastqs[0], fastqs[1])
    else:
        return "-i %s " % fastqs[0]


# def get_fastq4bwamem(wildcards):
#     """Get fastq files of given sample-unit."""
#     fastqs = samples.loc[wildcards.sample, "reads"].split(",")
#     if len(fastqs ) > 1:
#         return "%s %s " % (fastqs[0], fastqs[1])
#     else:
#         return "-p %s " % fastqs[0]


def get_batch_num(fastq):
    size = Path(fastq).stat().st_size
    return math.ceil(size/float(config['targetsize']))


def get_fastq_batch(wildcards):
    sample = wildcards.sample
    batch = wildcards.batch
    if is_interleaved(sample):
        return ["mapping/%s/fastq/%s.fastq.gz" % (sample, batch)]
    else:
        return ["mapping/%s/fastq/%s_R1.fastq.gz" % (sample, batch),
                "mapping/%s/fastq/%s_R2.fastq.gz" % (sample, batch)]


def get_fastqbatch4bwamem(wildcards):
    """Get fastq files of given sample-unit."""
    fastqs = get_fastq_batch(wildcards)
    if len(fastqs) == 1:
        return "-p %s " % fastqs[0]
    else:
        return " ".join(fastqs)


def get_read_group(wildcards):
    " Note sample name and platform in read group. "
    return r"'@RG\tID:{sample}\tLB:{sample}\tPL:ILLUMINA\tSM:{sample}'".format(
        sample=wildcards.sample)


def get_bwamem_params(config, default):
    if "bwamem_params" in config:
        return config["bwamem_params"]
    return default


def get_threads(rule, default):
    cluster_config = snakemake.workflow.cluster_config
    if rule in cluster_config and "threads" in cluster_config[rule]:
        return cluster_config[rule]["threads"]
    if "default" in cluster_config and "threads" in cluster_config["default"]:
        return cluster_config["default"]["threads"]
    return default


def get_mem(rule, default):
    cluster_config = snakemake.workflow.cluster_config
    if rule in cluster_config and "mem" in cluster_config[rule]:
        return cluster_config[rule]["mem"]
    if "default" in cluster_config and "mem" in cluster_config["default"]:
        return cluster_config["default"]["mem"]
    return default


def bwa_index_exists(genome):
    return os.path.isfile("%s.bwt" % genome)
