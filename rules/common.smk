
import os
import sys
from termcolor import cprint

# report: "../report/workflow.rst"

# Config file and sample sheets
configfile: "config.yaml"
validate(config, schema="../schemas/config.schema.yaml")


samples = pd.read_table(config["samples"],
                        dtype={"sample": str}).set_index("sample", drop=False)
validate(samples, schema="../schemas/samples.schema.yaml")


# Wildcard constraints
wildcard_constraints:
    sample="|".join(samples.index)


# Helper functions #####
def print_error_exit(message):
    """ Print soft error message in mangenta and exit """
    cprint("WARNING: " + message + ", exiting softly!", 'magenta',
           attrs=['bold'], file=sys.stderr)
    sys.exit(1)


def get_genome(wildcards):
    """ the genome file """
    return config["ref"]["genome"]

#
# def get_fai(wildcards):
#     """ the genome faid index"""
#     return config["ref"]["genome"] + ".fai"


def get_fastq(wildcards):
    """Get fastq files of given sample-unit."""
    fastqs = samples.loc[wildcards.sample, "reads"].split(",")
    return fastqs


def get_fastq4fastp(wildcards):
    """Get fastq files of given sample-unit."""
    fastqs = samples.loc[wildcards.sample, "reads"].split(",")
    if len(fastqs)>1:
        return "-i %s -I %s " % (fastqs[0], fastqs[1])
    else:
        return "-i %s " % fastqs[0]


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
