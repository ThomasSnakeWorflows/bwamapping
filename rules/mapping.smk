
import re

checkpoint splitfastq:
    input:
        get_fastq
    output:
        directory("mapping/{sample}/fastq")
    threads:
        12
    params:
        threads=12
    run:
        sample = wildcards.sample
        num_batches = get_batch_num(input[0])
        os.mkdir("mapping/%s/fastq" % sample)
        print("num batches %d" % num_batches)
        if is_interleaved(sample):
            command = "fastqsplitter -t %d -i %s " % (params.threads, input[0])
            for batch in range(num_batches):
                command += "-o mapping/%s/fastq/batch%s.fastq.gz " % (sample, batch)
            eprint(command)
            os.system(command)
        else:
            commandR1 = "fastqsplitter -t %d -i %s " % (params.threads, input[0])
            print("It is not interleaved")
            for batch in range(num_batches):
                commandR1 += "-o mapping/%s/fastq/batch%s_R1.fastq.gz " % (sample, batch)
            eprint(commandR1)
            os.system(commandR1)
            commandR2 = "fastqsplitter -t %d -i %s " % (params.threads, input[1])
            for batch in range(num_batches):
                commandR2 += "-o mapping/%s/fastq/batch%s_R2.fastq.gz " % (sample, batch)
            eprint(commandR2)
            os.system(commandR2)


rule bwamap:
    input:
        get_fastq_batch
    output:
        bam = temp("mapping/{sample}/{batch}.bam")
    params:
        reference = get_genome,
        rg = get_read_group,
        bwamem_params = get_bwamem_params(config, "-M -T 30"),
        fastqs = get_fastqbatch4bwamem
    threads:
        8
    log:
        stdout = "logs/bwa_{sample}_{batch}.o",
        stderr = "logs/bwa_{sample}_{batch}.e"
    shell:
        """
        rm -f {output.bam}*
        bwa mem -t {threads} {params.bwamem_params} -R {params.rg} \
           {params.reference} {params.fastqs} 2> {log.stderr} | \
           samtools view -bS - 2> {log.stderr} | \
           samtools sort -o {output.bam} 1> {log.stdout} 2> {log.stderr}
        """


def aggregate_input(wildcards):
    checkpoint_output = checkpoints.splitfastq.get(**wildcards).output[0]
    raw_batches = glob_wildcards(os.path.join(checkpoint_output,
                                      "{b}.fastq.gz")).b
    batches = list(set([re.sub("\_R1|\_R2","",batch) for batch in raw_batches]))
    return expand("mapping/{sample}/{batch}.bam",
                  sample=wildcards.sample,
                  batch=batches)

rule aggregate:
    input:
        aggregate_input
    output:
        "mapping/{sample}/{sample}.bam"
    shell:
        """
        rm -fr mapping/{wildcards.sample}/fastq
        samtools merge -O BAM {output} {input}
        """


rule bamindex:
    input:
        "mapping/{sample}/{sample}.bam"
    output:
        bai = "mapping/{sample}/{sample}.bam.bai",
        summary = "mapping/{sample}/{sample}.bam.flagstat"
    threads:
        1
    log:
        stdout = "logs/bamindex_{sample}.o",
        stderr = "logs/bamindex_{sample}.e"
    shell:
        "samtools index {input}; "
        "samtools flagstat {input} > {output.summary}"
