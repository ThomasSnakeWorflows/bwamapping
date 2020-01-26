



rule fastp:
    input:
        get_fastq
    output:
        html="stats/fastq/{sample}_fastp.html",
        json="stats/fastq/{sample}_fastp.json"
    params:
        fastqs=get_fastq4fastp
    log:
        stdout = "logs/map_{sample}.o",
        stderr = "logs/map_{sample}.e"
    shell:
        "./fastp -A -L -Q -h {output.html} -j {output.json} {params.fastqs} "
        " 1>{log.stdout} 2>{log.stderr}"

rule bwamap:
    input:
        get_fastq
    output:
        bam = "mapping/{sample}.bam"
    params:
        reference = get_genome,
        rg = get_read_group,
        bwamem_params = get_bwamem_params(config, "-M -T 30")
    threads:
        get_threads("bwamap", 8)
    log:
        stdout = "logs/fastp_{sample}.o",
        stderr = "logs/fastp_{sample}.e"
    shell:
        "rm -f {output.bam}*; "
        "bwa mem -t {threads} {params.bwamem_params} -R {params.rg} "
        "{params.reference} {input} 2> {log.stderr} | "
        "samtools view -bS - 2> {log.stderr} | "
        "samtools sort -o {output.bam} 1> {log.stdout} 2> {log.stderr}"

rule bamindex:
    input:
        "mapping/{sample}.bam"
    output:
        bai = "mapping/{sample}.bam.bai",
        summary = "mapping/{sample}.bam.flagstat"
    threads:
        1
    log:
        stdout = "logs/bamindex_{sample}.o",
        stderr = "logs/bamindex_{sample}.e"
    shell:
        "samtools index {input}; "
        "samtools flagstat {input} > {output.summary}"

rule qualimap:
    input:
        bam = "mapping/{sample}.bam",
        bai = "mapping/{sample}.bam.bai"
    output:
        report = "mapping/stats/{sample}/qualimapReport.html"
    threads:
        get_threads("qualimap", 8)
    params:
        mem = get_mem("qualimap", 4),
        outdir = "mapping/stats/{sample}"
    log:
        stdout = "logs/qualimap_{sample}.o",
        stderr = "logs/qualimap_{sample}.e"
    shell:
        "qualimap bamqc -bam {input.bam} -nt {threads} "
        " --java-mem-size={params.mem}G -c "
        " -outdir {params.outdir}"
