

rule fastp:
    input:
        get_fastq
    output:
        html = "stats/{sample}/fastp.html",
        json = "stats/{sample}/fastp.json"
    params:
        fastqs = get_fastq4fastp
    log:
        stdout = "logs/fastp_{sample}.o",
        stderr = "logs/fastp_{sample}.e"
    shell:
        "fastp -A -L -Q -h {output.html} -j {output.json} {params.fastqs} "
        " 1>{log.stdout} 2>{log.stderr}"


rule qualimap:
    input:
        bam = "mapping/{sample}/{sample}.bam",
        bai = "mapping/{sample}/{sample}.bam.bai"
    output:
        report = "mapping/{sample}/qualimap/qualimapReport.html"
    threads:
        8
    params:
        mem = 4,
        outdir = "mapping/{sample}/qualimap"
    log:
        stdout = "logs/qualimap_{sample}.o",
        stderr = "logs/qualimap_{sample}.e"
    shell:
        "qualimap bamqc -bam {input.bam} -nt {threads} "
        " --java-mem-size={params.mem}G -c "
        " -outdir {params.outdir}"

rule samtoolsstats:
    input:
        bam = "mapping/{sample}/{sample}.bam",
        bai = "mapping/{sample}/{sample}.bam.bai"
    output:
        report = "stats/{sample}_samtools.stats"
    shell:
        "samtools stats {input.bam} > {output}"

rule multiqc:
    input:
        expand("stats/{sample}_samtools.stats", sample=samples.index)
    output:
        "stats/multiqc/multiqc_report.html"
    shell:
        "multiqc stats -o stats/multiqc "
