

rule fastp:
    input:
        get_fastq
    output:
        html = "mapping/{sample}/fastp.html",
        json = "mapping/{sample}/fastp.json"
    params:
        fastqs = get_fastq4fastp
    log:
        stdout = "logs/map_{sample}.o",
        stderr = "logs/map_{sample}.e"
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
        get_threads("qualimap", 8)
    params:
        mem = get_mem("qualimap", 4),
        outdir = "mapping/{sample}/qualimap"
    log:
        stdout = "logs/qualimap_{sample}.o",
        stderr = "logs/qualimap_{sample}.e"
    shell:
        "qualimap bamqc -bam {input.bam} -nt {threads} "
        " --java-mem-size={params.mem}G -c "
        " -outdir {params.outdir}"
