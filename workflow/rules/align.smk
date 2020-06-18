#########
# ALIGN #
#########

rule index_fasta:
    input:
        fasta = "{name}.fasta",
    output:
        mmi = "{name}.mmi",
    conda:
        "../envs/minimap2.yaml"
    shell:
        "minimap2 -d {output.mmi} {input.fasta}"

# Generate alignment in SAM
rule align_sam:
    input:
        reference = config["reference"],
        fasta = config["assembly"]
    output:
        log = "pipeline/logs/H{hap}.sam.log",
        sam = temp("pipeline/alignments/H{hap}.sam")
    threads: 8
    conda:
        "../envs/minimap2.yaml"
    shell:
        "minimap2 -ax asm5 -r2k -t {threads} {input.reference} {input.fasta} 2> {output.log} > {output.sam}"

# Sort SAM
rule sort_sam:
    input:
        sam = "{name}.sam"
    output:
        bam = "{name}.sort.bam"
    conda:
        "../envs/samtools.yaml"
    shell:
        "samtools sort -m4G -@4 -o {output.bam} {input.sam}"

rule index_bam:
    input:
        "{name}.bam"
    output:
        "{name}.bam.bai"
    conda:
        "../envs/samtools.yaml"
    shell:
        "samtools index {input}"