########
# CALL #
########

def get_haplotype1(wildcards):
    return expand(config["assembly"][wildcards.assembly], hap=[1])

def get_haplotype2(wildcards):
    return expand(config["assembly"][wildcards.assembly], hap=[2])

# Call variants with dipcall
rule prepare_dipcall:
    input:
        reference = config["reference"],
        hap1 = get_haplotype1,
        hap2 = get_haplotype2,
        par = "bin/dipcall.kit/hs37d5.PAR.bed"
    output:
        "pipeline/{assembly}/calls/dipcall/HG002.mak"
    params:
        prefix = "pipeline/{assembly}/calls/dipcall/HG002"
    shell:
        "bin/dipcall.kit/run-dipcall -t 10 -x {input.par} {params.prefix} {input.reference} {input.hap1} {input.hap2} > {output}"

rule run_dipcall:
    input:
        reference = config["reference"],
        hap1 = get_haplotype1,
        hap2 = get_haplotype2,
        make = "pipeline/{assembly}/calls/dipcall/HG002.mak"
    output:
        "pipeline/{assembly}/calls/dipcall/HG002.dip.bed",
        "pipeline/{assembly}/calls/dipcall/HG002.dip.vcf.gz"
    threads: 40
    shell:
        "make -j {threads} -f {input.make}"

# Filter out small variants and set FILTER field to PASS because truvari requires it
rule filter_dipcall:
    input:
        "pipeline/{assembly}/calls/dipcall/HG002.dip.vcf.gz"
    output:
        "pipeline/{assembly}/calls/dipcall/variants.indel.vcf.gz"
    shell:
        "bcftools view -i 'STRLEN(REF)>19 | STRLEN(ALT)>19' {input} | awk 'OFS=\"\\t\" {{ if($1 ~ /^#/) {{ print $0 }} else {{ $7 = \"PASS\"; print $0 }} }}' | bcftools sort -Oz > {output}"

rule call_svim_diploid:
    input:
        reference = config["reference"],
        bam1 = "pipeline/{assembly}/alignments/H1.sort.bam",
        index1 = "pipeline/{assembly}/alignments/H1.sort.bam.bai",
        bam2 = "pipeline/{assembly}/alignments/H2.sort.bam",
        index2 = "pipeline/{assembly}/alignments/H2.sort.bam.bai"
    output:
        "pipeline/{assembly}/calls/svim/{rgt}_{rot}_{qgt}_{qot}_{med}/variants.vcf"
    params:
        working_dir = "pipeline/{assembly}/calls/svim/{rgt}_{rot}_{qgt}_{qot}_{med}"
    conda:
        "../envs/svimasm.yaml"
    shell:
        "svim-asm diploid {params.working_dir} {input.bam1} {input.bam2} {input.reference} --min_sv_size 20 \
        --tandem_duplications_as_insertions --interspersed_duplications_as_insertions --reference_gap_tolerance {wildcards.rgt} --reference_overlap_tolerance {wildcards.rot} \
        --query_gap_tolerance {wildcards.qgt} --query_overlap_tolerance {wildcards.qot} --max_edit_distance {wildcards.med} --sample HG002 --query_names"


rule format_svim_diploid:
    input:
        "pipeline/{assembly}/calls/svim/{rgt}_{rot}_{qgt}_{qot}_{med}/variants.vcf"
    output:
        "pipeline/{assembly}/calls/svim/{rgt}_{rot}_{qgt}_{qot}_{med}/variants.indel.vcf.gz"
    shell:
        "bcftools view -i \"SVTYPE = 'INS' | SVTYPE = 'DEL'\" {input} | bcftools sort -Oz > {output}"

rule tabix:
    input:
        "{name}.vcf.gz"
    output:
        "{name}.vcf.gz.tbi"
    shell:
        "tabix -p vcf {input}"
