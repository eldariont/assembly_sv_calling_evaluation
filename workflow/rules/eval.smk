########
# EVAL #
########


# Run truvari
rule truvari:
    input:
        reference = config["reference"],
        truth = config["truth"]["vcf"],
        regions = config["truth"]["bed"],
        calls = "pipeline/calls/{tool}/variants.indel.vcf.gz",
        index = "pipeline/calls/{tool}/variants.indel.vcf.gz.tbi"
    output:
        "pipeline/truvari/calls/{tool}/summary.txt",
        "pipeline/truvari/calls/{tool}/tp-call.vcf",
        "pipeline/truvari/calls/{tool}/tp-base.vcf",
        "pipeline/truvari/calls/{tool}/fp.vcf",
        "pipeline/truvari/calls/{tool}/fn.vcf"
    params:
        out_dir = "pipeline/truvari/calls/{tool}"
    conda:
        "../envs/truvari.yaml"
    shell:
        """rm -rf {params.out_dir} && \
        truvari bench -f {input.reference} \
                -b {input.truth} \
                --includebed {input.regions} \
                -o {params.out_dir} \
                --giabreport \
                --passonly \
                -r 1000 \
                -p 0 \
                -c {input.calls}"""

rule truvari_gtcomp:
    input:
        reference = config["reference"],
        truth = config["truth"]["vcf"],
        regions = config["truth"]["bed"],
        calls = "pipeline/calls/{tool}/variants.indel.vcf.gz",
        index = "pipeline/calls/{tool}/variants.indel.vcf.gz.tbi"
    output:
        "pipeline/truvari/gtcomp/{tool}/summary.txt",
        "pipeline/truvari/gtcomp/{tool}/tp-call.vcf",
        "pipeline/truvari/gtcomp/{tool}/tp-base.vcf",
        "pipeline/truvari/gtcomp/{tool}/fp.vcf",
        "pipeline/truvari/gtcomp/{tool}/fn.vcf"
    params:
        out_dir = "pipeline/truvari/gtcomp/{tool}"
    conda:
        "../envs/truvari.yaml"
    shell:
        """rm -rf {params.out_dir} && \
        truvari bench -f {input.reference} \
                -b {input.truth} \
                --includebed {input.regions} \
                -o {params.out_dir} \
                --giabreport \
                --passonly \
                --gtcomp \
                -r 1000 \
                -p 0 \
                -c {input.calls}"""

rule truvari_multimatch:
    input:
        reference = config["reference"],
        truth = config["truth"]["vcf"],
        regions = config["truth"]["bed"],
        calls = "pipeline/calls/{tool}/variants.indel.vcf.gz",
        index = "pipeline/calls/{tool}/variants.indel.vcf.gz.tbi"
    output:
        "pipeline/truvari_multimatch/calls/{tool}/summary.txt",
        "pipeline/truvari_multimatch/calls/{tool}/tp-call.vcf",
        "pipeline/truvari_multimatch/calls/{tool}/tp-base.vcf",
        "pipeline/truvari_multimatch/calls/{tool}/fp.vcf",
        "pipeline/truvari_multimatch/calls/{tool}/fn.vcf"
    params:
        out_dir = "pipeline/truvari_multimatch/calls/{tool}"
    conda:
        "../envs/truvari.yaml"
    shell:
        """rm -rf {params.out_dir} && \
        truvari bench -f {input.reference} \
                -b {input.truth} \
                --includebed {input.regions} \
                -o {params.out_dir} \
                --giabreport \
                --multimatch \
                --passonly \
                -r 1000 \
                -p 0 \
                -c {input.calls}"""

rule truvari_multimatch_gtcomp:
    input:
        reference = config["reference"],
        truth = config["truth"]["vcf"],
        regions = config["truth"]["bed"],
        calls = "pipeline/calls/{tool}/variants.indel.vcf.gz",
        index = "pipeline/calls/{tool}/variants.indel.vcf.gz.tbi"
    output:
        "pipeline/truvari_multimatch/gtcomp/{tool}/summary.txt",
        "pipeline/truvari_multimatch/gtcomp/{tool}/tp-call.vcf",
        "pipeline/truvari_multimatch/gtcomp/{tool}/tp-base.vcf",
        "pipeline/truvari_multimatch/gtcomp/{tool}/fp.vcf",
        "pipeline/truvari_multimatch/gtcomp/{tool}/fn.vcf"
    params:
        out_dir = "pipeline/truvari_multimatch/gtcomp/{tool}"
    conda:
        "../envs/truvari.yaml"
    shell:
        """rm -rf {params.out_dir} && \
        truvari bench -f {input.reference} \
                -b {input.truth} \
                --includebed {input.regions} \
                -o {params.out_dir} \
                --giabreport \
                --multimatch \
                --passonly \
                --gtcomp \
                -r 1000 \
                -p 0 \
                -c {input.calls}"""