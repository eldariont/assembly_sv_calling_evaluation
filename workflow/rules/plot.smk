#SV size evaluation

rule generate_calls_outside_regions:
    input:
        calls = "pipeline/calls/{tool}/variants.indel.vcf.gz",
        index = "pipeline/calls/{tool}/variants.indel.vcf.gz.tbi",
        regions = config["truth"]["bed"]
    output:
        "pipeline/{truvari,(truvari|truvari_multimatch)}/calls/{tool}_outside/outside-call.vcf"
    shell:
        "(zgrep \"#\" {input.calls} && bedtools intersect -v -f 1.0 -a {input.calls} -b {input.regions}) | bcftools view -i 'ABS(STRLEN(REF)-STRLEN(ALT))>=50 & ABS(STRLEN(REF)-STRLEN(ALT))<=50000 & FILTER==\"PASS\"' > {output}"


rule generate_base_outside_regions:
    input:
        base = config["truth"]["vcf"],
        regions = config["truth"]["bed"]
    output:
        "pipeline/{truvari,(truvari|truvari_multimatch)}/calls/{tool}_outside/outside-base.vcf"
    shell:
        "(zgrep \"#\" {input.base} && bedtools intersect -v -f 1.0 -a {input.base} -b {input.regions}) | bcftools view -i 'ABS(STRLEN(REF)-STRLEN(ALT))>=50 & ABS(STRLEN(REF)-STRLEN(ALT))<=50000 & FILTER==\"PASS\"' > {output}"

rule sort_vcf:
    input:
        "{name}.vcf"
    output:
        temp("{name}.sorted.vcf")
    shell:
        "bcftools sort {input} > {output}"

rule annotate_repeats:
    input:
        vcf = "{name}.sorted.vcf",
        conf = config["vcfanno_config"]
    output:
        "{name}.annotated.vcf"
    shell:
        "vcfanno_linux64 {input.conf} {input.vcf} > {output}"

rule plot:
    input:
        "{path}/tp-call.annotated.vcf",
        "{path}/tp-base.annotated.vcf",
        "{path}/fp.annotated.vcf",
        "{path}/fn.annotated.vcf",
        "{path}_outside/outside-call.annotated.vcf",
        "{path}_outside/outside-base.annotated.vcf" 
    output:
        "{path}/size_reps.png"
    shell:
        "python3 workflow/scripts/truvari_by_size.py {wildcards.path} {output}"

#Precision-Recall plots

rule reformat_truvari_results:
    input:
        "pipeline/{truvari}/{mode}/dipcall/summary.txt"
    output:
        "pipeline/{truvari,(truvari|truvari_multimatch)}/{mode,(calls|gtcomp)}/dipcall/pr_rec.txt"
    threads: 1
    shell:
        "cat {input} | grep 'precision\|recall' | grep -v 'gt' | tr -d ',' |sed 's/^[ \t]*//' | tr -d '\"' | tr -d ' ' | tr ':' '\t' | awk 'OFS=\"\\t\" {{ print \"dipcall\", \"{wildcards.truvari}\", \"{wildcards.mode}\", $1, $2 }}' > {output}"

rule reformat_truvari_results_svim:
    input:
        "pipeline/{truvari}/{mode}/svim/{parameters}/summary.txt"
    output:
        "pipeline/{truvari,(truvari|truvari_multimatch)}/{mode,(calls|gtcomp)}/svim/{parameters}/pr_rec.txt"
    threads: 1
    shell:
        "cat {input} | grep 'precision\|recall' | grep -v 'gt' | tr -d ',' |sed 's/^[ \t]*//' | tr -d '\"' | tr -d ' ' | tr ':' '\t' | awk 'OFS=\"\\t\" {{ print \"svim\", \"{wildcards.truvari}\", \"{wildcards.mode}\", $1, $2 }}' > {output}"


rule cat_truvari_results:
    input:
        svim = expand("pipeline/{truvari}/{mode}/svim/1000_1000_2000_2000_200/pr_rec.txt", 
                          truvari = ["truvari", "truvari_multimatch"],
                          mode = ["calls", "gtcomp"]),
        dipcall = expand("pipeline/{truvari}/{mode}/dipcall/pr_rec.txt", 
                          truvari = ["truvari", "truvari_multimatch"],
                          mode = ["calls", "gtcomp"])
    output:
        all = "pipeline/eval/all_results.txt"
    threads: 1
    shell:
        "cat {input.svim} {input.dipcall} > {output.all}"

rule cat_truvari_results_svim_parameters:
    input:
        svim = expand("pipeline/{truvari}/{mode}/svim/{parameters}/pr_rec.txt", 
                          truvari = ["truvari", "truvari_multimatch"],
                          mode = ["calls", "gtcomp"],
                          parameters = ["50_50_50_50_200", "100_100_2000_2000_200", "1000_1000_2000_2000_200"])
    output:
        all = "pipeline/eval/svim_parameter_results.txt"
    threads: 1
    run:
        with open(output.all, 'w') as output_file:
            for f in input.svim:
                parameters = f.split("/")[4]
                with open(f, 'r') as input_file:
                    for line in input_file:
                        print("%s\t%s" % (parameters, line.strip()), file=output_file)

rule plot_pr_tools:
    input:
        "pipeline/eval/all_results.txt"
    output:
        "pipeline/eval/results.tools.png"
    threads: 1
    log:
        "pipeline/logs/rplot.tools.log"
    shell:
        "Rscript --vanilla workflow/scripts/plot_tools.R {input} {output} > {log}"


rule plot_pr_svim_parameters:
    input:
        "pipeline/eval/svim_parameter_results.txt"
    output:
        "pipeline/eval/results.svim.parameters.png"
    threads: 1
    log:
        "pipeline/logs/rplot.svim.log"
    shell:
        "Rscript --vanilla workflow/scripts/plot_svim_parameters.R {input} {output} > {log}"
