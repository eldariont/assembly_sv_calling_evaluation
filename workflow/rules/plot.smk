#SV size evaluation

rule generate_calls_outside_regions:
    input:
        calls = "pipeline/{assembly}/calls/{tool}/variants.indel.vcf.gz",
        index = "pipeline/{assembly}/calls/{tool}/variants.indel.vcf.gz.tbi",
        regions = config["truth"]["bed"]
    output:
        "pipeline/{assembly}/{truvari,(truvari|truvari_multimatch)}/calls/{tool}_outside/outside-call.vcf"
    shell:
        "(zgrep \"#\" {input.calls} && bedtools intersect -v -f 1.0 -a {input.calls} -b {input.regions}) | bcftools view -i 'ABS(STRLEN(REF)-STRLEN(ALT))>=50 & ABS(STRLEN(REF)-STRLEN(ALT))<=50000 & FILTER==\"PASS\"' > {output}"


rule generate_base_outside_regions:
    input:
        base = config["truth"]["vcf"],
        regions = config["truth"]["bed"]
    output:
        "pipeline/{assembly}/{truvari,(truvari|truvari_multimatch)}/calls/{tool}_outside/outside-base.vcf"
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

rule plot_SV_size:
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

rule plot_performance_by_size:
    input:
        "{path}/tp-call.annotated.vcf",
        "{path}/tp-base.annotated.vcf",
        "{path}/fp.annotated.vcf",
        "{path}/fn.annotated.vcf",
    output:
        "{path}/performance_by_size.png"
    shell:
        "python3 workflow/scripts/performance_by_size.py {wildcards.path} {output}"


#Precision-Recall plots

rule reformat_truvari_results_error_types:
    input:
        tp_call = "pipeline/{assembly}/truvari/calls/dipcall/tp-call.vcf",
        tp_base = "pipeline/{assembly}/truvari/calls/dipcall/tp-base.vcf",
        fp = "pipeline/{assembly}/truvari/calls/dipcall/fp.vcf",
        fn = "pipeline/{assembly}/truvari/calls/dipcall/fn.vcf"
    output:
        "pipeline/{assembly,(wenger|garg)}/truvari/calls/dipcall/error_types.txt"
    threads: 1
    shell:
        "python workflow/scripts/compare_genotypes.py {input.tp_call} {input.tp_base} {input.fp} {input.fn} {wildcards.assembly} dipcall > {output}"

rule reformat_truvari_results_error_types_svim:
    input:
        tp_call = "pipeline/{assembly}/truvari/calls/svim/{parameters}/tp-call.vcf",
        tp_base = "pipeline/{assembly}/truvari/calls/svim/{parameters}/tp-base.vcf",
        fp = "pipeline/{assembly}/truvari/calls/svim/{parameters}/fp.vcf",
        fn = "pipeline/{assembly}/truvari/calls/svim/{parameters}/fn.vcf"
    output:
        "pipeline/{assembly,(wenger|garg)}/truvari/calls/svim/{parameters}/error_types.txt"
    threads: 1
    shell:
        "python workflow/scripts/compare_genotypes.py {input.tp_call} {input.tp_base} {input.fp} {input.fn} {wildcards.assembly} svim > {output}"


rule cat_truvari_results_error_types:
    input:
        svim = expand("pipeline/{assembly}/truvari/calls/svim/1000_1000_2000_2000_200/error_types.txt", 
                          assembly = ["wenger", "garg"]),
        dipcall = expand("pipeline/{assembly}/truvari/calls/dipcall/error_types.txt", 
                          assembly = ["wenger", "garg"])
    output:
        all = "pipeline/eval/all_results_error_types.txt"
    threads: 1
    shell:
        "cat {input.svim} {input.dipcall} > {output.all}"

rule cat_truvari_results_svim_parameters:
    input:
        svim = expand("pipeline/{assembly}/truvari/calls/svim/{parameters}/error_types.txt", 
                          assembly = ["wenger", "garg"],
                          parameters = ["50_50_50_50_200", "100_100_2000_2000_200", "1000_1000_2000_2000_200"])
    output:
        all = "pipeline/eval/svim_parameter_results.txt"
    threads: 1
    run:
        with open(output.all, 'w') as output_file:
            for f in input.svim:
                parameters = f.split("/")[5]
                with open(f, 'r') as input_file:
                    for line in input_file:
                        print("%s\t%s" % (parameters, line.strip()), file=output_file)

rule plot_pr_tools:
    input:
        "pipeline/eval/all_results_error_types.txt"
    output:
        png = "pipeline/eval/results.tools.png",
        tsv = "pipeline/eval/results.tools.tsv"
    threads: 1
    log:
        "pipeline/logs/rplot.tools.log"
    shell:
        "Rscript --vanilla workflow/scripts/plot_tools.R {input} {output.png} {output.tsv} > {log}"

rule plot_pr_tools_error_types:
    input:
        "pipeline/eval/all_results_error_types.txt"
    output:
        png = "pipeline/eval/results.tools.error_types.png",
        tsv = "pipeline/eval/results.tools.error_types.tsv"
    threads: 1
    log:
        "pipeline/logs/rplot.tools.error_types.log"
    shell:
        "Rscript --vanilla workflow/scripts/plot_tools_error_types.R {input} {output.png} {output.tsv} > {log}"



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


#SV counts from all classes

rule call_svim_diploid_all_types:
    input:
        reference = config["reference"],
        bam1 = "pipeline/{assembly}/alignments/H1.sort.bam",
        index1 = "pipeline/{assembly}/alignments/H1.sort.bam.bai",
        bam2 = "pipeline/{assembly}/alignments/H2.sort.bam",
        index2 = "pipeline/{assembly}/alignments/H2.sort.bam.bai"
    output:
        "pipeline/{assembly}/calls/svim/{rgt}_{rot}_{qgt}_{qot}_{med}_all_types/variants.vcf"
    params:
        working_dir = "pipeline/{assembly}/calls/svim/{rgt}_{rot}_{qgt}_{qot}_{med}_all_types"
    conda:
        "../envs/svimasm.yaml"
    shell:
        "svim-asm diploid {params.working_dir} {input.bam1} {input.bam2} {input.reference} --min_sv_size 20 \
        --reference_gap_tolerance {wildcards.rgt} --reference_overlap_tolerance {wildcards.rot} \
        --query_gap_tolerance {wildcards.qgt} --query_overlap_tolerance {wildcards.qot} --max_edit_distance {wildcards.med} --sample HG002 --query_names"

rule SV_length_plot_svim:
    input:
        "pipeline/{assembly}/calls/svim/{parameters}_all_types/variants.vcf"
    output:
        plot = "pipeline/{assembly}/SV-plots/SV-length_SVIM_{parameters}.png",
        counts = "pipeline/{assembly}/SV-plots/SV-counts_SVIM_{parameters}.txt",
    log:
        "logs/svplot/svlength_{assembly}_SVIM_{parameters}.log"
    conda:
        "../envs/cyvcf2.yaml"
    shell:
        "python workflow/scripts/SV-length-plot.py {input} --output {output.plot} --counts {output.counts} --filter 'GL000207.1,GL000226.1,GL000229.1,GL000231.1,GL000210.1,GL000239.1,GL000235.1,GL000201.1,GL000247.1,GL000245.1,GL000197.1,GL000203.1,GL000246.1,GL000249.1,GL000196.1,GL000248.1,GL000244.1,GL000238.1,GL000202.1,GL000234.1,GL000232.1,GL000206.1,GL000240.1,GL000236.1,GL000241.1,GL000243.1,GL000242.1,GL000230.1,GL000237.1,GL000233.1,GL000204.1,GL000198.1,GL000208.1,GL000191.1,GL000227.1,GL000228.1,GL000214.1,GL000221.1,GL000209.1,GL000218.1,GL000220.1,GL000213.1,GL000211.1,GL000199.1,GL000217.1,GL000216.1,GL000215.1,GL000205.1,GL000219.1,GL000224.1,GL000223.1,GL000195.1,GL000212.1,GL000222.1,GL000200.1,GL000193.1,GL000194.1,GL000225.1,GL000192.1,NC_007605,hs37d5'  --tool SVIM 2> {log}"

rule SV_length_plot_dipcall:
    input:
        "pipeline/{assembly}/calls/dipcall/HG002.dip.vcf.gz"
    output:
        plot = "pipeline/{assembly}/SV-plots/SV-length_dipcall.png",
        counts = "pipeline/{assembly}/SV-plots/SV-counts_dipcall.txt",
    log:
        "logs/svplot/svlength_{assembly}_dipcall.log"
    conda:
        "../envs/cyvcf2.yaml"
    shell:
        "python workflow/scripts/SV-length-plot.py {input} --output {output.plot} --counts {output.counts} --filter 'GL000207.1,GL000226.1,GL000229.1,GL000231.1,GL000210.1,GL000239.1,GL000235.1,GL000201.1,GL000247.1,GL000245.1,GL000197.1,GL000203.1,GL000246.1,GL000249.1,GL000196.1,GL000248.1,GL000244.1,GL000238.1,GL000202.1,GL000234.1,GL000232.1,GL000206.1,GL000240.1,GL000236.1,GL000241.1,GL000243.1,GL000242.1,GL000230.1,GL000237.1,GL000233.1,GL000204.1,GL000198.1,GL000208.1,GL000191.1,GL000227.1,GL000228.1,GL000214.1,GL000221.1,GL000209.1,GL000218.1,GL000220.1,GL000213.1,GL000211.1,GL000199.1,GL000217.1,GL000216.1,GL000215.1,GL000205.1,GL000219.1,GL000224.1,GL000223.1,GL000195.1,GL000212.1,GL000222.1,GL000200.1,GL000193.1,GL000194.1,GL000225.1,GL000192.1,NC_007605,hs37d5'  --tool DipCall 2> {log}"

rule merge_counts:
    input:
        svim = "pipeline/{assembly}/SV-plots/SV-counts_SVIM_1000_1000_2000_2000_200.txt",
        sniffles = "pipeline/{assembly}/SV-plots/SV-counts_dipcall.txt",
    output:
        "pipeline/{assembly}/SV-plots/SV-counts.merged.txt"
    shell:
        "cat {input} | grep -v '#' > {output}"

rule plot_counts:
    input:
        "pipeline/{assembly}/SV-plots/SV-counts.merged.txt"
    output:
        png = "pipeline/{assembly}/SV-plots/SV-counts.merged.png",
        tsv = "pipeline/{assembly}/SV-plots/SV-counts.merged.tsv"
    shell:
        "Rscript --vanilla workflow/scripts/plot_counts.R {input} {output.png} {output.tsv}"
