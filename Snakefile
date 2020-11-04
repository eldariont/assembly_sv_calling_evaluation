configfile: "config/config.yaml"

wildcard_constraints:
    chr="(\d+|genome)",
    rgt="\d+",
    rot="\d+",
    qgt="\d+",
    qot="\d+",
    med="\d+"

include: "workflow/rules/align.smk"
include: "workflow/rules/call.smk"
include: "workflow/rules/eval.smk"
include: "workflow/rules/plot.smk"

##### Target rules #####

rule all:
    input:
        # expand("pipeline/{eval}/calls/dipcall/summary.txt", eval=["truvari", "truvari_multimatch"]),
        # expand("pipeline/{eval}/calls/svim/{rgt_rot_qgt_qot}_{med}/summary.txt", eval=["truvari", "truvari_multimatch"], rgt_rot_qgt_qot=["50_50_50_50", "100_100_2000_2000", "1000_1000_2000_2000"], med=[200]),
        expand("pipeline/{assembly}/{eval}/calls/dipcall/size_reps.png", assembly=["garg", "wenger"], eval=["truvari", "truvari_multimatch"]),
        expand("pipeline/{assembly}/{eval}/calls/svim/{rgt_rot_qgt_qot}_{med}/size_reps.png", assembly=["garg", "wenger"], eval=["truvari", "truvari_multimatch"], rgt_rot_qgt_qot=["50_50_50_50", "100_100_2000_2000", "1000_1000_2000_2000"], med=[200]),
        expand("pipeline/{assembly}/{eval}/calls/dipcall/performance_by_size.png", assembly=["garg", "wenger"], eval=["truvari", "truvari_multimatch"]),
        expand("pipeline/{assembly}/{eval}/calls/svim/{rgt_rot_qgt_qot}_{med}/performance_by_size.png", assembly=["garg", "wenger"], eval=["truvari", "truvari_multimatch"], rgt_rot_qgt_qot=["50_50_50_50", "100_100_2000_2000", "1000_1000_2000_2000"], med=[200]),
        "pipeline/eval/results.tools.png",
        "pipeline/eval/results.tools.error_types.png",
        "pipeline/eval/results.tools_thesis.png",
        "pipeline/eval/results.svim.parameters.png",
        expand("pipeline/{assembly}/SV-plots/SV-counts.merged.png", assembly=["garg", "wenger"]),
        expand("pipeline/{assembly}/SV-plots/SV-length_SVIM_{rgt_rot_qgt_qot}_{med}.png", rgt_rot_qgt_qot=["50_50_50_50", "100_100_2000_2000", "1000_1000_2000_2000"], med=[200], assembly=["garg", "wenger"]),
        expand("pipeline/{assembly}/SV-plots/SV-length_dipcall.png", assembly=["garg", "wenger"]),
        # expand("coverage/assembly.H{hap}.cov.tsv", hap=["1.flt","1", "2.flt", "2"]),
        # "coverage/reads.cov.tsv"


###########
# GENERAL #
###########

# rule compute_coverage_assemblies:
#     input:
#         bam = "alignments/genome.H{hap}.sort.bam"
#     output:
#         "coverage/assembly.H{hap}.cov.tsv"
#     shell:
#         "samtools view -q 20 -b {input.bam} | bedtools genomecov -ibam stdin > {output}"


# rule compute_coverage_reads:
#     input:
#         bam = config["reads_bam"]
#     output:
#         "coverage/reads.cov.tsv"
#     shell:
#         "samtools view -q 20 -b {input.bam} | bedtools genomecov -ibam stdin > {output}"
