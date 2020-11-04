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

ASSEMBLIES = ["garg", "wenger"]

##### Target rules #####

rule all:
    input:
        "pipeline/eval/results.tools.png",
        "pipeline/eval/results.tools.error_types.png",
        "pipeline/eval/results.svim.parameters.png",
        expand("pipeline/{assembly}/SV-plots/SV-counts.merged.png", assembly=ASSEMBLIES),
        expand("pipeline/{assembly}/SV-plots/SV-length_SVIM_{rgt_rot_qgt_qot}_{med}.png", rgt_rot_qgt_qot=["50_50_50_50", "100_100_2000_2000", "1000_1000_2000_2000"], med=[200], assembly=ASSEMBLIES),
        expand("pipeline/{assembly}/SV-plots/SV-length_dipcall.png", assembly=ASSEMBLIES),
        expand("pipeline/{assembly}/{eval}/calls/dipcall/performance_by_size.png", assembly=ASSEMBLIES, eval=["truvari", "truvari_multimatch"]),
        expand("pipeline/{assembly}/{eval}/calls/svim/{rgt_rot_qgt_qot}_{med}/performance_by_size.png", assembly=ASSEMBLIES, eval=["truvari", "truvari_multimatch"], rgt_rot_qgt_qot=["50_50_50_50", "100_100_2000_2000", "1000_1000_2000_2000"], med=[200]),
        expand("pipeline/{assembly}/{eval}/calls/dipcall/size_reps.png", assembly=ASSEMBLIES, eval=["truvari", "truvari_multimatch"]),
        expand("pipeline/{assembly}/{eval}/calls/svim/{rgt_rot_qgt_qot}_{med}/size_reps.png", assembly=ASSEMBLIES, eval=["truvari", "truvari_multimatch"], rgt_rot_qgt_qot=["50_50_50_50", "100_100_2000_2000", "1000_1000_2000_2000"], med=[200])
