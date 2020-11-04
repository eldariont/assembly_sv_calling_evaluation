import sys
import collections

###############
# Analyze TPs #
###############

with open(sys.argv[1], 'r') as tps_call:
    call_lines = [line.strip() for line in tps_call.readlines() if not line.startswith("#")]

with open(sys.argv[2], 'r') as tps_base:
    base_lines = [line.strip() for line in tps_base.readlines() if not line.startswith("#")]

assert (len(call_lines) == len(base_lines)), "Line numbers different: {0} != {1}".format(len(call_lines), len(base_lines))

call_pairs = zip(call_lines, base_lines)
matching_gt = 0
different_gt = 0
genotypes = collections.Counter()
het_genotypes = ["1/0", "1|0", "1/.", "1|.", "0/1", "0|1", "./1", ".|1"]
hom_var_genotypes = ["1/1", "1|1"]
hom_ref_genotypes = ["./.", ".|.", "./0", ".|0", "0/.", "0|."]
for call_call, base_call in call_pairs:
    call_fields = call_call.split()
    base_fields = base_call.split()

    assert (call_fields[0] == base_fields[0]), "Chromosomes different: {0} != {1}".format(call_fields[0], base_fields[0])
    assert (abs(int(call_fields[1]) - int(base_fields[1])) < 100000), "Coordinates different: {0} != {1}".format(call_fields[1], base_fields[1])

    #Parse genotype of call
    call_gt_string = call_fields[9]
    call_gt = call_gt_string.strip().split(":")[0]
    call_gt_fields = call_gt.replace("|", "/").split("/")
    assert (len(call_gt_fields) == 2)

    if call_gt_fields[0] in [".", "0"] and call_gt_fields[1] in [".", "0"]:
        all_gt_option = "hom_ref"
    elif call_gt_fields[0] in [".", "0"] or call_gt_fields[1] in [".", "0"]:
        call_gt_option = "het"
    else:
        call_gt_option = "hom_var"

    #Parse genotype of base
    base_gt_string = base_fields[9]
    base_gt = base_gt_string.strip().split(":")[0]
    base_gt_fields = base_gt.replace("|", "/").split("/")
    assert (len(base_gt_fields) == 2)

    if base_gt_fields[0] in [".", "0"] and base_gt_fields[1] in [".", "0"]:
        all_gt_option = "hom_ref"
    elif base_gt_fields[0] in [".", "0"] or base_gt_fields[1] in [".", "0"]:
        base_gt_option = "het"
    else:
        base_gt_option = "hom_var"

    if call_gt_option == base_gt_option:
        matching_gt += 1
    else:
        different_gt += 1
        genotypes[call_gt_option+":"+base_gt_option] += 1

print(genotypes.most_common(), file=sys.stderr)
#######################
# Analyze FPs and FNs #
#######################

with open(sys.argv[3], 'r') as fps:
    fp_lines = [line.strip() for line in fps.readlines() if not line.startswith("#")]
with open(sys.argv[4], 'r') as fns:
    fn_lines = [line.strip() for line in fns.readlines() if not line.startswith("#")]

assembly = sys.argv[5]
caller = sys.argv[6]

print("\t".join([assembly, caller, str(matching_gt), str(different_gt), str(len(fp_lines)), str(len(fn_lines))]))
