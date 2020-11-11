import sys
from collections import Counter
import matplotlib.pyplot as plt
import numpy as np
import re

def read_variants(path, step_size_small = 50, step_size_large = 500, limit = 1000):
    ins_counter = Counter()
    del_counter = Counter()
    open_file = open(path, "r")
    for line in open_file:
        if not line.startswith("#"):
            fields = line.strip().split("\t")
            ref_length = len(fields[3])
            alt_length = len(fields[4].split(",")[0])
            match_len = re.search(r"SVLEN=(-?[0-9]+)", fields[7])
            if match_len:
                sv_len = int(match_len.group(1))
            else:
                sv_len = alt_length - ref_length
            if sv_len > 0:
                if sv_len >= limit:
                    bucket_lower_bound = sv_len // step_size_large * step_size_large #round down to next step_size_large
                    bucket_name = "[{0},{1})".format(bucket_lower_bound, bucket_lower_bound + step_size_large)
                    ins_counter[bucket_name] += 1
                else:
                    bucket_lower_bound = sv_len // step_size_small * step_size_small #round down to next step_size_small
                    bucket_name = "[{0},{1})".format(bucket_lower_bound, bucket_lower_bound + step_size_small)
                    ins_counter[bucket_name] += 1
            if sv_len < 0:
                abs_sv_len = abs(sv_len)
                if abs_sv_len >= limit:
                    bucket_lower_bound = abs_sv_len // step_size_large * step_size_large #round down to next step_size_large
                    bucket_name = "[{0},{1})".format(bucket_lower_bound, bucket_lower_bound + step_size_large)
                    del_counter[bucket_name] += 1
                else:
                    bucket_lower_bound = abs_sv_len // step_size_small * step_size_small #round down to next step_size_small
                    bucket_name = "[{0},{1})".format(bucket_lower_bound, bucket_lower_bound + step_size_small)
                    del_counter[bucket_name] += 1
    open_file.close()
    return (del_counter, ins_counter)

def forward(x):
    output = []
    for e in x:
        if e > 1000:
            output.append(1000 + (e-1000)/10)
        elif e < -1000:
            output.append(-1000 + (e+1000)/10)
        else:
            output.append(e)
    return np.array(output)

def inverse(x):
    output = []
    for e in x:
        if e > 1000:
            output.append(1000 + (e-1000)*10)
        elif e < -1000:
            output.append(-1000 + (e+1000)*10)
        else:
            output.append(e)
    return np.array(output)

def plot_length(ax, xs, precs, recs, calls_lt1000, calls_gt1000):
    ax.set_xlabel('SV length (bp)')
    ax.set_ylabel('Value')    
    ax.set_xlim(-10001, 10001)
    ax.set_ylim(0.0, 1.1)
    ax.set_xscale('function', functions=(forward, inverse))
    ax.set_xticks(np.concatenate([np.arange(-10000, -1000, 2500), np.arange(-1000, 0, 250), np.arange(0, 1001, 250), np.arange(2500, 10001, 2500)]))
    ax.set_yticks(np.arange(0.0, 1.05, 0.1))
    ax.tick_params(axis='x', labelrotation=320)

    ax2 = ax.twinx()
    ax2.set_ylabel('Variant calls')

    ax2.bar(xs, calls_lt1000, width = 40, label='Calls', color="grey")
    ax2.bar(xs, calls_gt1000, width = 400, label='Calls', color="grey")
    ax.plot(xs, precs, label='Precision')
    ax.plot(xs, recs, label='Recall')

    #ax.axvline(x=0)
    ax.grid()
    ax.legend()

def main():
    truvari_dir = sys.argv[1]
    tpcall_path = truvari_dir + "/tp-call.vcf"
    tpbase_path = truvari_dir + "/tp-base.vcf"
    fp_path = truvari_dir + "/fp.vcf"
    fn_path = truvari_dir + "/fn.vcf"

    step_size_small = 50
    step_size_large = 500
    limit = 1000
    tp_call_dels, tp_call_ins = read_variants(tpcall_path, step_size_small, step_size_large, limit)
    tp_base_dels, tp_base_ins = read_variants(tpbase_path, step_size_small, step_size_large, limit)
    fp_dels, fp_ins = read_variants(fp_path, step_size_small, step_size_large, limit)
    fn_dels, fn_ins = read_variants(fn_path, step_size_small, step_size_large, limit)

    buckets = ["[{0},{1})".format(i, i + step_size_small) for i in range(0, 1000, step_size_small)] + \
              ["[{0},{1})".format(i, i + step_size_large) for i in range(1000, 10000, step_size_large)]

    print("DELETIONS")
    xs = []
    precs = []
    recs = []
    calls_gt1000 = []
    calls_lt1000 = []
    print("prec\trec\trange\tTP_call\tTP_base\tFP\tFN")
    for b in reversed(buckets):
        if fp_dels[b] + tp_call_dels[b] > 3:
            prec = round(tp_call_dels[b] / (fp_dels[b] + tp_call_dels[b]), 2)
        else:
            prec = None
        if fn_dels[b] + tp_base_dels[b] > 3:
            rec = round(tp_base_dels[b] / (fn_dels[b] + tp_base_dels[b]), 2)
        else:
            rec = None
        limits = b[1:-1].split(",")
        xs.append((int(limits[1]) + int(limits[0])) / -2)
        precs.append(prec)
        recs.append(rec)
        if int(limits[0]) < 1000: 
            calls_lt1000.append(fp_dels[b] + tp_call_dels[b])
            calls_gt1000.append(0)
        else:
            calls_lt1000.append(0)
            calls_gt1000.append(fp_dels[b] + tp_call_dels[b])
        print("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}".format(prec, rec, b, tp_call_dels[b], tp_base_dels[b], fp_dels[b], fn_dels[b]))

    print("INSERTIONS")
    for b in buckets:
        if fp_ins[b] + tp_call_ins[b] > 3:
            prec = round(tp_call_ins[b] / (fp_ins[b] + tp_call_ins[b]), 2)
        else:
            prec = None
        if fn_ins[b] + tp_base_ins[b] > 3:
            rec = round(tp_base_ins[b] / (fn_ins[b] + tp_base_ins[b]), 2)
        else:
            rec = None
        limits = b[1:-1].split(",")
        xs.append((int(limits[1]) + int(limits[0])) / 2)
        precs.append(prec)
        recs.append(rec)
        if int(limits[0]) < 1000: 
            calls_lt1000.append(fp_ins[b] + tp_call_ins[b])
            calls_gt1000.append(0)
        else:
            calls_lt1000.append(0)
            calls_gt1000.append(fp_ins[b] + tp_call_ins[b])
        print("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}".format(prec, rec, b, tp_call_ins[b], tp_base_ins[b], fp_ins[b], fn_ins[b]))

    fig, (ax1) = plt.subplots(1, 1, figsize=(6,3))
    plot_length(ax1, xs, precs, recs, calls_lt1000, calls_gt1000)

    plt.tight_layout()
    output_path = sys.argv[2]
    fig.savefig(output_path)
    fig.clf()

if __name__ == "__main__":
    sys.exit(main())
