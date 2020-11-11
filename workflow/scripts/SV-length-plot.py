#! /usr/bin/env python
import sys
import re
from cyvcf2 import VCF
from collections import defaultdict
from argparse import ArgumentParser
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def main():
    args = get_args()
    len_dict = defaultdict(list)
    ignored_chroms = args.filter.split(",")
    for v in VCF(args.vcf):
        #filter by score
        if args.min_score != 0 and v.QUAL < args.min_score:
            continue
        #filter out HOM_REFs
        if len(v.gt_types) > 0 and v.gt_types[0] == 0:
            continue
        #filter by chromosome
        if v.CHROM in ignored_chroms:
            continue
        #filter out filtered records
        if v.FILTER in ["Decoy", "NearReferenceGap", "NearContigEnd", "hom_ref"]:
            continue
        #ignore TRAs and BNDs
        if v.INFO.get('SVTYPE') in ['TRA', 'BND']:
            if v.INFO.get('SVTYPE') == 'BND':
                alt_string = v.ALT[0]
                if alt_string.startswith('N'):
                    pos2_string = alt_string[2:-1]
                else:
                    pos2_string = alt_string[1:-2]
                chr2 = pos2_string.split(":")[0]
                if chr2 in ignored_chroms:
                    continue
            len_dict[v.INFO.get('SVTYPE')].append(0)
            continue
        try:
            #filter by SV length
            if abs(v.INFO.get('SVLEN')) >= 50:
                len_dict[v.INFO.get('SVTYPE')].append(abs(v.INFO.get('SVLEN')))
        except TypeError:
            #filter by inversion length
            if v.INFO.get('SVTYPE') == 'INV':
                if (v.end - v.start) >= 50:
                    len_dict[v.INFO.get('SVTYPE')].append(v.end - v.start)
                    sys.stderr.write("SVLEN field missing. Inferred SV length from END and POS:\n{}\n\n".format(v))
            elif re.match("^[ACGTN]+$", v.REF) != None and re.match("^[ACGTN]+$", v.ALT[0]) != None:
                svlen = len(v.ALT[0]) - len(v.REF)
                if svlen >= 50:
                    len_dict["INS"].append(svlen)
                elif svlen <= -50:
                    len_dict["DEL"].append(abs(svlen))
            else:
                sys.stderr.write("Exception when parsing variant:\n{}\n\n".format(v))
    with open(args.counts, 'w') as counts:
        counts.write("#Size class\tType\tCount\tCumulative length\tTool\n")
        for svtype, lengths in len_dict.items():
            if svtype in ['TRA', 'BND']:
                counts.write("{}\t{}\t{}\t{}\t{}\n".format("all", svtype, len(lengths), 0, args.tool))
            else:
                tiny = [l for l in lengths if l < 200]
                small = [l for l in lengths if l >= 200 and l < 1000]
                medium = [l for l in lengths if l >= 1000 and l < 10000]
                large = [l for l in lengths if l >= 10000 and l < 100000]
                huge = [l for l in lengths if l >= 100000]
                counts.write("{}\t{}\t{}\t{}\t{}\n".format("all", svtype, len(lengths), sum(lengths), args.tool))
                counts.write("{}\t{}\t{}\t{}\t{}\n".format("tiny", svtype, len(tiny), sum(tiny), args.tool))
                counts.write("{}\t{}\t{}\t{}\t{}\n".format("small", svtype, len(small), sum(small), args.tool))
                counts.write("{}\t{}\t{}\t{}\t{}\n".format("medium", svtype, len(medium), sum(medium), args.tool))
                counts.write("{}\t{}\t{}\t{}\t{}\n".format("large", svtype, len(large), sum(large), args.tool))
                counts.write("{}\t{}\t{}\t{}\t{}\n".format("huge", svtype, len(huge), sum(huge), args.tool))
    make_plot(dict_of_lengths=len_dict,
              output=args.output)


def make_plot(dict_of_lengths, output):
    """Makes two stacked bar charts
    Plotting two bar charts of number of SVs by length split by SV type
    Use a consistent colouring scheme for those in "standard_order" to
    make comparison reasonable

    First bar chart is up to 2kb with bins of 10bp
    Second bar chart is up to 20kb, with bins of 100bp
     and uses log scaling on the y-axis
    """
    try:
        del dict_of_lengths["TRA"]
    except KeyError:
        pass
    try:
        del dict_of_lengths["BND"]
    except KeyError:
        pass
    try:
        dict_of_lengths["DUP:TANDEM"] = dict_of_lengths.pop("DUP")
    except KeyError:
        pass
    try:
        dict_of_lengths["DUP:INT"] = dict_of_lengths.pop("cnv")
    except KeyError:
        pass
    try:
        dict_of_lengths["DUP:INT"] = dict_of_lengths.pop("DUP_INT")
    except KeyError:
        pass
    try:
        dict_of_lengths["Complex"] = dict_of_lengths.pop("DEL/INV")
    except KeyError:
        pass
    try:
        dict_of_lengths["Complex"] = dict_of_lengths.pop("DUP/INS")
    except KeyError:
        pass
    standard_order = ['DEL', 'INS', 'INV', 'DUP:TANDEM', 'DUP:INT', 'BND', 'Complex']
    if len(dict_of_lengths.keys()) > 0:
        spec_order = sorted([i for i in dict_of_lengths.keys() if i not in standard_order])
        sorter = standard_order + spec_order
        names, lengths = zip(
            *sorted([(svtype, lengths) for svtype, lengths in dict_of_lengths.items()],
                    key=lambda x: sorter.index(x[0])))
    else:
        names = []
        lengths = []
    plt.figure(num=None, figsize=(9, 3), dpi=300)
    plt.subplot(1, 2, 1)
    plt.hist(x=lengths,
             bins=[i for i in range(0, 2000, 10)],
             stacked=True,
             histtype='bar',
             label=names)
    plt.ylim(1, 3200)
    plt.xlabel('Length of structural variant')
    plt.ylabel('Number of variants')
    plt.legend(frameon=False,
               fontsize="small")

    plt.subplot(1, 2, 2)
    plt.hist(x=lengths,
             bins=[i for i in range(0, 20000, 100)],
             stacked=True,
             histtype='bar',
             label=names,
             log=True)
    plt.ylim(1, 10200)
    plt.xlabel('Length of structural variant')
    plt.ylabel('Number of variants')
    plt.legend(frameon=False,
               fontsize="small")
    plt.tight_layout()
    plt.savefig(output)


def get_args():
    parser = ArgumentParser(description="create stacked bar plot of the SV lengths split by type")
    parser.add_argument("vcf", help="vcf file to parse")
    parser.add_argument("-m", "--min_score", 
                        type=int,
                        default=0,
                        help="filter out calls with lower score (default: 0, do not filter)")
    parser.add_argument("-o", "--output",
                        help="output file to write figure to",
                        default="SV-length.png")
    parser.add_argument("-c", "--counts",
                        help="output file to write counts to",
                        default="SV-length.txt")
    parser.add_argument("-t", "--tool",
                        help="tool name",
                        default="Unknown")
    parser.add_argument("-f", "--filter",
                        type=str,
                        help="comma-separated list of chromosomes to filter out (e.g. 'hs37d5')",
                        default="")
    return parser.parse_args()


if __name__ == '__main__':
    main()
