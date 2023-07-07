#!/usr/bin/env python3

import sys
import statistics
from pysam import VariantFile
from matplotlib import pyplot as plt
import seaborn as sns
import pandas as pd

import warnings

warnings.simplefilter(action="ignore", category=FutureWarning)


def main():
    in_prefix = sys.argv[1]
    out_prefix = sys.argv[2]
    tau_t = 0.7

    TRUTH = {}
    CALL = {}
    REGIONS = {}

    TREGIONS = set()
    CREGIONS = set()
    for line in open(in_prefix + ".regions.bed"):
        chrom, s, e, t, l = line.strip("\n").split("\t")
        region = f"{chrom}:{s}-{e}"
        if t == "T":
            TREGIONS.add(region)
        elif t == "C":
            CREGIONS.add(region)

    for line in open(in_prefix + ".results.txt"):
        t, idx, region, score = line.strip("\n").split(" ")
        score = float(score)
        if score < 0:
            score = 0
            # continue
        if t == "T":
            TRUTH[idx] = float(score)
        elif t == "C":
            CALL[idx] = float(score)
        else:
            REGIONS[region] = float(score)

    df = []
    for tau in range(0, 101, 1):
        tau = tau / 100
        R = sum([x >= tau for x in TRUTH.values()]) / len(TRUTH)
        P = sum([x >= tau for x in CALL.values()]) / len(CALL)
        regP = sum([x >= tau for x in REGIONS.values()]) / len(CREGIONS)
        regR = sum([x >= tau for x in REGIONS.values()]) / len(TREGIONS)
        df.append([tau, "Precision", P])  # "Correcteness"
        df.append([tau, "Recall", R])  # "Completeness"
        df.append([tau, "Precision (Regions)", regP])
        df.append([tau, "Recall (Regions)", regR])
    df = pd.DataFrame(df, columns=["tau", "Metric", "Value"])
    print(df[df["tau"] == tau_t])

    sns.lineplot(df, x="tau", y="Value", hue="Metric")

    plt.xlim(0.9, 1.05)
    plt.ylim(-0.05, 1.05)

    plt.tight_layout()
    plt.savefig(out_prefix + ".plot.png")


if __name__ == "__main__":
    main()
