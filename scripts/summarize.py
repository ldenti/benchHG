#!/usr/bin/env python3

import sys
import statistics
from pysam import VariantFile
from matplotlib import pyplot as plt
import seaborn as sns
import pandas as pd


def main():
    prefix_path = sys.argv[1]

    TRUTH = {}
    CALL = {}
    REGIONS = {}

    TREGIONS = set()
    CREGIONS = set()
    for line in open(prefix_path + ".regions.bed"):
        chrom, s, e, t, l = line.strip("\n").split("\t")
        region = f"{chrom}:{s}-{e}"
        if t == "T":
            TREGIONS.add(region)
        elif t == "C":
            CREGIONS.add(region)

    for line in open(prefix_path + ".results.txt"):
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
        df.append([tau, "Correcteness", P])
        df.append([tau, "Completeness", R])
        df.append([tau, "Precision (Regions)", regP])
        df.append([tau, "Recall (Regions)", regR])
    df = pd.DataFrame(df, columns=["tau", "Metric", "Value"])
    # print(df)

    sns.lineplot(df, x="tau", y="Value", hue="Metric")

    plt.tight_layout()
    plt.savefig("plot.png")


if __name__ == "__main__":
    main()
