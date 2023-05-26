#!/usr/bin/env python3

import sys
import statistics
from pysam import VariantFile
from matplotlib import pyplot as plt
import seaborn as sns
import pandas as pd


def main():
    truth_vcf_path = sys.argv[1]
    call_vcf_path = sys.argv[2]
    scores_path = sys.argv[3]
    # tau = float(sys.argv[4])

    TRUTH = {}
    TRUTH_L = {}
    for rec in VariantFile(truth_vcf_path):
        idx = rec.id
        TRUTH[idx] = 0
        TRUTH_L[idx] = (
            rec.info["SVLEN"]
            if type(rec.info["SVLEN"]) == int
            else statistics.fmean(rec.info["SVLEN"])
        )
    CALL = {}
    CALL_L = {}
    for rec in VariantFile(call_vcf_path):
        idx = rec.id
        CALL[idx] = 0
        CALL_L[idx] = (
            rec.info["SVLEN"]
            if type(rec.info["SVLEN"]) == int
            else statistics.fmean(rec.info["SVLEN"])
        )

    for line in open(scores_path):
        t, idx, score = line.strip("\n").split(" ")
        if t == "T":
            TRUTH[idx] = float(score)
        else:
            CALL[idx] = float(score)

    df = []
    for tau in range(0, 101, 10):
        tau = tau / 100
        R = sum([x >= tau for x in TRUTH.values()]) / len(TRUTH)
        P = sum([x >= tau for x in CALL.values()]) / len(CALL)
        print(sum([x >= tau for x in CALL.values()]), len(CALL))
        F1 = 2 * P * R / (P + R)
        df.append([tau, P, R, F1])
    df = pd.DataFrame(df, columns=["tau", "P", "R", "F1"])
    print(df)

    # Rdf = pd.DataFrame([[TRUTH[k], TRUTH_L[k]] for k in TRUTH], columns=["Score", "L"])
    # Pdf = pd.DataFrame([[CALL[k], CALL_L[k]] for k in CALL], columns=["Score", "L"])

    # fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2)
    # sns.histplot(
    #     data=Rdf, x="Score", ax=ax1, element="step", fill=False, cumulative=True
    # )
    # sns.histplot(
    #     data=Pdf, x="Score", ax=ax2, element="step", fill=False, cumulative=True
    # )

    # sns.histplot(Rdf, x="Score", y="L", ax=ax3)
    # sns.histplot(Pdf, x="Score", y="L", ax=ax4)

    # ax1.set_title("R")
    # ax2.set_title("P")
    # plt.tight_layout()
    # # plt.savefig("plot.png")
    # plt.savefig(scores_path + ".png")


if __name__ == "__main__":
    main()
