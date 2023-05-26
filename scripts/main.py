#!/usr/bin/env python3

import sys
import os
import argparse
import logging

from Bio import SeqIO
import networkx as nx
import pandas as pd
from matplotlib import pyplot as plt
import seaborn as sns
from pysam import VariantFile
import Levenshtein as lev
from intervaltree import IntervalTree
from rich.console import Console
from rich.logging import RichHandler
from rich.progress import track

console = Console(stderr=True)

logging.basicConfig(
    format="%(message)s",
    level=logging.INFO,
    datefmt="[%x %X]",
    handlers=[RichHandler(console=console)],
)


def find_unique_kmer(seq, k, from_start):
    KMERS = {}
    for i in range(len(seq) - k + 1):
        kmer = seq[i : i + k]
        if kmer not in KMERS:
            KMERS[kmer] = i
        else:
            KMERS[kmer] = -1
    for kmer in KMERS if from_start else reversed(KMERS):
        if KMERS[kmer] != -1:
            return KMERS[kmer] if from_start else KMERS[kmer] + k
    return 0 if from_start else len(seq)


class SubHaplo:
    def __init__(self):
        self.region = ""
        self.seq = ""
        self.records = []
        self.genotypes = []


class HaploBlock:
    def __init__(self):
        self.records = []
        self.chrom = ""
        self.start = -1  # 0-based inclusive
        self.end = -1  # 0-based exclusive
        self.ustart = -1  # 0-based inclusive
        self.uend = -1  # 0-based exclusive
        self.candidates_haplotypes = []
        self.candidates_haplotypes_sequences = []

    def append(self, record):
        self.records.append(record)
        self.chrom = record.chrom  # assuming parent is checking for same chromosome
        self.start = min([r.start for r in self.records])
        self.end = max([r.stop for r in self.records])

    def region(self):
        return f"{self.chrom}:{self.start}-{self.end}"

    def uregion(self):
        return f"{self.chrom}:{self.ustart}-{self.uend}"

    def __len__(self):
        return len(self.records)

    def __iter__(self):
        return self.records

    def set_extension(self, us, ue):
        self.ustart = us
        self.uend = ue

    def extend(self, seq, klen, wsize):
        # Extend the region to unique kmers
        unique_kmer_s = find_unique_kmer(
            seq[self.start - wsize - 1 : self.start], klen, True
        )
        self.ustart = self.start - (wsize + unique_kmer_s)
        unique_kmer_e = find_unique_kmer(seq[self.end : self.end + wsize], klen, False)
        self.uend = self.end + unique_kmer_e

    def get_candidate(self, i, j):
        return [
            x for x in self.candidates_haplotypes[i][j] if x[0] != -1
        ], self.candidates_haplotypes_sequences[i][j]

    def build_candidates(self, seq, truth=False):
        # To understand if two SVs are compatible (no overlaps), we build a graph where each node is a SV and two SVs are linked if they do not overlap on the reference. Then we remove transitive edges and finally each path in the graph (with potentially multiple weakly connected components) is a potential haplotype
        G = nx.DiGraph()
        for i, record in enumerate(self.records):
            G.add_node(i)
            for node in G:
                if node == i:
                    continue
                if self.records[node].stop <= record.start:
                    G.add_edge(node, i)
        G = nx.transitive_reduction(G)

        # Any isolated node is a potential haplotype
        self.haplos = [
            [n] for n in G.nodes() if G.in_degree(n) == 0 and G.out_degree(n) == 0
        ]
        sources = [n for n in G.nodes() if G.in_degree(n) == 0 and G.out_degree(n) > 0]
        targets = [n for n in G.nodes() if G.in_degree(n) > 0 and G.out_degree(n) == 0]
        # Any path from a source to a target is a potential haplotype
        for s in sources:
            for t in targets:
                for path in nx.all_simple_paths(G, source=s, target=t):
                    self.haplos.append(path)

        self.graphs = []
        if len(self.haplos) > 2:
            print("More than 2 components in", self.uregion())
            return

        for haplo in self.haplos:
            # print("Haplotype:", haplo)
            G = nx.DiGraph()
            n = 0
            last_p = self.ustart
            G.add_node(n, v=-1, seq=seq[last_p : self.records[0].start].upper(), gt=0)
            prepre_nodes = [0]
            pre_nodes = [0]
            n += 1
            for i in haplo:
                record = self.records[i]
                curr_nodes = []
                if len(pre_nodes) == 1:
                    # previous node is reference
                    # print("A: Predecessor of node", i, "is reference")
                    for a, all_seq in enumerate(record.alleles):
                        G.add_node(n, v=i, seq=all_seq.upper(), gt=a)
                        for pn in pre_nodes:
                            G.add_edge(pn, n)
                        curr_nodes.append(n)
                        n += 1
                    prepre_nodes = pre_nodes
                    pre_nodes = curr_nodes
                else:
                    # previous nodes are alleles
                    last_v = G.nodes[pre_nodes[0]]["v"]
                    if self.records[last_v].stop > record.start:
                        # overlap, so no reference
                        # print("B1: We have an overlap with node", i)
                        for a, all_seq in enumerate(record.alleles):
                            G.add_node(n, v=i, seq=all_seq.upper(), gt=a)
                            for pn in prepre_nodes:
                                G.add_edge(pn, n)
                            curr_nodes.append(n)
                            n += 1
                        prepre_nodes = pre_nodes
                        pre_nodes.extend(curr_nodes)
                    else:
                        # no overlap, so reference, then variation
                        # print("B2: There is no overlap with node", i)
                        G.add_node(n, v=-1, seq=seq[last_p : record.start].upper(), gt=0)
                        for pn in pre_nodes:
                            G.add_edge(pn, n)
                        pre_nodes = [n]
                        curr_nodes = []
                        n += 1
                        for a, all_seq in enumerate(record.alleles):
                            G.add_node(n, v=i, seq=all_seq.upper(), gt=a)
                            for pn in pre_nodes:
                                G.add_edge(pn, n)
                            curr_nodes.append(n)
                            n += 1
                        prepre_nodes = pre_nodes
                        pre_nodes = curr_nodes
                last_p = record.stop
            G.add_node(n, v=-1, seq=seq[last_p : self.uend+1].upper(), gt=0)
            for pn in pre_nodes:
                G.add_edge(pn, n)
            n+=1
            self.graphs.append(G)


def main(args):
    logging.getLogger().setLevel(logging.DEBUG if args.verbose else logging.INFO)

    logging.info("Parsing reference..")
    reference = {}
    for record in SeqIO.parse(args.FA, "fasta"):
        reference[record.id] = str(record.seq).upper()

    CONF_TREES = {}
    if args.conf != "":
        logging.info("Parsing confidence regions..")
        for line in open(args.conf):
            chrom, s, e = line.strip("\n").split("\t")
            s, e = int(s), int(e)
            if chrom not in CONF_TREES:
                CONF_TREES[chrom] = IntervalTree()
            CONF_TREES[chrom][s : e + 1] = 0  # CHECKME: w?
        for tree in CONF_TREES.values():
            tree.merge_overlaps()

    TRF_TREES = {}
    if args.trf != "":
        logging.info("Parsing TRF regions..")
        for line in open(args.trf):
            chrom, s, e = line.strip("\n").split("\t")
            s, e = int(s), int(e)
            if chrom not in TRF_TREES:
                TRF_TREES[chrom] = IntervalTree()
            TRF_TREES[chrom][s - 100 : e + 100 + 1] = 0
        for tree in TRF_TREES.values():
            tree.merge_overlaps()

    logging.info("Parsing truth VCF..")
    vcf1 = VariantFile(args.VCF1)
    haploblocks_1 = []
    last_trf = (-1, -1)
    for record in vcf1 if args.region == "" else vcf1.fetch(region=args.region):
        if record.chrom in CONF_TREES and not (
            CONF_TREES[record.chrom][record.start]
            and CONF_TREES[record.chrom][record.stop]
        ):
            continue
        if record.info["SVLEN"] < args.l:
            continue

        trf = (-1, -1)
        if record.chrom in TRF_TREES:
            trf_overlaps = list(
                TRF_TREES[record.chrom].overlap(record.start, record.stop)
            )
            if len(trf_overlaps) != 0:
                trf = (trf_overlaps[0].begin, trf_overlaps[0].end)
            if last_trf != (-1, -1) and last_trf == trf:
                haploblocks_1[-1].append(record)
                continue
            else:
                last_trf = trf

        if (
            len(haploblocks_1) == 0
            or haploblocks_1[-1].chrom != record.chrom
            or (
                len(haploblocks_1[-1]) != 0
                and record.start - (haploblocks_1[-1].end + 1) > args.w
            )
        ):
            haploblocks_1.append(HaploBlock())
        haploblocks_1[-1].append(record)
    logging.info(f"Extracted {len(haploblocks_1)} HaploBlocks from truthset")
    # truth_trees = {}
    TRUTH_VAR = []
    TRUTH_NODES = []
    TRUTH_IDGR = []
    TRUTH_ODGR = []
    i = 0  # tracking an enumerate has strange behaviour
    for h in track(haploblocks_1):
        h.extend(reference[h.chrom], args.k, args.w)
        h.build_candidates(reference[h.chrom])
        TRUTH_VAR.append(len(h))
        TRUTH_NODES.append([0,0,0])
        # pflag = True
        for G in h.graphs:
            for node in G.nodes():
                ind, outd = G.in_degree(node), G.out_degree(node)
                TRUTH_IDGR.append(ind)
                TRUTH_ODGR.append(outd)
                i = 0 if G.nodes[node]["v"] == -1 else 1
                TRUTH_NODES[-1][i]+=1
                TRUTH_NODES[-1][2]+=1
        # if pflag:
        #     print(h.uregion())
        #     for i, G in enumerate(h.graphs):
        #         nx.draw(G, with_labels=True, pos=nx.nx_agraph.graphviz_layout(G))
        #         plt.savefig("graphs/" + h.uregion() + "." + str(i) + ".png")
        #         plt.clf()
        # if h.chrom not in truth_trees:
        #     truth_trees[h.chrom] = IntervalTree()
        # truth_trees[h.chrom][h.ustart : h.uend] = i
        # i += 1
    # print(TRUTH_VAR, TRUTH_NODES, TRUTH_IDGR, TRUTH_ODGR)

    logging.info("Parsing call VCF..")
    vcf2 = VariantFile(args.VCF2)
    haploblocks_2 = []
    last_trf = (-1, -1)
    for record in vcf2 if args.region == "" else vcf2.fetch(region=args.region):
        if record.chrom in CONF_TREES and not (
            CONF_TREES[record.chrom][record.start]
            and CONF_TREES[record.chrom][record.stop]
        ):
            continue
        if abs(record.info["SVLEN"]) < args.l:
            continue

        trf = (-1, -1)
        if record.chrom in TRF_TREES:
            trf_overlaps = list(
                TRF_TREES[record.chrom].overlap(record.start, record.stop)
            )
            if len(trf_overlaps) != 0:
                trf = (trf_overlaps[0].begin, trf_overlaps[0].end)
            if last_trf != (-1, -1) and trf != (-1, -1) and last_trf == trf:
                haploblocks_2[-1].append(record)
                continue
            else:
                last_trf = trf

        if (
            len(haploblocks_2) == 0
            or haploblocks_2[-1].chrom != record.chrom
            or (
                len(haploblocks_2[-1]) != 0
                and record.start - (haploblocks_2[-1].end + 1) > args.w
            )
        ):
            haploblocks_2.append(HaploBlock())
        haploblocks_2[-1].append(record)
    logging.info(f"Extracted {len(haploblocks_2)} HaploBlocks from callset")

    CALL_VAR = []
    CALL_NODES = []
    CALL_IDGR = []
    CALL_ODGR = []
    for h in track(haploblocks_2):
        h.extend(reference[h.chrom], args.k, args.w)
        h.build_candidates(reference[h.chrom])
        # print(h.uregion())
        CALL_VAR.append(len(h))
        CALL_NODES.append([0,0,0])
        for G in h.graphs:
            for node in G.nodes():
                ind, outd = G.in_degree(node), G.out_degree(node)
                CALL_IDGR.append(ind)
                CALL_ODGR.append(outd)
                i = 0 if G.nodes[node]["v"] == -1 else 1
                CALL_NODES[-1][i]+=1
                CALL_NODES[-1][2]+=1
    # print(CALL_VAR, CALL_NODES, CALL_IDGR, CALL_ODGR)

    VAR = pd.DataFrame([[i, "Truth"] for i in TRUTH_VAR] + [[i, "Call"] for i in CALL_VAR], columns = ["V", "Graph"])
    NODES = pd.DataFrame([[i,j,z, "Truth"] for i,j,z in TRUTH_NODES] + [[i,j,z, "Call"] for i,j,z in CALL_NODES], columns = ["R", "B", "N", "Graph"])
    IDGR = pd.DataFrame([[i, "Truth"] for i in TRUTH_IDGR] + [[i, "Call"] for i in CALL_IDGR], columns = ["D", "Graph"])
    ODGR = pd.DataFrame([[i, "Truth"] for i in TRUTH_ODGR] + [[i, "Call"] for i in CALL_ODGR], columns = ["D", "Graph"])
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2,2, figsize=(11,11))
    sns.boxplot(data=VAR, x="V", y="Graph", ax=ax1)
    sns.boxplot(data=NODES, x="N", y="Graph", ax=ax2)
    sns.boxplot(data=IDGR, x="D", y="Graph", ax=ax3)
    sns.boxplot(data=ODGR, x="D", y="Graph", ax=ax4)
    plt.tight_layout()
    plt.savefig("plot.png")
    

    # CALLS = {}
    # TRUES = {}
    # for call_hb in track(haploblocks_2):
    #     hits = truth_trees[h.chrom].overlap(call_hb.ustart, call_hb.uend)
    #     logging.debug("H", len(hits))
    #     if len(hits) == 0:
    #         for record in call_hb.records:
    #             record_idx = record.id
    #             CALLS[record_idx] == min(
    #                 CALLS[record_idx], float("inf")
    #             ) if record_idx in CALLS else float("inf")
    #     else:
    #         if len(hits) > 1:
    #             # TODO
    #             logging.info("Skipped call_hb since >1 hits")
    #             continue
    #         truth_hb = haploblocks_1[list(hits)[0].data]
    #         if len(truth_hb.records) > 5 or len(call_hb.records) > 5:
    #             # TODO
    #             logging.info("Skipped hit since too complex")
    #             continue

    #         ustart = min(call_hb.ustart, truth_hb.ustart)
    #         uend = max(call_hb.uend, truth_hb.uend)

    #         truth_hb.set_extension(ustart, uend)
    #         call_hb.set_extension(ustart, uend)

    #         # print(truth_hb.uregion(), call_hb.uregion())
    #         call_hb.build_candidates(reference[call_hb.chrom])
    #         truth_hb.build_candidates(reference[truth_hb.chrom], truth=True)

    #         assert len(truth_hb.candidates_haplotypes) <= 2
    #         assert len(call_hb.candidates_haplotypes) <= 2
    #         ALL_ED = [None, None, None, None]
    #         for ti1, truth_cand in enumerate(truth_hb.candidates_haplotypes_sequences):
    #             # print(ti1, len(truth_cand))
    #             ED = []
    #             for ti2, truth_seq in enumerate(truth_cand):
    #                 print(truth_hb.get_candidate(ti1, ti2)[0])
    #                 for ci1, call_cand in enumerate(
    #                     call_hb.candidates_haplotypes_sequences
    #                 ):
    #                     for ci2, call_seq in enumerate(call_cand):
    #                         print(call_hb.get_candidate(ci1, ci2)[0])
    #                         ed = lev.distance(call_seq, truth_seq)
    #                         ED.append((ed, (ci1, ci2), (ti1, ti2)))
    #             ED = sorted(ED, key=lambda x: x[0])
    #             ALL_ED[2 * ti1] = ED[0]
    #             if len(ED) > 1:
    #                 ALL_ED[2 * ti1 + 1] = ED[1]

    #         HAP1 = ALL_ED[0]
    #         HAP2 = ALL_ED[1]
    #         if ALL_ED[2] != None:
    #             HAP2 = ALL_ED[2]

    #         print(HAP1, HAP2)
    #         for hap in [HAP1, HAP2]:
    #             if hap == None:
    #                 continue
    #             ed, (c1, c2), (t1, t2) = hap

    #             # callTP/FPs
    #             call_hap, call_hap_seq = call_hb.get_candidate(c1, c2)
    #             for vidx, gt in call_hap:
    #                 record_idx = call_hb.records[vidx].id
    #                 if gt == 0:
    #                     CALLS[record_idx] == min(
    #                         CALLS[record_idx], float("inf")
    #                     ) if record_idx in CALLS else float("inf")
    #                 else:
    #                     CALLS[record_idx] = (
    #                         min(CALLS[record_idx], ed) if record_idx in CALLS else ed
    #                     )

    #             # trueTP/FN
    #             true_hap, true_hap_seq = truth_hb.get_candidate(t1, t2)
    #             print(true_hap)
    #             for vidx, gt in true_hap:
    #                 record_idx = truth_hb.records[vidx].id
    #                 print(record_idx, gt)
    #                 TRUES[record_idx] = (
    #                     min(TRUES[record_idx], ed) if record_idx in TRUES else ed
    #                 )
    # print(CALLS)
    # print(TRUES)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        prog="ProgramName",
        description="What the program does",
        epilog="Text at the bottom of help",
    )
    parser.add_argument("FA", help=".fa for reference")
    parser.add_argument("VCF1", help=".vcf for truth calls")
    parser.add_argument("VCF2", help=".vcdf for calls")
    parser.add_argument(
        "--region",
        dest="region",
        default="",
        type=str,
        help="region of interest",
    )
    parser.add_argument(
        "-k",
        dest="k",
        default=13,
        type=int,
        help="k-mer size (default: 13)",
    )
    parser.add_argument(
        "-w",
        dest="w",
        default=150,
        type=int,
        help="Window size (default: 150)",
    )
    # CHECKME: do we want this? What if two 30bp SVs can be merged into a single 60bp SV?
    parser.add_argument(
        "--minl",
        dest="l",
        default=50,
        type=int,
        help="Minimum SV length to consider (default: 50)",
    )
    parser.add_argument(
        "--trf",
        dest="trf",
        default="",
        type=str,
        help="Tandem repeats .bed file",
    )
    parser.add_argument(
        "--conf",
        dest="conf",
        default="",
        type=str,
        help="Confident regions .bed file",
    )
    parser.add_argument(
        "-v",
        dest="verbose",
        default=False,
        action="store_true",
        help="verbose mode",
    )
    parser.add_argument(
        "-o",
        dest="odir",
        default=".",
        type=str,
        help="Output directory",
    )
    main(parser.parse_args())
