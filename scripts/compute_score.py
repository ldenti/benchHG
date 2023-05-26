#!/usr/bin/env python3

import sys
import re
from itertools import product

from pysam import VariantFile


def is_subpath(P, p):
    for i in range(len(P)):
        if P[i : i + len(p)] == p:
            return True
    return False

def check_path(all_path, is_ins, true_path, in_edges, out_edges):
    correct = False
    if is_ins:
        if is_subpath(true_path, all_path):
            correct = True
    else:
        pre_nodes = in_edges[all_path[0]]
        post_nodes = out_edges[all_path[-1]]
        for pair in product(pre_nodes, post_nodes):
            if is_subpath(true_path, list(pair)):
                correct = True
                break # CHECKME: is this correct?
    return correct

def main():
    gfa_path = sys.argv[1]
    tvcf_path = sys.argv[2]
    cvcf_path = sys.argv[3]
    gaf_path = sys.argv[4]

    region = gaf_path.split("/")[-2]
    s, e = region.split(":")[1].split("-")
    l = int(e) - int(s) + 1

    # Computing scores from GAF
    SCORES = []
    for line in open(gaf_path):
        qname, ql, qs, qe, strand, gs, ge, score, cigar, path = line.strip("\n").split(
            "\t"
        )
        ql = int(ql)
        OPs = {"X": 0, "=": 0, "I": 0, "D": 0}
        for match in re.finditer("[0-9]+[=XID]", cigar):
            match = match.group(0)
            n, op = int(match[:-1]), match[-1]
            OPs[op] += n
        il = OPs["X"] + OPs["="] + OPs["I"]
        score = 1 - (OPs["I"] + OPs["D"] + OPs["X"] + abs(ql - il)) / ql # CHECKME
        path = path.split("-")
        SCORES.append([path, score])

    # Extracting paths from GFA and intersecting them with call VCF
    # to assign each path to the correct variant
    nodes = {}
    for line in open(gfa_path):
        if line.startswith("S"):
            _, idx, seq = line.strip("\n").split("\t")
            nodes[idx] = seq
    paths = []
    for line in open(gfa_path):
        if line.startswith("P"):
            _, idx, path, _ = line.strip("\n").split("\t")
            if idx[0] != "_":
                continue
            path = [x[:-1] for x in path.split(",")]
            path_seq = "".join([nodes[x] for x in path])
            paths.append((path, path_seq))
    in_edges = {path[0]: [] for path, _ in paths}
    out_edges = {path[-1]: [] for path, _ in paths}
    for line in open(gfa_path):
        if line.startswith("L"):
            _, id1, _, id2, _, _ = line.strip("\n").split("\t")
            if id1 in out_edges:
                out_edges[id1].append(id2)
            if id2 in in_edges:
                in_edges[id2].append(id1)
    # Assign scores to truth calls
    tvcf = VariantFile(tvcf_path)
    TRUTHS = {}
    for rec in tvcf:
        idx = rec.id
        gt = rec.samples[0]["GT"]

        if gt == (1, 0):
            score = SCORES[0][1]
        elif gt == (0, 1):
            score = SCORES[1][1]
        elif gt == (0, 0):
            score = -1
        elif gt[0] == gt[1]:
            score = max(SCORES[0][1], SCORES[1][1])
        else:
            score = (SCORES[0][1] + SCORES[1][1]) / 2
        TRUTHS[idx] = score

    # Assign scores to calls
    CALLS = {}
    cvcf = VariantFile(cvcf_path)
    for rec in cvcf:
        refall = rec.ref
        altalls = rec.alts
        assert len(altalls) <= 2  # TODO: assuming diploid
        # retrieve path by matching the sequences
        all_paths = [None for _ in altalls]
        for a, altall in enumerate(altalls):
            is_ins = len(refall) < len(altall)
            for p, s in paths:
                if is_ins and s == altall[1:].upper():
                    all_paths[a] = (p, is_ins)
                    break
                elif not is_ins and s == refall[1:].upper():
                    all_paths[a] = (p, is_ins)
                    break
        # We assume that we must have found the paths since graph is built from same VCF. If not, issue could be in graph construction
        score = 0
        if len(all_paths) == 1:
            assert all_paths[0] != None
            # we assign the score if the allele is covered by at least one path
            if check_path(all_paths[0][0], all_paths[0][1], SCORES[0][0], in_edges, out_edges) or check_path(all_paths[0][0], all_paths[0][1], SCORES[1][0], in_edges, out_edges):
                score = (SCORES[0][1] +  SCORES[1][1])/2
        else:
            assert all_paths[0] != None and all_paths[1] != None
            score1, score2 = 0, 0
            is_covered_1_by_1 = check_path(all_paths[0][0], all_paths[0][1], SCORES[0][0], in_edges, out_edges)
            is_covered_1_by_2 = check_path(all_paths[0][0], all_paths[0][1], SCORES[1][0], in_edges, out_edges)
            is_covered_2_by_1 = check_path(all_paths[1][0], all_paths[1][1], SCORES[0][0], in_edges, out_edges)
            is_covered_2_by_2 = check_path(all_paths[1][0], all_paths[1][1], SCORES[1][0], in_edges, out_edges)
            # CHECKME: assuming that the same path cannot cover both alleles
            if (is_covered_1_by_1 or is_covered_1_by_2) and (is_covered_2_by_1 or is_covered_2_by_2):
                # arbitrary
                score1 = SCORES[0][1]/2
                score2 = SCORES[1][1]/2
            else:
                if is_covered_1_by_1 or is_covered_1_by_2:
                    if is_covered_1_by_1:
                        score1 = SCORES[0][1]/2
                    else:
                        score1 = SCORES[1][1]/2
                elif is_covered_2_by_1 or is_covered_2_by_2:
                    if is_covered_2_by_1:
                        score2 = SCORES[0][1]/2
                    else:
                        score2 = SCORES[1][1]/2
            score = score1 + score2
        CALLS[rec.id] = score

    for idx, score in TRUTHS.items():
        print("T", idx, score)
    for idx, score in CALLS.items():
        print("C", idx, score)


if __name__ == "__main__":
    main()
