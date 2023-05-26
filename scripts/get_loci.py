#!/usr/bin/env python3

import sys
import os
import argparse
import logging

from Bio import SeqIO
from pysam import VariantFile
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

    def extend(self, seq, klen, wsize):
        # Extend the region to unique kmers
        unique_kmer_s = find_unique_kmer(
            seq[self.start - wsize - 1 : self.start], klen, True
        )
        self.ustart = self.start - (wsize + unique_kmer_s)
        unique_kmer_e = find_unique_kmer(seq[self.end : self.end + wsize], klen, False)
        self.uend = self.end + unique_kmer_e


def main(args):
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
    ovcf1 = VariantFile(os.path.join(args.outd, "truth.vcf"), 'w', header=vcf1.header)

    haploblocks_1 = []
    last_trf = (-1, -1)
    for record in vcf1 if args.region == "" else vcf1.fetch(region=args.region):
        if record.chrom in CONF_TREES and not (
            CONF_TREES[record.chrom][record.start]
            and CONF_TREES[record.chrom][record.stop]
        ):
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
        ovcf1.write(record)
    ovcf1.close()
    logging.info(f"Extracted {len(haploblocks_1)} HaploBlocks from truthset")

    truth_trees = {}
    i = 0  # tracking an enumerate has strange behaviour
    for h in track(haploblocks_1, console=console):
        h.extend(reference[h.chrom], args.k, args.w)
        if h.chrom not in truth_trees:
            truth_trees[h.chrom] = IntervalTree()
        truth_trees[h.chrom][h.ustart : h.uend] = i
        i += 1

    logging.info("Parsing call VCF..")
    vcf2 = VariantFile(args.VCF2)
    ovcf2 = VariantFile(os.path.join(args.outd, "call.vcf"), 'w', header=vcf2.header)
    haploblocks_2 = []
    last_trf = (-1, -1)
    for record in vcf2 if args.region == "" else vcf2.fetch(region=args.region):
        if record.chrom in CONF_TREES and not (
            CONF_TREES[record.chrom][record.start]
            and CONF_TREES[record.chrom][record.stop]
        ):
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
        ovcf2.write(record)
    ovcf2.close()
    logging.info(f"Extracted {len(haploblocks_2)} HaploBlocks from callset")

    loci_trees = {}
    bad_loci_f = open(os.path.join(args.outd, "loci.bad.txt"), "w")

    for h in track(haploblocks_2, console=console):
        h.extend(reference[h.chrom], args.k, args.w)
        hits = list(truth_trees[h.chrom].overlap(h.ustart, h.uend))
        if len(hits) == 0:
            logging.info(f"Skipped {h.uregion()} since 0 hits")
            print(f"{h.chrom}:{h.ustart}-{h.uend}", file=bad_loci_f)
            continue
        else:
            ustart = min(
                h.ustart, min([haploblocks_1[hit.data].ustart for hit in hits])
            )
            uend = max(h.uend, max([haploblocks_1[hit.data].uend for hit in hits]))
            if h.chrom not in loci_trees:
                loci_trees[h.chrom] = IntervalTree()
            loci_trees[h.chrom][ustart:uend] = 0
    bad_loci_f.close()

    good_loci_f = open(os.path.join(args.outd, "loci.good.txt"), "w")
    for chrom in loci_trees:
        loci_trees[chrom].merge_overlaps()
        for interval in loci_trees[chrom]:
            print(f"{chrom}:{interval.begin}-{interval.end}", file=good_loci_f)
    good_loci_f.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        # prog="ProgramName",
        # description="What the program does",
        # epilog="Text at the bottom of help",
    )
    parser.add_argument("FA", help=".fa for reference")
    parser.add_argument("VCF1", help=".vcf for truth calls")
    parser.add_argument("VCF2", help=".vcdf for calls")
    parser.add_argument(
        "--region", dest="region", default="", type=str, help="region of interest"
    )
    parser.add_argument(
        "-k", dest="k", default=13, type=int, help="k-mer size (default: 13)"
    )
    parser.add_argument(
        "-w", dest="w", default=250, type=int, help="Window size (default: 250)"
    )
    parser.add_argument(
        "--trf", dest="trf", default="", type=str, help="Tandem repeats .bed file"
    )
    parser.add_argument(
        "--conf", dest="conf", default="", type=str, help="Confident regions .bed file"
    )
    parser.add_argument(
        "-o",
        dest="outd",
        default=".",
        type=str,
        help="Output directory (default: .)",
    )
    main(parser.parse_args())
