#!/usr/bin/env python3

import sys
from pysam import VariantFile


def main():
    vcf_path = sys.argv[1]
    min_l = int(sys.argv[2])
    max_l = 30000

    vcf = VariantFile(vcf_path)
    if "SVLEN" not in vcf.header.info:
        vcf.header.add_line(
            '##INFO=<ID=SVLEN,Number=A,Type=Integer,Description="Difference in length between REF and ALT alleles">'
        )
    if "SVTYPE" not in vcf.header.info:
        vcf.header.add_line(
            '##INFO=<ID=SVTYPE,Number=A,Type=String,Description="Type of structural variant">'
        )
    print(vcf.header, end="")
    i = 1
    for rec in vcf:
        if "SVTYPE" in rec.info:
            if rec.info["SVTYPE"] in ["DUP", "INV", "BND"]:
                continue
        if "SVLEN" in rec.info:
            assert len(rec.alts) == 1
            l = rec.info["SVLEN"]
            if type(l) == tuple:
                l = l[0]
            if abs(l) < min_l or abs(l) > max_l:
                continue
        else:
            # here we need to take only "valid" alleles, add info to INFO and fix genotype
            # assuming single allele thanks to bcftools norm -m-
            alt = rec.alts[0]
            if alt == "*":
                continue
            l = len(alt) - len(rec.ref)
            if abs(l) < min_l or abs(l) > max_l:
                continue
            if len(alt) > max_l and len(rec.ref) > max_l:
                continue

            gt1 = rec.samples[0]["GT"][0]
            gt1 = 0 if gt1 == None else gt1
            gt2 = rec.samples[0]["GT"][1]
            gt2 = 0 if gt2 == None else gt2
            is_phased = rec.samples[0].phased

            rec.samples[0]["GT"] = (gt1, gt2)
            rec.samples[0].phased = is_phased
            rec.info["SVTYPE"] = ["INS" if l > 0 else "DEL"]
            rec.info["SVLEN"] = [l]

            if rec.id is None:
                rec.id = f"v{i}"
                i += 1
        print(rec, end="")


if __name__ == "__main__":
    main()
