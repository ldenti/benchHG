#!/usr/bin/env python3

import sys
from pysam import VariantFile


def main():
    vcf_path = sys.argv[1]
    min_l = 20

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
            l = rec.info["SVLEN"]
            if abs(l) < min_l:
                continue
        else:
            # here we need to take only "valid" alleles, add to INFO and fix genotype
            alts = []
            Ls = []
            types = []
            toremove = []
            for a, aseq in enumerate(rec.alts):
                if aseq == "*":
                    toremove.append(a + 1)
                    continue
                l = len(aseq) - len(rec.ref)
                if abs(l) < min_l:
                    toremove.append(a + 1)
                    continue
                alts.append(aseq)
                Ls.append(l)
                if l > 0:
                    types.append("INS")
                else:
                    types.append("DEL")
            if len(Ls) == 0 or len(Ls) > 2:
                continue
            alts = tuple(alts)
            Ls = tuple(Ls)
            types = tuple(types)
            if len(toremove) > 0:
                gt1 = rec.samples[0]["GT"][0]
                gt2 = rec.samples[0]["GT"][1]
                gt1 = 0 if gt1 == None or gt1 in toremove else gt1
                for tr in toremove:
                    if gt1 > tr:
                        gt1 -= 1
                gt2 = 0 if gt2 == None or gt2 in toremove else gt2
                for tr in toremove:
                    if gt2 > tr:
                        gt2 -= 1
                is_phased = rec.samples[0].phased
                rec.samples[0]["GT"] = (gt1, gt2)
                rec.samples[0].phased = is_phased

            rec.alts = alts
            rec.info["SVTYPE"] = types
            rec.info["SVLEN"] = Ls
            if "AD" in rec.samples[0]:
                rec.samples[0].pop("AD")

            if rec.id is None:
                rec.id = f"v{i}"
                i += 1
        print(rec, end="")


if __name__ == "__main__":
    main()
