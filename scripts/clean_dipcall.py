#!/usr/bin/env python3

import sys
from pysam import VariantFile


def main():
    vcf_path = sys.argv[1]

    vcf = VariantFile(vcf_path)
    vcf.header.add_line(
        '##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Difference in length between REF and ALT alleles">'
    )
    print(vcf.header, end="")
    last_pos = -1
    last_len = -1
    for rec in vcf.fetch():
        assert len(rec.alts) == 1
        l = len(rec.ref) - len(rec.alts[0])
        if l < 50:
            continue
        if l == last_len and rec.pos == last_pos:
            continue
        last_len = l
        last_pos = rec.pos
        rec.info["SVLEN"] = l
        print(rec, end="")


if __name__ == "__main__":
    main()
