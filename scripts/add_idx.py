#!/usr/bin/env python3

import sys
from pysam import VariantFile


def main():
    vcf_path = sys.argv[1]

    vcf = VariantFile(vcf_path)
    ovcf = VariantFile("-", "w", header=vcf.header)
    i = 1
    for rec in vcf.fetch():
        if rec.id is None:
            rec.id = f"v{i}"
            i += 1
        ovcf.write(rec)


if __name__ == "__main__":
    main()
