# sveval

### Compilation
``` sh
git clone --recursive https://github.com/ldenti/sveval.git

cd sveval/deps/vg
. ./source_me.sh && make -j8

cd ../htslib
autoreconf -i
./configure
make -j8

cd ../bcftools
git apply ../../patches/bcftools.diff
make -j8

cd ../PaSGAL
git apply ../../patches/pasgal.diff

cd ../..
make -j8
```

### Usage
``` sh
# preprocess your VCF
bcftools norm -m- --threads 8 [truth.vcf.gz] | python3 scripts/clean_vcf.py | bcftools norm -m+ --threads 8 -Oz > [truth.sv.vcf.gz]
tabix -p vcf [truth.sv.vcf.gz]
./sveval [reference.fa] [truth.sv.vcf.gz] [predictions.sv.vcf.gz] [--conf <regions.bed>] [--trf <trf.bed>] [-@ <threads>] > [output]
```

## TODO
- [ ] more chromosomes
- [ ] plots and analysis
- [ ] cmake
- [ ] static binary
- [ ] conda
