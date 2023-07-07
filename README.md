# sveval

### Compilation
``` sh
git clone --recursive https://github.com/ldenti/sveval.git

cd sveval/deps/vg
git apply ../../patches/vg.diff
cd deps/tabixpp
git apply ../../../../patches/tabixpp.diff
cd ../..
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

cd ../ksw2
make

cd ../..
make -j8
```

### Usage
``` sh
# preprocess your VCF
bcftools norm -m- --threads 8 [truth.vcf.gz] | python3 scripts/clean_vcf.py | bcftools norm -m+ --threads 8 -Oz > [truth.sv.vcf.gz]
tabix -p vcf [truth.sv.vcf.gz]
# run evaluation on preprocessed VCFs
./benchHG [reference.fa] [truth.sv.vcf.gz] [predictions.sv.vcf.gz] [--conf <regions.bed>] [--trf <trf.bed>] [-@ <threads>] > [output]
```
