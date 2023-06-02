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
./sveval [reference.fa] [truth.vcf.gz] [predictions.vcf.gz] [--conf <regions.bed>] [--trf <trf.bed>] [-@ <threads>] > [output]
```

## TODO
- [ ] more chromosomes
- [ ] plots and analysis
- [ ] cmake
- [ ] static binary
- [ ] conda
