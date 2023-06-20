##
# sveval
#
# @file
# @version 0.0.1

CXXFLAGS=-Wall -Wno-sign-compare -std=c++14 -I. -I./deps -I./deps/spdlog/include -I./deps/interval-tree/include -I./deps/PaSGAL/src/include -I./deps/vg/src -I./deps/vg/include -fopenmp -pthread
LIBS=-L./deps/vg/lib -lvg -lvgio -lhts -ldeflate -lbdsg -lhandlegraph -lvcflib -lprotobuf -lxg -lgbwtgraph -lgbwt -lgcsa2 -ltabixpp -lsdsl -ldivsufsort -ldivsufsort64 -ldwfl -ldw -ldl -ldwelf -lelf -lebl -lzstd -lz -lbz2 -llzma -latomic
OPT=-O2

.PHONY: all

all: sveval

sveval: ./deps/bcftools/smpl_ilist.o ./deps/bcftools/filter.o ./deps/bcftools/smpl_ilist.o ./deps/bcftools/consensus.o locator.o consenser.o graph.o aligner.o tscorer.o cscorer.o main.o
	@echo "* Linking $@"
	$(CXX) $(OPT) $(CXXFLAGS) -o $@ $^ $(LIBS)

%.o: %.cpp
	@echo '* Compiling $<'
	$(CXX) $(OPT) $(CXXFLAGS) -o $@ -c $<

debug: OPT=-O0 -g -DPDEBUG
debug: sveval

clean:
	rm -rf *.o

# end
