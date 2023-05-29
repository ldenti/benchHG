##
# Project Title
#
# @file
# @version 0.1

CXXFLAGS=-g -Wall -Wno-sign-compare -O0 -std=c++14 -I. -I./spdlog/include -I/home/ld/code/PaSGAL/src/include -I./vg/src -I./vg/include -fopenmp -pthread
LIBS=-L./vg/lib -lvg -lvgio -lhts -lbdsg -lhandlegraph -lvcflib -lprotobuf -lxg -lgbwtgraph -lgbwt -lgcsa2 -ltabixpp -lsdsl -ldivsufsort -ldivsufsort64 -ldwfl -ldw -ldwelf -lelf -lebl -lzstd -lz -lbz2 -llzma -latomic

# LD_LIB_FLAGS := -lvcflib -ltabixpp -lgssw -lssw -lsublinearLS -lpthread -lncurses -lgcsa2 -lgbwtgraph -lgbwt -lkff -ldivsufsort -ldivsufsort64 -lvcfh -lraptor2 -lpinchesandcacti -l3edgeconnected -lsonlib -lfml -lstructures -lbdsg -lxg -lsdsl -lzstd -lhandlegraph

.PHONY: all

all: sveval

# vg/obj/region.o
sveval: bcftools/smpl_ilist.o bcftools/filter.o bcftools/smpl_ilist.o bcftools/consensus.o consenser.o main.o
	@echo "* Linking $@"
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LIBS)

%.o: %.cpp
	@echo '* Compiling $<'
	$(CXX) $(CXXFLAGS) -o $@ -c $<

clean:
	rm -rf *.o

# end
