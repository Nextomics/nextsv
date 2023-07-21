CXX = g++
CXXFLAGS = -g -Og -std=c++11
LIBS = -lz

all: longreadqc

longreadqc: longreadqc.cpp filter_fq.cpp  qc_bam.cpp  qc_fq.cpp qc_fa.cpp qc_paf.cpp tk.cpp
	$(CXX) $(CXXFLAGS) $^ -o $@ $(LIBS)

clean:
	rm -f longreadqc
