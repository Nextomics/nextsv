all: longreadqc pigz minimap2 ngmlr winnowmap seqtk

longreadqc:
	cd external_src/longreadqc && make && cp longreadqc ../../bin/ && make clean

pigz:
	cd external_src/pigz && make && cp pigz ../../bin/ && cp unpigz ../../bin/ && make clean

minimap2:
	cd external_src/minimap2 && make && cp minimap2 ../../bin/ && make clean

ngmlr:
	cd external_src/ngmlr/ && mkdir -p build && cd build && cmake .. && make && cp ../bin/ngmlr-0.2.8/ngmlr ../../../bin/ && cd ..  && rm -r build bin

winnowmap:
	cd external_src/Winnowmap/ && make && mv bin/* ../../bin/

seqtk:
	cd external_src/seqtk/ && make && mv seqtk ../../bin/

