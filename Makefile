all: longreadqc pigz

longreadqc:
	cd external_src/longreadqc && make && cp longreadqc ../../bin/ && make clean

pigz:
	cd external_src/pigz && make && cp pigz ../../bin/ && cp unpigz ../../bin/ && make clean
