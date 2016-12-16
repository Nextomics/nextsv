CC = gcc   -O0 -g -I ./include -L ./lib/

SOBJ = main.o 

BAMSTAT: $(SOBJ)
	$(CC) -o bamstat $(SOBJ) -l hts -l z -l m -l pthread

clean:
	rm -f *.o bamstat 
