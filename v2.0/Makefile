GCC= g++ -O3
LIB= -L. -lgzstream -lz

fqclean: main.cc main.h fq_read.h UID_type.h Makefile
	$(GCC) $< -o $@ $(LIB)

.PHONY: clean

clean:
	-rm fqclean
