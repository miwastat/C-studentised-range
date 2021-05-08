#!sh

OBJ=strng_lq.o strng_lp.o rng_lp.o nrml_p.o
CC=gcc

# Strip *.exe files in Windows_NT
ifeq ($(OS),Windows_NT)
	EXE=.exe
endif

strng_tbl: strng_tbl.o $(OBJ)
	$(CC) strng_tbl.o $(OBJ) -o strng_tbl -lm
	strip strng_tbl$(EXE)

strng_tbl.o: strng_tbl.c
	$(CC) -c strng_tbl.c

strng_lq_tst: strng_lq_tst.o $(OBJ)
	$(CC) strng_lq_tst.o $(OBJ) -o strng_lq_tst -lm
	strip strng_lq_tst$(EXE)

strng_lq_tst.o: strng_lq_tst.c
	$(CC) -c strng_lq_tst.c

strng_lq.o: strng_lq.c
	$(CC) -c strng_lq.c

strng_lp.o: strng_lp.c
	$(CC) -c strng_lp.c

rng_lp_tst: rng_lp_tst.o rng_lp.o nrml_p.o
	$(CC) rng_lp_tst.o rng_lp.o nrml_p.o -o rng_lp_tst -lm
	strip rng_lp_tst$(EXE)

rng_lp_tst.o: rng_lp_tst.c
	$(CC) -c rng_lp_tst.c

rng_lp.o: rng_lp.c 
	$(CC) -c rng_lp.c

nrml_p.o: nrml_p.c
	$(CC) -c nrml_p.c
