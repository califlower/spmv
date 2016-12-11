#/*
# *********************************************
# *  314 Principles of Programming Languages  *
# *  Fall 2016                                *
# *********************************************
# */

CCFLAGS = -Wall -fopenmp

all: sequential parallelS parallelD parallelN

sequential: spmv_seq.c mmio.c mmio.h utils.c utils.h 
	gcc $(CCFLAGS) spmv_seq.c mmio.c utils.c utils.h -o spmv_seq

parallelS: spmv_N_thread_static.c mmio.c mmio.h utils.c utils.h
	gcc $(CCFLAGS) spmv_N_thread_static.c mmio.c utils.c -o pstatic

parallelD: spmv_N_thread_dynamic.c mmio.c mmio.h utils.c utils.h
	gcc $(CCFLAGS) spmv_N_thread_dynamic.c mmio.c utils.c -o pdynamic

parallelN: spmv_N_thread_newalg.c mmio.c mmio.h utils.c utils.h
	gcc $(CCFLAGS) spmv_N_thread_newalg.c mmio.c utils.c -o pnewalg 

clean:
	rm -f -r spmv_seq pstatic pdynamic pnewalg




