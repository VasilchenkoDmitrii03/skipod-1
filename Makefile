prog_omp: poisson_omp.c problem.h
	xlc_r -O1 -qsmp=omp -qarch=pwr8 poisson_omp.c -o prog_omp -lm

upload:
	scp poisson_omp.c problem.h Makefile omp_*.lsf polus:~/skmodel/

local:
	gcc -O1 -Wno-unknown-pragmas poisson_omp.c -o prog -lm
	./prog 20 20
	./prog 40 40
	./prog 80 80

run400:
	bsub < omp_1_400x600.lsf
	bsub < omp_2_400x600.lsf
	bsub < omp_4_400x600.lsf
	bsub < omp_8_400x600.lsf
	bsub < omp_16_400x600.lsf

run800:
	bsub < omp_1_800x1200.lsf
	bsub < omp_4_800x1200.lsf
	bsub < omp_8_800x1200.lsf
	bsub < omp_16_800x1200.lsf
	bsub < omp_32_800x1200.lsf

results:
	grep -h "^M=" omp_*.out | sort

clean:
	rm -f *.out *.err

.PHONY: upload local run400 run800 results clean
