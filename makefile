graph.pdf time.pdf : atpos.dat graph.py times.dat timer.py
    #python3 cluster requirement
	python3 graph.py
	python3 timer.py

atpos.dat times.dat : a.out
    # real, user, sys
	/usr/bin/time -o times.dat -f "%e \t%U \t%S" ./a.out 1
	/usr/bin/time -a -o times.dat -f "%e \t%U \t%S" ./a.out 2
	/usr/bin/time -a -o times.dat -f "%e \t%U \t%S" ./a.out 4
	/usr/bin/time -a -o times.dat -f "%e \t%U \t%S" ./a.out 8
	/usr/bin/time -a -o times.dat -f "%e \t%U \t%S" ./a.out 16

a.out : atoms.c
	gcc -lm -fopenmp -O3 atoms.c

clean :
	rm a.out atpos.dat graph.pdf time.pdf times.dat
