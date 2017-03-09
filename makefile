graph.pdf time.pdf : Energies.dat graph.py times.dat timer.py
    #python3 cluster requirement
	python3 graph.py
	python3 timer.py

Energies.dat times.dat : a.out
    # real, user, sys
	/usr/bin/time -o times.dat -f "%e \t%U \t%S" ./a.out 1
#	/usr/bin/time -a -o times.dat -f "%e \t%U \t%S" ./a.out 2
#	/usr/bin/time -a -o times.dat -f "%e \t%U \t%S" ./a.out 4

a.out : atoms.c
	gcc -lm -fopenmp -O3 atoms.c

clean :
	rm a.out Energies.dat graph.pdf time.pdf times.dat
