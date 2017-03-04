graph.pdf time.pdf : atpos.dat graph.py times.dat timer.py
    #python3 cluster requirement
	python3 graph.py
	python3 timer.py

atpos.dat times.dat : a.out
    # real, user, sys
	time -o times.dat -f "%e \t%U \t%S" ./a.out 1
	time -a -o times.dat -f "%e \t%U \t%S" ./a.out 2
	time -a -o times.dat -f "%e \t%U \t%S" ./a.out 4

a.out : atoms.c
	gcc -lm -fopenmp atoms.c

clean :
	rm a.out atpos.dat graph.pdf time.pdf times.dat
