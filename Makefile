graph.pdf : atpos.dat graph.py
	python graph.py

atpos.dat : atoms.x
	mpiexec -n 1 atoms.x

atoms.x : atoms.c
	mpicc atoms.c -o atoms.x -lm

clean :
	rm atoms.x atpos.dat graph.pdf
