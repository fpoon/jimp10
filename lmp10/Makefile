aprox: main.o splines.o points.o aproksymator_na_bazie.o gaus/libge.a intrp legendre.o
	$(CC) -o aprox  main.o splines.o points.o interpolator.o legendre.o -L gaus -l ge -lm

intrp: main.o splines.o points.o interpolator.o gaus/libge.a legendre.o
	$(CC) -o intrp  main.o splines.o points.o interpolator.o legendre.o -L gaus -l ge -lm

prosta: main.o splines.o points.o prosta.o
	$(CC) -o prosta  main.o splines.o points.o prosta.o	

aproksymator_na_bazie.o: makespl.h points.h gaus/piv_ge_solver.h
	$(CC) -I gaus -c aproksymator_na_bazie.c -o aproksymator_na_bazie.o

interpolator.o: makespl.h points.h gaus/piv_ge_solver.h 
	$(CC) -I gaus -c interpolator.c

legendre.o:
	$(CC) -c legendre.c -o legendre.o
.PHONY: clean

clean:
	-rm *.o aprox intrp prosta
