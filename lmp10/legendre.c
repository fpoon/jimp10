#include "makespl.h"
#include "gaus/piv_ge_solver.h"

#include <stdio.h>
#include <stdlib.h>
#include <float.h>

#include <stdio.h>
#include <math.h>

#define DERIV_H 0.001

long newton(int n, int k)
{
  if( n < k || n < 0 || k < 0) return -1;
  if( k== 0 || n == k) return 1;
  return newton(n-1,k-1)+newton(n-1,k);
}

double legendre(double x, int n)
{
  if (n == 0) return 1;
  
  int i=0;
  double ret=0;
  
  double factor = pow(2,-n);
  
  //Dodane ze względów optymalizacyjnych
  switch ( n )
  {
	case 0:
		return 1.0;
	case 1:
		return x;
	case 2:
		return (3.0*pow(x,2.0) - 1.0) / 2.0;
	case 3:
		return (5.0*pow(x,3.0) - 3.0*x) / 2.0;
	case 4:
		return (35.0*pow(x,4.0) - 30.0*pow(x,2.0) + 3.0) / 8.0;
	case 5:
		return (63.0*pow(x,5.0) - 70.0*pow(x,3.0) + 15.0*x) / 8.0;
	case 6:
		return (231.0*pow(x,6.0) - 315*pow(x,4.0) + 105.0*pow(x,2.0) - 5.0) / 16.0;
	case 7:
		return (429.0*pow(x,7.0) - 693.0*pow(x,5.0) + 315.0*pow(x,3.0) - 35.0*x) / 16.0;
	case 8:
		return (6435.0*pow(x,8.0) - 12012.0*pow(x,6.0) + 6930.0*pow(x,4.0) - 1260.0*pow(x,2.0) + 35.0) / 128.0;
	case 9:
		return (12155*pow(x,9.0) - 25740*pow(x,7.0) + 18018*pow(x,5.0) - 4620*pow(x,3.0) + 315*x) / 128.0;
	case 10:
		return (46189*pow(x,10.0) - 109395*pow(x,8.0) + 90090*pow(x,6.0) - 30030*pow(x,4.0) + 3465*pow(x,2.0) - 63.0) / 256.0;
	case 11:
		return (88179*pow(x,11.0) - 230945*pow(x,9.0) + 218790*pow(x,7.0) - 90090*pow(x,5.0) + 15015*pow(x, 3.0) - 693.0*x) / 256.0;
  }
  
  //Na wypadek gdyby trzeba było policzyć wielomian wyższego stopnia
  for(i = 0; i <= n/2; i++)
  {
    ret += newton(n,i)*newton(2*n-2*i,n)*pow(x,n-2*i)*((i%2)==0?1:-1);
  }
  
  return ret*factor;
}

double derivative(double y1, double y2)
{
   double h = 2*DERIV_H;
   return (y2-y1)/h;
}


/* UWAGA: liczbę używanych f. bazowych można ustawić przez wartość
          zmiennej środowiskowej APPROX_BASE_SIZE
*/

/*
 * Funkcje bazowe: n - liczba funkcji a,b - granice przedzialu aproksymacji i
 * - numer funkcji x - wspolrzedna dla ktorej obliczana jest wartosc funkcji
 */
double fi(double a, double b, int n, int i, double x)
{	
	return legendre(x,i);
}

/* Pierwsza pochodna fi */
double dfi(double a, double b, int n, int i, double x)
{
	return derivative(legendre(x-DERIV_H,i),legendre(x+DERIV_H,i));
}

/* Druga pochodna fi */
double d2fi(double a, double b, int n, int i, double x)
{	
	return derivative(
	  derivative(legendre(x-2*DERIV_H,i),legendre(x,i)),
	  derivative(legendre(x,i),legendre(x+2*DERIV_H,i))
 			);
}

/* Trzecia pochodna fi */
double d3fi(double a, double b, int n, int i, double x)
{	
	return derivative(
	      derivative(
			      derivative(legendre(x-4*DERIV_H,i),legendre(x-2*DERIV_H,i)),
			      derivative(legendre(x-2*DERIV_H,i),legendre(x,i))
		    ),
	      derivative(
			      derivative(legendre(x,i),legendre(x+2*DERIV_H,i)),
			      derivative(legendre(x+2*DERIV_H,i),legendre(x+4*DERIV_H,i))
 			)
 		);
}

/* Pomocnicza f. do rysowania bazy */
double
xfi(double a, double b, int n, int i, FILE *out)
{
	double		h = (b - a) / (n - 1);
	double		h3 = h * h * h;
	int		hi         [5] = {i - 2, i - 1, i, i + 1, i + 2};
	double		hx      [5];
	int		j;

	for (j = 0; j < 5; j++)
		hx[j] = a + h * hi[j];

	fprintf( out, "# nb=%d, i=%d: hi=[", n, i );
	for( j= 0; j < 5; j++ )
		fprintf( out, " %d", hi[j] );
	fprintf( out, "] hx=[" );
	for( j= 0; j < 5; j++ )
		fprintf( out, " %g", hx[j] );
	fprintf( out, "]\n" );
}

void
make_spl(points_t * pts, spline_t * spl)
{

	matrix_t       *eqs= NULL;
	double         *x = pts->x;
	double         *y = pts->y;
	double		a = x[0];
	double		b = x[pts->n - 1];
	int		i, j, k;
	int		nb = pts->n - 3 > 10 ? 10 : pts->n - 3;
  char *nbEnv= getenv( "APPROX_BASE_SIZE" );

	if( nbEnv != NULL && atoi( nbEnv ) > 0 )
		nb = atoi( nbEnv );

	eqs = make_matrix(nb, nb + 1);

#ifdef DEBUG
#define TESTBASE 500
	{
		FILE           *tst = fopen("debug_base_plot.txt", "w");
		double		dx = (b - a) / (TESTBASE - 1);
		for( j= 0; j < nb; j++ )
			xfi( a, b, nb, j, tst );
		for (i = 0; i < TESTBASE; i++) {
			fprintf(tst, "%g", a + i * dx);
			for (j = 0; j < nb; j++) {
				fprintf(tst, " %g", fi  (a, b, nb, j, a + i * dx));
				fprintf(tst, " %g", dfi (a, b, nb, j, a + i * dx));
				fprintf(tst, " %g", d2fi(a, b, nb, j, a + i * dx));
				fprintf(tst, " %g", d3fi(a, b, nb, j, a + i * dx));
			}
			fprintf(tst, "\n");
		}
		fclose(tst);
	}
#endif

	for (j = 0; j < nb; j++) {
		for (i = 0; i < nb; i++)
			for (k = 0; k < pts->n; k++)
				add_to_entry_matrix(eqs, j, i, fi(a, b, nb, i, x[k]) * fi(a, b, nb, j, x[k]));

		for (k = 0; k < pts->n; k++)
			add_to_entry_matrix(eqs, j, nb, y[k] * fi(a, b, nb, j, x[k]));
	}

#ifdef DEBUG
	write_matrix(eqs, stdout);
#endif

	if (piv_ge_solver(eqs)) {
		spl->n = 0;
		return;
	}
#ifdef DEBUG
	write_matrix(eqs, stdout);
#endif

	if (alloc_spl(spl, nb) == 0) {
		for (i = 0; i < spl->n; i++) {
			double xx = spl->x[i] = a + i*(b-a)/(spl->n-1);
			xx+= 10.0*DBL_EPSILON;  // zabezpieczenie przed ulokowaniem punktu w poprzednim przedziale
			spl->f[i] = 0;
			spl->f1[i] = 0;
			spl->f2[i] = 0;
			spl->f3[i] = 0;
			for (k = 0; k < nb; k++) {
				double		ck = get_entry_matrix(eqs, k, nb);
				spl->f[i]  += ck * fi  (a, b, nb, k, xx);
				spl->f1[i] += ck * dfi (a, b, nb, k, xx);
				spl->f2[i] += ck * d2fi(a, b, nb, k, xx);
				spl->f3[i] += ck * d3fi(a, b, nb, k, xx);
			}
		}
	}

#ifdef DEBUG
	{
		FILE           *tst = fopen("debug_spline_plot.txt", "w");
		double		dx = (b - a) / (TESTBASE - 1);
		for (i = 0; i < TESTBASE; i++) {
			double yi= 0;
			double dyi= 0;
			double d2yi= 0;
			double d3yi= 0;
			double xi= a + i * dx;
			for( k= 0; k < nb; k++ ) {
							yi += get_entry_matrix(eqs, k, nb) * fi(a, b, nb, k, xi);
							dyi += get_entry_matrix(eqs, k, nb) * dfi(a, b, nb, k, xi);
							d2yi += get_entry_matrix(eqs, k, nb) * d2fi(a, b, nb, k, xi);
							d3yi += get_entry_matrix(eqs, k, nb) * d3fi(a, b, nb, k, xi);
			}
			fprintf(tst, "%g %g %g %g %g\n", xi, yi, dyi, d2yi, d3yi );
		}
		fclose(tst);
	}
#endif

}
