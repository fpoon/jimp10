#include <stdio.h>
#include <math.h>

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
  for(i = 0; i <= n/2; i++)
  {
    ret += newton(n,i)*newton(2*n-2*i,n)*pow(x,n-2*i)*((i%2)==0?1:-1);
  }
  
  return ret*factor;
}
/*
int main()
{
  int n,k;
  double x;
  scanf("%lf %d",&x, &n);
  printf("%g\n",legendre(x,n));
  return 0;
}*/