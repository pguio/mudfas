#ifndef _BMUTIL_HXX
#define _BMUTIL_HXX

#define FREE_ARG char*
#define NR_END 1

#include <stdlib.h>
#include <limits.h>
#include "Types.hxx"

// Constants for random generator
#ifdef INT_MAX
long const a = 16807, c = 0, m = INT_MAX;
#elif RAND_MAX
long const a = 16807, c = 0, m = RAND_MAX;
#else
long const a = 16807, c = 0, m = 2147483647;
#endif
int const erfk = 4, ritmax = 15,  kkmax = 20000;

#define POW2(x) ((x)*(x)) // x^2
#define M_SQRT2PI (2.0*M_SQRT2/M_2_SQRTPI) // sqrt(2 pi)
#define M_SQRTPI (2.0/M_2_SQRTPI) // sqrt(pi)
#define M_SQRT6 2.4494897427831779  // sqrt(6)

inline double maxwl(double x, double vm, double vt)
{
  return fabs(x)*exp(-POW2((x-vm)/vt)/2.0);
}
inline double pow2(double x)
{
  return x*x;
}
inline int Minus1PowN(int n)
{
  return (((n+1)%2)<<1)-1;
}
inline long randU(long seed)
{
  return (a*seed+c)&m;
}
inline double scale(long seed, double lo, double hi)
{
  return (double)seed*(hi-lo)/m + lo;
}
inline double Fz(double z, double v, double s)
{
  return 0.5*v*erf(z) - s*exp(-POW2(z))/M_SQRT2PI;
}
inline double dFz(double z, double v, double s)
{
  return (v+M_SQRT2*s*z)*exp(-POW2(z))/M_SQRTPI;
}

int ArgParsing(int argc, char *argv[], const char *Option, double *Value);

int ipow(int k, int n);
int **imatrix(long nrl, long nrh, long ncl, long nch);
int ***itensor(int nts[]);
double ***dtensor(int nts[]);
double **dmatrix(long nrl, long nrh, long ncl, long nch);
void free_dmatrix(double **m, long nrl, long nrh, long ncl, long nch);
void free_imatrix(int **m, long nrl, long nrh, long ncl, long nch);
void free_itensor(int ***ts);
void free_dtensor(double ***ts);
#ifdef TWO_D
void ifill0(int **u, int nx[]);
void fill0(double **u, int nx[]);
#endif
#ifdef THREE_D
void ifill0(int ***u, int nx[]);
void fill0(double ***u, int nx[]);
#endif
void ifill30(int ***u, int nx[]);
void ErrorHandler(const char FunctName[], const char FileName[], int line);
void copy4(double **aout, double **ain, int n);
void copy3D(double ***aout, double ***ain, int na[]);
void copy(double **aout, double **ain, int na[]);
void box_mueller(PhaseSpaceType s[], long &seed, double smean[], double sstd[]);
double integrate(double u, double vth, double vmx);
double integrate1(double u, double vth, double vmx);
double rayleigh0(long &seed, double sstd);
double rayleigh(long &seed, double u, double uth, double c, double Au, double z0);
void uniform(PhaseSpaceType s[], long &seed, double **sb);
void symmetrize(double **ain, double **aout, int nin[]);
void unsymmetrize(double **ain, double **aout, int nout[]);
double **boundary3D(int *ll);
void copy_boundary3D(double **b,double **bnd,int ll[]);
void free_boundary3D(double **bnd);
double ***HDFformat0_1(double ***m, int nts[]);
void tilfil3D(double ***m, int nts[],char tekst[]);
void tilfil3D_byte(double ***m, int nts[],char tekst[]);
void boundary3D_tester(double **b,int ll[]);
void from_rho_grid(double ***a,double ***s,int ll[]);
double  ***to_rho_grid(double ***b,int ll[]);
void graf(int m,char tekst[],int step);
#endif


