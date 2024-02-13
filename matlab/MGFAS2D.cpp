
/********************************************************************
 *
 * $Id: MGFAS2D.cpp,v 1.7 2011/03/26 12:56:39 patrick Exp $
 *
********************************************************************/

#include <stdio.h>
#include <math.h>
#include <float.h>
#include <errno.h>
#include <string.h>
#include <iostream>
#include <fstream>
using namespace std;
#include "mex.h"
#include "bmpoisson.hxx"

#define FREE_ARG char*
#define NR_END 1

double  ll[2];
double Temin, Temax; // Electron temperatures
double sigmaT; // Standard deviation of Te

int pw(int x, int i)
{
  int fx;
  fx = x;
  if (i==0) return 1;
  for (int j = 1; j < i; j++)
    fx*= x;
  return fx;
}

int ipow(int k, int n)
{

  if (n==0) return 1;
  if (n==1) return k;
  return k*ipow(k,n-1) ;
}

void copy(double **aout, double **ain, int na[])
{
  int i,j;
  for (i=0; i<=na[1]; i++)
    for (j=0; j<=na[0]; j++)
      aout[j][i] = ain[j][i];
}

void copy4(double **aout, double **ain, int n)
{
  int i,j;
  for (i=0; i<=n; i++)
    for (j=0; j<=3; j++)
      aout[j][i]=ain[j][i];
}

void ErrorHandler(const char FunctName[], const char FileName[], int line)
{
  fprintf(stderr,"%s\n\tDetected in %s() [%s line %d]\n", strerror(errno),
          FunctName, FileName, line);
}

double **dmatrix(long nrl, long nrh, long ncl, long nch)
{
  static char func[] = "dmatrix";
  long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
  double **m;

  m=(double **) mxCalloc(nrow+NR_END,sizeof(double *));

  if (!m) ErrorHandler(func,__FILE__,__LINE__);
  m += NR_END;
  m -= nrl;
  m[nrl]=(double *) mxCalloc(nrow*ncol+NR_END,sizeof(double));
  if (!m[nrl]) ErrorHandler(func,__FILE__,__LINE__);
  m[nrl] += NR_END;
  m[nrl] -= ncl;
  for (i=nrl+1; i<=nrh; i++) m[i]=m[i-1]+ncol;
  return m;
}

void free_dmatrix(double **m, long nrl, long nrh, long ncl, long nch)
{
  mxFree((FREE_ARG) (m[nrl]+ncl-NR_END));
  mxFree((FREE_ARG) (m+nrl-NR_END));
}

double Te(int j, int i, int jmax, int imax)
{
  double x;
  x=2.0*(double)i/(double)imax-1.0; // Grid [0,imax] <==> x [-1,1]
  x /= sigmaT; // x scaling
  return (Temax-Temin)*exp(-x*x/2.0)+Temin;
}

double boundary(int s, int t)
{
  double fx;
  fx = 0.0;
  return fx;
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{

  double **phi, **rho, **bnd;
  int nr, nc, n2[2], nbnd;
  int i, j;
  double *ptr1, *ptr2;
  int mgfasflag;

  n2[0]=nr=mxGetM(prhs[0])-1;
  n2[1]=nc=mxGetN(prhs[0])-1;
  phi=dmatrix(0,n2[0],0,n2[1]);
  rho=dmatrix(0,n2[0],0,n2[1]);

  if (n2[0] >= n2[1]) nbnd = n2[0];
  else nbnd = n2[1];
  bnd = dmatrix(0,3,0,nbnd);
  for (i = 0; i < 4; i++)
    for (j = 0; j <= nbnd; j++)
      bnd[i][j] = boundary(i,j);

  ptr1=mxGetPr(prhs[0]);
  for (j=0; j<=nc; j++)
    for (i=0; i<=nr; i++)
      rho[i][j] = *ptr1++;

  ll[0]=*(mxGetPr(prhs[1])+1);
  ll[1]=*(mxGetPr(prhs[1])+3);

  Temin=*mxGetPr(prhs[2]);
  Temax=*(mxGetPr(prhs[2])+1);
  sigmaT=*(mxGetPr(prhs[2])+2);

  cout<<"Temin "<<Temin<<" Temax "<<Temax<<" sigmaT "<<sigmaT<<endl;

  mgfasflag=(int)*mxGetPr(prhs[3]);
  switch (mgfasflag) {
  case 0:
    mgfas(phi,rho,bnd,&Te,ll,n2);
    break;
  case 1:
    s_mgfas(phi,rho,bnd,&Te,ll,n2);
    break;
  case 2:
    s1_mgfas(phi,rho,bnd,&Te,ll,n2);
    break;
  }

  plhs[0] = mxCreateDoubleMatrix(nr+1, nc+1, mxREAL);
  plhs[1] = mxCreateDoubleMatrix(nr+1, nc+1, mxREAL);
  ptr1=mxGetPr(plhs[0]);
  ptr2=mxGetPr(plhs[1]);
  for (j=0; j<=nc; j++)
    for (i=0; i<=nr; i++) {
      *ptr1++=phi[i][j];
      *ptr2++=Te(i,j,nr,nc);
    }

  free_dmatrix(rho,0,n2[0],0,n2[1]);
  free_dmatrix(phi,0,n2[0],0,n2[1]);
  free_dmatrix(bnd,0,3,0,nbnd);

}
