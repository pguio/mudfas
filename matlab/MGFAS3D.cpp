
/********************************************************************
 *
 * $Id: MGFAS3D.cpp,v 1.5 2011/03/26 12:56:39 patrick Exp $
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
#include "bmpoisson3D.hxx"

#define FREE_ARG char*
#define NR_END 1

double  ll[3];
double Temin, Temax; // Electron temperatures
double sigmaT; // Standard deviation of Te

int ipow(int k, int n)
{

  if (n==0) return 1;
  if (n==1) return k;
  return k*ipow(k,n-1) ;
}

void copy3D(double ***aout, double ***ain, int na[])
{
  int i,j,k;

  for (i=0; i<=na[0]; i++)
    for (j=0; j<=na[1]; j++)
      for (k=0; k<=na[2]; k++)
        aout[i][j][k] = ain[i][j][k];

}

void ErrorHandler(const char FunctName[], const char FileName[], int line)
{
  fprintf(stderr,"%s\n\tDetected in %s() [%s line %d]\n", strerror(errno),
          FunctName, FileName, line);
}

double ***dtensor(int nts[])
{
//long i,j,nrow=nrh-nrl+1,ncol=nch-ncl+1,ndep=ndh-ndl+1;
  long nrl=0,nrh=nts[0],ncl=0,nch=nts[1],ndl=0,ndh=nts[2];
  long i,j,k,nrow=nrh-nrl+1,ncol=nch-ncl+1,ndep=ndh-ndl+1;
  double ***t;
  int nr_end = 0;

  /* allocate pointers to pointers to rows */
  t=(double ***) mxCalloc(nrow+nr_end,sizeof(double**));
  if (!t) cout<<"allocation failure 1 in f3tensor()"<<endl;
  t += nr_end;
  t -= nrl;

  /* allocate pointers to rows and set pointers to them */
  t[nrl]=(double **) mxCalloc(nrow*ncol+nr_end,sizeof(double*));
  if (!t[nrl]) cout<<"allocation failure 2 in f3tensor()"<<endl;
  t[nrl] += nr_end;
  t[nrl] -= ncl;

  /* allocate rows and set pointers to them */
  t[nrl][ncl]=(double *) mxCalloc(nrow*ncol*ndep+nr_end,sizeof(double));
  if (!t[nrl][ncl]) cout<<"allocation failure 3 in f3tensor()"<<endl;
  t[nrl][ncl] += nr_end;
  t[nrl][ncl] -= ndl;

  for (j=ncl+1; j<=nch; j++) t[nrl][j]=t[nrl][j-1]+ndep;
  for (i=nrl+1; i<=nrh; i++) {
    t[i]=t[i-1]+ncol;
    t[i][ncl]=t[i-1][ncl]+ncol*ndep;
    for (j=ncl+1; j<=nch; j++) t[i][j]=t[i][j-1]+ndep;
  }

  for (i=0; i <= nts[0]; i++)
    for (j=0; j <= nts[1]; j++)
      for (k=0; k <= nts[2]; k++) {
        t[i][j][k] = 0.0;
      }
  return t;
}

void free_dtensor(double ***ts)
{
  mxFree(ts[0][0]);
  mxFree(ts[0]);
  mxFree(ts);
}

double **boundary3D(int ll[])
{
  int i;
//  boudary comndition:
  double **b;
//cout<<"double **boundary3D  DIM :"<<ll[0]<<"  "<<ll[1]<<"  "<<ll[2]<<endl;

  b = (double**)mxCalloc(6,sizeof(double*));

  b[0] = (double*)mxCalloc((ll[0]+1)*(ll[1]+1),sizeof(double)); //(x,y) z = 0
  b[1] = (double*)mxCalloc((ll[0]+1)*(ll[1]+1),sizeof(double)); //(x,y) z = zmax
  b[2] = (double*)mxCalloc((ll[2]+1)*(ll[1]+1),sizeof(double)); //(y,z) x = 0
  b[3] = (double*)mxCalloc((ll[2]+1)*(ll[1]+1),sizeof(double)); //(y,z) x = xmax
  b[4] = (double*)mxCalloc((ll[0]+1)*(ll[2]+1),sizeof(double)); //(x,z) y = 0
  b[5] = (double*)mxCalloc((ll[0]+1)*(ll[2]+1),sizeof(double)); //(x,y) y = ymax

  for ( i =0; i< (ll[0]+1)*(ll[1]+1); i++) {
    *(b[0]+i) = 0.0;
    *(b[1]+i) = 0.0;
  }

  for ( i =0; i< (ll[2]+1)*(ll[1]+1); i++) {
    *(b[2]+i) = 0.0;
    *(b[3]+i) = 0.0;
  }

  for ( i =0; i< (ll[0]+1)*(ll[2]+1); i++) {
    *(b[4]+i) = 0.0;
    *(b[5]+i) = 0.0;
  }

  return b;
}

void copy_boundary3D(double **b,double **bnd,int ll[])
{
  int i;
//  boudary comndition:

  for ( i =0; i< (ll[0]+1)*(ll[1]+1); i++) {
    *(b[0]+i) = *(bnd[0]+i);
    *(b[1]+i) = *(bnd[1]+i);
  }

  for ( i =0; i< (ll[2]+1)*(ll[1]+1); i++) {
    *(b[2]+i) = *(bnd[2]+i);
    *(b[3]+i) = *(bnd[3]+i);
  }

  for ( i =0; i< (ll[0]+1)*(ll[2]+1); i++) {
    *(b[4]+i) = *(bnd[4]+i);
    *(b[5]+i) = *(bnd[5]+i);
  }

}

void free_boundary3D(double **bnd)
{
  for ( int i = 0; i <6; i++)
    mxFree(bnd[i]);

  mxFree(bnd);
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{

  double ***phi, ***rho, ***rho0, **bnd;
  int n2[3];
  const int *n;
  int i, j, k;
  double *ptr1, *ptr2;

  n=mxGetDimensions(prhs[0]);
  for (i=0; i<3; i++) {
    n2[i]=n[i]-1;
    //cout<<"n2["<<i<<"]="<<n2[i]<<" ";
  }
  //cout<<endl;

  phi=dtensor(n2);
  rho=dtensor(n2);
  rho0=dtensor(n2);

  bnd = boundary3D(n2);

  ptr1=mxGetPr(prhs[0]);
  for (k=0; k<=n2[2]; k++)
    for (j=0; j<=n2[1]; j++)
      for (i=0; i<=n2[0]; i++) {
        rho[i][j][k] = *ptr1++;
        rho0[i][j][k] = 1.0;
        phi[i][j][k] = 0.0;
      }


//		cout<<"rho[0][0][0]= "<<rho[0][0][0]<<endl;
//		cout<<"rho[1][0][0]= "<<rho[1][0][0]<<
//			    "\trho[end][0][0]= "<<rho[n2[0]][0][0]<<endl;
//		cout<<"rho[0][1][0]= "<<rho[0][1][0]<<
//					"\trho[0][end][0]= "<<rho[0][n2[1]][0]<<endl;
//		cout<<"rho[0][0][1]= "<<rho[0][0][1]<<
//					"\trho[0][0][end]= "<<rho[0][0][n2[2]]<<endl;

  ll[0]=*(mxGetPr(prhs[1])+1);
  ll[1]=*(mxGetPr(prhs[1])+3);
  ll[2]=*(mxGetPr(prhs[1])+5);

//		cout<<"ll = "<<ll[0]<<" "<<ll[1]<<" "<<ll[2]<<endl;

  Temin=*mxGetPr(prhs[2]);
  Temax=*(mxGetPr(prhs[2])+1);
  sigmaT=*(mxGetPr(prhs[2])+2);

//		cout<<"Temin "<<Temin<<" Temax "<<Temax<<" sigmaT "<<sigmaT<<endl;
  mgfas(phi,rho,bnd,rho0,Temin,ll,n2);

  plhs[0] = mxCreateNumericArray(3,n,mxDOUBLE_CLASS,mxREAL);
//		plhs[1] = mxCreateNumericArray(3,n,mxDOUBLE_CLASS,mxREAL);
  ptr1=mxGetPr(plhs[0]);
//		ptr2=mxGetPr(plhs[1]);
  for (k=0; k<=n2[2]; k++)
    for (j=0; j<=n2[1]; j++)
      for (i=0; i<=n2[0]; i++) {
        *ptr1++=phi[i][j][k];
//					*ptr2++=rho[i][j][k];
      }

  free_dtensor(rho);
  free_dtensor(rho0);
  free_dtensor(phi);
  free_boundary3D(bnd);

}
