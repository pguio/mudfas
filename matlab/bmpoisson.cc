#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
using namespace std;
#include "bmpoisson.hxx"
#include "bmutil.hxx"

double **irho[NGMAX+1], **irhs[NGMAX+1], **iu[NGMAX+1],
**irnd[NGMAX+1], **ibnd[NGMAX+1];

void interp(double **uf, double **uc, int nf[])
{
  //
  //Coarse to fine prolongation by linear interpolation
  // **uf fine-grid solution returned
  // **uc coarse-grid input
  // nf[] dimensions of the fine-grid solution
  //
  int ic,iif,jc,jf,nc[2];

  for (ic=0;ic<2;ic++) nc[ic]=nf[ic]/2; // Coarse-grid dimension

  for (jc=0,jf=0;jc<=nc[1];jc++,jf+=2) // Do elements that are copies
    for (ic=0;ic<=nc[0];ic++)
      uf[2*ic][jf] = uc[ic][jc];

  for (jf=0;jf<=nf[1];jf+=2) // Do odd-numbered columns,
    for (iif=1;iif<nf[0];iif+=2) // Interpolating vertically
      uf[iif][jf] = 0.5*(uf[iif+1][jf]+uf[iif-1][jf]);

  for (jf=1;jf<nf[1];jf+=2) // Do even-numbered columns,
    for (iif=0;iif <= nf[0];iif++) // Interpolating horizontally
      uf[iif][jf] = 0.5*(uf[iif][jf+1]+uf[iif][jf-1]);
}

void rstrct(double **uc, double **uf, int nc[])
{
  //
  // Half-weighting restriction
  //
  // **uc coarse-grid solution
  // **uf fine-grid input
  // nc[] coarse grid dimension
  //
  int ic,iif,jc,jf,ncc[2];

  for (ic=0;ic<2;ic++) ncc[ic]=2*nc[ic]; // Fine-grid dimension

  // Interior points
  for (jf=2,jc=1;jc<nc[1];jc++,jf+=2) {  // Loop over y
    for (iif=2,ic=1; ic<nc[0]; ic++,iif+=2) { // Loop over x
      uc[ic][jc]=0.5*uf[iif][jf]
                 +0.125*(uf[iif+1][jf]+uf[iif-1][jf]
                         +uf[iif][jf+1]+uf[iif][jf-1]);
    }
  }
  for (jc=0,ic=0; ic<=nc[0]; ic++,jc+=2) { // Boundary points
    uc[ic][0]=uf[jc][0]; // ymin
    uc[ic][nc[1]]=uf[jc][ncc[1]]; // ymax
  }


  for (jc=2,ic=1; ic<nc[1]; ic++,jc+=2) { // Boundary points
    uc[0][ic]=uf[0][jc]; // xmin
    uc[nc[0]][ic]=uf[ncc[0]][jc]; // xmax
  }


}

void s_rstrct(double **uc, double **uf, int nc[])
{
  //
  // Half-weighting restriction
  //
  // **uc coarse-grid solution
  // **uf fine-grid input
  // nc[] coarse grid dimension
  //
  int ic,iif,jc,jf,ncc[2];

  for (ic=0; ic<2; ic++) ncc[ic]=2*nc[ic]; // Fine-grid dimension

  // Interior points
  for (jf=2,jc=1; jc<nc[1]; jc++,jf+=2) { // Loop over y
    for (iif=2,ic=1; ic<nc[0]; ic++,iif+=2) { // Loop over x
      uc[ic][jc]=0.5*uf[iif][jf]
                 +0.125*(uf[iif+1][jf]+uf[iif-1][jf]
                         +uf[iif][jf+1]+uf[iif][jf-1]);
    }
  }
  for (jc=0,ic=0; ic<=nc[0]; ic++,jc+=2) { // Boundary points
    uc[ic][0]=uf[jc][0]; // ymin
    uc[ic][nc[1]]=uf[jc][ncc[1]]; // ymax
  }


  for (jc=2,ic=1; ic<nc[1]; ic++,jc+=2) { // Boundary points
    uc[0][ic]=uf[0][jc]; // xmin
    // xmax
    uc[nc[0]][ic]= 0.5*uf[ncc[0]][jc]+0.25*uf[ncc[0]-1][jc]
                   +0.125*(uf[ncc[0]][jc+1]+uf[ncc[0]][jc-1]);
  }
}

void s1_rstrct(double **uc, double **uf, int nc[])
{
  //
  // Half-weighting restriction
  //
  // **uc coarse-grid solution
  // **uf fine-grid input
  // nc[] coarse grid dimension
  //
  int ic,iif,jc,jf,ncc[2];

  for (ic=0; ic<2; ic++) ncc[ic]=2*nc[ic]; // Fine-grid dimension

  // Interior points
  for (jf=2,jc=1; jc<nc[1]; jc++,jf+=2) { // Loop over y
    for (iif=2,ic=1; ic<nc[0]; ic++,iif+=2) { // Loop over x
      uc[ic][jc]=0.5*uf[iif][jf]
                 +0.125*(uf[iif+1][jf]+uf[iif-1][jf]
                         +uf[iif][jf+1]+uf[iif][jf-1]);
    }
  }
  for (jc=0,ic=0; ic<=nc[0]; ic++,jc+=2) { // Boundary points
    uc[ic][0] = uf[jc][0]; // ymin
    uc[ic][nc[1]] = uf[jc][ncc[1]]; // ymax
  }


  for (jc=0,ic=0; ic<=nc[1]; ic++,jc+=2) { // Boundary points
    uc[0][ic] = uf[0][jc]; // xmin
    uc[nc[0]][ic] = uf[ncc[0]][jc]; // xmax
  }


}

void smooth(double **uc, double **uf, int nc)
{
  //
  // Smooth the boundary from fine to coarse grid
  // **uc Coarse-grid boundary solution
  // **uf Fine-grid boundary input
  // nc Coarse-grid dimension
  //
  int i, j, jj, nnc;
  double vkt = 0.25;

  nnc = 2*nc; // Fine-grid dimension

  // Corner points
  uc[0][0]=uc[1][0]  = (1.-2*vkt)*uf[0][0]+vkt*(uf[0][1]+uf[1][1]);
  uc[1][nc]=uc[2][0] = (1.-2*vkt)*uf[1][nnc]+vkt*(uf[1][nnc-1]+uf[2][1]);
  uc[2][nc]=uc[3][nc]= (1.-2*vkt)*uf[2][nnc]+vkt*(uf[2][nnc-1]+uf[3][nnc-1]);
  uc[0][nc]=uc[3][0] = (1.-2*vkt)*uf[0][nnc]+vkt*(uf[0][nnc-1]+uf[3][1]);

  for (i=0;i<4; i++) // Loop over the 4  boundaries
    for (j=1;j<nc;j++) { // Inside
      jj = 2*j;
      uc[i][j] = (1.-2*vkt)*uf[i][jj]+vkt*(uf[i][jj-1]+uf[i][jj+1]);
    }
}

int inklud(int dim, int til)
{
  //
  // return TRUE if (-1)^dim+til == 1 or 0
  // return FALSE otherwise
  //    if ((ipow(-1,dim) + til)==1) return 1;
  //    if ((ipow(-1,dim) + til)==0) return 1;
  //		return 0;
  //
  if ( (Minus1PowN(dim)+til)==1 || (Minus1PowN(dim)+til)==0 ) return 1;
  return 0;
}

void GausS_DvN(double **a, double **b, double **rh,
               double (*Te)(int i, int j, int imax, int jmax),
               double ll[], int nl[])
{
  //
  // Gauss-Seidel relaxation scheme
  //
  // **a potential input and solution
  // **b boundary
  // **rh rho
  // (*Te) temperature
  // ll size of the grid
  // nl number of grid points
  //
  static const char *func = "GausS_DvN";
  int s, t, k, nytst;
  double bst, err, kk[2], k2, ggi[2], gg;
  int svts;
  int it=0;

  cout<<func<<": nl=("<<nl[0]<<","<<nl[1]<<")"<<endl;

  for (s=0; s<2; s++) kk[s] = ll[s]/nl[s]; // Grid increment dx, dy

  // Calculate the weights of the Laplacian
  k2 = 1.0/(2.0*(pow2(kk[0])+pow2(kk[1])));
  for (s=0; s<2; s++) ggi[s] = pow2(kk[1-s])*k2;
  gg  = pow2(kk[0]*kk[1])*k2;

  for (t=0; t<=nl[1]; t++) { // Initialise boundary
    a[0][t] = b[1][t]; // xmin
    a[nl[0]][t] = b[2][t]; // xmax
  }


  do { // Loop over Gauss-Seidel iteration
    err=0;
    for (svts=0; svts <= 1; svts++) { // Odd-even ordering
      // Calculate the new potential at the boundary ymin
      for (t=1; (2*t-svts)<nl[0]; t++) { // Loop over x
        k = 2*t-svts; // Calculate the index along y axis
        bst = a[k][0]; // save the old potential
        a[k][0] =
          2*ggi[1]*(a[k][1]-kk[1]*b[0][k])
          +ggi[0]*(a[k+1][0]+a[k-1][0])
          +gg*(rh[k][0]
               -exp(a[k][0]/Te(k,0,nl[0],nl[1])));
        // Residual at grid point [k][0]
        err += fabs(a[k][0]-bst);
      }
      // Calculate the new potential inside the boundary
      for (t=1; t<nl[1]; t++) { // Loop over y
        // if (t-svts odd) nytst=0
        // else nytst=-1
        nytst = - (t-svts) % 2;
        for (s=1; (2*s+nytst)<nl[0]; s++) { // Loop over x
          k = 2*s+nytst; // Calculate the index along x axis
          bst=a[k][t]; // save the old potential
          a[k][t] =
            ggi[0]*(a[k+1][t]+a[k-1][t])
            +ggi[1]*(a[k][t+1]+a[k][t-1])
            +gg*(rh[k][t]
                 -exp(a[k][t]/Te(k,t,nl[0],nl[1])));
          // Residual at grid point [k][t]
          err += fabs(a[k][t]-bst);
        }
      }
      nytst = !inklud(nl[1],svts);
      // Calculate the new potential at the boundary ymax
      for (s=1; (2*s-nytst)<nl[0]; s++) { // Loop over x
        k = 2*s-nytst;
        bst=a[k][nl[1]];
        a[k][nl[1]] =
          ggi[0]*(a[k+1][nl[1]]+a[k-1][nl[1]])
          +2*ggi[1]*(a[k][nl[1]-1]+kk[1]*b[3][k])
          +gg*(rh[k][nl[1]]
               -exp(a[k][nl[1]]/Te(k,nl[1],nl[0],nl[1])));
        // Residual at grid point [k][nl[1]]
        err += fabs(a[k][nl[1]]-bst);
      }
    }
    it++;
    err /= (nl[0]-1)*(nl[1]+1); // Average the residual
  } while ( it <= maxit && err > minerr);
  cout<<"it "<<it<<" (maxit "<<maxit<<")"<<
  " err "<<err<<" (minerr "<<minerr<<")"<<endl;
}

void s_GausS_DvN(double **a, double **b, double **rh,
                 double (*Te)(int i, int j, int imax, int jmax),
                 double ll[], int nl[])
{
  //
  // Gauss-Seidel relaxation scheme
  //
  // **a potential input and solution
  // **b boundary
  // **rh rho
  // (*Te) temperature
  // ll size of the grid
  // nl number of grid points
  //
  static const char *func="s_GausS_DvN";
  int s, t, k, nytst;
  double bst, err, kk[2], k2, ggi[2], gg;
  int svts;
  int it=0;

  cout<<func<<": nl=("<<nl[0]<<","<<nl[1]<<")"<<endl;
  for (s=0; s<2; s++) kk[s] = ll[s]/nl[s]; // Grid increment dx, dy

  // Calculate the weights for the Laplacian
  k2 = 1/(2.*(pow2(kk[0])+pow2(kk[1])));
  for (s=0; s<2; s++) ggi[s] = pow2(kk[1-s])*k2;
  gg  = pow2(kk[0]*kk[1])*k2;

  for (t=0; t<=nl[0]; t++) a[t][0] = b[0][t]; // Set up boundary at ymin
  for (t=0; t<=nl[1]; t++) a[0][t] = b[1][t]; // Set up boundary at xmin
  do { // Gauss-Seidel iteration loop
    err=0;
    for (svts = 0; svts<=1; svts ++) {
      for (t=1; t<nl[1]; t++) {
        // nytst = ipow(-1,(t-svts));
        // if (nytst == 1) nytst = 0;
        nytst = - (t-svts)%2;
        for (s=1; (2*s+nytst)<nl[0]; s++) {
          k = 2*s+nytst;
          bst=a[k][t];
          a[k][t] =
            ggi[0]*(a[k+1][t]+a[k-1][t])
            +ggi[1]*(a[k][t+1]+a[k][t-1])
            +gg*(rh[k][t]
                 -exp(a[k][t]/Te(k,t,nl[0],nl[1])));
          err += fabs(a[k][t]-bst);
        }
      }
      nytst = !inklud(nl[0],svts);
      for (s=1; (2*s-nytst)<nl[1]; s++) {
        k = 2*s-nytst;
        bst=a[nl[0]][k];
        a[nl[0]][k] =
          2.*ggi[0]*a[nl[0]-1][k]
          +ggi[1]*(a[nl[0]][k-1]+a[nl[0]][k+1])
          +gg*(rh[nl[0]][k]
               -exp(a[nl[0]][k]/Te(nl[0],k,nl[0],nl[1])));
        err += fabs(a[nl[0]][k]-bst);
      }
      nytst = !inklud(nl[1],svts);
      for (s=1; (2*s-nytst)<nl[0]; s++) {
        k = 2*s-nytst;
        bst=a[k][nl[1]];
        a[k][nl[1]] =
          ggi[0]*(a[k+1][nl[1]]+a[k-1][nl[1]])
          +2.*ggi[1]*(a[k][nl[1]-1]+kk[1]*b[3][k])
          +gg*(rh[k][nl[1]]
               -exp(a[k][nl[1]]/Te(k,nl[1],nl[0],nl[1])));
        err += fabs(a[k][nl[1]]-bst);
      }
      if (inklud(nl[0]+nl[1],svts)) {
        bst=a[nl[0]][nl[1]];
        a[nl[0]][nl[1]] =
          2.*ggi[0]*a[nl[0]-1][nl[1]]
          +2.*ggi[1]*(a[nl[0]][nl[1]-1]+kk[1]*b[3][nl[0]])
          +gg*(rh[nl[0]][nl[1]]
               -exp(a[nl[0]][nl[1]]/Te(nl[0],nl[1],nl[0],nl[1])));
        err += fabs(a[nl[0]][nl[1]]-bst);
      }
    }
    it++;
    err /= nl[0]*nl[1];
  } while ( it<=maxit && err>minerr );
  cout<<"it "<<it<<" (maxit "<<maxit<<")"<<
  " err "<<err<<" (minerr "<<minerr<<")"<<endl;
}

void s1_GausS_DvN(double **a, double **b, double **rh,
                  double (*Te)(int i, int j, int imax, int jmax),
                  double ll[], int nl[])
{
  //
  // Gauss-Seidel relaxation scheme
  //
  // **a potential input and solution
  // **b boundary
  // **rh rho
  // (*Te) temperature
  // ll size of the grid
  // nl number of grid points

  static const char *func = "s1_GausS_DvN";
  int s, t, k, d, nytst, kx[2]={0,nl[0]}, sign[2]={1,-1};
  double bst, err, kk[2], k2, ggi[2], gg, min1err = .1*minerr;
  int svts;
  int it=0;

  cout<<func<<": nl=("<<nl[0]<<","<<nl[1]<<")"<<endl;
  for (s=0; s<2; s++) kk[s] = ll[s]/nl[s]; // grid increment dx, dy

  // Calculate the weights for the Laplacian
  k2 = 1/(2.*(pow2(kk[0])+pow2(kk[1])));
  for (s=0; s<2; s++) ggi[s] = pow2(kk[1-s])*k2;
  gg  = pow2(kk[0]*kk[1])*k2;

  for (t=0; t<=nl[0]; t++) a[t][0] = b[0][t]; // Set up boundary at ymin
  d = 1; // Select y direction
  do { // Gauss-Seidel loop
    err=0;
    d = 1-d; // Switch direction
    for (svts=0; svts<=1; svts ++) {
      nytst = -!inklud(kx[d],svts);
      for (s=1; (2*s+svts)<nl[1]; s++) {
        k = 2*s+nytst;
        bst=a[kx[d]][k];
        a[kx[d]][k] =
          2.*ggi[0]*a[kx[d]+sign[d]][k]
          +ggi[1]*(a[kx[d]][k-1]+a[kx[d]][k+1])
          +gg*(rh[kx[d]][k]
               -exp(a[kx[d]][k]/Te(kx[d],k,nl[0],nl[1])));
        err += fabs(a[kx[d]][k]-bst);
      }
      if (inklud(kx[d]+nl[1],svts)) {
        bst=a[kx[d]][nl[1]];
        a[kx[d]][nl[1]] =
          2.*ggi[0]*a[kx[d]+sign[d]][nl[1]]
          +2.*ggi[1]*(a[kx[d]][nl[1]-1]+kk[1]*b[3][kx[d]])
          +gg*(rh[kx[d]][nl[1]]
               -exp(a[kx[d]][nl[1]]/Te(kx[d],nl[1],nl[0],nl[1])));
        err += fabs(a[kx[d]][nl[1]]-bst);
      }
      for (t=1; t<nl[1]; t++) {
        // nytst = ipow(-1,(t-svts));
        // if (nytst == 1) nytst = 0;
        nytst = - (t-svts)%2;
        for (s=1; (2*s+nytst)<nl[0]; s++) {
          k = 2*s+nytst;
          bst=a[k][t];
          a[k][t] =
            ggi[0]*(a[k+1][t]+a[k-1][t])
            +ggi[1]*(a[k][t+1]+a[k][t-1])
            +gg*(rh[k][t]
                 -exp(a[k][t]/Te(k,t,nl[0],nl[1])));
          err += fabs(a[k][t]-bst);
        }
      }
      nytst = -!inklud(kx[1-d],svts);
      for (s=1; (2*s+svts)<nl[1]; s++) {
        k = 2*s+nytst;
        bst=a[kx[1-d]][k];
        a[kx[1-d]][k] =
          2.*ggi[0]*a[kx[1-d]+sign[1-d]][k]
          +ggi[1]*(a[kx[1-d]][k-1]+a[kx[1-d]][k+1])
          +gg*(rh[kx[1-d]][k]
               -exp(a[kx[1-d]][k]/Te(kx[1-d],k,nl[0],nl[1])));
        err += fabs(a[kx[1-d]][k]-bst);
      }
      nytst = !inklud(nl[1],svts);
      for (s=1; (2*s-nytst)<nl[0]; s++) {
        k = 2*s-nytst;
        k = kx[d] + sign[d]*k;
        bst=a[k][nl[1]];
        a[k][nl[1]] =
          ggi[0]*(a[k+1][nl[1]]+a[k-1][nl[1]])
          +2.*ggi[1]*(a[k][nl[1]-1]+kk[1]*b[3][k])
          +gg*(rh[k][nl[1]]
               -exp(a[k][nl[1]]/Te(k,nl[1],nl[0],nl[1])));
        err += fabs(a[k][nl[1]]-bst);
      }
      if (inklud(kx[1-d]+nl[1],svts)) {
        bst=a[kx[1-d]][nl[1]];
        a[kx[1-d]][nl[1]] =
          2.*ggi[0]*a[kx[1-d]+sign[1-d]][nl[1]]
          +2.*ggi[1]*(a[kx[1-d]][nl[1]-1]+kk[1]*b[3][kx[1-d]])
          +gg*(rh[kx[1-d]][nl[1]]
               -exp(a[kx[1-d]][nl[1]]/Te(kx[1-d],nl[1],nl[0],nl[1])));
        err += fabs(a[kx[1-d]][nl[1]]-bst);
      }
    }
    it++;
    err /= (nl[0]+1)*nl[1];
  } while ( it<=maxit && err>min1err );
  cout<<"it "<<it<<" (maxit "<<maxit<<")"<<
  " err "<<err<<" (min1err "<<min1err<<")"<<endl;
}

void mgfas(double **u, double **r, double **b,
           double (*T)(int i, int j, int imax, int jmax),
           double lx[], int n[])
{
  //
  // Multi Grid Full Approximation Storage Algorithm
  // Solves the nonlinear elliptic equation
  // d^2u/dx^2+d^2u/dy^2+exp(u(x,y)/Te)=r(x,y)
  // with the following boundary conditions
  // u(xmin)=0
  // and
  // u(xmax)=0
  // Number of grids set to cn
  //
  // Variables:
  // **u estimated potential
  // **r given density
  // **b boundary
  // (*T) temperature
  // lx size of the grid
  // n number of grid points
  //
  static const char *func = "mgfas";
  int i, j, nn[2], ngmin[2], ib, ngrid = cn;

  cout<<func<<": n=("<<n[0]<<","<<n[1]<<")"<<endl;

  if (n[0] >= n[1]) ib = 0; // Boundary index set to x
  else ib = 1; // Boundary index set to y
  for (i=0; i<2; i++) { // Loop over dimension
    nn[i] = n[i]/2;  // Coarsest grid dimension=finest/2
    ngmin[i] = n[i]; // Finest grid dimension
  }
  for (i=0; i<ngrid; i++) // Loop over multigrids
    for (j=0; j<2; j++) // Loop over dimension
      ngmin[j]/= 2; // Coarsest grid dimension /2
  cout<<"Initialise coarser grid 0"<<endl;
  irnd[ngrid]=dmatrix(0,3,0,nn[ib]); // Boundary array
  //cout<<"---> irnd["<<ngrid<<"]=ALLOC("<<3<<","<<nn[ib]<<")"<<endl;
  smooth(irnd[ngrid],b,nn[ib]); // Smooth boundary
  irho[ngrid]=dmatrix(0,nn[0],0,nn[1]); //Allocate storage for density
  //cout<<"---> irho["<<ngrid<<"]=ALLOC("<<nn[0]<<","<<nn[1]<<")"<<endl;
  rstrct(irho[ngrid],r,nn); // Fine (r) to coarse grid (irho[ngrid])
  iu[ngrid]=dmatrix(0,nn[0],0,nn[1]); // Allocate storage for potential
  //cout<<"---> iu["<<ngrid<<"]=ALLOC("<<nn[0]<<","<<nn[1]<<")"<<endl;
  rstrct(iu[ngrid],u,nn); // Fine (u) to coarse grid (iu[ngrid])
  for (i=1; i<ngrid; i++) { // Loop over grid
    //cout<<"Initialise coarser grid "<<i<<endl;
    for (i=0; i<2; i++) nn[i]/= 2; // Coarsest grid dimension /2
    irnd[--ngrid]=dmatrix(0,3,0,nn[ib]); //  Boundary array
    //cout<<"---> irnd["<<ngrid<<"]=ALLOC("<<3<<","<<nn[ib]<<")"<<endl;
    smooth(irnd[ngrid],irnd[ngrid+1],nn[ib]); // Smooth boundary
    irho[ngrid]=dmatrix(0,nn[0],0,nn[1]); // Allocate storage for density
    //cout<<"---> irho["<<ngrid<<"]=ALLOC("<<nn[0]<<","<<nn[1]<<")"<<endl;
    rstrct(irho[ngrid],irho[ngrid+1],nn); // Fine to coarse grid
    iu[ngrid]=dmatrix(0,nn[0],0,nn[1]); // Allocate storage for potential
    //cout<<"---> iu["<<ngrid<<"]=ALLOC("<<nn[0]<<","<<nn[1]<<")"<<endl;
    rstrct(iu[ngrid],iu[ngrid+1],nn); // Fine to coarse grid
  }
  for (i=0; i<2; i++) nn[i] = ngmin[i]; //
  ngrid = cn+1;
  iu[ngrid]=dmatrix(0,n[0],0,n[1]); //
  //cout<<"---> iu["<<ngrid<<"]=ALLOC("<<n[0]<<","<<n[1]<<")"<<endl;
  for (j=1;j<=ngrid;j++) { // From coarsest to finest grid
    irhs[j]=dmatrix(0,nn[0],0,nn[1]);
    //cout<<"---> irhs["<<j<<"]=ALLOC("<<nn[0]<<","<<nn[1]<<")"<<endl;
    ibnd[j]=dmatrix(0,3,0,nn[ib]);
    //cout<<"---> ibnd["<<j<<"]=ALLOC("<<3<<","<<nn[ib]<<")"<<endl;
    // WARNING!!!!
    // WARNING!!!! ERROR IF CN>=3 ==> copy NULL pointer
    // WARNING!!!! (same in s_mgfas and s1_mgfas)
    // WARNING!!!!
    // if ( (j != ngrid ? irnd[j] : b)==NULL )
    //	ErrorHandler(func,__FILE__,__LINE__);
    copy4(ibnd[j],(j != ngrid ? irnd[j] : b),nn[ib]); // Set up boundary
    // Interpolate from coarse grid to next finer grid
    if (j > 1) interp(iu[j],iu[j-1],nn);
    //if ( (j != ngrid ? irho[j] : r)==NULL )
    //	ErrorHandler(func,__FILE__,__LINE__);
    copy(irhs[j],(j != ngrid ? irho[j] : r),nn); // Set up r
    GausS_DvN(iu[j],ibnd[j],irhs[j],T,lx,nn); // Gauss-Seidel relaxation
    for (i=0; i<2; i++) nn[i]*= 2;
  }
  copy(u,iu[ngrid],n); // Copy final solution
  for (i=0; i<2; i++) nn[i] = n[i];
  for (j=cn+1; j>=1; j--) { // Loop to deallocate memory
    free_dmatrix(irhs[j],0,nn[0],0,nn[1]);
    //cout<<"---> irhs["<<j<<"]=FREE("<<nn[0]<<","<<nn[1]<<")"<<endl;
    free_dmatrix(ibnd[j],0,3,0,nn[ib]);
    //cout<<"---> ibnd["<<j<<"]=FREE("<<3<<","<<nn[ib]<<")"<<endl;
    free_dmatrix(iu[j],0,nn[0],0,nn[1]);
    //cout<<"---> iu["<<j<<"]=FREE("<<nn[0]<<","<<nn[1]<<")"<<endl;
    if (j != cn+1) {
      free_dmatrix(irho[j],0,nn[0],0,nn[1]);
      //cout<<"---> irho["<<j<<"]=FREE("<<nn[0]<<","<<nn[1]<<")"<<endl;
    }
    if (j != cn+1) {
      free_dmatrix(irnd[j],0,3,0,nn[ib]);
      //cout<<"---> irnd["<<j<<"]=FREE("<<3<<","<<nn[ib]<<")"<<endl;
    }
    for (i=0; i<2; i++) nn[i]/= 2;
  }
}

void s_mgfas(double **u, double **r, double **b,
             double (*T)(int i, int j, int imax, int jmax),
             double lx[], int n[])
{
  //
  // Multi Grid Full Approximation Storage Algorithm
  // Solves the nonlinear elliptic equation
  // d^2u/dx^2+d^2u/dy^2+exp(u(x,y)/Te)=r(x,y)
  // with the following boundary conditions
  // u(xmin)=0
  // and
  // u(ymin)=0
  // Number of grids set to cn
  //
  // Variables:
  // **u estimated potential
  // **r given density
  // **b boundary
  // (*T) temperature
  // lx size of the grid
  // n number of grid points
  //
  static const char *func = "s_mgfas";
  int i, j, nn[2], ngmin[2], ib, ngrid = cn;

  cout<<func<<": n=("<<n[0]<<","<<n[1]<<")"<<endl;
  if (n[0] >= n[1]) ib = 0;
  else ib = 1;
  for (i=0; i<2; i++) {
    nn[i] = n[i]/2;
    ngmin[i] = n[i];
  }
  for (i=0; i<ngrid; i++)
    for (j=0; j<2; j++) ngmin[j]/= 2;
  irnd[ngrid]=dmatrix(0,3,0,nn[ib]);
  smooth(irnd[ngrid],b,nn[ib]);
  irho[ngrid]=dmatrix(0,nn[0],0,nn[1]);
  s_rstrct(irho[ngrid],r,nn);
  iu[ngrid]=dmatrix(0,nn[0],0,nn[1]);
  s_rstrct(iu[ngrid],u,nn);
  for (i=1; i<ngrid; i++) {
    for (i=0; i<2; i++) nn[i]/= 2;
    irnd[--ngrid]=dmatrix(0,3,0,nn[ib]);
    smooth(irnd[ngrid],irnd[ngrid+1],nn[ib]);
    irho[ngrid]=dmatrix(0,nn[0],0,nn[1]);
    s_rstrct(irho[ngrid],irho[ngrid+1],nn);
    iu[ngrid]=dmatrix(0,nn[0],0,nn[1]);
    s_rstrct(iu[ngrid],iu[ngrid+1],nn);
  }
  for (i=0; i<2; i++) nn[i] = ngmin[i];
  ngrid = cn+1;
  iu[ngrid]=dmatrix(0,n[0],0,n[1]);
  for (j=1;j<=ngrid;j++) {
    irhs[j]=dmatrix(0,nn[0],0,nn[1]);
    ibnd[j]=dmatrix(0,3,0,nn[ib]);
    copy4(ibnd[j],(j != ngrid ? irnd[j] : b),nn[ib]);
    if (j > 1) interp(iu[j],iu[j-1],nn);
    copy(irhs[j],(j != ngrid ? irho[j] : r),nn);
    s_GausS_DvN(iu[j],ibnd[j],irhs[j],T,lx,nn);
    for (i=0; i<2; i++) nn[i]*= 2;
  }
  copy(u,iu[ngrid],n);
  for (i=0; i<2; i++) nn[i] = n[i];
  for (j=cn+1; j>=1; j--) {
    free_dmatrix(irhs[j],0,nn[0],0,nn[1]);
    free_dmatrix(ibnd[j],0,3,0,nn[ib]);
    free_dmatrix(iu[j],0,nn[0],0,nn[1]);
    if (j != cn+1) free_dmatrix(irho[j],0,nn[0],0,nn[1]);
    if (j != cn+1) free_dmatrix(irnd[j],0,3,0,nn[ib]);
    for (i=0; i<2; i++) nn[i]/= 2;
  }
}


void s1_mgfas(double **u, double **r, double **b,
              double (*T)(int i, int j, int imax, int jmax),
              double lx[], int n[])
{
  //
  // Multi Grid Full Approximation Storage Algorithm
  // Solves the nonlinear elliptic equation
  // d^2u/dx^2+d^2u/dy^2+exp(u(x,y)/Te)=r(x,y)
  // with the following boundary conditions
  // u(ymin)=0
  // Number of grids set to cn
  //
  // Variables:
  // **u estimated potential
  // **r given density
  // **b boundary
  // (*T) temperature
  // lx size of the grid
  // n number of grid points
  //
  static const char *func = "s1_mgfas";
  int i, j, nn[2], ngmin[2], ib, ngrid = cn;

  cout<<func<<": n=("<<n[0]<<","<<n[1]<<")"<<endl;
  if (n[0] >= n[1]) ib = 0;
  else ib = 1;
  for (i=0; i<2; i++) {
    nn[i] = n[i]/2;
    ngmin[i] = n[i];
  }
  for (i=0; i<ngrid; i++)
    for (j=0; j<2; j++) ngmin[j]/= 2;
  irnd[ngrid]=dmatrix(0,3,0,nn[ib]);
  smooth(irnd[ngrid],b,nn[ib]);
  irho[ngrid]=dmatrix(0,nn[0],0,nn[1]);
  s1_rstrct(irho[ngrid],r,nn);
  iu[ngrid]=dmatrix(0,nn[0],0,nn[1]);
  s1_rstrct(iu[ngrid],u,nn);
  for (i=1; i<ngrid; i++) {
    for (i=0; i<2; i++) nn[i]/= 2;
    irnd[--ngrid]=dmatrix(0,3,0,nn[ib]);
    smooth(irnd[ngrid],irnd[ngrid+1],nn[ib]);
    irho[ngrid]=dmatrix(0,nn[0],0,nn[1]);
    s1_rstrct(irho[ngrid],irho[ngrid+1],nn);
    iu[ngrid]=dmatrix(0,nn[0],0,nn[1]);
    s1_rstrct(iu[ngrid],iu[ngrid+1],nn);
  }
  for (i=0; i<2; i++) nn[i] = ngmin[i];
  ngrid = cn+1;
  iu[ngrid]=dmatrix(0,n[0],0,n[1]);
  for (j=1;j<=ngrid;j++) {
    irhs[j]=dmatrix(0,nn[0],0,nn[1]);
    ibnd[j]=dmatrix(0,3,0,nn[ib]);
    copy4(ibnd[j],(j != ngrid ? irnd[j] : b),nn[ib]);
    if (j > 1) interp(iu[j],iu[j-1],nn);
    copy(irhs[j],(j != ngrid ? irho[j] : r),nn);
    s1_GausS_DvN(iu[j],ibnd[j],irhs[j],T,lx,nn);
    for (i=0; i<2; i++) nn[i]*= 2;
  }
  copy(u,iu[ngrid],n);
  for (i=0; i<2; i++) nn[i] = n[i];
  for (j=cn+1; j>=1; j--) {
    free_dmatrix(irhs[j],0,nn[0],0,nn[1]);
    free_dmatrix(ibnd[j],0,3,0,nn[ib]);
    free_dmatrix(iu[j],0,nn[0],0,nn[1]);
    if (j != cn+1) free_dmatrix(irho[j],0,nn[0],0,nn[1]);
    if (j != cn+1) free_dmatrix(irnd[j],0,3,0,nn[ib]);
    for (i=0; i<2; i++) nn[i]/= 2;
  }
}
