#ifdef THREE_D
#include <math.h>
#include <stdlib.h>
#include <iostream>
using namespace std;
#include "bmpoisson3D.hxx"
#include "bmutil.hxx"


//-----------------------------------------------------------
void interp(double ***uf, double ***uc, int nf[])
//-----------------------------------------------------------
{
  // coarse => fine grid
  cout<<"void interp(double ***uf, double ***uc, int nf[])"<<endl;


  int ic,iif,jjf,jc,jf,kf,kc,nc[3];

  for (ic=0;ic<3;ic++) {
    nc[ic]=nf[ic]/2;
    //cout<<"W "<<nc[ic]<<"  "<<nf[ic]<<endl;
  }

  for (kc=0,kf=0;kc<=nc[2];kc++,kf+=2)
    for (jc=0,jf=0;jc<=nc[1];jc++,jf+=2)
      for (ic=0;ic<=nc[0];ic++) {
        //cout<<"W "<<ic<<"  "<<jf<<"  "<<kf<<"  "<<jc<<"  "<<kc<<endl;
        uf[2*ic][jf][kf]=uc[ic][jc][kc];
        //if( uf[2*ic][jf][kf] < 0.5 ) cout<<" inetr A "<<ic<<"  "<<jf<<"  "<<kf<<endl;
      }
  for (kc=0,kf=0;kc<=nc[2];kc++,kf+=2) {
    for (jf=0;jf<=nf[1];jf+=2)
      for (iif=1;iif<nf[0];iif+=2) {
        uf[iif][jf][kf]=0.5*(uf[iif+1][jf][kf]+uf[iif-1][jf][kf]);
        //if( uf[iif][jf][kf] < 0.5 ) cout<<" inetr B "<<iif<<"  "<<jf<<"  "<<kf<<endl;
      }
    for (jf=1;jf<nf[1];jf+=2)
      for (iif=0;iif <=nf[0];iif++) {
        uf[iif][jf][kf]= 0.5 * ( uf[iif][jf+1][kf]+uf[iif][jf-1][kf]);
        //if( uf[iif][jf][kf] < 0.5 ) cout<<" inetr C "<<iif<<"  "<<jf<<"  "<<kf<<endl;
      }


  }

  for (kf=1;kf<nf[2];kf+=2)
    for (iif=0;iif <=nf[0];iif++)
      for (jjf=0;jjf <=nf[1];jjf++) {
        uf[iif][jjf][kf]=0.5*(uf[iif][jjf][kf+1]+uf[iif][jjf][kf-1]);
        //if( uf[iif][jjf][kf] < 0.5 ) cout<<" inetr D "<<iif<<"  "<<jjf<<"  "<<kf<<endl;
      }

  //tilfil3D(uf,nf,"uf.dta");
  //tilfil3D(uc,nc,"uc.dta");

}
//-----------------------------------------------------------
void rstrct(double ***uc, double ***uf, int nc[])
//-----------------------------------------------------------
{
  //  fine => coarse grid

  int ic,iif,jc,jf,kf,kc,ncc[3];

  for (ic=0;ic<3;ic++)
    ncc[ic]=2*nc[ic];

  for (kf=2,kc=1;kc<nc[2];kc++,kf+=2) {
    for (jf=2,jc=1;jc<nc[1];jc++,jf+=2) {
      for (iif=2,ic=1; ic<nc[0]; ic++,iif+=2) {
        //cout<<"void rstrct(  "<<ic<<" "<<jc<<" "<<kc<<endl;
        uc[ic][jc][kc] =0.5*uf[iif][jf][kf]
                        + 0.08333*(uf[iif+1][jf][kf]+uf[iif-1][jf][kf])
                        + 0.08333*(uf[iif][jf+1][kf]+uf[iif][jf-1][kf])
                        + 0.08333*(uf[iif][jf+1][kf]+uf[iif][jf][kf-1]);
      }
    }
  }
  //tilfil(uf,ncc,"uf.dta");

  for (kc=0,kf=0;kc<=nc[2];kc++,kf+=2) {

    for (jc=0,ic=0; ic<=nc[0]; ic++,jc+=2) {
      uc[ic][0][kc]= uf[jc][0][kf];
      uc[ic][nc[1]][kc]=uf[jc][ncc[1]][kf];
//  	  cout <<"ncc[0]="<<nc[0]
//  	       <<"  ncc[1]="<<nc[1]
//  	       <<"  ncc[2]="<<nc[2]
//  	       <<"  kf:"<<kf
//  	       <<"  ic: "<<ic<<"    jc: "<<jc<<"  "<<uf[jc][ncc[1]][kf]<<endl;

    }

    for (jc=0,ic=0; ic<=nc[1]; ic++,jc+=2) {
      uc[0][ic][kc]=uf[0][jc][kf];
      uc[nc[0]][ic][kc]=uf[ncc[0]][jc][kf];
    }

  }


  for (kc=0,kf=0;kc<=nc[0];kc++,kf+=2)
    for (jc=0,ic=0; ic<=nc[1]; ic++,jc+=2) {
      uc[kc][ic][0]= uf[kf][jc][0];
      uc[kc][ic][nc[2]]=uf[kf][jc][ncc[2]];
    }
  //      for (jc=2,ic=1; ic<nc[1]; ic++,jc+=2)
  //      for (jc=0,ic=0; ic<=nc[0]; ic++,jc+=2)

}


//-----------------------------------------------------------
void smooth(double **uc, double **uf, int nc[])
//-----------------------------------------------------------
{
  int i, j, jj,k,kk, nnc[3],lx,ly,lz;
  double vkt = 1.0/6.0;

  for ( int t = 0; t < 3; t++)
    nnc[t] = 2*nc[t];

  // corners whith  y =0

  uc[0][0]=uc[2][0]=uc[4][0] =
                      (1.0- 3.0*vkt)* uf[0][0]
                      + vkt*(uf[0][1]+
                             uf[0][nnc[0]+1]+
                             uf[4][nnc[0]+1]);

  uc[0][nc[0]]=uc[4][nc[0]]=uc[3][0] =
                              (1.0- 3.0*vkt)* uf[0][nnc[0]]
                              + vkt*(uf[0][nnc[0]-1]+
                                     uf[3][1]+
                                     uf[3][nnc[2]+1]);

  uc[4][(nc[0]+1)*(nc[2])-1]=uc[2][(nc[1]+1)*(nc[2])-1]=uc[1][0]=
                               (1.0- 3.0*vkt)* uf[1][0]
                               +vkt*( uf[1][1]+
                                      uf[1][nnc[0]+1]+
                                      uf[4][(nnc[0]+1)*nnc[2]]);

  uc[1][nc[0]] = uc[4][(nc[0]+1)*(nc[2]+1)-1] = uc[3][(nc[1]+1)*(nc[2])-1]=
                   (1.0- 3.0*vkt)* uf[1][nnc[0]]
                   +vkt*(uf[1][nnc[0]-1]+
                         uf[1][2*(nnc[0]+1)-1]+
                         uf[4][(nnc[0]+1)*(nnc[2]+1)-1]);

  if ( uc[1][nc[0]] != 0.0) {
    cout<<" ER DET HER A"<<uc[1][nc[0]]<<endl;
//        //cout<<k<<"  "<<j<<"  "<<nc[2]<<endl;
//        exit(1);
  }


  // corners whith y = ymax


  uc[0][(nc[0]+1)*(nc[1])-1]=uc[2][nc[1]+1]=uc[5][0] =
                               (1.0- 3.0*vkt)* uf[5][0]
                               + vkt*(uf[5][1] +
                                      uf[0][(nnc[0]+1)*nnc[1] - 1] +
                                      uf[5][nnc[0]+1]);

  uc[0][(nc[0]+1)*(nc[1]+1)-1]=uc[5][(nc[0]+1)]=uc[3][(nc[1]+1)] =
                                 (1.0- 3.0*vkt)* uf[5][nnc[0]]
                                 + vkt*(uf[0][(nnc[1]+1)*(nnc[0]+1)-2]+
                                        uf[3][nnc[1]]+
                                        uf[3][2*(nnc[1]+1)-1]);
  if ( uc[5][(nc[0]+1)] != 0.0) {
    cout<<" ER DET HER R"<<uc[5][(nc[0]+1)]<<endl;
//        //cout<<k<<"  "<<j<<"  "<<nc[2]<<endl;
    exit(1);
  }

  uc[5][(nc[0]+1)*(nc[2]-1)-1]=uc[2][(nc[1]+1) * (nc[2]+1)-1]=uc[1][(nc[0]+1)*nc[1]]=
                                 (1.0- 3.0*vkt)* uf[1][(nnc[0]+1)*nnc[1]]
                                 +vkt*( uf[1][(nnc[0]+1)*(nnc[1]-1)]+
                                        uf[1][(nnc[0]+1)*(nnc[1])+1]+
                                        uf[5][(nnc[0]+1)*(nnc[2]-1)]);

  if ( uc[5][(nc[0]+1)*(nc[2]-1)-1] != 0.0) {
    cout<<" ER DET HER Q"<<uc[5][(nc[0]+1)*(nc[2]-1)-1]<<endl;
    cout<<(nc[0]+1)*(nc[2]-1)-1<<endl;
    exit(1);
  }

  uc[1][(nc[1]+1)*(nc[0]+1)-1]=uc[5][(nc[0]+1) * (nc[2]+1)-1]=uc[3][(nc[1]+1) * (nc[2]+1)-1]=
                                 (1.0- 3.0*vkt)* uf[1][(nnc[1]+1)*(nnc[0]+1)-1]
                                 +vkt*(uf[1][(nnc[1]+1)*(nnc[0]+1)-2]+
                                       uf[1][(nnc[1]+1)*nnc[0]-1]+
                                       uf[5][(nnc[0]+1)*nnc[2]-1]);

  if ( uc[5][(nc[0]+1)*(nc[2]-1)-1] != 0.0) {
    cout<<" ER DET HER Z"<<uc[5][(nc[0]+1)*(nc[2]-1)-1]<<endl;
    exit(1);
  }

  vkt = 1.0/6.0;


  for ( j = 1; j < nc[0]; j++) {
    jj = 2*j;

    uc[4][j] = uc[0][j] =
                 (1.-4.0*vkt)*uf[4][jj]
                 +vkt*( uf[4][jj+1] +
                        uf[4][jj-1] +
                        uf[4][nnc[0]+1 +jj] +
                        uf[0][nnc[0]+jj]);

    uc[4][j+(nc[0]+1)*nc[2]] = uc[1][j] =
                                 (1.-4.0*vkt)*uf[4][jj+(nc[0]+1)*nc[2]]
                                 +vkt*( uf[4][jj+1+(nc[0]+1)*nc[2]] +
                                        uf[4][jj-1+(nc[0]+1)*nc[2]]+
                                        uf[4][nnc[0] +jj]+
                                        uf[1][nnc[0]+jj]);
    if ( uc[4][j] != 0.0 ) {
      cout<<" ER DET HER H  "<<uc[4][j]<<endl;
//  	  cout<<";aksjdflkjsadlhf"<< uc[4][j+(nc[0]+1)*nc[2]]<<endl;
      exit(1);
    }

    uc[5][j] = uc[0][j+(nc[0]+1)*nc[1]] =
                 (1.-4.0*vkt)*uf[5][jj]
                 +vkt*( uf[5][jj+1] +
                        uf[5][jj-1] +
                        uf[5][nnc[0]+1 +jj] +
                        uf[0][(nc[0]+1)*(nc[1]-1)+jj]);

    if ( uc[5][j] != 0.0 ) {
      cout<<" ER DET HER 0  "<<uc[5][j]<<endl;
      exit(1);
    }

    uc[5][j+(nc[0]+1)*nc[2]] = uc[0][j+(nc[0]+1)*nc[1]] =
                                 (1.-4.0*vkt)*uf[5][jj+(nc[0]+1)*nc[2]]
                                 +vkt*( uf[5][jj+1+(nc[0]+1)*nc[2]] +
                                        uf[5][jj-1+(nc[0]+1)*nc[2]]+
                                        uf[5][(nnc[0]+1)*(nnc[2]-1) +jj] +
                                        uf[0][(nc[0]+1)*(nc[1]-1)+jj]);

    if ( uc[5][j+(nc[0]+1)*nc[2]] != 0.0 ) {
      cout<<" ER DET HER P  "<<uc[5][j+(nc[0]+1)*nc[2]]<<endl;
      exit(1);
    }

  }

  for (j= 1; j<(nc[1]); j++) {

    jj= j*2;

    uc[0][j*(nc[0]+1)]=uc[2][j]=
                         (1.-4.0*vkt)*uf[2][jj]
                         +vkt*( uf[2][jj-1]+
                                uf[2][jj+1]+
                                uf[2][jj+nnc[2]+1] +
                                uf[0][(jj+1)*(nnc[0]+1)]);

    if (  uc[2][j] != 0.0) {
      cout<<" ER DET HER M "<< uc[2][j]<<endl;
    }
    uc[3][j] = uc[0][j*(nc[0]+1)-1]=
                 (1.-4.0*vkt)*uf[3][jj]
                 +vkt*( uf[2][jj-1]+
                        uf[2][jj+1]+
                        uf[2][jj+nnc[2]+1] +
                        uf[0][(jj+1)*(nnc[0]+1)-1]);

    uc[1][j*(nc[0]+1)]=uc[2][j+(nc[2])*(nc[1]+1)]=
                         (1.-4.0*vkt)*uf[2][jj+(nc[2])*(nc[1]+1)]
                         +vkt*( uf[2][jj-1+(nnc[2])*(nnc[1]+1)]+

                                uf[2][jj+1+(nnc[2])*(nnc[1]+1)]+

                                uf[2][jj+(nnc[2]-1)*(nnc[1]+1)] +
                                uf[1][(jj+1)*(nnc[0]+1)]);


    if (uc[1][j*(nc[0]+1)]  != 0.0) {
      cout<<" ER DET HER B "<< uc[1][j*(nc[0]+1)]<<endl;
//  	    cout<<uf[2][jj+(nc[2])*(nc[1]+1)]<<endl;
//  	    cout<< uf[2][jj-1+(nnc[2]+1)*nnc[1]]<<endl;
//  	    cout<<uf[2][jj+1+(nnc[2]+1)*nnc[1]]<<endl; //<===================
//  	    cout<< uf[2][jj+1+nnc[2]*nnc[1]]<<endl;
//  	    cout<< uf[1][(jj+1)*(nnc[0]+1)]<<endl;
//  	    cout<<jj<<"  "<<j<<"  "<<nc[1]<<"  "<<nnc[1]<<" "<<jj+1+(nnc[2]+1)*nnc[1]<<endl;
//  	    exit(1);
    }

    uc[1][j*(nc[0]+1)-1]=uc[3][j+(nc[2])*(nc[1]+1)]=
                           (1.-4.0*vkt)*uf[3][jj+(nc[2])*(nc[1]+1)]
                           +vkt*( uf[3][jj-1+(nnc[2])*(nnc[1]+1)]+
                                  uf[3][jj+1+(nnc[2])*(nnc[1]+1)]+
                                  uf[3][jj+1+(nnc[2]-1)*(nnc[1]+1)] +
                                  uf[1][jj*(nnc[0]+1) -1]);

    if (uc[1][j*(nc[0]+1)-1]  != 0.0) {
      cout<<" ER DET HER C"<<uc[1][j*(nc[0]+1)-1]<<endl;

//  	  cout<<"  "<<j<<"  "<<nc[2]<<endl;
      //  	  exit(1);
    }
  }

  //cout<<" mitten void smooth "<<endl;
  for (j= 1; j < nc[2]; j++) {
    jj= j*2;

    uc[4][(nc[0]+1)] = uc[2][(nc[1]+1)*j]=
                         (1.-4.0*vkt)*uf[4][(nnc[0]+1)*jj]
                         +vkt*( uf[4][(nnc[0]+1)*(jj-1)]+
                                uf[4][(nnc[0]+1)*(jj+1)]+
                                uf[4][1+(nnc[0]+1)*jj]+
                                uf[2][1+(nnc[1]+1)*jj]);

    if (  uc[4][(nc[0]+1)]  != 0.0) {
      cout<<" ER DET HER D"<<uc[4][(nc[0]+1)]<<endl;
    }

    uc[4][nc[0]+(nc[0]+1)*j] = uc[3][(nc[1]+1)*j]=
                                 (1.-4.0*vkt)*uf[0][nnc[0]+(nnc[0]+1)*jj]
                                 +vkt*( uf[4][nnc[0]+(nnc[0]+1)*jj]+
                                        uf[4][nnc[0] +(nnc[0]+1)*(jj+1)]+
                                        uf[4][nnc[0]-1 +(nnc[0]+1)*jj]+
                                        uf[3][1 + (nnc[1]+1)*jj] );

    if ( uc[4][nc[0]+(nc[0]+1)*j]  != 0.0) {
      cout<<" ER DET HER E"<<endl;
    }

    uc[5][(nc[0]+1)*j] = uc[2][nc[1]+(nc[1]+1)*j]=
                           (1.-4.0*vkt)*uf[5][(nnc[0]+1)*jj]
                           +vkt*( uf[5][(nnc[0]+1)*(jj+1)]+
                                  uf[5][(nnc[0]+1)*(jj-1)]+
                                  uf[5][1+(nnc[0]+1)*jj]+
                                  uf[2][nnc[1]-1 + (nnc[1]+1)*jj]);

    if (  uc[5][(nc[0]+1)*j] != 0.0) {
      cout<<" ER DET HER K "<<  uc[5][(nc[0]+1)*j]<<endl;
      exit(1);
    }

    uc[5][nc[0] + (nc[0]+1)*j]= uc[3][nc[1]+(nc[1]+1)*j]=
                                  (1.-4.0*vkt)*uf[5][nnc[0]+(nnc[0]+1)*jj]
                                  +vkt*( uf[5][nnc[0]+(nnc[0]+1)*(jj-1)]+
                                         uf[5][nnc[0]+(nnc[0]+1)*(jj+1)]+
                                         uf[5][nnc[0]-1 +(nnc[0]+1)*jj]+
                                         uf[3][nnc[1]-1 +(nnc[1]+1)*jj]);

    if ( uc[5][nc[0] + (nc[0]+1)*j] != 0.0) {
      cout<<" ER DET HER K "<<uc[5][nc[0] + (nc[0]+1)*j]<<endl;
      exit(1);
    }

  }



  //cout<<" mitten void smooth "<<endl;

  lx = nnc[0]+1;
  ly = nnc[1]+1;
  lz = nnc[2]+1;

  vkt = 1.0/6.0;


  for (j=1;j<nc[1];j++)   /* (x,y) z = 0,zmax */
    for (k=1;k<nc[0];k++) {
      jj = 2*j;
      kk = 2*k;

      uc[0][k+j*(nc[0]+1)] = (1.-4.0*vkt)*uf[0][kk+jj*lx]
                             + vkt*(uf[0][kk+(jj-1)*lx]+uf[0][kk+(jj+1)*lx ])
                             + vkt*(uf[0][kk+1+(jj)*lx]+uf[0][kk+1+(jj)*lx ]);
      if ( uc[0][k+j*(nc[0]+1)] != 0.0) {
        cout<<" ER DET HER F  "<< uc[0][k+j*(nc[1]+1)]<<endl;
        exit(1);
      }

      uc[1][k+j*(nc[0]+1)] = (1.-4.0*vkt)*uf[1][kk+jj*lx]
                             + vkt*(uf[1][ kk+(jj-1)*lx ]+uf[1][kk+(jj+1)*lx ])
                             + vkt*(uf[1][ kk+1+(jj)*lx ]+uf[1][kk+1+(jj)*lx ]);

      if (uc[1][k+j*(nc[0] +1)] != 0.0) {
        cout<<" ER DET HER G  "<<uc[1][k+j*(nc[1] +1)]<<endl;
//  	    cout<<" ER DET HER II IGJEN"<<uc[4][k+j*nc[2]]<<endl;
//  	    cout<<k<<"  "<<j<<"  "<<nc[2]<<endl;
        exit(1);
      }


    }



  for (j=1;j<nc[2];j++)  /* (Y,Z) X = 0,Xmax */
    for (k=1;k<nc[1];k++) {
      jj = 2*j;
      kk = 2*k;
      uc[2][k+j*(nc[1]+1)] = (1.-4.0*vkt)*uf[2][kk+jj*ly]
                             + vkt*(uf[2][kk+(jj-1)*ly]+uf[2][kk+(jj+1)*ly ])
                             + vkt*(uf[2][kk+1+(jj)*ly]+uf[2][kk+1+(jj)*ly ]);

      if ( uc[2][k+j*(nc[1]+1)] != 0.0 ) {
        cout<<" ER DET HER L  "<<uc[2][k+j*(nc[1]+1)]<<endl;
        exit(1);
      }

      uc[3][k+j*(nc[1]+1)] = (1.-4.0*vkt)*uf[3][kk+jj*ly]
                             + vkt*(uf[3][kk+(jj-1)*ly]+uf[3][kk+(jj+1)*ly ])
                             + vkt*(uf[3][kk+1+(jj)*ly]+uf[3][kk+1+(jj)*ly ]);
    }


  for (j=1;j<nc[2];j++)  /* (X,Z) Y = 0,Ymax */
    for (k=1;k<nc[0];k++) {
      jj = 2*j;
      kk = 2*k;
      uc[4][k+j*(nc[0]+1)] = (1.-4.0*vkt)*uf[4][kk+jj*lx]
                             + vkt*(uf[4][kk+(jj-1)*lx]+uf[4][kk+(jj+1)*lx ])
                             + vkt*(uf[4][kk+1+(jj)*lx]+uf[4][kk+1+(jj)*lx ]);

      if (uc[4][k+j*(nc[0]+1)] != 0.0) {
        cout<<" ER DET HER I "<<uc[4][k+j*nc[0]]<<endl;
//  	    cout<<k<<"  "<<j<<"  "<<nc[2]<<endl;
        exit(1);
      }

      uc[5][k+j*(nc[0]+1)] = (1.-4.0*vkt)*uf[5][kk+jj*lx]
                             + vkt*(uf[5][kk+(jj-1)*lx]+uf[5][kk+(jj+1)*lx ])
                             + vkt*(uf[5][kk+1+(jj)*lx]+uf[5][kk+1+(jj)*lx ]);

      if (uc[5][k+j*(nc[0] +1)] != 0.0) {
        cout<<" ER DET HER J  "<<uc[5][k+j*nc[0]]<<endl;
//  	    cout<<k<<"  "<<j<<"  "<<nc[2]<<endl;
        exit(1);
      }
    }


  //  boundary3D_tester(uc,nc);


}

//-----------------------------------------------------------
int inklud(int dim, int til)
//-----------------------------------------------------------
{
  if ( (Minus1PowN(dim)+til)==1 || (Minus1PowN(dim)+til)==0 ) return 1;
  return 0;
}

//-----------------------------------------------------------
void GausS_DvN(double ***a, double **b, double ***rh, double ***rh0,  const double Te,
               double ll[], int nl[])
//-----------------------------------------------------------
{
  /*
     a  : phi
     b  : baundary condition
     rh : rho
  */
  int s, t, k,p, nytst;
  double bst,svts, err, kk[3], k2, ggi[3], gg, res;
  unsigned int it;

  err = 1.0e-2;
  //tilfil(rh,nl,"rh.dta");

  for (s=0; s<3; s++)
    kk[s] = ll[s]/nl[s];

  k2 = 1/(2.*( POW2(kk[0])*POW2(kk[1])
               +POW2(kk[1])*POW2(kk[2])
               +POW2(kk[2])*POW2(kk[0]) ));

  ggi[0] = ( POW2(kk[1])*POW2(kk[2]))*k2;
  ggi[1] = ( POW2(kk[0])*POW2(kk[2]))*k2;
  ggi[2] = ( POW2(kk[0])*POW2(kk[1]))*k2;

  gg  = POW2(kk[0]*kk[1]*kk[2])*k2;
  cout<<"gg/ggi  "<<ggi[0]+ggi[1]+ggi[2]<<"  "<<gg;

  err=minerr+1.0;
  it=0;
  nytst = 0;


  for (p=0; p<=nl[2]; p++)       //Dirichlet (X=0,Xmax) YZ plan
    for (t = 0; t <= nl[1]; t++) {

#ifndef SYMFIELD_THREED_X0
      a[0][t][p]     = b[2][ t + p*(nl[1]+1) ];
#endif
      a[nl[0]][t][p] = b[3][ t + p*(nl[1]+1) ];
//  	if( a[0][t][p]  != 0.0  || a[nl[0]][t][p] != 0.0  )
//  	  cout<<"###t "<<t<<"  p"<<p<<"  "<<t + p*nl[1] <<"  "<<nl[1]<<endl;
    }

  //tilfil3D(a,nl,"a.dta");

  for (p=0; p<=nl[1]; p++)       //Dirichlet (Z=0,Zmax) XY plan
    for (t = 0; t <= nl[0]; t++) {
#ifndef SYMFIELD_THREED_Z0
      a[t][p][0]     = b[0][ t + p*(nl[0]+1) ];
#endif
      a[t][p][nl[2]] = b[1][ t + p*(nl[0]+1) ];
//  	if(/* a[t][p][0] != 0.0 ||*/ a[t][p][nl[2]] != 0.0 ){
//  	  cout<<"%%%t "<<t<<"  p"<<p<<"  "<<t + p*nl[1] <<"  "<<nl[1]<<endl;
//  	  cout<< b[1][ t + p*nl[1] ]<<endl;
//  	}
    }


  //tilfil3D(a,nl,"a.dta");
  while (( it<=maxit) && (err>minerr )) {
    err=0;
    for ( svts = 0; svts <=  1; svts ++) {
      for (p=1; p<nl[2]; p++) {
#ifdef SYMFIELD_THREED_X0
        // von Neumann x =0 , zy
        a[0][0][p] = 2.0*ggi[0]*(a[1][0][p]-kk[0]*b[2][p*(nl[1]+1)]);
        a[0][0][p]+= 2.0*ggi[1]*(a[0][1][p]-kk[1]*b[4][p*(nl[0]+1)]);
        a[0][0][p]+= ggi[2]*(a[0][0][p+1]+a[0][0][p-1]);
        a[0][0][p]+= gg*(rh[0][0][p] - /*rh0[0][0][p]*/exp(a[0][0][p]/Te));
        //	     a[0][0][p] =3.0;

        a[0][nl[1]][p] = 2.0*ggi[0]*(a[1][nl[1]][p]-kk[0]*b[2][nl[1]+p*(nl[1]+1)]);
        a[0][nl[1]][p]+= 2.0*ggi[1]*(a[0][1][p] - kk[1]*b[5][p*(nl[0]+1)]);
        a[0][nl[1]][p]+= ggi[2]*(a[0][nl[1]][p+1]+a[0][nl[1]][p-1]);
        a[0][nl[1]][p]+= gg*(rh[0][nl[1]][p] - /*rh0[0][nl[1]][p]*/exp(a[0][nl[1]][p]/Te));
        //a[0][nl[1]][p] =3.0;

        for (t=1; k = 2*t - svts,k <= nl[1]; t++)    //x =0 , zy
        {                                 // von Neumann
          //k = 2*t - svts;
          bst = a[0][k][p];

          a[0][k][p] = 2.0*ggi[0]*(a[1][k][p]-kk[0]*b[2][k+p*(nl[1]+1)]);
          a[0][k][p]+= ggi[1]*(a[0][k+1][p]+a[0][k-1][p]);
          a[0][k][p]+= ggi[2]*(a[0][k][p+1]+a[0][k][p-1]);
          a[0][k][p]+= gg*(rh[0][k][p] - /*rh0[0][k][p]*/exp(a[0][k][p]/Te));

          err+=fabs(a[0][k][p] - bst);
          // a[0][k][p] = 3.0;
        }
#endif

        for (t=1; k = 2*t - (int)svts,k <nl[0]; t++)    //y =0 , zx
        {                                 // von Neumann
          //k = 2*t - svts;
          bst = a[k][0][p];
          a[k][0][p] = 2.0*ggi[1]*(a[k][1][p]-kk[1]*b[4][k+p*(nl[0]+1)]);
          a[k][0][p]+= ggi[0]*(a[k+1][0][p]+a[k-1][0][p]);
          a[k][0][p]+= ggi[2]*(a[k][0][p+1]+a[k][0][p-1]);
          a[k][0][p]+= gg*(rh[k][0][p] - /*rh0[k][0][p]*/exp(a[k][0][p]/Te));
          err+=fabs(a[k][0][p] - bst);
          //a[k][0][p]=0.0;
          //cout<<"y=0 k ="<<k<<"  , p ="<<p<<endl;
        }



        for (t=1; t<nl[1]; t++) { // indre punkter
          nytst = ipow(-1,(t+p-(int)svts));

          if (nytst == 1)
            nytst = 0;

          for (s=1; k = 2*s+nytst,k<nl[0]; s++) {
            //k = 2*s+nytst;
            bst=a[k][t][p];
            a[k][t][p] = ggi[0]*(a[k+1][t][p]+a[k-1][t][p]);
            a[k][t][p]+= ggi[1]*(a[k][t+1][p]+a[k][t-1][p]);
            a[k][t][p]+= ggi[2]*(a[k][t][p+1]+a[k][t][p-1]);
            a[k][t][p]+= gg*(rh[k][t][p] -/*rh0[k][t][p]*/exp(a[k][t][p]/Te));
            //if(  (rh0[k][t][p] < 0.999) || (rh0[k][t][p] > 1.001) ) {
            //  cout<<"her er noe galt rh0"<< rh0[k][t][p]<<"  "<<k<<"  "<<t<<"  "<<p<<endl;
            //}
            //cout<<"k= "<<k<<"  t= "<<t<<"  p= "<<p;
            //cout<<"    a[k][t][p] = "<<a[k][t][p]<<endl;
            err+=fabs(a[k][t][p] - bst);
            //a[k][t][p] =0.0;

            //k = 2*s+nytst;
//  		     //cout<<"k= "<<k<<"  t= "<<t<<"  p= "<<p;
//  		     bst=a[k][t][p];

//  		     res = ( ggi[0]*(a[k+1][t][p]+a[k-1][t][p])
//  			     + ggi[1]*(a[k][t+1][p]+a[k][t-1][p])
//  			     + ggi[2]*(a[k][t][p+1]+a[k][t][p-1])
//  			     - 6.0*gg*a[k][t][p]);
//  		     res -=  exp(a[k][t][p]/Te) + rh[k][t][p];
//  		     a[k][t][p] -=  res/((-6.0)*gg - 1.0/Te * exp(a[k][t][p]/Te));

//  		     //cout<<"    a[k][t][p] = "<<a[k][t][p]<<endl;
//  		     err+=fabs(a[k][t][p] - bst);
            //a[k][t][p] =0.0;
          }
        }


        nytst = !inklud(nl[1],(int)svts);

        for (s=1; k = 2*s-nytst,k<nl[0]; s++)   //iy =ny , ix varierer
        {                                 // von Neumann
          //k = 2*s-nytst;
          bst=a[k][nl[1]][p];
          //cout<<"y = max  k ="<<k<<"  , p ="<<p<<"           "<<a[k][nl[1]][p]<<endl;
          a[k][nl[1]][p] = ggi[0]*(a[k+1][nl[1]][p] + a[k-1][nl[1]][p]);
          a[k][nl[1]][p]+= 2.0*ggi[1]*(a[k][nl[1]-1][p] - kk[1]*b[5][k+p*(nl[0]+1)]);
          a[k][nl[1]][p]+= ggi[2]*(a[k][nl[1]][p+1] + a[k][nl[1]][p-1]);
          a[k][nl[1]][p]+= gg*(rh[k][nl[1]][p] - /*rh0[k][nl[1]][p]*/exp(a[k][nl[1]][p]/Te));

          //a[k][nl[1]][p]=10.0;
          if (b[5][k+p*(nl[0]+1)] != 0.0 ) cout<<"b[5][k+p*(nl[0]+1)]) ="<<b[5][k+p*(nl[0]+1)]<<endl;

          err+=fabs(a[k][nl[1]][p] - bst);

          //a[k][nl[1]][p]=0.5;
          //cout<<"y = max  k ="<<k<<"  , p ="<<p<<endl;
        }
      } //end p: nl[2] loop
      //exit(1);

#ifdef SYMFIELD_THREED_Z0
#ifdef SYMFIELD_THREED_X0

      a[0][0][0] = 2.0*ggi[2]*(a[0][0][1]-kk[2]*b[0][0]);
      a[0][0][0]+= 2.0*ggi[1]*(a[0][1][0]-kk[1]*b[4][0]);
      a[0][0][0]+= 2.0*ggi[0]*(a[1][0][0]-kk[0]*b[2][0]);
      a[0][0][0]+= gg*(rh[0][0][0] - /*rh0[0][0][0]*/exp(a[0][0][0]/Te));

      a[0][nl[1]][0] = 2.0*ggi[2]*(a[0][nl[1]][1]-kk[2]*b[0][nl[0]+1]);
      a[0][nl[1]][0]+= 2.0*ggi[1]*(a[0][1][0]-kk[1]*b[4][nl[0]+1]);
      a[0][nl[1]][0]+= 2.0*ggi[0]*(a[1][nl[1]][0]-kk[0]*b[2][nl[1]]);
      a[0][nl[1]][0]+= gg*(rh[0][nl[1]][0] -/*rh0[0][nl[1]][0]*/exp(a[0][nl[1]][0]/Te));

      for (s=1; k = 2*s - svts,k <= nl[1]; s++) {
        bst = a[0][k][0];

        a[0][k][0] = 2.0*ggi[2]*(a[0][k][1]-kk[2]*b[0][k*(nl[0]+1)]);
        a[0][k][0]+= ggi[1]*(a[0][k+1][0]+a[0][k-1][0]);
        a[0][k][0]+= 2.0*ggi[0]*(a[1][k][0]-kk[0]*b[2][k*(nl[2]+1)]);
        a[0][k][0]+= gg*(rh[0][k][0] - /*rh0[0][k][0]*/exp(a[0][k][0]/Te));

        err+=fabs(a[0][k][0] - bst);

        // a[t][k][0] =0.0;

      }




#endif
#endif




#ifdef SYMFIELD_THREED_Z0


      for (t=1; t<nl[0]; t++) { // von Neuman z=0 ,xy
        nytst = ipow(-1,(t-svts));

        //cout<<"void GausS_DvN(:a[t][0][0] t="<<t<<endl;
        //cout<<"void GausS_DvN(:a[t][0][0]="<<a[0][0][0]<<endl;
        a[t][0][0] = 2.0*ggi[2]*(a[t][0][1]-kk[2]*b[0][t]);
        a[t][0][0]+= 2.0*ggi[1]*(a[t][1][0]-kk[1]*b[4][t]);
        a[t][0][0]+= ggi[0]*(a[t+1][0][0]+a[t-1][0][0]);
        a[t][0][0]+= gg*(rh[t][0][0] - /*rh0[t][0][0]*/exp(a[t][0][0]/Te));

        a[t][nl[1]][0] = 2.0*ggi[2]*(a[t][nl[1]][1]-kk[2]*b[0][t]);
        a[t][nl[1]][0]+= 2.0*ggi[1]*(a[t][nl[1]-1][0]-kk[1]*b[5][t]);
        a[t][nl[1]][0]+= ggi[0]*(a[t+1][nl[1]][0]+a[t-1][nl[1]][0]);
        a[t][nl[1]][0]+= gg*(rh[t][nl[1]][0] - /*rh0[t][nl[1]][0]*/exp(a[t][nl[1]][0]/Te));

        for (s=1; k = 2*s - svts,k <= nl[1]; s++) {
          bst = a[t][k][0];

          a[t][k][0] = 2.0*ggi[2]*(a[t][k][1]-kk[2]*b[0][t+k*(nl[0]+1)]);
          a[t][k][0]+= ggi[1]*(a[t][k+1][0]+a[t][k-1][0]);
          a[t][k][0]+= ggi[0]*(a[t+1][k][0]+a[t-1][k][0]);
          a[t][k][0]+= gg*(rh[t][k][0] - /*rh0[t][k][0]*/exp(a[t][k][0]/Te));

          err+=fabs(a[t][k][0] - bst);

          // a[t][k][0] =0.0;

        }

      }
#endif
    } // end svtgs loop

    //tilfil(a,nl,"a.dta");
    //tilfil(rh,nl,"rh.dta");

    //for (int tid=1; tid < 10e10; tid++);

    it=++it;
    //cout<<err/((nl[0]-1)*(nl[1]+1)*(nl[2]-1))<<"  "<<it<<"  "<<(nl[0]-1)<<"  "<<(nl[1]+1)<<endl;


    err/= (nl[0]-1)*(nl[1]+1)*(nl[2]-1);
    //graf(err,"err.dta",int(it), nl);
    //cout<<"err: "<<err<<"    minerr: "<<minerr<<"   it: "<<it<<"  "<<(nl[0]-1)<<"  "<<(nl[1]+1)<<endl;
    //for (int tid=1; tid < 10e10; tid++);
    if ( it >= maxit) {
      cout<<" Ikke funnet stabil loesning, err: "<<err<<endl;
    }

  }
  //cout<< "break"<<endl;
  //for (int tid=1; tid < 10e20; tid++);
  cout<<"   err = "<<err<<"  it= "<<it<<endl;

}

//-------------------------------------------------------------------------
void mgfas(double ***u, double ***r, double **b, double ***rho0,const double T, double lx[],
           int n[])
//-------------------------------------------------------------------------
{
  int i, j,k, ng = 0, nn[3], ngmin[3], ib, ngrid = cn,**bnd;
  double ***irho[NGMAX+1],***irhs[NGMAX+1],***iu[NGMAX+1],***ru[NGMAX+1],
  **irnd[NGMAX+1],**ibnd[NGMAX+1],***irho0[NGMAX+1];

  // cout<<"void mgfas::::  ngrid = "<< ngrid<<endl;
//    //  boudary comndition:
//    bnd = new double*[6];

//    bnd[0] = new double[ ll[0] * ll[1] ]; //(x,y) z = 0
//    bnd[1] = new double[ ll[0] * ll[1] ]; //(x,y) z = zmax
//    bnd[2] = new double[ ll[2] * ll[1] ]; //(y,z) x = 0
//    bnd[3] = new double[ ll[2] * ll[1] ]; //(y,z) x = xmax
//    bnd[4] = new double[ ll[0] * ll[2] ]; //(x,z) y = 0
//    bnd[5] = new double[ ll[0] * ll[2] ]; //(x,y) y = ymax

  //     cout<<" kommer dette paa skjermen"<<endl;
  if (n[0] >= n[1])
    ib = 0;
  else
    ib = 1;

  for (i=0; i<3; i++) {
    nn[i] = n[i]/2;
    ngmin[i] = n[i];
  }

  for (i=0; i<ngrid; i++)
    for (j=0; j<3; j++)
      ngmin[j]/= 2;   // dimetion on coarstest grid

  irnd[ngrid]=boundary3D(nn);  //boundery condition
  smooth(irnd[ngrid],b,nn);

  irho[ngrid]=dtensor(nn);
  rstrct(irho[ngrid],r,nn);

  iu[ngrid]=dtensor(nn);
  rstrct(iu[ngrid],u,nn);

  irho0[ngrid]=dtensor(nn);
  rstrct(irho0[ngrid],rho0,nn);

  for (i=1; i<ngrid; i++) {
    for (j=0; j<3; j++)
      nn[j]/= 2;

    cout<<" kommer dette paa skjermen II i:"<<i<<endl;
    irnd[--ngrid]=boundary3D(nn);
    smooth(irnd[ngrid],irnd[ngrid+1],nn);
    irho[ngrid]=dtensor(nn);
    rstrct(irho[ngrid],irho[ngrid+1],nn);
    iu[ngrid]=dtensor(nn);
    rstrct(iu[ngrid],iu[ngrid+1],nn);
    irho0[ngrid]=dtensor(nn);
    rstrct(irho0[ngrid],irho0[ngrid+1],nn);
  }

  for (i=0; i<3; i++)
    nn[i] = ngmin[i];

  ngrid = cn+1;
  iu[ngrid]=dtensor(n);
  irho0[ngrid]=dtensor(n);

  for (j=1;j<=ngrid;j++) {
    irhs[j]=dtensor(nn);
    ibnd[j]=boundary3D(nn);
    copy_boundary3D(ibnd[j],(j != ngrid ? irnd[j] : b),nn);
    if (j > 1) {
      interp(iu[j],iu[j-1],nn);
      interp(irho0[j],irho0[j-1],nn);

    } else {
      //tilfil(iu[j],nn,"iu.dta");
    }

    //cout<<nn[0]<<"  "<<nn[1]<<"  "<<nn[2]<<"  j="<<j<<"   ngrid "<<ngrid<<endl;
    copy3D(irhs[j],(j != ngrid ? irho[j] : r),nn);

    for (int rr =0;rr<=nn[0];rr++)
      for (int ss =0;ss<=nn[1];ss++)
        for (int tt =0;tt<=nn[2];tt++)
          if (irho0[j][rr][ss][tt] <0.5)
            cout<<" Feilen er foer dette "<<irho0[j][rr][ss][tt]<<"  "<<rr<<"  "<<ss<<"  "<<tt<<endl;


    GausS_DvN(iu[j],ibnd[j],irhs[j],irho0[j],T,lx,nn);





    for (i=0; i<3; i++)
      nn[i]*= 2;
  }


  copy3D(u,iu[ngrid],n);
  copy3D(rho0,irho0[ngrid],n);

  for (i=0; i<3; i++)
    nn[i] = n[i];

  for (j=cn+1; j>=1; j--) {
    free_dtensor(irhs[j]);
    //free_dtensor(ibnd[j],0,3,0,nn[ib],0,nn[2]);
    free_boundary3D(ibnd[j]);
    free_dtensor(iu[j]);
    free_dtensor(irho0[j]);
    if (j != cn+1)
      free_dtensor(irho[j]);
    //        if (j != cn+1)
    //  	{
    //  	  //free_dtensor(irnd[j],0,3,0,nn[ib],0,nn[2]);
    //  	  free_dtensor(irho[j]);
    //  	}
    for (i=0; i<3; i++)
      nn[i]/= 2;
  }
}

//  //-------------------------------------------------------------------------
//  double **boundary3D(int ll[])
//  //-------------------------------------------------------------------------
//  {
//    int i;
//    //  boudary comndition:
//    b = new double*[6];

//    b[0] = new double[ (ll[0]+1) * (ll[1]+1) ]; //(x,y) z = 0
//    b[1] = new double[ (ll[0]+1) * (ll[1]+1) ]; //(x,y) z = zmax
//    b[2] = new double[ (ll[2]+1) * (ll[1]+1) ]; //(y,z) x = 0
//    b[3] = new double[ (ll[2]+1) * (ll[1]+1) ]; //(y,z) x = xmax
//    b[4] = new double[ (ll[0]+1) * (ll[2]+1) ]; //(x,z) y = 0
//    b[5] = new double[ (ll[0]+1) * (ll[2]+1) ]; //(x,y) y = ymax

//    for( i =0; i< (ll[0]+1)*(ll[1]+1); i++)
//      {
//        *(b[0]+i) = 0.0;
//        *(b[1]+i) = 0.0;
//      }

//    for( i =0; i< (ll[2]+1)*(ll[1]+1); i++)
//      {
//        *(b[2]+i) = 0.0;
//        *(b[3]+i) = 0.0;
//      }

//    for( i =0; i< (ll[0]+1)*(ll[2]+1); i++)
//      {
//        *(b[4]+i) = 0.0;
//        *(b[5]+i) = 0.0;
//      }

//    return b;

//  }


//  //-------------------------------------------------------------------------
//  double **copy_boundary3D(double **bnd,int ll[])
//  //-------------------------------------------------------------------------
//  {
//    int i;
//    //  boudary comndition:
//    b = new double*[6];

//    b[0] = new double[ (ll[0]+1) * (ll[1]+1) ]; //(x,y) z = 0
//    b[1] = new double[ (ll[0]+1) * (ll[1]+1) ]; //(x,y) z = zmax
//    b[2] = new double[ (ll[2]+1) * (ll[1]+1) ]; //(y,z) x = 0
//    b[3] = new double[ (ll[2]+1) * (ll[1]+1) ]; //(y,z) x = xmax
//    b[4] = new double[ (ll[0]+1) * (ll[2]+1) ]; //(x,z) y = 0
//    b[5] = new double[ (ll[0]+1) * (ll[2]+1) ]; //(x,y) y = ymax

//    for( i =0; i< (ll[0]+1)*(ll[1]+1); i++)
//      {
//        *(b[0]+i) = *(bnd[0]+i);
//        *(b[1]+i) = *(bnd[1]+i);
//      }


//    for( i =0; i< (ll[2]+1)*(ll[1]+1); i++)
//      {
//        *(b[2]+i) = *(bnd[2]+i);
//        *(b[3]+i) = *(bnd[3]+i);
//      }

//    for( i =0; i< (ll[0]+1)*(ll[2]+1); i++)
//      {
//        *(b[4]+i) = *(bnd[4]+i);
//        *(b[5]+i) = *(bnd[5]+i);
//      }

//    return b;

//  }


//  //-------------------------------------------------------------------------
//  void free_boundary3D(double **bnd)
//  //-------------------------------------------------------------------------
//  {

//    for ( int i = 0; i <6; i++)
//      delete[] bnd[i];

//    delete[] bnd;

//  }


#endif


















