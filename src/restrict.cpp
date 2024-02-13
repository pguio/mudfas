/**************************************************************************
 *
 * $Id: restrict.cpp,v 1.29 2011/03/26 12:56:40 patrick Exp $
 *
 * Copyright (c) 2000-2011 Patrick Guio <patrick.guio@gmail.com>
 * All Rights Reserved.
 *
 * This program is free software; you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the
 * Free Software Foundation; either version 2.  of the License, or (at your
 * option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 *
 **************************************************************************/

/*
 * Restriction operators in both 2D and 3D adapted from
 * MUDPACK routines mud2sp.f and mudcom.f
 * (http://www.scd.ucar.edu/css/software/mudpack/)
 */

#include <mudfas.h>

namespace mudfas {

#if (DIM==2)
  template <class Ti, class Tj>
  void MultiGridSolver::FullWeighting(
    Field& fine, Ti &IM1, Ti &I, Ti &IP1, Tj &JM1, Tj &J, Tj &JP1,
    Field& coarse, Ti &IC, Tj &JC)
  {
    coarse(IC,JC) =
      (fine(IM1,JM1) + fine(IP1,JM1) +
       fine(IM1,JP1) + fine(IP1,JP1) +
       2.0*(
         fine(IM1,J) + fine(IP1,J) +
         fine(I,JM1) + fine(I,JP1)
       ) +
       4.0*fine(I,J)
      )*0.0625;
  }
#elif (DIM==3)
  template <class Ti, class Tj, class Tk>
  void MultiGridSolver::FullWeighting(
    Field &fine, Ti &IM1, Ti &I, Ti &IP1, Tj &JM1, Tj &J, Tj &JP1,
    Tk &KM1, Tk &K, Tk &KP1, Field &coarse, Ti &IC, Tj &JC, Tk &KC)
  {
    coarse(IC,JC,KC) =
      ((fine(IM1,JM1,KM1) + fine(IP1,JM1,KM1) +
        fine(IM1,JP1,KM1) + fine(IP1,JP1,KM1) +
        2.0*(
          fine(IM1,J,KM1) + fine(IP1,J,KM1) +
          fine(I,JM1,KM1) + fine(I,JP1,KM1)
        ) +
        4.0*fine(I,J,KM1)
       )*0.0625 +
       2.0*(
         fine(IM1,JM1,K) + fine(IP1,JM1,K) +
         fine(IM1,JP1,K) + fine(IP1,JP1,K) +
         2.0*(fine(IM1,J,K) + fine(IP1,J,K) +
              fine(I,JM1,K) + fine(I,JP1,K)
             ) +
         4.0*fine(I,J,K)
       )*0.0625 +
       (fine(IM1,JM1,KP1) + fine(IP1,JM1,KP1) +
        fine(IM1,JP1,KP1) + fine(IP1,JP1,KP1) +
        2.0*(
          fine(IM1,J,KP1) + fine(IP1,J,KP1) +
          fine(I,JM1,KP1) + fine(I,JP1,KP1)
        ) +
        4.0*fine(I,J,KP1)
       )*0.0625
      )*0.25;
  }
#endif

#if (DIM==2)
  void MultiGridSolver::restrict_op(Field &phif, Field &phic, int kgrid)
  {
    // restrict phi grid residual in phif to coarse grid in phic
    // using full weighting for all double precision 2d codes

    int nx  = mgGridSize(kgrid+1)(0);
    int ny  = mgGridSize(kgrid+1)(1);
    int ncx = mgGridSize(kgrid)(0);
    int ncy = mgGridSize(kgrid)(1);

    // set x,y coarsening integer subscript scales
    int ix  = (ncx==nx ? 0 : 1);
    int jy  = (ncy==ny ? 0 : 1);

    Range IM1, I, IP1, IC;
    Range JM1, J, JP1, JC;

    // restrict on interior
    IC.setRange(2,ncx-1);
    I.setRange(2+ix,adjfn(2+ix,ncx-1+ix*(ncx-2),1+ix),1+ix);
    IM1=I-1;
    IP1=I+1;

    JC.setRange(2,ncy-1);
    J.setRange(2+jy,adjfn(2+jy,ncy-1+jy*(ncy-2),1+jy),1+jy);
    JM1=J-1;
    JP1=J+1;

    FullWeighting(phif,IM1,I,IP1,JM1,J,JP1,phic,IC,JC);

    // set residual on boundaries

    Vector2di im1, i, ip1, ic;
    Vector2di jm1, j, jp1, jc;

    // y=yc,yd boundaries
    jc  = 1, ncy;
    j   = 1, ny;
    jm1 = (bdyc.min(1) == periodic ? ny-1 : 2), ny-1;
    jp1 = 2, (bdyc.max(1) == periodic ? 2 : ny-1);

    // y=yc,yd and x=xa,xb corners
    ic  = 1, ncx;
    i   = 1, nx;
    im1 = (bdyc.min(0) == periodic ? nx-1 : 2), nx-1;
    ip1 = 2, (bdyc.max(0) == periodic ? 2 : nx-1);

    for (int bx=0; bx<2; ++bx) {
      for (int by=0; by<2; ++by) {
        FullWeighting(phif,im1(bx),i(bx),ip1(bx),jm1(by),j(by),jp1(by),
                      phic,ic(bx),jc(by));
      }
    }

    // set y=yc,yd interior edges
    IC.setRange(2,ncx-1);
    I.setRange(2+ix,adjfn(2+ix,ncx-1+ix*(ncx-2),1+ix),1+ix);
    IM1 = I-1;
    IP1 = I+1;

    for (int by=0; by<2; ++by) {
      FullWeighting(phif,IM1,I,IP1,jm1(by),j(by),jp1(by),phic,IC,jc(by));
    }

    // set x=xa,xb interior edges
    JC.setRange(2,ncy-1);
    J.setRange(2+jy,adjfn(2+jy,ncy-1+jy*(ncy-2),1+jy),1+jy);
    JM1 = J-1;
    JP1 = J+1;

    for (int bx=0; bx<2; ++bx) {
      FullWeighting(phif,im1(bx),i(bx),ip1(bx),JM1,J,JP1,phic,ic(bx),JC);
    }
  }
#elif (DIM==3)
  void MultiGridSolver::restrict_op(Field &phif, Field &phic, int kgrid)
  {
    // restrict phi grid residual in phif to coarse grid in phic
    // using full weighting for all double precision 2d codes

    int nx  = mgGridSize(kgrid+1)(0);
    int ny  = mgGridSize(kgrid+1)(1);
    int nz  = mgGridSize(kgrid+1)(2);
    int ncx = mgGridSize(kgrid)(0);
    int ncy = mgGridSize(kgrid)(1);
    int ncz = mgGridSize(kgrid)(2);

    // set x,y coarsening integer subscript scales
    int ix  = (ncx==nx ? 0 : 1);
    int jy  = (ncy==ny ? 0 : 1);
    int kz  = (ncz==nz ? 0 : 1);

    Range IM1, I, IP1, IC;
    Range JM1, J, JP1, JC;
    Range KM1, K, KP1, KC;

    // restrict on interior

    IC.setRange(2,ncx-1);
    I.setRange(2+ix,adjfn(2+ix,ncx-1+ix*(ncx-2),1+ix),1+ix);
    IM1 = I-1;
    IP1 = I+1;

    JC.setRange(2,ncy-1);
    J.setRange(2+jy,adjfn(2+jy,ncy-1+jy*(ncy-2),1+jy),1+jy);
    JM1 = J-1;
    JP1 = J+1;

    KC.setRange(2,ncz-1);
    K.setRange(2+kz,adjfn(2+kz,ncz-1+kz*(ncz-2),1+kz),1+kz);
    KM1 = K-1;
    KP1 = K+1;

    FullWeighting(phif,IM1,I,IP1,JM1,J,JP1,KM1,K,KP1,phic,IC,JC,KC);

    // set residual on boundaries

    Vector2di im1, i, ip1, ic;
    Vector2di jm1, j, jp1, jc;
    Vector2di km1, k, kp1, kc;

    // x=xa and x=xb
    ic  = 1, ncx;
    i   = 1, nx;
    im1 = (bdyc.min(0) == periodic ? nx-1 : 2), nx-1;
    ip1 = 2, (bdyc.max(0) == periodic ? 2 : nx-1);

    // (y,z) interior
    JC.setRange(2,ncy-1);
    J.setRange(2+jy,adjfn(2+jy,ncy-1+jy*(ncy-2),1+jy),1+jy);
    JM1 = J-1;
    JP1 = J+1;

    KC.setRange(2,ncz-1);
    K.setRange(2+kz,adjfn(2+kz,ncz-1+kz*(ncz-2),1+kz),1+kz);
    KM1 = K-1;
    KP1 = K+1;

    for (int bx=0; bx<2; ++bx) {
      FullWeighting(phif,im1(bx),i(bx),ip1(bx),JM1,J,JP1,KM1,K,KP1,
                    phic,ic(bx),JC,KC);
    }

    // x=xa,xb and y=yc,yd interior edges
    jc  = 1, ncy;
    j   = 1, ny;
    jm1 = (bdyc.min(1) == periodic ? ny-1 : 2), ny-1;
    jp1 = 2, (bdyc.max(1) == periodic ? 2 : ny-1);

    for (int bx=0; bx<2; ++bx) {
      for (int by=0; by<2; ++by) {
        FullWeighting(phif,im1(bx),i(bx),ip1(bx),jm1(by),j(by),jp1(by),KM1,K,KP1,
                      phic,ic(bx),jc(by),KC);
      }
    }

    // x=xa,xb; y=yc,yd; z=ze,zf corners
    kc  = 1, ncz;
    k   = 1, nz;
    km1 = (bdyc.min(2) == periodic ? nz-1 : 2), nz-1;
    kp1 = 2, (bdyc.max(2) == periodic ? 2 : nz-1);

    for (int bx=0; bx<2; ++bx) {
      for (int by=0; by<2; ++by) {
        for (int bz=0; bz<2; ++bz) {
          FullWeighting(phif,im1(bx),i(bx),ip1(bx),jm1(by),j(by),jp1(by),
                        km1(bz),k(bz),kp1(bz),phic,ic(bx),jc(by),kc(bz));
        }
      }
    }

    // x=xa,xb and z=ze,zf interior edges
    for (int bx=0; bx<2; ++bx) {
      for (int bz=0; bz<2; ++bz) {
        FullWeighting(phif,im1(bx),i(bx),ip1(bx),JM1,J,JP1,
                      km1(bz),k(bz),kp1(bz),phic,ic(bx),JC,kc(bz));
      }
    }

    // (x,z) interior
    IC.setRange(2,ncx-1);
    I.setRange(2+ix,adjfn(2+ix,ncx-1+ix*(ncx-2),1+ix),1+ix);
    IM1 = I-1;
    IP1 = I+1;

    for (int by=0; by<2; ++by) {
      FullWeighting(phif,IM1,I,IP1,jm1(by),j(by),jp1(by),KM1,K,KP1,
                    phic,IC,jc(by),KC);
    }

    // y=yc,yd and z=ze,zf interior edges
    for (int by=0; by<2; ++by) {
      for (int bz=0; bz<2; ++bz) {
        FullWeighting(phif,IM1,I,IP1,jm1(by),j(by),jp1(by),km1(bz),k(bz),kp1(bz),
                      phic,IC,jc(by),kc(bz));
      }
    }

    // (x,y) interior
    for (int bz=0; bz<2; ++bz) {
      FullWeighting(phif,IM1,I,IP1,JM1,J,JP1,km1(bz),k(bz),kp1(bz),
                    phic,IC,JC,kc(bz));
    }
  }

#endif

}
