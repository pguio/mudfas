/**************************************************************************
 *
 * $Id: prolon.cpp,v 1.21 2011/03/26 12:56:40 patrick Exp $
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
 *
 * Prolongator operators in both 2D and 3D adapted from
 * MUDPACK routines mud2sp.f and mudcom.f
 * (http://www.scd.ucar.edu/css/software/mudpack/)
 *
 * Prolongator are of 2 types linear and cubic interpolations
 *
 * Code uses the Blitz C++ library (http://oonumerics.org/blitz)
 *
 */

#include <mudfas.h>

namespace mudfas {

  void MultiGridSolver::prolon1(Field &p, int kgrid, Field &q)
  {
    // 1D prolongator

    int ncx    = mgGridSize(kgrid-1)(0);

    int nx     = mgGridSize(kgrid)(0);

    int ist    = (bdyc.min(0) == dirichlet ? 2 : 1);
    int ioddst = (bdyc.min(0) == dirichlet ? 3 : 1);
    int ifn    = (bdyc.max(0) == dirichlet ? nx-1 : nx);
    int ioddfn = (bdyc.max(0) == dirichlet ?  nx-2 : nx);

    Range I;
    Range IC;
    Range all(Range::all());

    if (intpol == linear || ncx < 4) { // linearly interpolate in x

      if (ncx<nx) { // every other point of nx grid is ncx grid
        I.setRange(ioddst,adjfn(ioddst,ioddfn,2),2);
        IC.setRange((ioddst>>1)+1, (ioddfn>>1)+1);

#if (DIM==2)

        q(I,all)=p(IC,all);
#elif (DIM==3)

        q(I,all,all)=p(IC,all,all);
#endif

        I.setRange(2,adjfn(2,ifn,2),2);
#if (DIM==2)

        q(I,all)=0.5*(q(I-1,all)+q(I+1,all));
#elif (DIM==3)

        q(I,all,all)=0.5*(q(I-1,all,all)+q(I+1,all,all));
#endif

      } else {
        // nx grid equals ncx grid
        I.setRange(ist,ifn);
#if (DIM==2)

        q(I,all)=p(I,all);
#elif (DIM==3)

        q(I,all,all)=p(I,all,all);
#endif

      }

      // set virtual end points if periodic
      if (bdyc.min(0) == periodic) {
#if (DIM==2)
        q(0,all)=q(nx-1,all);
        q(nx+1,all)=q(2,all);
#elif (DIM==3)

        q(0,all,all)=q(nx-1,all,all);
        q(nx+1,all,all)=q(2,all,all);
#endif

      }
      return ;
    } else { // cubically interpolate in x

      if (ncx < nx) {
        I.setRange(ioddst,adjfn(ioddst,ioddfn,2),2);
        IC.setRange((ioddst>>1)+1,(ioddfn>>1)+1);
#if (DIM==2)

        q(I,all)=p(IC,all);
#elif (DIM==3)

        q(I,all,all)=p(IC,all,all);
#endif

        // set deep interior with symmetric formula

        I.setRange(4,adjfn(4,nx-3,2),2);
#if (DIM==2)

        q(I,all)=(-q(I-3,all)+9.0*(q(I-1,all)+q(I+1,all))-q(I+3,all))*0.0625;
#elif (DIM==3)

        q(I,all,all)=(-q(I-3,all,all)+9.0*(q(I-1,all,all)+q(I+1,all,all))-
                      q(I+3,all,all))*0.0625;
#endif

        // interpolate from q at i=2 and i=nx-1
        if (bdyc.min(0) != periodic) {
          // asymmetric formula near nonperiodic bndys
#if (DIM==2)
          q(2,all)=(5.0*q(1,all)+15.0*q(3,all)-5.0*q(5,all)+q(7,all))*0.0625;
          q(nx-1,all)=(5.*q(nx,all)+15.0*q(nx-2,all)-
                       5.0*q(nx-4,all)+q(nx-6,all))*0.0625;
#elif (DIM==3)

          q(2,all,all)=(5.0*q(1,all,all)+15.0*q(3,all,all)-
                        5.0*q(5,all,all)+q(7,all,all))*0.0625;
          q(nx-1,all,all)=(5.0*q(nx,all,all)+15.0*q(nx-2,all,all)-
                           5.0*q(nx-4,all,all)+q(nx-6,all,all))*0.0625;
#endif

        } else {
          // periodicity in x alows symmetric formula near bndys
#if (DIM==2)
          q(2,all)=(-q(nx-2,all)+9.0*(q(1,all)+q(3,all))-q(5,all))*0.0625;
          q(nx-1,all)=(-q(nx-4,all)+9.0*(q(nx-2,all)+q(nx,all))-q(3,all))*0.0625;
          q(nx+1,all)=q(2,all);
          q(0,all)=q(nx-1,all);
#elif (DIM==3)

          q(2,all,all)=(-q(nx-2,all,all)+9.0*(q(1,all,all)+q(3,all,all))-
                        q(5,all,all))*0.0625;
          q(nx-1,all,all)=(-q(nx-4,all,all)+9.0*(q(nx-2,all,all)+
                                                 q(nx,all,all))-q(3,all,all))*0.0625;
          q(nx+1,all,all)=q(2,all,all);
          q(0,all,all)=q(nx-1,all,all);
#endif

        }
        return ;
      } else {
        // ncx grid equals nx grid
        I.setRange(ist,ifn);
#if (DIM==2)

        q(I,all)=p(I,all);
#elif (DIM==3)

        q(I,all,all)=p(I,all,all);
#endif

        if(bdyc.min(0) == periodic) {
#if (DIM==2)
          q(0,all)=q(nx-1,all);
          q(nx+1,all)=q(2,all);
#elif (DIM==3)

          q(0,all,all)=q(nx-1,all,all);
          q(nx+1,all,all)=q(2,all,all);
#endif

        }
        return ;
      }
    }
  }

#if (DIM==2)
  void MultiGridSolver::prolon(Field &p, int kgrid, Field &q)
#elif (DIM==3)
  void MultiGridSolver::prolon2(Field &p, int kgrid, Field &q)
#endif
  {
    // 2D prolongator

    int ncy    = mgGridSize(kgrid-1)(1);

    int nx     = mgGridSize(kgrid)(0);
    int ny     = mgGridSize(kgrid)(1);

    int ist    = (bdyc.min(0) == dirichlet ? 2 : 1);
    int ifn    = (bdyc.max(0) == dirichlet ? nx-1 : nx);
    int jst    = (bdyc.min(1) == dirichlet ? 2 : 1);
    int jfn    = (bdyc.max(1) == dirichlet ? ny-1 : ny);
    int joddst = (bdyc.min(1) == dirichlet ? 3 : 1);
    int joddfn = (bdyc.max(1) == dirichlet ? ny-2 : ny);

    Range I, J;
    Range JC;
    Range all = Range::all();
    Field P, Q;
    if (intpol == linear || ncy < 4) {
      // linearly interpolate in y

      if (ncy < ny) {
        // ncy grid is an every other point subset of ny grid
        // set odd j lines interpolating in x and then set even
        // j lines by averaging odd j lines
        J.setRange(joddst,adjfn(joddst,joddfn,2),2);
        JC.setRange((joddst>>1)+1, (joddfn>>1)+1);
#if (DIM==2)

        P.resize(p(all,JC).shape());
        P=p(all,JC);
        Q.resize(q.rows(),P.cols());
        Q=0.0;
#elif (DIM==3)

        P.resize(p(all,JC,all).shape());
        P=p(all,JC,all);
        Q.resize(q.rows(),P.cols(),P.depth());
        Q=0.0;
#endif

        prolon1(P, kgrid, Q);
#if (DIM==2)

        q(all,J)=Q(all,JC-JC.first());
#elif (DIM==3)

        q(all,J,all)=Q(all,JC-JC.first(),all);
#endif

        J.setRange(2,adjfn(2,jfn,2),2);
        I.setRange(ist,ifn);
#if (DIM==2)

        q(I,J)=0.5*(q(I,J-1)+q(I,J+1));
#elif (DIM==3)

        q(I,J,all)=0.5*(q(I,J-1,all)+q(I,J+1,all));
#endif

        // set periodic virtual boundaries if necessary
        if (bdyc.min(1) == periodic) {
          I.setRange(ist,ifn);
#if (DIM==2)

          q(I,0)=q(I,ny-1);
          q(I,ny+1)=q(I,2);
#elif (DIM==3)

          q(I,0,all)=q(I,ny-1,all);
          q(I,ny+1,all)=q(I,2,all);
#endif

        }
        return ;
      } else {
        // ncy grid equals ny grid so interpolate in x only
        J.setRange(jst,jfn);
        JC=J;
#if (DIM==2)

        P.resize(p(all,JC).shape());
        P=p(all,JC);
        Q.resize(q.rows(),P.cols());
        Q=0.0;
#elif (DIM==3)

        P.resize(p(all,JC,all).shape());
        P=p(all,JC,all);
        Q.resize(q.rows(),P.cols(),P.depth());
        Q=0.0;
#endif

        prolon1(P, kgrid, Q);
#if (DIM==2)

        q(all,J)=Q(all,JC-JC.first());
#elif (DIM==3)

        q(all,J,all)=Q(all,JC-JC.first(),all);
#endif

        // set periodic virtual boundaries if necessary
        if (bdyc.min(1) == periodic) {
          I.setRange(ist,ifn);
#if (DIM==2)

          q(I,0)=q(I,ny-1);
          q(I,ny+1)=q(I,2);
#elif (DIM==3)

          q(I,0,all)=q(I,ny-1,all);
          q(I,ny+1,all)=q(I,2,all);
#endif

        }
        return ;
      }
    } else {
      // cubically interpolate in y

      if (ncy < ny) {
        // set every other point of ny grid by interpolating in x
        J.setRange(joddst,adjfn(joddst,joddfn,2),2);
        JC.setRange((joddst>>1)+1,(joddfn>>1)+1);
#if (DIM==2)

        P.resize(p(all,JC).shape());
        P=p(all,JC);
        Q.resize(q.rows(),P.cols());
        Q=0.0;
#elif (DIM==3)

        P.resize(p(all,JC,all).shape());
        P=p(all,JC,all);
        Q.resize(q.rows(),P.cols(),P.depth());
        Q=0.0;
#endif

        prolon1(P, kgrid, Q);
#if (DIM==2)

        q(all,J)=Q(all,JC-JC.first());
#elif (DIM==3)

        q(all,J,all)=Q(all,JC-JC.first(),all);
#endif

        // set deep interior of ny grid using values just
        // generated and symmetric cubic interpolation in y

        J.setRange(4,adjfn(4,ny-3,2),2);
        I.setRange(ist,ifn);
#if (DIM==2)

        q(I,J)=(-q(I,J-3)+9.0*(q(I,J-1)+q(I,J+1))-q(I,J+3))*0.0625;
#elif (DIM==3)
        q(I,J,all)=(-q(I,J-3,all)+9.0*(q(I,J-1,all)+q(I,J+1,all))-
                    q(I,J+3,all))*0.0625;
#endif

        // interpolate from q at j=2 and j=ny-1
        if (bdyc.min(1) != periodic) {
          // asymmetric formula near nonperiodic y boundaries
          I.setRange(ist,ifn);
#if (DIM==2)

          q(I,2)=(5.*q(I,1)+15.0*q(I,3)-5.*q(I,5)+q(I,7))*0.0625;
          q(I,ny-1)=(5.0*q(I,ny)+15.0*q(I,ny-2)-5.*q(I,ny-4)+q(I,ny-6))*0.0625;
#elif (DIM==3)

          q(I,2,all)=(5.0*q(I,1,all)+15.0*q(I,3,all)-
                      5.0*q(I,5,all)+q(I,7,all))*0.0625;
          q(I,ny-1,all)=(5.0*q(I,ny,all)+15.0*q(I,ny-2,all)-
                         5.0*q(I,ny-4,all)+q(I,ny-6,all))*0.0625;
#endif

        } else {
          // periodicity in y alows symmetric formula near bndys
          I.setRange(ist,ifn);
#if (DIM==2)

          q(I,2)=(-q(I,ny-2)+9.0*(q(I,1)+q(I,3))-q(I,5))*0.0625;
          q(I,ny-1)=(-q(I,ny-4)+9.0*(q(I,ny-2)+q(I,ny))-q(I,3))*0.0625;
          q(I,ny+1)=q(I,2);
          q(I,0)=q(I,ny-1);
#elif (DIM==3)

          q(I,2,all)=(-q(I,ny-2,all)+9.0*(q(I,1,all)+
                                          q(I,3,all))-q(I,5,all))*0.0625;
          q(I,ny-1,all)=(-q(I,ny-4,all)+9.0*(q(I,ny-2,all)+
                                             q(I,ny,all))-q(I,3,all))*0.0625;
          q(I,ny+1,all)=q(I,2,all);
          q(I,0,all)=q(I,ny-1,all);
#endif

        }
        return ;
      } else {
        // ncy grid is equals ny grid so interpolate in x only
        J.setRange(jst,jfn);
        JC=J;
#if (DIM==2)

        P.resize(p(all,JC).shape());
        P=p(all,JC);
        Q.resize(q.rows(),P.cols());
        Q=0.0;
#elif (DIM==3)

        P.resize(p(all,JC,all).shape());
        P=p(all,JC,all);
        Q.resize(q.rows(),P.cols(),P.depth());
        Q=0.0;
#endif

        prolon1(P, kgrid, Q);
#if (DIM==2)

        q(all,J)=Q(all,JC-JC.first());
#elif (DIM==3)

        q(all,J,all)=Q(all,JC-JC.first(),all);
#endif

        // set periodic virtual boundaries if necessary
        if (bdyc.min(1) == periodic) {
          I.setRange(ist,ifn);
#if (DIM==2)

          q(I,0)=q(I,ny-1);
          q(I,ny+1)=q(I,2);
#elif (DIM==3)

          q(I,0,all)=q(I,ny-1,all);
          q(I,ny+1,all)=q(I,2,all);
#endif

        }
        return ;
      }
    }
  }

#if (DIM==3)
  void MultiGridSolver::prolon(Field &p, int kgrid, Field &q)
  {
    // 3D prolongator

    int ncz    = mgGridSize(kgrid-1)(2);

    int nx     = mgGridSize(kgrid)(0);
    int ny     = mgGridSize(kgrid)(1);
    int nz     = mgGridSize(kgrid)(2);

    int ist    = (bdyc.min(0) == dirichlet ? 2 : 1);
    int ifn    = (bdyc.max(0) == dirichlet ? nx-1 : nx);
    int jst    = (bdyc.min(1) == dirichlet ? 2 : 1);
    int jfn    = (bdyc.max(1) == dirichlet ? ny-1 : ny);
    int kst    = (bdyc.min(2) == dirichlet ? 2 : 1);
    int kfn    = (bdyc.max(2) == dirichlet ? nz-1 : nz);
    int koddst = (bdyc.min(2) == dirichlet ? 3 : 1);
    int koddfn = (bdyc.max(2) == dirichlet ? nz-2 : nz);

    Range I, J, K;
    Range KC;
    Range all = Range::all();
    Field P, Q;
    if (intpol == linear || ncz < 4) { // linearly interpolate in z

      if (ncz < nz) {
        // ncz grid is an every other point subset of nz grid
        // set odd k planes interpolating in x&y and then set even
        // k planes by averaging odd k planes
        K.setRange(koddst,adjfn(koddst,koddfn,2),2);
        KC.setRange((koddst>>1)+1,(koddfn>>1)+1);
        P.resize(p(all,all,KC).shape());
        P=p(all,all,KC);
        Q.resize(q.rows(),q.cols(),P.depth());
        Q=0.0;
        prolon2(P, kgrid, Q);
        q(all,all,K)=Q(all,all,KC-KC.first());

        I.setRange(ist,ifn);
        J.setRange(jst,jfn);
        K.setRange(2,adjfn(2,kfn,2),2);
        q(I,J,K)=.5*(q(I,J,K-1)+q(I,J,K+1));

        // set periodic virtual boundaries if necessary
        if (bdyc.min(2) == periodic) {
          I.setRange(ist,ifn);
          J.setRange(jst,jfn);
          q(I,J,0)=q(I,J,nz-1);
          q(I,J,nz+1)=q(I,J,2);
        }
      } else {
        // ncz grid is equals nz grid so interpolate in x&y only
        K.setRange(kst,kfn);
        KC=K;
        P.resize(p(all,all,KC).shape());
        P=p(all,all,KC);
        Q.resize(q.rows(),q.cols(),P.depth());
        Q=0.0;
        prolon2(P, kgrid, Q);
        q(all,all,K)=Q(all,all,KC-KC.first());

        // set periodic virtual boundaries if necessary
        if (bdyc.min(2) == periodic) {
          J.setRange(jst,jfn);
          I.setRange(ist,ifn);
          q(I,J,0)=q(I,J,nz-1);
          q(I,J,nz+1)=q(I,J,2);
        }
      }
    } else { // cubically interpolate in z

      if (ncz < nz) {
        // set every other point of nz grid by interpolating in x&y
        K.setRange(koddst,adjfn(koddst,koddfn,2),2);
        KC.setRange((koddst>>1)+1,(koddfn>>1)+1);
        P.resize(p(all,all,KC).shape());
        P=p(all,all,KC);
        Q.resize(q.rows(),q.cols(),P.depth());
        Q=0.0;
        prolon2(P, kgrid, Q);
        q(all,all,K)=Q(all,all,KC-KC.first());

        // set deep interior of nz grid using values just
        // generated and symmetric cubic interpolation in z

        K.setRange(4,adjfn(4,nz-3,2),2);
        J.setRange(jst,jfn);
        I.setRange(ist,ifn);
        q(I,J,K)=(-q(I,J,K-3)+9.0*(q(I,J,K-1)+q(I,J,K+1))-q(I,J,K+3))*0.0625;

        // interpolate from q at k=2 and k=nz-1
        if (bdyc.min(2) != periodic) {
          // asymmetric formula near nonperiodic z boundaries
          J.setRange(jst,jfn);
          I.setRange(ist,ifn);
          q(I,J,2)=(5.0*q(I,J,1)+15.0*q(I,J,3)-5.0*q(I,J,5)+q(I,J,7))*0.0625;
          q(I,J,nz-1)=(5.0*q(I,J,nz)+15.0*q(I,J,nz-2)-
                       5.0*q(I,J,nz-4)+q(I,J,nz-6))*0.0625;
        } else {
          // periodicity in y alows symmetric formula near bndys
          J.setRange(jst,jfn);
          I.setRange(ist,ifn);
          q(I,J,2)=(-q(I,J,nz-2)+9.0*(q(I,J,1)+q(I,J,3))-q(I,J,5))*0.0625;
          q(I,J,nz-1)=(-q(I,J,nz-4)+9.0*(q(I,J,nz-2)+q(I,J,nz))-q(I,J,3))*0.0625;
          q(I,J,nz+1)=q(I,J,2);
          q(I,J,0)=q(I,J,nz-1);
        }
      } else {
        // ncz grid is equals nx grid so interpolate in x&y only
        K.setRange(kst,kfn);
        KC=K;
        P.resize(p(all,all,KC).shape());
        P=p(all,all,KC);
        Q.resize(q.rows(),q.cols(),P.depth());
        Q=0.0;
        prolon2(P, kgrid, Q);
        q(all,all,K)=Q(all,all,KC-KC.first());

        // set periodic virtual boundaries if necessary
        if (bdyc.min(2) == periodic) {
          J.setRange(jst,jfn);
          I.setRange(ist,ifn);
          q(I,J,0)=q(I,J,nz-1);
          q(I,J,nz+1)=q(I,J,2);
        }
      }
    }
  }

#endif

}
