/**************************************************************************
 *
 * $Id: init.cpp,v 1.54 2011/11/07 18:39:29 patrick Exp $
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
 * Initialisation routines for the mudfas nonlinear multigrid method to
 * solve the nonlinear 2D/3D problems
 *
 * nabla^2(phi) - exp(phi/Te) = -rho
 *
 * largely inspired from the linear multigrid method used in the MUDPACK
 * fortran package (http://www.scd.ucar.edu/css/software/mudpack/)
 *
 * Code uses the Blitz C++ library (http://oonumerics.org/blitz)
 *
 */

#include <mudfas.h>
#if !defined(BLITZ2)
#include <blitz/tinyvec-et.h>
#endif

namespace mudfas {

  using blitz::abs;
  using blitz::where;

  void MultiGridSolver::getMultiGridParameters(int grid_size,
      int &coarsestGridSize, int &gridNumber) const
  {
    // give a try with coarsestGridSize set to 3
    for (int ng=1; ng<=grid_size/2+1; ++ng) {
      const int cgs = 2;
      int n = cgs*(1<<(ng-1))+1;
      if (n > grid_size)
        break;
      if (n == grid_size) {
        coarsestGridSize = cgs;
        gridNumber = ng;
        return;
      }
    }
    // try with coarsestGridSize larger than 3
    for (int cgs=3; cgs<=grid_size/2; cgs+=2) {
      for (int ng=1; ng<=grid_size/2+1; ++ng) {
        int n = cgs*(1<<(ng-1))+1;
        if (n > grid_size)
          break;
        if (n == grid_size) {
          coarsestGridSize = cgs;
          gridNumber = ng;
          return;
        }
      }
    }
    // if here couple (coarsestGridSize,gridNumber) not found
    ostringstream os;
    os << "triplet (n,p,q) not found for n = " << grid_size ;
    throw ClassException("MultiGridSolver", os.str());
  }

  void MultiGridSolver::setGrids(int kgrid, FieldVeci &coarsestGridSize,
                                 FieldVeci &gridNumber)
  {
    for (int d=0; d<DIM; ++d) {
      int q = (kgrid+1+gridNumber(d)-numGrid > 0 ?
               kgrid+1+gridNumber(d)- numGrid : 1);
      mgGridSize(kgrid)(d) = coarsestGridSize(d)*(1<<(q-1))+1;
    }

    // mgPhi array of size n_i+2 in dimension i
    // with indices varying from 0 to n_i+1
    mgPhi(kgrid).resize(mgGridSize(kgrid)+2);
    mgPhi(kgrid).reindexSelf(FieldVeci(0));
    mgPhi(kgrid) = 0.0;

    // mgRhs array of size n_i in dimension i
    // with indices varying from 1 to n_i
    mgRhs(kgrid).resize(mgGridSize(kgrid));
    mgRhs(kgrid).reindexSelf(FieldVeci(1));
    mgRhs(kgrid) = 0.0;
  }

  void MultiGridSolver::setSk(int kgrid)
  {
    mgS(kgrid) = 1.0;
  }

  void MultiGridSolver::setCycle(int kgrid)
  {
    switch (kcycle) {
    case v:
      mgCycle(kgrid) = 1;
      break;
    case w:
      mgCycle(kgrid) = 2;
      break;
    case f:
      mgCycle(kgrid) = (kgrid == 0 ? 1 : 2);
      break;
    }
  }


  void MultiGridSolver::setCofk(int kgrid)
  {
    using blitz::tensor::i;
    using blitz::tensor::j;

    FieldVeci gridsize(mgGridSize(kgrid));
    FieldVecr dr((domain.max-domain.min)/(gridsize-1));
    FieldVecr dr2(dr*dr);
    MultiDimCof cof;

    for (int id=0; id<DIM; ++id) {
      // set subscript limits to bypass specified boundaries
      int nr = gridsize(id);
      Range size(1, nr);
      Array1dr r(size), crr(size), cr(size), cer(size);
      r = domain.min(id)+(i-1)*dr(id);
      setCof(id, r, crr, cr, cer);
      crr = where(crr > abs(cr)*dr(id)*0.5, crr , abs(cr)*dr(id)*0.5);
      for (int ic=0; ic<3; ++ic) {
        cof(id,ic).resize(Range(1,nr));
        //			cof(id,ic).reindexSelf(blitz::TinyVector<int,1>(1));
        cof(id,ic) = 0.0;
      }
      int start = (bdyc.min(id) == dirichlet ? 2 : 1);
      int end   = (bdyc.max(id) == dirichlet ? nr-1 : nr);
      Range I(start, end);
      cof(id,0)(I) = crr(I)/dr2(id)-cr(I)/(2.0*dr(id));
      cof(id,1)(I) = crr(I)/dr2(id)+cr(I)/(2.0*dr(id));
      cof(id,2)(I) = cer(I)-(cof(id,0)(I)+cof(id,1)(I));

      // adjust discretization for mixed derivative b.c.
      if (bdyc.min(id) == mixed) {
        real c = cof(id,0)(1);
        cof(id,0)(1)  = 0.0;
        cof(id,1)(1) += c;
        cof(id,2)(1) += 2.0*dr(id)*alfa.min(id)*c;
      }
      if (bdyc.max(id) == mixed) {
        real c = cof(id,1)(nr);
        cof(id,0)(nr) += c;
        cof(id,1)(nr)  = 0.0;
        cof(id,2)(nr) -= 2.0*dr(id)*alfa.max(id)*c;
      }
    }

    for (int ic=0; ic<3; ++ic) {
      for (int id=0; id<DIM; ++id) {
        mgCof(id,ic)(kgrid).resize(gridsize);
        mgCof(id,ic)(kgrid).reindexSelf(FieldVeci(1));
      }
      Range all(Range::all());
      Array1dr cx(cof(0,ic)(all));
      Array1dr cy(cof(1,ic)(all));
#if (DIM==2)

      mgCof(0,ic)(kgrid) = cx(i)+0.0*j;
      mgCof(1,ic)(kgrid) = 0.0*i+cy(j);
#elif (DIM==3)

      using blitz::tensor::k;

      Array1dr cz(cof(2,ic)(all));
      mgCof(0,ic)(kgrid) = cx(i)+0.0*j+0.0*k;
      mgCof(1,ic)(kgrid) = 0.0*i+cy(j)+0.0*k;
      mgCof(2,ic)(kgrid) = 0.0*i+0.0*j+cz(k);
#endif

    }
  }

  void MultiGridSolver::swk()
  {
    // set phi, rhs input in arrays which include
    // virtual boundaries for phi
    int nx = gridSize(0);
    int ny = gridSize(1);
    Range I(1, nx), J(1, ny);
#if (DIM==2)

    mgPhi(numGrid-1)(I,J)   = Phi;
#elif (DIM==3)

    int nz = gridSize(2);
    Range K(1, nz);
    mgPhi(numGrid-1)(I,J,K) = Phi;
#endif

    mgRhs(numGrid-1) = Rhs;
  }

  void MultiGridSolver::trsfc(int kc)
  {
    // transfer fine grid to coarse grid
    int nx  = mgGridSize(kc+1)(0);
    int ncx = mgGridSize(kc)(0);
    int ix  = (ncx == nx ? 0 : 1);
    Range IC(1,ncx);
    Range I(1,adjfn(1,ncx+ix*(ncx-1),1+ix),1+ix);

    int ny  = mgGridSize(kc+1)(1);
    int ncy = mgGridSize(kc)(1);
    int jy  = (ncy == ny ? 0 : 1);
    Range JC(1,ncy);
    Range J(1,adjfn(1,ncy+jy*(ncy-1),1+jy),1+jy);

#if (DIM==2)

    mgPhi(kc)(IC,JC) = mgPhi(kc+1)(I,J);
    mgRhs(kc)(IC,JC) = mgRhs(kc+1)(I,J);

#elif (DIM==3)

    int nz  = mgGridSize(kc+1)(2);
    int ncz = mgGridSize(kc)(2);
    int kz  = (ncz == nz ? 0 : 1);
    Range KC(1,ncz);
    Range K(1,adjfn(1,ncz+kz*(ncz-1),1+kz),1+kz);

    mgPhi(kc)(IC,JC,KC) = mgPhi(kc+1)(I,J,K);
    mgRhs(kc)(IC,JC,KC) = mgRhs(kc+1)(I,J,K);

#endif
  }


  void MultiGridSolver::adjmd(int kgrid)
  {
    using blitz::firstDim;
    using blitz::secondDim;
#if (DIM==3)
    using blitz::thirdDim;
#endif

    FieldVeci gridsize(mgGridSize(kgrid));
    FieldVecr dr((domain.max-domain.min)/(gridsize-1));
    FieldVecr dr2(dr*dr);

    FieldVeci st, fn;
    for (int id=0; id<DIM; ++id) {
      st(id) = (bdyc.min(id) == dirichlet ? 2 : 1);
      fn(id) = (bdyc.max(id) == dirichlet ? gridsize(id)-1 : gridsize(id));
    }

    // adjust right hand side at derivative boundaries
    for (int id=0; id<DIM; ++id) {
      Range J(st((id+1)%DIM), fn((id+1)%DIM));
#if (DIM==3)

      Range K(st((id+2)%DIM), fn((id+2)%DIM));
#endif

      if (bdyc.min(id) == mixed) {
        real crr  = 1.0;
        real cr   = 0.0;
        real cr2  = cr*dr(id)*0.5;
        crr = std::max(crr, std::abs(cr2));
        real c    = crr/dr2(id)-cr/(2.0*dr(id));
#if (DIM==2)

        mgRhs(kgrid)(1,J)   += 2.0*dr(id)*c*gbdy.min(id);
#elif (DIM==3)

        mgRhs(kgrid)(1,J,K) += 2.0*dr(id)*c*gbdy.min(id);
#endif

      }
      if (bdyc.max(id) == mixed) {
        int nr  = gridsize(id);
        real crr  = 1.0;
        real cr   = 0.0;
        real cr2  = cr*dr(id)*0.5;
        crr = std::max(crr, std::abs(cr2));
        real c    = crr/dr2(id)+cr/(2.0*dr(id));
#if (DIM==2)

        mgRhs(kgrid)(nr,J)   -= 2.0*dr(id)*c*gbdy.max(id);
#elif (DIM==3)

        mgRhs(kgrid)(nr,J,K) -= 2.0*dr(id)*c*gbdy.max(id);
#endif

      }
#if (DIM==2)
      mgRhs(kgrid).transposeSelf(secondDim, firstDim);
#elif (DIM==3)

      mgRhs(kgrid).transposeSelf(secondDim,thirdDim,firstDim);
#endif

    }

    // set specified b.c.
    for (int id=0; id<DIM; ++id) {
      int nr = gridsize(id);
      Range I(1,nr), J(1,gridsize((id+1)%DIM));
#if (DIM==3)

      Range K(1, gridsize((id+2)%DIM));
#endif

      if (bdyc.min(id) == dirichlet) {
#if (DIM==2)
        mgRhs(kgrid)(1,J)   = mgPhi(kgrid)(1,J);
#elif (DIM==3)

        mgRhs(kgrid)(1,J,K) = mgPhi(kgrid)(1,J,K);
#endif

      }
      if (bdyc.max(id) == dirichlet) {
#if (DIM==2)
        mgRhs(kgrid)(nr,J)   = mgPhi(kgrid)(nr,J);
#elif (DIM==3)

        mgRhs(kgrid)(nr,J,K) = mgPhi(kgrid)(nr,J,K);
#endif

      }
#if (DIM==2)
      mgRhs(kgrid).transposeSelf(secondDim, firstDim);
#elif (DIM==3)

      mgRhs(kgrid).transposeSelf(secondDim,thirdDim,firstDim);
#endif

    }
  }


#if (DIM==2)
  void MultiGridSolver::pervb(Field& phi, int kgrid)
  {
    // set virtual periodic boundaries from interior values
    // in two dimensions (for all 2-d solvers)
    int nx = mgGridSize(kgrid)(0);
    int ny = mgGridSize(kgrid)(1);
    Range I(1,nx), J(1,ny);

    if (bdyc.min(0) == periodic) {
      phi(0,J)    = phi(nx-1,J);
      phi(nx,J)   = phi(1,J);
      phi(nx+1,J) = phi(2,J);
    }
    if (bdyc.min(1) == periodic) {
      phi(I,0)    = phi(I,ny-1);
      phi(I,ny)   = phi(I,1);
      phi(I,ny+1) = phi(I,2);
    }
  }
#elif (DIM==3)
  void MultiGridSolver::pervb(Field& phi, int kgrid)
  {
    // set virtual periodic boundaries from interior values
    // in three dimensions (for all 3-d solvers)
    int nx = mgGridSize(kgrid)(0);
    int ny = mgGridSize(kgrid)(1);
    int nz = mgGridSize(kgrid)(2);
    Range I(1,nx), J(1,ny), K(1,nz);

    if (bdyc.min(0) == periodic) {
      phi(0,J,K)    = phi(nx-1,J,K);
      phi(nx,J,K)   = phi(1,J,K);
      phi(nx+1,J,K) = phi(2,J,K);
    }
    if (bdyc.min(1) == periodic) {
      phi(I,0,K)    = phi(I,ny-1,K);
      phi(I,ny,K)   = phi(I,1,K);
      phi(I,ny+1,K) = phi(I,2,K);
    }
    if (bdyc.min(2) == periodic) {
      phi(I,J,0)    = phi(I,J,nz-1);
      phi(I,J,nz)   = phi(I,J,1);
      phi(I,J,nz+1) = phi(I,J,2);
    }
  }
#endif

}
