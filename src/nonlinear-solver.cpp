/**************************************************************************
 *
 * $Id: nonlinear-solver.cpp,v 1.25 2011/03/26 12:56:40 patrick Exp $
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
 * Class NonLinearSolver to solve the nonlinear 2D/3D problems using the
 * multigrid method
 *
 * L(phi) + F(phi) = rhs
 *
 * Short description of the method
 *
 * -> Relaxation algorithm: Gauss-Seidel-Newton method
 *
 * -> Prolongation algorithm: Linear or cubic interpolation
 *
 * -> Scheduling:
 *
 * ---> Adaptive full multigrid (FMG) schedule (nested iterations)
 *
 * ---> followed by (if needed):  nmax adaptive multi grid iterations
 *
 * -> Argument parsing of the configuration parameters such as: number of
 * pre/post smoothing, type of prolongator, V/W cycles, tolerance
 * for the solution
 *
 * Code uses the Blitz C++ library (http://oonumerics.org/blitz)
 *
 */

#include <nonlinear-solver.h>

namespace mudfas {

  using blitz::product;

  NonLinearSolver::NonLinearSolver(int nargs, char* args[]) :
    MultiGridSolver(nargs, args)
  {}

  NonLinearSolver::~NonLinearSolver()
  {}


#if (DIM==2)
  void NonLinearSolver::diffop(Field &phi, Field &lphi, int kgrid)
  // 2D non linear operator
  {
    int nx = mgGridSize(kgrid)(0);
    int ny = mgGridSize(kgrid)(1);

    Field F(Range(1,nx), Range(1,ny));
    getNonLinearFunction(phi, F, kgrid);

    int ist = (bdyc.min(0) == dirichlet ? 2 : 1);
    int ifn = (bdyc.max(0) == dirichlet ? nx-1 : nx);
    int jst = (bdyc.min(1) == dirichlet ? 2 : 1);
    int jfn = (bdyc.max(1) == dirichlet ? ny-1 : ny);

    Range I(ist,ifn), J(jst,jfn);

    // compute fine grid residual
    // L(phi) = d2(phi)/dx2 + d2(phi)/dy2 + F(phi)

    Field cx1(mgCof(0,0)(kgrid)), cx2(mgCof(0,1)(kgrid)), cx3(mgCof(0,2)(kgrid)),
          cy1(mgCof(1,0)(kgrid)), cy2(mgCof(1,1)(kgrid)), cy3(mgCof(1,2)(kgrid));

    lphi(I,J) =
      cx1(I,J)*phi(I-1,J) + cx2(I,J)*phi(I+1,J) +
      cy1(I,J)*phi(I,J-1) + cy2(I,J)*phi(I,J+1) +
      ( cx3(I,J) + cy3(I,J) )*phi(I,J) + F(I,J);

  }
#elif (DIM==3)
  void NonLinearSolver::diffop(Field &phi, Field &lphi, int kgrid)
  // 3D non linear operator
  {
    int nx = mgGridSize(kgrid)(0);
    int ny = mgGridSize(kgrid)(1);
    int nz = mgGridSize(kgrid)(2);

    Field F(Range(1,nx), Range(1,ny), Range(1,nz));
    getNonLinearFunction(phi, F, kgrid);

    // set loop limits
    int ist = (bdyc.min(0) == dirichlet ? 2 : 1);
    int ifn = (bdyc.max(0) == dirichlet ? nx-1 : nx);
    int jst = (bdyc.min(1) == dirichlet ? 2 : 1);
    int jfn = (bdyc.max(1) == dirichlet ? ny-1 : ny);
    int kst = (bdyc.min(2) == dirichlet ? 2 : 1);
    int kfn = (bdyc.max(2) == dirichlet ? nz-1 : nz);

    Range I(ist,ifn), J(jst,jfn), K(kst,kfn);

    // compute fine grid residual
    // L(phi) = d2(phi)/dx2 + d2(phi)/dy2 + d2(phi)/dz2 + F(phi)

    Field cx1(mgCof(0,0)(kgrid)), cx2(mgCof(0,1)(kgrid)), cx3(mgCof(0,2)(kgrid)),
          cy1(mgCof(1,0)(kgrid)), cy2(mgCof(1,1)(kgrid)), cy3(mgCof(1,2)(kgrid)),
          cz1(mgCof(2,0)(kgrid)), cz2(mgCof(2,1)(kgrid)), cz3(mgCof(2,2)(kgrid));

    lphi(I,J,K) =
      cx1(I,J,K)*phi(I-1,J,K)+cx2(I,J,K)*phi(I+1,J,K)+
      cy1(I,J,K)*phi(I,J-1,K)+cy2(I,J,K)*phi(I,J+1,K)+
      cz1(I,J,K)*phi(I,J,K-1)+cz2(I,J,K)*phi(I,J,K+1)+
      (cx3(I,J,K)+cy3(I,J,K)+cz3(I,J,K))*phi(I,J,K)+F(I,J,K);
  }
#endif

#if (DIM==2)
  void NonLinearSolver::gaussSeidelNewton(Field& phi, Field& rhs,
                                          Field& F, Field& dF, int kgrid,
                                          Range& I, Range& J)
  {
    // L(phi) = d2(phi)/dx2 + d2(phi)/dy2 + F(phi)
    // phi = phi - ( L(phi) - rhs ) / d(L(phi))/dphi

    Field cx1(mgCof(0,0)(kgrid)), cx2(mgCof(0,1)(kgrid)), cx3(mgCof(0,2)(kgrid)),
          cy1(mgCof(1,0)(kgrid)), cy2(mgCof(1,1)(kgrid)), cy3(mgCof(1,2)(kgrid));

    phi(I,J) -=
      (
        cx1(I,J)*phi(I-1,J) + cx2(I,J)*phi(I+1,J) +
        cy1(I,J)*phi(I,J-1) + cy2(I,J)*phi(I,J+1) +
        (cx3(I,J) + cy3(I,J))*phi(I,J) + F(I,J) - rhs(I,J)
      ) /
      (cx3(I,J) + cy3(I,J) + dF(I,J));
  }
#elif (DIM==3)
  void NonLinearSolver::gaussSeidelNewton(Field& phi, Field& rhs,
                                          Field& F, Field& dF, int kgrid,
                                          Range &I, Range &J, Range &K)
  {
    // L(phi)=d2(phi)/dx2 + d2(phi)/dy2 + d2(phi)/dz2 + F(phi)
    // phi = phi - ( L(phi)-rhs ) / d(L(phi))/dphi

    Field cx1(mgCof(0,0)(kgrid)), cx2(mgCof(0,1)(kgrid)), cx3(mgCof(0,2)(kgrid)),
          cy1(mgCof(1,0)(kgrid)), cy2(mgCof(1,1)(kgrid)), cy3(mgCof(1,2)(kgrid)),
          cz1(mgCof(2,0)(kgrid)), cz2(mgCof(2,1)(kgrid)), cz3(mgCof(2,2)(kgrid));

    phi(I,J,K) -=
      (
        cx1(I,J,K)*phi(I-1,J,K) + cx2(I,J,K)*phi(I+1,J,K) +
        cy1(I,J,K)*phi(I,J-1,K) + cy2(I,J,K)*phi(I,J+1,K) +
        cz1(I,J,K)*phi(I,J,K-1) + cz2(I,J,K)*phi(I,J,K+1) +
        (cx3(I,J,K) + cy3(I,J,K) + cz3(I,J,K))*phi(I,J,K) +
        F(I,J,K) - rhs(I,J,K)
      ) /
      (cx3(I,J,K) + cy3(I,J,K) + cz3(I,J,K) + dF(I,J,K));
  }
#endif


#if (DIM==2)
  void NonLinearSolver::smooth(Field& phi, Field& rhs, int kgrid)
  {
    // gauss-seidel-newton red/black point relaxation

    // L(phi)=d2(phi)/dx2 + d2(phi)/dy2 -exp(phi/Te)
    // phi=phi - ( L(phi)-rhs ) / d(L(phi))/dphi

    int nx = mgGridSize(kgrid)(0);
    int ny = mgGridSize(kgrid)(1);

    Range I(1,nx), J(1,ny);
    Field F(I,J), dF(I,J);
    getNonLinearFunctionAndDerivative(phi, F, dF, kgrid);

    // set periodic b.c. indicator
    int nper = product(bdyc.min);

    int ist = (bdyc.min(0) == dirichlet ? 3 : 1);
    int ifn = (bdyc.max(0) == dirichlet ? nx-1 : nx);
    int jst = (bdyc.min(1) == dirichlet ? 3 : 1);
    int jfn = (bdyc.max(1) == dirichlet ? ny-1 : ny);

    // set periodic boundaries if necessary
    if (nper == 0)
      pervb(phi,kgrid);

    // L(phi)=d2(phi)/dx2 + d2(phi)/dy2 -exp(phi/Te)
    // phi=phi - ( L(phi)-rhs ) / d(L(phi))/dphi

    // relax on red grid points
    I.setRange(ist,adjfn(ist,ifn,2),2);
    J.setRange(jst,adjfn(jst,jfn,2),2);
    gaussSeidelNewton(phi,rhs,F,dF,kgrid,I,J);

    I.setRange(2,adjfn(2,ifn,2),2);
    J.setRange(2,adjfn(2,jfn,2),2);
    gaussSeidelNewton(phi,rhs,F,dF,kgrid,I,J);

    if (nper == 0)
      pervb(phi,kgrid);

    // relax on black grid points
    I.setRange(ist,adjfn(ist,ifn,2),2);
    J.setRange(2,adjfn(2,jfn,2),2);
    gaussSeidelNewton(phi,rhs,F,dF,kgrid,I,J);

    I.setRange(2,adjfn(2,ifn,2),2);
    J.setRange(jst,adjfn(jst,jfn,2),2);
    gaussSeidelNewton(phi,rhs,F,dF,kgrid,I,J);

    // set periodic boundaries
    if (nper == 0)
      pervb(phi,kgrid);
  }

#elif (DIM==3)

  void NonLinearSolver::smooth(Field& phi, Field& rhs, int kgrid)
  {
    // gauss-seidel point relaxation with red/black ordering
    // in three dimensions for nonseparable pde
    // relax in order:
    // (1) red (x,y) on odd z planes
    // (2) black (x,y) on even z planes
    // (3) black (x,y) on odd z planes
    // (4) red (x,y) on even z planes

    int nx = mgGridSize(kgrid)(0);
    int ny = mgGridSize(kgrid)(1);
    int nz = mgGridSize(kgrid)(2);

    Range I(1,nx), J(1,ny), K(1,nz);
    Field F(I,J,K), dF(I,J,K);
    getNonLinearFunctionAndDerivative(phi, F, dF, kgrid);

    // set periodic b.c. indicator
    int nper = product(bdyc.min);

    // set loop limits to avoid specified boundaries
    // in red/black sweeps
    int ist = (bdyc.min(0) == dirichlet ? 3 : 1);
    int ifn = (bdyc.max(0) == dirichlet ? nx-1 : nx);
    int jst = (bdyc.min(1) == dirichlet ? 3 : 1);
    int jfn = (bdyc.max(1) == dirichlet ? ny-1 : ny);
    int kst = (bdyc.min(2) == dirichlet ? 3 : 1);
    int kfn = (bdyc.max(2) == dirichlet ? nz-1 : nz);

    // set periodic boundaries if necessary
    if (nper == 0)
      pervb(phi,kgrid);

    // red (x,y) on odd z planes
    I.setRange(ist,adjfn(ist,ifn,2),2);
    J.setRange(jst,adjfn(jst,jfn,2),2);
    K.setRange(kst,adjfn(kst,kfn,2),2);
    gaussSeidelNewton(phi,rhs,F,dF,kgrid,I,J,K);

    I.setRange(2,adjfn(2,ifn,2),2);
    J.setRange(2,adjfn(2,jfn,2),2);
    gaussSeidelNewton(phi,rhs,F,dF,kgrid,I,J,K);

    if (nper == 0)
      pervb(phi,kgrid);

    // black (x,y) or even z planes
    I.setRange(ist,adjfn(ist,ifn,2),2);
    K.setRange(2,adjfn(2,kfn,2),2);
    gaussSeidelNewton(phi,rhs,F,dF,kgrid,I,J,K);

    I.setRange(2,adjfn(2,ifn,2),2);
    J.setRange(jst,adjfn(jst,jfn,2),2);
    gaussSeidelNewton(phi,rhs,F,dF,kgrid,I,J,K);

    if (nper == 0)
      pervb(phi,kgrid);

    // black (x,y) on odd z planes
    I.setRange(ist,adjfn(ist,ifn,2),2);
    J.setRange(2,adjfn(2,jfn,2),2);
    K.setRange(kst,adjfn(kst,kfn,2),2);
    gaussSeidelNewton(phi,rhs,F,dF,kgrid,I,J,K);

    I.setRange(2,adjfn(2,ifn,2),2);
    J.setRange(jst,adjfn(jst,jfn,2),2);
    gaussSeidelNewton(phi,rhs,F,dF,kgrid,I,J,K);

    if (nper == 0)
      pervb(phi,kgrid);

    // red(x,y) on even z planes
    I.setRange(ist,adjfn(ist,ifn,2),2);
    J.setRange(jst,adjfn(jst,jfn,2),2);
    K.setRange(2,adjfn(2,kfn,2),2);
    gaussSeidelNewton(phi,rhs,F,dF,kgrid,I,J,K);

    I.setRange(2,adjfn(2,ifn,2),2);
    J.setRange(2,adjfn(2,jfn,2),2);
    gaussSeidelNewton(phi,rhs,F,dF,kgrid,I,J,K);

    // set periodic boundaries if necessary
    if (nper == 0)
      pervb(phi,kgrid);
  }

#endif

}
