/**************************************************************************
 *
 * $Id: mudfas.h,v 1.103 2017/12/23 18:42:24 patrick Exp $
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
 * Abstract Base Class MultiGridSolver to discretize and solve with the
 * multigrid iteration method two and three dimensional, real, linear and
 * nonlinear elliptic partial differential equations on rectangular domains
 * with any combination of periodic, specified (Dirichlet), or mixed
 * derivative boundary conditions
 *
 * Inspired from the linear multigrid method used in the MUDPACK fortran
 * package  (http://www.scd.ucar.edu/css/software/mudpack/)
 *
 * and the book "An introduction to multigrid methods" by P. Wesseling J.
 * Wiley & Sons, 1992 (ISBN 0 471 93083 0)
 *
 * Short description of the method
 *
 * -> Relaxation algorithm:
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

#ifndef MUDFAS_H
#define MUDFAS_H

#include <iostream>
#include <string>

#include <mudfas-defs.h>
#include <boundary.h>
#include <classexception.h>

namespace mudfas {

  // Definition of the Abstract Base Class MultiGridSolver
  class MultiGridSolver : public parser::Parser {
  public:

    typedef std::string string;
    typedef std::ostream ostream;

    typedef blitz::Range Range;

    typedef blitz::Array<Field,1> MultiGridField;
    typedef blitz::Array<FieldVeci,1> MultiGridFieldVeci;

    typedef blitz::TinyMatrix<Array1dr,DIM,3> MultiDimCof;
    typedef blitz::TinyMatrix<MultiGridField,DIM,3> MultiGridMultiDimCof;

    typedef Parser::LUT LUT;
    typedef Parser::LUTPair LUTPair;

    enum bc_enum     { periodic=0, dirichlet, mixed };
    enum kcycle_enum { v=1, w, f };
    enum intpol_enum { linear=1, cubic=3 };
    enum axis_enum   { x=0, y, z };

    friend ostream& operator<<(ostream& os, const MultiGridSolver &s);

    MultiGridSolver(int nargs, char* args[]);
    virtual ~MultiGridSolver();

    virtual void initialise();

    virtual void setPhiAndRhs(const Field &phi, const Field &rhs) = 0;

    void getPhi(Field &phi) const;
    void getRhs(Field &rhs) const;

    bool isInitGuess() const;

    FieldVeci getGridSize() const {
      return gridSize;
    }
    void getBoundaryCondition(Boundary<int,DIM> &_bdyc) const {
      _bdyc = bdyc;
    }
    void setGbdr(Boundary<real,DIM> &_gbdy) {
      gbdy = _gbdy;
    }

    virtual void solve();

  protected:

    enum parser_enum {
      _gridsize=1, _bcmin, _bcmax, _rmin, _rmax, 
      _alfamin, _alfamax, _gbdymin, _gbdymax,
      _tolmax, _delta, _maxcy, _iguess, _kcycle, _ipre, _ipost, _intpol, next
    };

    FieldVeci gridSize;
    // Boundary condition setup
    Boundary<int,DIM> bdyc;
    // Parameters for boundary condition
    // For example at x=xa
    // d(phi(x, y))/dx + alfa(y) * phi(x, y) = gbdy(y)
    Boundary<real,DIM> gbdy, alfa;
    // Range of the domain
    Boundary<real,DIM> domain;

    real tolmax;    // Tolerance for convergence
    real delta;     // Adaptive parameter
    int   maxcy;     // Maximum number of cycles
    bool  iguess;    // initial guess flag
    int   kcycle;
    int   ipre;
    int   ipost;
    int   intpol;

    static LUT bcMap;
    static LUT kcycleMap;
    static LUT intpolMap;

    Field Phi, Rhs;

    // Multigrid stuff
    int numGrid;
    MultiGridFieldVeci mgGridSize;
    Array1dr mgCycle;            // Param for cycling
    Array1dr mgS;                // Param for nonlinear multigrid for all grids
    MultiGridMultiDimCof mgCof;  // Coefs of the discretised PDE for all grids
    MultiGridField mgPhi, mgRhs;

    virtual void printOn(ostream& os) const;

    void pervb(Field& phi, int kgrid);

    void getDefect(int kgrid, real &meanDefect, real &stdDefect);
    // recursive (non)linear multigrid algorithm
    virtual void mgd(MultiGridField &Phik, MultiGridField &Rhsk, int kgrid);

    // (non)linear differential operator
    virtual void diffop(Field &phi, Field &lphi, int kgrid) = 0;
    // smoothing operator
    virtual void smooth(Field& phi, Field& rhs, int kgrid) = 0;

    // Restriction operator (in prolon.cpp)
    void restrict_op(Field &resf, Field &rhsc, int kgrid);
#if (DIM==2)

    template <class Ti, class Tj>
    void FullWeighting(Field& fine, Ti &IM1, Ti &I, Ti &IP1,
                       Tj &JM1, Tj &J, Tj &JP1,
                       Field& coarse, Ti &IC, Tj &JC);
#elif (DIM==3)

    template <class Ti, class Tj, class Tk>
    void FullWeighting(Field& fine, Ti &IM1, Ti &I, Ti &IP1,
                       Tj &JM1, Tj &J, Tj &JP1, Tk &KM1, Tk &K, Tk &KP1,
                       Field& coarse, Ti &IC, Tj &JC, Tk &KC);
#endif

    // Prolongator operator (in restrict.cpp)
    void prolon1(Field &p, int kgrid, Field &q);
#if (DIM==3)

    void prolon2(Field &p, int kgrid, Field &q);
#endif

    void prolon(Field &p, int kgrid, Field &q);

  private:

    void swk();
    void trsfc(int kc);
    void adjmd(int kgrid);

    void setGrids(int kgrid, FieldVeci &coarsestGridSize, FieldVeci &gridNumber);
    void setCofk(int kgrid);
    void setSk(int kgrid);
    void setCycle(int kgrid);


    void getMultiGridParameters(int grid_size, int &coarsest_grid_size,
                                int &grid_number) const;

    virtual void setCof(int dim, const Array1dr &r, Array1dr &crr, Array1dr &cr,
                        Array1dr &cer) const = 0;

    void initParsing(int nargs, char *args[]);
    void paramParsing();
    void checkParam() const;
  };


  inline
  int adjfn(int st, int fn, int stride)
  {
    return fn - (fn-st) % stride;
  }

}

#endif // MUDFAS_H
