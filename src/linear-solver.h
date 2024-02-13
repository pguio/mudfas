/**************************************************************************
 *
 * $Id: linear-solver.h,v 1.15 2011/03/26 12:56:40 patrick Exp $
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
 * Class LinearSolver to solve the nonlinear 2D/3D problems using the
 * multigrid method
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
 * pre/post smoothin, type of prolongator, V/W cycles, tolerance
 * for the solution
 *
 * Code uses the Blitz C++ library (http://oonumerics.org/blitz)
 *
 */

#ifndef LINEAR_SOLVER_H
#define LINEAR_SOLVER_H

#include <mudfas.h>

namespace mudfas {

  class LinearSolver : public MultiGridSolver {
  public:

    LinearSolver(int nargs, char* args[]);
    virtual ~LinearSolver();

  protected:

#if 0

    virtual void mgd(MultiGridField &Phik, MultiGridField &Rhsk, int kgrid);
#endif

    virtual void diffop(Field &phi, Field &lphi, int kgrid);
    virtual void smooth(Field& phi, Field& rhs, int kgrid);

  private:

#if (DIM==2)

    void gaussSeidel(Field &phi, Field &rhs, int kgrid, Range &I, Range &J);
#elif (DIM==3)

    void gaussSeidel(Field &phi, Field &rhs, int kgrid, Range &I, Range &J,
                     Range &K);
#endif
  };


}

#endif // LINEAR_SOLVER_H
