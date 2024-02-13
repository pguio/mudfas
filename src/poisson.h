/**************************************************************************
 *
 * $Id: poisson.h,v 1.26 2011/11/07 18:39:29 patrick Exp $
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
 * Class PoissonSolver to solve the linear 2D/3D problems using the
 * multigrid method
 *
 * nabla^2(phi) = rho
 *
 */

#ifndef POISSON_SOLVER_H
#define POISSON_SOLVER_H

#include <linear-solver.h>
#include <integrate.h>
#if !defined(BLITZ2)
#include <blitz/tinyvec-et.h>
#endif 

namespace mudfas {

  class PoissonSolver : public LinearSolver {
  public:

    friend ostream& operator<<(ostream& os, const PoissonSolver &s);

    PoissonSolver(int nargs, char* args[]);
    virtual ~PoissonSolver();

    virtual void setPhiAndRhs(const Field &phi, const Field &rho);

    template <class T_numtype>
    void getElectricFieldEnergy(T_numtype &ESE) const;

  protected:

  private:

    void initParsing(int nargs, char *args[]);
    void paramParsing();

    virtual void setCof(int dim, const Array1dr &r, Array1dr &crr, Array1dr &cr,
                        Array1dr &cer) const;
  };

}

#endif // POISSON_SOLVER_H
