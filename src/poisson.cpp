/**************************************************************************
 *
 * $Id: poisson.cpp,v 1.29 2011/03/26 12:56:40 patrick Exp $
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

#include <poisson.h>

namespace mudfas {

  using std::cout;
  using std::endl;
  using std::ostream;

#define ID "\n$Id: poisson.cpp,v 1.29 2011/03/26 12:56:40 patrick Exp $\n"

  ostream& operator<<(ostream& os, const PoissonSolver &s)
  {
    s.printOn(os);
    return os;
  }


  PoissonSolver::PoissonSolver(int nargs, char* args[])
    : LinearSolver(nargs, args)
  {
    initParsing(nargs, args);
    paramParsing();
  }

  PoissonSolver::~PoissonSolver()
  {}

  void PoissonSolver::initParsing(int nargs, char *args[])
  {
    registerClass("PoissonSolver");
    registerPackage(PACKAGE, VERSION "\n" ID "\n", MUDFAS_COPYRIGHT);
  }

  void PoissonSolver::paramParsing()
  {}

  void PoissonSolver::setPhiAndRhs(const Field &phi, const Field &rho)
  {
    Phi.resize(gridSize);
    Phi = phi;

    Rhs.resize(gridSize);
    // \nabla^2\phi = -n_i + n_e
    Rhs = -rho;
  }

  template <class T_numtype>
  void PoissonSolver::getElectricFieldEnergy(T_numtype &ESE) const
  {
    Range I(1, gridSize(0));
    Range J(1, gridSize(1));

#if (DIM==2)

    Field phi(mgPhi(numGrid-1)(I,J));
#elif (DIM==3)

    Range K(1, gridSize(2));
    Field phi(mgPhi(numGrid-1)(I,J,K));
#endif

    Field rho(phi.shape());
    rho = -mgRhs(numGrid-1);

    // ESE from Plasma physics via computer simulation
    // Birsdall and Langdon, p.74
    Field Energy(rho*phi);

#if 0

    cout
        << "rho.shape()    = " << rho.shape() << endl
        << "phi.shape()    = " << phi.shape() << endl
        << "Energy.shape() = " << Energy.shape() << endl
        << "mean(rho)      = " << mean(rho) << endl
        << "mean(phi)      = " << mean(phi) << endl
        << "mean(Energy)   = " << Energy.shape() << endl;
#endif

    FieldVecr dr((domain.max-domain.min)/(gridSize-1));
    mudfas::real dI;
    ESE = 0.5*integrate(Energy, dr, dI);
  }

  template
  void PoissonSolver::getElectricFieldEnergy<float>
  (float &ESE) const;

  template
  void PoissonSolver::getElectricFieldEnergy<double>
  (double &ESE) const;

  void PoissonSolver::setCof(int dim, const Array1dr &r, Array1dr &crr,
                             Array1dr &cr, Array1dr &cer) const
  {
    crr = 1.0;
    cr  = 0.0;
    cer = 0.0;
  }

#undef ID

}
