/**************************************************************************
 *
 * $Id: TestSolver.cpp,v 1.11 2019/05/10 16:38:06 patrick Exp $
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
 * Test driver for the nonlinear Poisson solver class PoissonSolver to
 * solve thee following problem
 *
 * nabla^2(phi) - exp(phi/Te) = -rho
 *
 * by the multigrid method
 *
 * Code uses the Blitz C++ library (http://oonumerics.org/blitz)
 *
 *
 */

#if defined(HAVE_CONFIG_H)
#include <mudfas-config.h>
#endif

#include <new> // to get std::set_new_handler
#include <iostream>
#include <iomanip>
#include <fstream>
using namespace std;

#include <test-solver.h>

void my_new_handler()
{
  cerr << "Out of memory" << endl;
  abort();
}

int main(int nargs, char *args[])
{
  try {
    std::set_new_handler(my_new_handler);
    cout.precision(4);

    enum { _saveFields=1 };
    bool saveFields = false;

#if defined(HAVE_MPI)
    MPI_Init(&nargs, &args);
#endif

    parser::Parser parser(nargs, args);
    parser.registerProgram(args[0]);
    parser.registerPackage(PACKAGE, VERSION, MUDFAS_COPYRIGHT);

    parser.insertOption(_saveFields, "savefields", parser::types::boolean,
                        "Save phi and rhs", Any(saveFields));

    parser.parseOption(_saveFields, saveFields);

    mudfas::TestSolver solver(nargs, args);

#define PARSE(Fun)  \
if (parser.Fun()) { \
	solver.Fun();     \
	return 0;         \
}                   \
 
    PARSE(parseHelp)
    PARSE(parseVersion)
    PARSE(parseTemplate)

#undef PARSE

    solver.initialise();
    cout << solver << endl;
    solver.solve();

    if (saveFields) solver.savePhiAndRhs();

    mudfas::real ESE;
    solver.getElectricFieldEnergy(ESE);
    cout << "Electric field energy = " << ESE << endl;

#if defined(POISSON_BOLTZMANN)
    mudfas::real KT, KP;
    solver.getFluidEnergy(KT, KP);
    cout << "Kinetic energy = " << KT-KP
         << " (" << KT << ", " << KP << ")" << endl;
#endif

#if defined(HAVE_MPI)
    MPI_Finalize();
#endif

    return 0;

  } catch(int status) {
    return status;
  } catch(ClassException& c) {
    cerr << c.what() << endl;
    return !0;
  }
}
