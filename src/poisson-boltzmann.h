/**************************************************************************
 *
 * $Id: poisson-boltzmann.h,v 1.46 2011/11/07 18:39:29 patrick Exp $
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
 * Class PoissonBoltzmannSolver to solve the nonlinear 2D/3D problems using the
 * multigrid method
 *
 * nabla^2(phi) - exp(phi/Te) = -rho
 *
 * -> Possibility to include new temperature models
 *
 */

#ifndef POISSON_BOLTZMANN_H
#define POISSON_BOLTZMANN_H

#include <nonlinear-solver.h>
#include <integrate.h>
#if !defined(BLITZ2)
#include <blitz/tinyvec-et.h>
#endif

namespace mudfas {

  class PoissonBoltzmannSolver : public NonLinearSolver {
  public:

    friend ostream& operator<<(ostream& os,
                               const PoissonBoltzmannSolver &s);

#if (DIM==2)

    enum teshape_enum { constant=0, ygauss, xgauss, xdep, xerf, xatan,
                        last_temp=xatan+1
                      };
#elif (DIM==3)

    enum teshape_enum { constant=0, xygauss, xzgauss, last_temp=xzgauss+1 };
#endif

    PoissonBoltzmannSolver(int nargs, char* args[]);
    virtual ~PoissonBoltzmannSolver();

    virtual void initialise();

    virtual void setPhiAndRhs(const Field &phi, const Field &rho);

    void getTemp(Field &Te) const;
    void getDens(Field &Te) const;

    void getElectronDensity(Field &Ne) const;

    template <class T_numtype>
    void getElectricFieldEnergy(T_numtype &ESE) const;

    template <class T_numtype>
    void getFluidEnergy(T_numtype &KT, T_numtype &KP) const;

    real getMinTemp() const {
      return minTe;
    }
    real getMaxTemp() const {
      return maxTe;
    }

  protected:

    enum { _minTe=MultiGridSolver::next, _maxTe, _refTe, _stdTe, _teShape, _alpha,
           _minNe, _maxNe, _stdNe, _neShape, _nhmfilename, next
         };


    virtual void printOn(ostream& os) const;

  private:

    real minTe, maxTe, refTe, stdTe;
    int teShape;
    real alpha;
    //real kappa;

    real minNe, maxNe, stdNe;
    int neShape;
    // non-homogeneous model primary variables
    std::string nhmFilename;
    //FieldVeci gridSize;
    Field relDens;

    static LUT teShapeMap;
    static LUT neShapeMap;

    real beta;
    real dT, Tmin, s, s2;

    MultiGridField mgTemp;
    MultiGridField mgDens;

    void constMin(Field &Te) const;
#if (DIM==2)

    void yGauss(Field &Te) const;
    void xGauss(Field &Te) const;
    void xDep(Field &Te) const;
    void xErf(Field &Te) const;
    void xAtan(Field &Te) const;
#elif (DIM==3)

    void xyGauss(Field &Te) const;
    void xzGauss(Field &Te) const;
#endif

    void (PoissonBoltzmannSolver::*Temp[last_temp])(Field &Te) const;

    virtual void getNonLinearFunction(Field &phi, Field &F, int kgrid) const;
    virtual void getNonLinearFunctionAndDerivative(Field &phi, Field &F,
        Field &dF, int kgrid) const;

    virtual void setCof(int dim, const Array1dr &r, Array1dr &crr, Array1dr &cr,
                        Array1dr &cer) const;

    void QuasiNeutralityCorrection(Field &phi, Field &rhs);
    void PhiCorrection(Field &phi, Field &rhs);

    void interpolate(const Field & srcArray, Field & dstArray);

    void initParsing(int nargs, char *args[]);
    void paramParsing();
    void checkParam() const;
  };

}

#endif // POISSON_BOLTZMANN_H
