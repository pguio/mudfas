/**************************************************************************
 *
 * $Id: poisson-boltzmann.cpp,v 1.76 2017/11/26 18:18:54 patrick Exp $
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

#include <poisson-boltzmann.h>

namespace mudfas {

  using std::cout;
  using std::cerr;
  using std::endl;
  using std::ios;

#define ID "$Id: poisson-boltzmann.cpp,v 1.76 2017/11/26 18:18:54 patrick Exp $"

#if (DIM==2)

  const PoissonBoltzmannSolver::LUTPair teShapePair[] = {
    PoissonBoltzmannSolver::LUTPair(PoissonBoltzmannSolver::constant, "constant"),
    PoissonBoltzmannSolver::LUTPair(PoissonBoltzmannSolver::ygauss  , "y-Gaussian"),
    PoissonBoltzmannSolver::LUTPair(PoissonBoltzmannSolver::xgauss  , "x-Gaussian"),
    PoissonBoltzmannSolver::LUTPair(PoissonBoltzmannSolver::xdep    , "x-depletion"),
    PoissonBoltzmannSolver::LUTPair(PoissonBoltzmannSolver::xerf    , "x-erf"),
    PoissonBoltzmannSolver::LUTPair(PoissonBoltzmannSolver::xatan   , "x-atan")
  };

  const PoissonBoltzmannSolver::LUTPair neShapePair[] = {
    PoissonBoltzmannSolver::LUTPair(PoissonBoltzmannSolver::constant, "constant"),
    PoissonBoltzmannSolver::LUTPair(PoissonBoltzmannSolver::xgauss  , "x-Gaussian"),
  };

#elif (DIM==3)

  const PoissonBoltzmannSolver::LUTPair teShapePair[] = {
    PoissonBoltzmannSolver::LUTPair(PoissonBoltzmannSolver::constant, "constant"),
    PoissonBoltzmannSolver::LUTPair(PoissonBoltzmannSolver::xygauss , "xy-Gaussian"),
    PoissonBoltzmannSolver::LUTPair(PoissonBoltzmannSolver::xzgauss , "xz-Gaussian")
  };

  const PoissonBoltzmannSolver::LUTPair neShapePair[] = {
    PoissonBoltzmannSolver::LUTPair(PoissonBoltzmannSolver::constant, "constant"),
    PoissonBoltzmannSolver::LUTPair(PoissonBoltzmannSolver::xygauss , "xy-Gaussian"),
  };

#endif

  const PoissonBoltzmannSolver::LUTPair *teShapePairEnd(teShapePair +
      sizeof teShapePair / sizeof teShapePair[0]);

  PoissonBoltzmannSolver::LUT PoissonBoltzmannSolver::teShapeMap(teShapePair, teShapePairEnd);

  const PoissonBoltzmannSolver::LUTPair *neShapePairEnd(neShapePair +
      sizeof neShapePair / sizeof neShapePair[0]);

  PoissonBoltzmannSolver::LUT PoissonBoltzmannSolver::neShapeMap(neShapePair, neShapePairEnd);

  const double M_1_SQRT2PI = 0.5*M_2_SQRTPI*M_SQRT1_2;

  std::ostream& operator<<(std::ostream& os, const PoissonBoltzmannSolver &s)
  {
    s.printOn(os);
    return os;
  }

  PoissonBoltzmannSolver::PoissonBoltzmannSolver(int nargs, char* args[])
    : NonLinearSolver(nargs, args),
      minTe(2.0), maxTe(12.0), refTe(1.0),
      stdTe(0.25), teShape(constant), alpha(0.0),
      minNe(1.0), maxNe(1.0), stdNe(0.25), neShape(constant)
  {
    Temp[constant] = &PoissonBoltzmannSolver::constMin;

#if (DIM==2)

    Temp[ygauss]   = &PoissonBoltzmannSolver::yGauss;
    Temp[xgauss]   = &PoissonBoltzmannSolver::xGauss;
    Temp[xdep]     = &PoissonBoltzmannSolver::xDep;
    Temp[xerf]     = &PoissonBoltzmannSolver::xErf;
    Temp[xatan]    = &PoissonBoltzmannSolver::xAtan;
#elif (DIM==3)

    Temp[xygauss]  = &PoissonBoltzmannSolver::xyGauss;
    Temp[xzgauss]  = &PoissonBoltzmannSolver::xzGauss;
#endif

    initParsing(nargs, args);
    paramParsing();
    checkParam();
  }

  PoissonBoltzmannSolver::~PoissonBoltzmannSolver()
  {}

  void PoissonBoltzmannSolver::getTemp(Field &Te) const
  {
    (this->*Temp[teShape])(Te);
  }

  template <typename T_numtype>
  T_numtype leftInt(T_numtype x, T_numtype x0,
                    T_numtype m, T_numtype s, T_numtype dx)
  {
    using blitz::pow2;
    // \int 1/dx*( (x-x0)/dx+1)/sqrt(2*pi)/s.*exp(-(x-m).^2/2/s^2)
    return 1.0/pow2(dx)*(-s*M_1_SQRT2PI*std::exp(-0.5*pow2((x-m)/s))+
                         0.5*(m-x0+dx)*erf((x-m)*M_SQRT1_2/s));
  }

  template <typename T_numtype>
  T_numtype rightInt(T_numtype x, T_numtype x0,
                     T_numtype m, T_numtype s, T_numtype dx)
  {
    using blitz::pow2;
    // \int 1/dx*(-(x-x0)/dx+1)/sqrt(2*pi)/s.*exp(-(x-m).^2/2/s^2)
    return 1.0/pow2(dx)*( s*M_1_SQRT2PI*std::exp(-0.5*pow2((x-m)/s))-
                          0.5*(m-x0-dx)*erf((x-m)*M_SQRT1_2/s));
  }

  template <typename T_numtype>
  void oneDimGaussProbInterp(T_numtype min, T_numtype max, int n,
                             T_numtype m, T_numtype s, Array1dr & Prob)
  {
    double dx = (max-min)/(n-1);
    double xl = min;
    double xh = xl+dx;
    double mx = m;
    double sx = s;

    Prob(0) = rightInt(xh, xl, mx, sx, dx)-rightInt(xl, xl, mx, sx, dx);

    for (int i=1; i<n-1; ++i) {
      Prob(i) = leftInt(xh, xh, mx, sx, dx)-leftInt(xl, xh, mx, sx, dx);
      xl = xh;
      xh = xl+dx;
      Prob(i) += rightInt(xh, xl, mx, sx, dx)-rightInt(xl, xl, mx, sx, dx);
    }

    Prob(n-1) = leftInt(xh, xh, mx, sx, dx)-leftInt(xl, xh, mx, sx, dx);
  }

#if (DIM==2)

  void PoissonBoltzmannSolver::getDens(Field &Ne) const
  {
    int nx = Ne.rows();
    Array1dr xProb(nx);
    real mx = (domain.min(0)+domain.max(0))/2.0;
    real sx =  stdNe;

    oneDimGaussProbInterp(domain.min(0), domain.max(0), nx, mx, sx, xProb);

    int ny = Ne.cols();
    Array1dr yProb(ny);
    real my = (domain.min(1)+domain.max(1))/2.0;
    real sy =  stdNe;

#if 0

    yProb = 1.0;
#else

    oneDimGaussProbInterp(domain.min(1), domain.max(1), ny, my, sy, yProb);
#endif

    Ne = xProb(blitz::tensor::i)*yProb(blitz::tensor::j);

    Ne = ((maxNe-minNe)/blitz::max(Ne))*Ne + minNe;
  }

#elif (DIM==3)

  void PoissonBoltzmannSolver::getDens(Field &Ne) const
  {
    int nx = Ne.rows();
    Array1dr xProb(nx);
    real mx = (domain.min(0)+domain.max(0))/2.0;
    real sx =  stdNe;
    oneDimGaussProbInterp(domain.min(0), domain.max(0), nx, mx, sx, xProb);

    int ny = Ne.cols();
    Array1dr yProb(ny);
    yProb = 1.0;

    int nz = Ne.depth();
    Array1dr zProb(nz);
    real mz = (domain.min(2)+domain.max(2))/2.0;
    real sz =  stdNe;
    oneDimGaussProbInterp(domain.min(2), domain.max(2), nz, mz, sz, zProb);

    Ne = xProb(blitz::tensor::i)*yProb(blitz::tensor::j)*zProb(blitz::tensor::k);

    Ne = ((maxNe-minNe)/blitz::max(Ne))*Ne + minNe;

#if 0

    Range all(Range::all());

    int nx = Ne.rows();
    int nz = Ne.depth();

    real dNe = maxNe-minNe;

    real dx = (domain.max(0)-domain.min(0))/(nx-1);
    real xl = domain.min(0);
    real xh = xl+dx;
    real mx = (domain.min(0)+domain.max(0))/2.0;
    real sx = stdNe;

    real dz = (domain.max(2)-domain.min(2))/(nz-1);
    real zl = domain.min(2);
    real zh = zl+dz;
    real mz = (domain.min(2)+domain.max(2))/2.0;
    real sz = stdNe;

    real probx = rightInt(xh, xl, mx, sx, dx)-rightInt(xl, xl, mx, sx, dx);
    real probz = rightInt(zh, zl, mz, sz, dz)-rightInt(zl, zl, mz, sz, dz);
    real prob = probx*probz;
    Ne(0,all,0) = prob;
    real max = prob;
    for (int i=1; i<nx-1; ++i) {
      probx = leftInt(xh, xh, mx, sx, dx)-leftInt(xl, xh, mx, sx, dx);
      xl = xh;
      xh = xl+dx;
      probx += rightInt(xh, xl, mx, sx, dx)-rightInt(xl, xl, mx, sx, dx);

      zl = domain.min(2);
      zh = zl+dz;
      for (int k=1; k<nz-1; ++k) {
        probz = leftInt(zh, zh, mz, sz, dz)-leftInt(zl, zh, mz, sz, dz);
        zl = zh;
        zh = zl+dz;
        probz += rightInt(zh, zl, mz, sz, dz)-rightInt(zl, zl, mz, sz, dz);
        prob = probx*probz;
        Ne(i, all, k) = prob;
        max = std::max(max, prob);
      }
    }
    probx = leftInt(xh, xh, mx, sx, dx)-leftInt(xl, xh, mx, sx, dx);
    probz = leftInt(zh, zh, mz, sz, dz)-leftInt(zl, zh, mz, sz, dz);
    prob = probx*probz;
    Ne(nx-1,all,nz-1) = prob;
    max = std::max(max, prob);

    Ne = (dNe/max)*Ne + minNe;
#endif

  }

#endif

  void PoissonBoltzmannSolver::initialise()
  {

    NonLinearSolver::initialise();

    // From Mamun and Cairns, Stability of solitary waves in a magnetized
    // non-thermal plasma, J. Plasma Physics, 56, 175-185, 1996
    beta = 4.0*alpha/(1.0+3.0*alpha);

    Tmin = minTe/refTe;
    dT   = (maxTe-minTe)/refTe;
    s    = stdTe;
    s2   = blitz::pow2(stdTe);

    mgTemp.resize(numGrid);
    for (int kgrid=0; kgrid<numGrid; ++kgrid) {
      mgTemp(kgrid).resize(mgGridSize(kgrid));
      getTemp(mgTemp(kgrid));
    }
#if 0
    {
      // save the electron temperature field on finest grid
      std::ofstream id("Te.out");
      id << "Te = " << mgTemp(numGrid-1) << endl;
    }
#endif
    if ( !nhmFilename.empty() ) {
      std::ifstream ifs(nhmFilename.c_str());
      if (nhmFilename.empty() || ifs.fail()) {
        ostringstream os;
        os << "Unable to open file:  " << nhmFilename;
        throw ClassException("PoissonBoltzmannSolver", os.str());
      }
      FieldVeci modelSize;
      // file generated using matlab function write_nhbgrd_model
      for (int d=0; d<DIM; d++) {
        ifs >> modelSize(d);
        if (ifs.fail()) {
          cerr << "modelSize(" <<  d << ")=" << modelSize(d) << endl;
          throw ClassException("PoissonBoltzmannSolver", "ifs.bad()");
        }
      }
#if 0
      {
      // read first line that contains comment
      std::string comment;
      std::getline(ifs,comment);
      }
#endif
      relDens.resize(modelSize);
#if (DIM==3)
      for (int k=0; k<modelSize(2); ++k)
#endif
        for (int j=0; j<modelSize(1); ++j)
          for (int i=0; i<modelSize(0); ++i) { // Fastest changing variable saved
            if (ifs.eof() || ifs.bad()) {
              ostringstream os;
              os << "Premature end of file on " << nhmFilename;
              throw ClassException("NonHomBackgrd", os.str());
            }
#if (DIM==2)
            ifs >> relDens(i,j);
#elif (DIM==3)
            ifs >> relDens(i,j,k);
#endif
          }
      mgDens.resize(numGrid);
      for (int kgrid=0; kgrid<numGrid; ++kgrid) {
        mgDens(kgrid).resize(mgGridSize(kgrid));
        interpolate(relDens, mgDens(kgrid));
      }
    } else {
      if (neShape != constant) {
        mgDens.resize(numGrid);
        for (int kgrid=0; kgrid<numGrid; ++kgrid) {
          mgDens(kgrid).resize(mgGridSize(kgrid));
          getDens(mgDens(kgrid));
        }

      }
    }
		if (debugLevel() >= 5) {
      // save the electron temperature field on finest grid
		  saveMatlab("mudfas.m",ios::out|ios::trunc, "Te", mgTemp(numGrid-1));
      if (neShape != constant || !nhmFilename.empty() ) {
        // save the unperturbed electron density field on finest grid
        //std::ofstream id("ne.out");
        //id << "ne = " << mgDens(numGrid-1) << endl;
        saveMatlab("mudfas.m",ios::out|ios::app, "ne", mgDens(numGrid-1));
      }
    }
  }


  void PoissonBoltzmannSolver::setPhiAndRhs(const Field &phi, const Field &rho)
  {
    Phi.resize(gridSize);
    Phi = phi;

    Rhs.resize(gridSize);
    // \nabla^2\phi -n_e\exp(\phi/T_e) = -n_i
    Rhs = -rho;
  }

  void PoissonBoltzmannSolver::getElectronDensity(Field &Ne) const
  {
    using blitz::pow2;

    Range I(1, gridSize(0));
    Range J(1, gridSize(1));
#if (DIM==2)

    Field phi(mgPhi(numGrid-1)(I,J));
#elif (DIM==3)

    Range K(1, gridSize(2));
    Field phi(mgPhi(numGrid-1)(I,J,K));
#endif

    Field Te(mgTemp(numGrid-1));

    Ne.resize(phi.shape());
    Ne = (1.0-beta*phi/Te+beta*pow2(phi/Te))*exp(phi/Te);

    if (neShape != constant || !nhmFilename.empty()) {
      Ne *= mgDens(numGrid-1);
    }
  }

  template <class T_numtype>
  void PoissonBoltzmannSolver::getElectricFieldEnergy(T_numtype &ESE) const
  {
    using blitz::pow2;

    Range I(1, gridSize(0));
    Range J(1, gridSize(1));
#if (DIM==2)

    Field phi(mgPhi(numGrid-1)(I,J));
#elif (DIM==3)

    Range K(1, gridSize(2));
    Field phi(mgPhi(numGrid-1)(I,J,K));
#endif

    Field Te(mgTemp(numGrid-1));

    Field rho(phi.shape());

    if (neShape != constant || !nhmFilename.empty()) {
      rho = -mgRhs(numGrid-1)-mgDens(numGrid-1)*(1.0-beta*phi/Te+beta*pow2(phi/Te))*exp(phi/Te);
    } else
      rho = -mgRhs(numGrid-1)-(1.0-beta*phi/Te+beta*pow2(phi/Te))*exp(phi/Te);

#if 0

    cout << "mean rho=" << blitz::mean(rho) << std::endl;
#endif
    // ESE from Plasma physics via computer simulation
    // Birsdall and Langdon, p.74
    Field Energy(rho*phi);
    FieldVecr dr((domain.max-domain.min)/(gridSize-1));
    mudfas::real dI;
    ESE = static_cast<T_numtype>(0.5*integrate(Energy, dr, dI));
  }

  template
  void PoissonBoltzmannSolver::getElectricFieldEnergy<float>
  (float &ESE) const;

  template
  void PoissonBoltzmannSolver::getElectricFieldEnergy<double>
  (double &ESE) const;

  template <class T_numtype>
  void PoissonBoltzmannSolver::getFluidEnergy(T_numtype &KT, T_numtype &KP) const
  {
    using blitz::pow2;

    Range I(1, gridSize(0));
    Range J(1, gridSize(1));
#if (DIM==2)

    Field phi(mgPhi(numGrid-1)(I,J));
#elif (DIM==3)

    Range K(1, gridSize(2));
    Field phi(mgPhi(numGrid-1)(I,J,K));
#endif

    Field Te(mgTemp(numGrid-1));

    FieldVecr dr((domain.max-domain.min)/(gridSize-1));
    mudfas::real dI;

    if (neShape != constant || !nhmFilename.empty()) {
      // first term kinetic energy (see chapter Energy in ~/doc/picsim.tex)
      Field Energy(mgDens(numGrid-1)*Te*(1.0-beta*phi/Te+beta*pow2(phi/Te))*exp(phi/Te));
      KT = static_cast<T_numtype>(integrate(Energy, dr, dI));
      // second term kinetic energy (see chapter Energy in ~/doc/picsim.tex)
      Energy = -phi*mgDens(numGrid-1)*(1.0-beta*phi/Te+beta*pow2(phi/Te))*exp(phi/Te);
      KP = static_cast<T_numtype>(integrate(Energy, dr, dI));
    } else {

      // first term kinetic energy (see chapter Energy in ~/doc/picsim.tex)
      Field Energy(Te*(1.0-beta*phi/Te+beta*pow2(phi/Te))*exp(phi/Te));
      KT = static_cast<T_numtype>(integrate(Energy, dr, dI));
      // second term kinetic energy (see chapter Energy in ~/doc/picsim.tex)
      Energy = -phi*(1.0-beta*phi/Te+beta*pow2(phi/Te))*exp(phi/Te);
      KP = static_cast<T_numtype>(integrate(Energy, dr, dI));
    }
#if 0
    // Normalisation
    Energy = (1.0-beta*phi/Te+beta*pow2(phi/Te))*exp(phi/Te);
    real C(integrate(Energy, dr, dI));
    KT /= C;
    KP /= C;
#endif

  }

  template
  void PoissonBoltzmannSolver::getFluidEnergy<float>
  (float &KT,
   float &KP) const;

  template
  void PoissonBoltzmannSolver::getFluidEnergy<double>
  (double &KT,
   double &KP) const;

  void PoissonBoltzmannSolver::printOn(std::ostream& os) const
  {
    MultiGridSolver::printOn(os);

    using parser::header;
    using parser::map_elt;

    os << header(DIMSTR "-d Poisson-Boltzmann setup")
       << header(DIMSTR "-d electron temperature model setup")
       << "Te min              = " << minTe << '\n'
       << "Te max              = " << maxTe << '\n'
       << "Ref Temp            = " << refTe << '\n'
       << "Te spatial std      = " << stdTe << '\n'
       << "Temp shape          = " << map_elt(teShapeMap, teShape) << '\n'
       << header(DIMSTR "-d electron density model setup")
       << "Non thermal alpha   = " << alpha << '\n';
    if (nhmFilename.empty()) {
      os
          << "Ne min              = " << minNe << '\n'
          << "Ne max              = " << maxNe << '\n'
          << "Ne spatial std      = " << stdNe << '\n'
          << "Ne shape            = " << map_elt(neShapeMap, neShape) << std::endl;
    } else {
      os
          << "NH Model Filename   = " << nhmFilename << '\n';
    }
  }


  void PoissonBoltzmannSolver::constMin(Field &Te) const
  {
    Te = Tmin;
  }

#if (DIM==2)

  void PoissonBoltzmannSolver::yGauss(Field &Te) const
  {
    using blitz::secondDim;
    using blitz::tensor::j;

    int ny1 = Te.cols()-1;
    int jst = Te.lbound(secondDim);
    // dT*exp(-y^2/(2stdTe^2)) where y \in [-1,1]
    // is equivalent to
    // dT*exp(-Y^2/\sigma^2) where Y \in [Ymin,Ymax]
    // and with \sigma = \sqrt{2}(Ymin-Ymax) stdTe
    Te = dT*exp(-pow2(2.0*(j-jst)/ny1-1.0)/(2.0*s2))+Tmin;
  }

  void PoissonBoltzmannSolver::xGauss(Field &Te) const
  {
    using blitz::firstDim;
    using blitz::tensor::i;

    int nx1 = Te.rows()-1;
    int ist = Te.lbound(firstDim);
    Te = dT*exp(-pow2(2.0*(i-ist)/nx1-1.0)/(2.0*s2))+Tmin;
  }

  void PoissonBoltzmannSolver::xDep(Field &Te) const
  {
    using blitz::firstDim;
    using blitz::tensor::i;

    int nx1 = Te.rows()-1;
    int ist = Te.lbound(blitz::firstDim);
    Te = dT*(1.0-exp(-pow2(2.0*(i-ist)/nx1-1.0)/(2.0*s2)))+Tmin;
  }

  void PoissonBoltzmannSolver::xErf(Field &Te) const
  {
    using blitz::firstDim;
    using blitz::tensor::i;

    int nx1 = Te.rows()-1;
    int ist = Te.lbound(firstDim);
    Te = 0.5*dT*(1.0+blitz::erf((2.0*(i-ist)/nx1-1.0)/s))+Tmin;
  }

  void PoissonBoltzmannSolver::xAtan(Field &Te) const
  {
    using blitz::firstDim;
    using blitz::tensor::i;

    int nx1 = Te.rows()-1;
    int ist = Te.lbound(firstDim);
    Te = 0.5*dT*(2.0*blitz::atan((2.0*(i-ist)/nx1-1.0)/s)/M_PI+1.0)+Tmin;
  }

#elif (DIM==3)

  void PoissonBoltzmannSolver::xyGauss(Field &Te) const
  {
    using blitz::firstDim;
    using blitz::tensor::i;
    using blitz::secondDim;
    using blitz::tensor::j;

    int nx1 = Te.rows()-1;
    int ny1 = Te.cols()-1;
    int ist = Te.lbound(firstDim);
    int jst = Te.lbound(secondDim);
    Te = dT*exp(-(pow2(2.0*(i-ist)/nx1-1.0) +
                  pow2(2.0*(j-jst)/ny1-1.0))/(2.0*s2))+Tmin;
  }

  void PoissonBoltzmannSolver::xzGauss(Field &Te) const
  {
    using blitz::firstDim;
    using blitz::tensor::i;
    using blitz::thirdDim;
    using blitz::tensor::k;

    int nx1 = Te.rows()-1;
    int nz1 = Te.depth()-1;
    int ist = Te.lbound(firstDim);
    int kst = Te.lbound(thirdDim);
    Te = dT*exp(-(pow2(2.0*(i-ist)/nx1-1.0) +
                  pow2(2.0*(k-kst)/nz1-1.0))/(2.0*s2))+Tmin;
  }

#endif

  void PoissonBoltzmannSolver::getNonLinearFunction(Field &phi,
      Field &F, int kgrid) const
  {

    int nx = mgGridSize(kgrid)(0);
    int ny = mgGridSize(kgrid)(1);
    Range I(1,nx), J(1,ny);
#if (DIM==3)

    int nz = mgGridSize(kgrid)(2);
    Range K(1,nz);
#endif

    Field Te(mgTemp(kgrid));

    // non-linear function = -\exp(\phi/T_e)
    if (alpha == 0.0) {
#if (DIM==2)
      F = -exp(phi(I, J)/Te);
#elif (DIM==3)

      F = -exp(phi(I, J, K)/Te);
#endif

    } else {
      // From Mamun and Cairns, Stability of solitary waves in a magnetized
      // non-thermal plasma, J. Plasma Physics, 56, 175-185, 1996
#if (DIM==2)
      F = phi(I, J)/Te;
#elif (DIM==3)

      F = phi(I, J, K)/Te;
#endif

      F = -(1.0-beta*F+beta*pow2(F))*exp(F);
    }

    if (neShape != constant || !nhmFilename.empty()) {
      F *= mgDens(kgrid);
    }
  }

  void PoissonBoltzmannSolver::getNonLinearFunctionAndDerivative(Field &phi,
      Field &F, Field &dF, int kgrid) const
  {
    int nx = mgGridSize(kgrid)(0);
    int ny = mgGridSize(kgrid)(1);
    Range I(1,nx), J(1,ny);
#if (DIM==3)

    int nz = mgGridSize(kgrid)(2);
    Range K(1,nz);
#endif

    Field Te(mgTemp(kgrid));

    // Nonlinear function = -\exp(\phi/T_e)
    if (alpha == 0.0) {
#if (DIM==2)
      F = -exp(phi(I, J)/Te);
#elif (DIM==3)

      F = -exp(phi(I, J, K)/Te);
#endif

      dF = 1.0/Te*F;
    } else {
      // From Mamun and Cairns, Stability of solitary waves in a magnetized
      // non-thermal plasma, J. Plasma Physics, 56, 175-185, 1996
#if (DIM==2)
      F = phi(I, J)/Te;
#elif (DIM==3)

      F = phi(I, J, K)/Te;
#endif

      dF = -(-beta/Te+beta*2.0*F/Te)*exp(F);
      F  = -(1.0-beta*F+beta*pow2(F))*exp(F);
      dF = dF+1.0/Te*F;
    }

    if (neShape != constant || !nhmFilename.empty()) {
      F  *= mgDens(kgrid);
      dF *= mgDens(kgrid);
    }
  }

  void PoissonBoltzmannSolver::setCof(int dim, const Array1dr &r, Array1dr &crr,
                                      Array1dr &cr, Array1dr &cer) const
  {
    crr = 1.0;
    cr  = 0.0;
    cer = 0.0;
  }

  void PoissonBoltzmannSolver::QuasiNeutralityCorrection(Field &Phi, Field &rhs)
  {
    using blitz::exp;
    using std::endl;

    FieldVecr dr(domain.max-domain.min);
    real V = blitz::product(dr);

    dr = dr / (gridSize-1);

    real ni = -integrate(rhs, dr)/V;
    cout << "ni= " << ni << endl;

    Field Te(Phi.shape());
    getTemp(Te);

    Field Ne((1.0-beta*Phi/Te+beta*pow2(Phi/Te))*exp(Phi/Te));

    real ne = integrate(Ne, dr)/V;
    cout << "ne " << ne << endl;
    cout << "ne-ni= " << std::abs(ne-ni) << endl;

    real MeanPhi = blitz::mean(blitz::abs(Phi));
#if 0

    cout << "MeanPhi= " << MeanPhi << endl;
#endif

    real Cst1 = -2.0*MeanPhi;
    Ne = (1.0-beta*(Phi-Cst1)/Te+beta*pow2((Phi-Cst1)/Te))*exp((Phi-Cst1)/Te);
    real dn1 = integrate(Ne, dr)/V-ni;
    real Cst2 = 2.0*MeanPhi;
    Ne = (1.0-beta*(Phi-Cst2)/Te+beta*pow2((Phi-Cst2)/Te))*exp((Phi-Cst2)/Te);
    real dn2 = integrate(Ne, dr)/V-ni;

    if ( dn1*dn2 > 0 ) {
      std::cerr << "*** Exception in QuasiNeutralityCorrection: dn1*dn2 > 0"
                << dn1*dn2<< endl;
      throw(EXIT_FAILURE);
    }

    real Newdn, eps = tolmax, NewCst;
    do {
      NewCst = 0.5*(Cst1+Cst2);
      Ne = (1.0-beta*(Phi-NewCst)/Te+beta*pow2((Phi-NewCst)/Te))*
           exp((Phi-NewCst)/Te);
      Newdn = integrate(Ne, dr)/V-ni;
      if ( Newdn*dn2 < 0 ) {
        Cst1 = NewCst;
        dn1 = Newdn;
      } else {
        Cst2 = NewCst;
        dn1 = dn2;
        dn2 = Newdn;
      }
    } while ( std::abs(Cst1-Cst2) > eps );
    cout << "Quasi Neutrality Correction= " << 0.5*(Cst1+Cst2) << endl;

    Phi = Phi-.5*(Cst1+Cst2);
    Ne = (1.0-beta*Phi/Te+beta*pow2(Phi/Te))*exp(Phi/Te);
    ne = integrate(Ne, dr)/V;
    cout << "ne " << ne << endl;
  }

  void PoissonBoltzmannSolver::PhiCorrection(Field &phi, Field &rhs)
  {
    Field phi_corrected(rhs.shape());

    Range I(1,rhs.rows()), J(1,rhs.cols());

    phi_corrected = phi(I,J);
    QuasiNeutralityCorrection(phi_corrected, rhs);
    phi(I,J) = phi_corrected;

    // Calculate new defect
    real mdk, sdk;
    getDefect(numGrid-1, mdk, sdk);
    std::cout << "mean |defect|= " << mdk << " std |defect|= "<< sdk << std::endl;
  }

  void PoissonBoltzmannSolver::interpolate(const Field & src, Field & dst)
  {
    using blitz::TinyVector;

    FieldVeci srcSize(src.shape()), dstSize(dst.shape());
    FieldVecr srcDr(1.0/(srcSize-1))  , dstDr(1.0/(dstSize-1));
    Field src1(FieldVeci(srcSize+1));
    src1 = 0.0;
#if (DIM==2)
    src1(Range(0,srcSize(0)-1),Range(0,srcSize(1)-1)) = src;
#elif (DIM==3)
    src1(Range(0,srcSize(0)-1),Range(0,srcSize(1)-1),Range(0,srcSize(2)-1)) = src;
#endif

#if (DIM==3)
    for (int k=0; k<dstSize(2); ++k)
#endif
      for (int j=0; j<dstSize(1); ++j)
        for (int i=0; i<dstSize(0); ++i) {
#if (DIM==2)
          FieldVecr rDest(dstDr(0)/srcDr(0)*i,dstDr(1)/srcDr(1)*j);
#elif (DIM==3)
          FieldVecr rDest(dstDr(0)/srcDr(0)*i,dstDr(1)/srcDr(1)*j,dstDr(2)/srcDr(2)*k);
#endif
          FieldVeci I(blitz::floor(rDest));
          FieldVecr Wl(rDest-I);
          FieldVecr Wu(1.0-Wl);
          //cout << "i,j="<<i<<","<<j<<
          //" rDest="<<rDest<<" I="<<I<<" Wl="<<Wl<<" Wu="<<Wu<<std::endl;
#if (DIM==2)
          TinyVector<double,4> a(src1(I(0), I(1)),
                                 src1(I(0), I(1)+1),
                                 src1(I(0)+1, I(1)),
                                 src1(I(0)+1, I(1)+1));
          TinyVector<double,4> b(Wu(0)*Wu(1), Wu(0)*Wl(1),
                                 Wl(0)*Wu(1), Wl(0)*Wl(1));
          dst(i,j) = static_cast<real>(blitz::dot(a,b));
#elif (DIM==3)
          TinyVector<double,8> a(src1(I(0),I(1),I(2)),
                                 src1(I(0),I(1),I(2)+1),
                                 src1(I(0),I(1)+1,I(2)),
                                 src1(I(0),I(1)+1,I(2)+1),
                                 src1(I(0)+1,I(1),I(2)),
                                 src1(I(0)+1,I(1),I(2)+1),
                                 src1(I(0)+1,I(1)+1,I(2)),
                                 src1(I(0)+1,I(1)+1,I(2)+1));
          TinyVector<double,8> b(Wu(0)*Wu(1)*Wu(2), Wu(0)*Wu(1)*Wl(2),
                                 Wu(0)* Wl(1)*Wu(2), Wu(0)* Wl(1)*Wl(2),
                                 Wl(0)*Wu(1)*Wu(2),  Wl(0)*Wu(1)*Wl(2),
                                 Wl(0)* Wl(1)*Wu(2),  Wl(0)* Wl(1)*Wl(2));
          dst(i,j,k) = static_cast<real>(blitz::dot(a,b));;
#endif
        }

  }

  void PoissonBoltzmannSolver::initParsing(int nargs, char *args[])
  {
    registerClass("PoissonBoltzmannSolver");
    registerPackage(PACKAGE, VERSION "\n" ID "\n", MUDFAS_COPYRIGHT);

    using parser::types::integer;
    using parser::types::real;
    using parser::types::charStr;

    insertOption(_minTe  , "temin"  , real   , "Minimum temperature"     , Any(minTe));
    insertOption(_maxTe  , "temax"  , real   , "Maximum temperature"     , Any(maxTe));
    insertOption(_refTe  , "reftemp", real   , "Reference temperature"   , Any(refTe));
    insertOption(_stdTe  , "stdte"  , real   , "Spatial temperature STD" , Any(stdTe));
    insertOption(_teShape, "teshape", integer, "Temperature shape"       , Any(teShape));
    insertOption(_alpha  , "alpha"  , real   , "Near-Boltzmann expansion", Any(alpha));

    insertOption(_minNe  , "nemin"  , real   , "Minimum density"         , Any(minNe));
    insertOption(_maxNe  , "nemax"  , real   , "Maximum density"         , Any(maxNe));
    insertOption(_stdNe  , "stdne"  , real   , "Spatial density STD"     , Any(stdNe));
    insertOption(_neShape, "neshape", integer, "Density shape"           , Any(neShape));

    insertOption(_nhmfilename, "nhmfilename", charStr,
                 "Filename for non-homogeneous model", Any(nhmFilename));
  }

  void PoissonBoltzmannSolver::paramParsing()
  {
    parseOption(_minTe  , minTe);
    parseOption(_maxTe  , maxTe);
    parseOption(_refTe  , refTe);
    parseOption(_stdTe  , stdTe);
    parseOption(_teShape, teShape);

    parseOption(_nhmfilename, nhmFilename);

    parseOption(_alpha  , alpha);

    if (nhmFilename.empty()) {
      parseOption(_minNe  , minNe);
      parseOption(_maxNe  , maxNe);
      parseOption(_stdNe  , stdNe);
      parseOption(_neShape, neShape);
    }
  }

  void PoissonBoltzmannSolver::checkParam() const
  {
    checkMap(_teShape, teShapeMap, teShape);
    if (nhmFilename.empty()) {
      checkMap(_neShape, neShapeMap, neShape);
    }
  }

#undef ID

}
