/**************************************************************************
 *
 * $Id: mudfas.cpp,v 1.102 2017/12/23 18:42:24 patrick Exp $
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
 * package (http://www.scd.ucar.edu/css/software/mudpack/)
 *
 * and the book "An introduction to multigrid methods" by P. Wesseling J.
 * Wiley & Sons, 1992, ISBN 0 471 93083 0
 *
 * Short description of the method:
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
 * -> Possibility to include new temperature models
 *
 * Code uses the Blitz C++ library (http://oonumerics.org/blitz)
 *
 */

#include <mudfas.h>

namespace mudfas {

  using blitz::abs;
  using blitz::max;
  using blitz::mean;
  using blitz::sqr;
  using blitz::sqrt;
  using blitz::sum;

  using parser::header;
  using parser::map_elt;
  using parser::yesno;

  using std::cout;
  using std::endl;
  using std::flush;
  using std::ios;
  using std::ostream;


#define ID "$Id: mudfas.cpp,v 1.102 2017/12/23 18:42:24 patrick Exp $"


  const MultiGridSolver::LUTPair bcPairs[] = {
    MultiGridSolver::LUTPair(MultiGridSolver::periodic , " periodic"),
    MultiGridSolver::LUTPair(MultiGridSolver::dirichlet, "Dirichlet"),
    MultiGridSolver::LUTPair(MultiGridSolver::mixed    , "    mixed"),
  };

  MultiGridSolver::LUT MultiGridSolver::bcMap(bcPairs, bcPairs+3);

  const MultiGridSolver::LUTPair kcyclePairs[] = {
    MultiGridSolver::LUTPair(MultiGridSolver::v, "V cycle"),
    MultiGridSolver::LUTPair(MultiGridSolver::w, "W cycle"),
    MultiGridSolver::LUTPair(MultiGridSolver::f, "F cycle"),
  };

  MultiGridSolver::LUT MultiGridSolver::kcycleMap(kcyclePairs, kcyclePairs+3);

  const MultiGridSolver::LUTPair intpolPairs[] = {
    MultiGridSolver::LUTPair(MultiGridSolver::linear, "linear"),
    MultiGridSolver::LUTPair(MultiGridSolver::cubic , " cubic"),
  };

  MultiGridSolver::LUT MultiGridSolver::intpolMap(intpolPairs, intpolPairs+2);

  ostream& operator<<(ostream& os, const MultiGridSolver &s)
  {
    s.printOn(os);
    return os;
  }


  MultiGridSolver::MultiGridSolver(int nargs, char* args[])
    : Parser(nargs, args), gridSize(DEFAULT_GRIDSIZE),
      bdyc(mixed,mixed), gbdy(0.0,0.0), alfa(0.0,0.0),
      domain(DEFAULT_MIN,DEFAULT_MAX),
      tolmax(5e-4), delta(5e-1), maxcy(10), iguess(false),
      kcycle(w), ipre(2), ipost(1), intpol(cubic)
  {
    initParsing(nargs, args);
    paramParsing();
    checkParam();
  }

  MultiGridSolver::~MultiGridSolver()
  {}



  void MultiGridSolver::initialise()
  {
    FieldVeci coarsestGridSize, gridNumber;
    for (int d=0; d<DIM; ++d) {
      getMultiGridParameters(gridSize(d), coarsestGridSize(d), gridNumber(d));
    }
    numGrid = max(gridNumber);

    // Allocate arrays
    mgGridSize.resize(numGrid);

    mgS.resize(numGrid);
    mgCycle.resize(numGrid);

    for (int dim=0; dim<DIM; ++dim) {
      for (int k=0; k<3; ++k) {
        mgCof(dim,k).resize(numGrid);
      }
    }

    mgPhi.resize(numGrid);
    mgRhs.resize(numGrid);

    for (int kg=0; kg<numGrid; ++kg) {
      setGrids(kg, coarsestGridSize, gridNumber);
      setSk(kg);
      setCycle(kg);
      setCofk(kg);
    }
  }

  void MultiGridSolver::getPhi(Field &phi) const
  {
    phi.resize(gridSize);

    Range I(1, gridSize(0));
    Range J(1, gridSize(1));
#if (DIM==2)

    phi = mgPhi(numGrid-1)(I,J);
#elif (DIM==3)

    Range K(1, gridSize(2));
    phi = mgPhi(numGrid-1)(I,J,K);
#endif

  }

  void MultiGridSolver::getRhs(Field &rhs) const
  {
    rhs.resize(gridSize);
    rhs = mgRhs(numGrid-1);
  }

  bool MultiGridSolver::isInitGuess() const
  {
    return iguess;
  }

  void MultiGridSolver::solve()
  {
    string sp(40u,' ');

    swk();  // set Phi, Rhs in mgPhi, mgRhs and adjust right hand side
    if (! iguess ) { // no initial guess at finest grid level!
      // transfer down to all grid levels
      for (int kg=numGrid-1; kg>=1; --kg)
        trsfc(kg-1);
      // adjust right hand side at all grid levels in case
      // rhs or specified b.c. in phi or gbdy changed
      for (int kg=0; kg<=numGrid-1; ++kg)
        adjmd(kg);
      // choose an initial guess
      smooth(mgPhi(0), mgRhs(0), 0);
      // execute one full multigrid cycle
      for (int kg=1; kg<=numGrid-1; ++kg) {
        // lift or prolong approximation mgPhi from kg-1 to kg
        prolon(mgPhi(kg-1), kg, mgPhi(kg));
        if (debugLevel() >= 5) {
          cout <<"Entering mgd with k=" << kg+1 << sp << endl;
        }
        mgd(mgPhi, mgRhs, kg);
      }
    } else { // initial guess at finest grid level!
      // transfer down to all grid levels
      for (int kg=numGrid-1; kg>=1; --kg)
        trsfc(kg-1);
      // adjust rhs at finest grid level only
      adjmd(numGrid-1);
    }
    // See discussion p.176 in Wesseling
    int kg = numGrid-1;
    real epsk = tolmax*mean(abs(mgRhs(kg)));
    real mdk, sdk;
    getDefect(kg, mdk, sdk);
    real tk = mdk-epsk;
    if (tk < 0.0) {
      std::ios::fmtflags f = cout.flags() & std::ios::floatfield;
      int p = cout.precision();
#if defined(HAVE_MPI)
#if 1 // display only masterProc
      if (rankProc == masterProc) {
#else // display diagnostic for all CPUS
      int flag;
      MPI_Status status;
      if (rankProc > 0)
        MPI_Recv(&flag, 1, MPI_INT,
                 rankProc-1, 100, MPI_COMM_WORLD, &status);
#endif
        cout << std::setprecision(2) << std::scientific
             << "tk= " << tk << " mean |defect|= " << mdk << " std |defect|= " << sdk
             << std::resetiosflags(std::ios::scientific)
             << std::setiosflags(f) << std::setprecision(p) << endl;
#if 1 // display only masterProc
      }
#else // display diagnostic for all CPUS
        if (rankProc < nbProc-1)
          MPI_Send(&rankProc, 1, MPI_INT,
                   rankProc+1, 100, MPI_COMM_WORLD);
#endif
#else // ! defined(HAVE_MPI)
      cout << std::setprecision(2) << std::scientific
           << "tk= " << tk << " mean |defect|= " << mdk << " std |defect|= " << sdk
           << std::resetiosflags(std::ios::scientific)
           << std::setiosflags(f) << std::setprecision(p) << endl;
#endif
    } else
    {
      // execute maxcy more multigrid k cycles from finest level
      for (int i=1; i<=maxcy; ++i) {
        if (debugLevel() >= 5) {
          cout << "Entering iter# " << i << " of mgd with k="
               << kg+1 << sp << endl;
        }
        mgd(mgPhi, mgRhs, kg);
        getDefect(kg, mdk, sdk);
        tk = mdk-epsk;
        if (tk < 0.0 || i == maxcy) {
          std::ios::fmtflags f = cout.flags() & std::ios::floatfield;
          int p = cout.precision();
          cout << std::setprecision(2) << std::scientific
               << "tk= " << tk << " mean |defect|= " << mdk << " std |defect|= " << sdk
               << std::resetiosflags(std::ios::scientific)
               << std::setiosflags(f) << std::setprecision(p) << endl;
          break;
        }
      }
    }
  }

  void MultiGridSolver::getDefect(int kg, real &meanDefect, real &stdDefect)
  {
    Field dk(mgRhs(kg).lbound(), mgRhs(kg).shape());
    dk = 0.0;
    diffop(mgPhi(kg), dk, kg);
    dk -= mgRhs(kg);
    meanDefect = mean(abs(dk));
    stdDefect = std::sqrt(sum(sqr(abs(dk)-meanDefect))/(dk.size()-1.0));
  }

#if 0
  void MultiGridSolver::mgd(MultiGridField &mgU, MultiGridField &mgF, int kg)
  {

    string ln(2*(kg+1),'-'),sp(2*(numGrid-kg+1),' '),bk(17+2*(numGrid+2),'\b')

    if (debugLevel() >= 10) {
      cout << (ln+'>'+sp+"In  mgd with k=") << kg+1 << bk << flush;
    }

    if (kg == 0) {
      smooth(mgU(0), mgF(0), 0);
    } else {
      for (int i=1; i<=ipre; ++i) { // pre-smoothing
        smooth(mgU(kg), mgF(kg), kg);
      }
      // allocate array for the correction to propagate
      MultiGridField mg1F(kg);
      for (int i=0; i<=kg-1; ++i) {
        mg1F(i).resize(mgF(i).shape());
        mg1F(i).reindexSelf(FieldVeci(1));
        mg1F(i) = mgF(i);
      }
      // Calculate the residual rk (-defect)
      Field rk(mgF(kg).lbound(),mgF(kg).shape());
      rk = 0.0;
      diffop(mgU(kg), rk, kg);
      rk = -rk + mgF(kg);

      Field rkc(mgF(kg-1).lbound(),mgF(kg-1).shape());
      restrict_op(rk, rkc, kg-1);

      // Initialise adaptive stuff p.175 in Wesseling
      real tkc  = mean(abs(rk));
      real ekc  = delta*mgS(kg-1)*tkc;
      tkc -= mgS(kg)*mean(abs(mgF(kg)));

      // choose mgU(kg) p.171 in Wesseling
      restrict_op(mgU(kg), mgU(kg-1), kg-1);

      Field Lukc(mgF(kg-1).lbound(),mgF(kg-1).shape());
      Lukc = 0.0;
      diffop(mgU(kg-1), Lukc, kg-1);
      mg1F(kg-1) = Lukc + mgS(kg-1)*rkc;

      // keep the coarse Phi for fine correction
      Field ukc(mgU(kg-1).lbound(),mgU(kg-1).shape());
      ukc = mgU(kg-1);

      for (int i=1; i<=mgCycle(kg); ++i) {
        mgd(mgU, mg1F, kg-1);
        // Adaptive stuff p.175 in Wesseling
        Field t(mg1F(kg-1).lbound(), mg1F(kg-1).shape());
        t = 0.0;
        diffop(mgU(kg-1), t, kg-1);
        t -= mg1F(kg-1);
        tkc = mean(abs(t))-ekc;
        if (tkc < 0.0)
          break;
        // End Adaptive stuff
      }
      // calculate an estimate for the fine correction:
      // build the coarse correction, prolongate it and inject it
      Field ckc(mgU(kg-1).lbound(),mgU(kg-1).shape());
      ckc = mgU(kg-1)-ukc;
      Field ck(mgU(kg).lbound(),mgU(kg).shape());
      ck = 0.0;
      prolon(ckc, kg, ck);
      mgU(kg) += 1.0/mgS(kg-1)*ck;

      for (int i=1; i<=ipost; ++i) { // post-smoothing
        smooth(mgU(kg), mgF(kg), kg);
      }
    }
    if (debugLevel() >= 10) {
      cout << ('<'+ln+sp+"Out mgd with k=") << kg+1 << bk << flush;
    }
  }
#endif

  void MultiGridSolver::mgd(MultiGridField &mgU, MultiGridField &mgF, int kg)
  {
    string ln(2*(kg+1),'-'),sp(2*(numGrid-kg+1),' '),bk(17+2*(numGrid+2),'\b');

    if (debugLevel() >= 10) {
      cout << (ln+'>'+sp+"In  mgd with k=") << kg+1 << bk << flush;
    }

    if (kg == 0) {
      smooth(mgU(0), mgF(0), 0);
    } else {
      for (int i=1; i<=ipre; ++i) { // pre-smoothing
        smooth(mgU(kg), mgF(kg), kg);
      }
      // Calculate the residual (-defect) rk (2) p.169 in Wesseling
      Field rk(mgF(kg).lbound(),mgF(kg).shape());
      rk = 0.0;
      diffop(mgU(kg), rk, kg);
      rk = -rk + mgF(kg);

      // Initialise adaptive stuff p.175 in Wesseling
      real tkc  = mean(abs(rk));
      real ekc  = delta*mgS(kg-1)*tkc;
      tkc -= mgS(kg)*mean(abs(mgF(kg)));

      // choose mgU(kg) p. 171 in Wesseling
      restrict_op(mgU(kg), mgU(kg-1), kg-1);

      // (4) p.169 in Wesseling
      Field rkc(mgF(kg-1).lbound(),mgF(kg-1).shape());
      restrict_op(rk, rkc, kg-1);
      Field Lukc(mgF(kg-1).lbound(),mgF(kg-1).shape());
      Lukc = 0.0;
      diffop(mgU(kg-1), Lukc, kg-1);
      mgF(kg-1) = Lukc + mgS(kg-1)*rkc;

      // keep a copy of the coarse u_{k-1} for fine grid correction correction
      Field ukc(mgU(kg-1).lbound(),mgU(kg-1).shape());
      ukc = mgU(kg-1);

      for (int i=1; i<=mgCycle(kg); ++i) {
        mgd(mgU, mgF, kg-1);
        // Adaptive stuff p.176 in Wesseling
        Field t(mgF(kg-1).lbound(), mgF(kg-1).shape());
        t = 0.0;
        diffop(mgU(kg-1), t, kg-1);
        t -= mgF(kg-1);
        tkc = mean(abs(t))-ekc;
        if (tkc < 0.0)
          break;
      }
      // calculate an estimate for the fine correction (6)  p.169 in Wesseling
      // build the coarse correction, prolongate it and inject it
      Field ckc(mgU(kg-1).lbound(),mgU(kg-1).shape());
      ckc = mgU(kg-1)-ukc;
      Field ck(mgU(kg).lbound(),mgU(kg).shape());
      ck = 0.0;
      prolon(ckc, kg, ck);
      mgU(kg) += 1.0/mgS(kg-1)*ck;
      // post-smoothing
      for (int i=1; i<=ipost; ++i) {
        smooth(mgU(kg), mgF(kg), kg);
      }
    }
    if (debugLevel() >= 10) {
      cout << ('<'+ln+sp+"Out mgd with k=") << kg+1 << bk << flush;
    }
  }

  void MultiGridSolver::printOn(ostream& os) const
  {
    os << header(DIMSTR "-d multigrid solver setup")

       << "gridsize            = "<< gridSize << '\n'
       << "## Number of Grids  = " << numGrid << '\n';
    for (int k=0; k<numGrid; ++k)
      os << "## Grid " << k+1 << " = " << mgGridSize(k) << '\n';

    os
        << "rmin                = " << domain.min << '\n'
        << "rmax                = " << domain.max << '\n'
        << "bcmin               = " << map_elt(bcMap, bdyc.min) << '\n'
        << "bcmax               = " << map_elt(bcMap, bdyc.max) << '\n'
        << "alfamin             = " << alfa.min << '\n'
        << "alfamax             = " << alfa.max << '\n'
        << "gbdymin             = " << gbdy.min << '\n'
        << "gbdymax             = " << gbdy.max << '\n'

        << "tolmax              = " << tolmax << '\n'
        << "delta               = " << delta << '\n'
        << "maxcy               = " << maxcy << '\n'
        << "iguess              = " << yesno(iguess) << '\n'
        << "kcycle              = " << map_elt(kcycleMap, kcycle) << '\n'
        << "ipre                = " << ipre << '\n'
        << "ipost               = " << ipost << '\n'
        << "intpol              = " << map_elt(intpolMap, intpol) << endl;
  }

  void MultiGridSolver::initParsing(int nargs, char *args[])
  {
    registerClass(DIMSTR "-D MultiGridSolver");
    registerPackage(PACKAGE, VERSION "\n" ID "\n", MUDFAS_COPYRIGHT);

    parseLevelDebugOption("MultiGridSolver::dl");

    using parser::types::boolean;
    using parser::types::integer;
    using parser::types::real;
    using parser::types::intVect;
    using parser::types::realVect;

    insertOption(_gridsize, "gridsize", intVect , "Grid size, dim=" DIMSTR                  , Any(gridSize));
    insertOption(_bcmin   , "bcmin"   , intVect , "Condition at lower boundary, dim=" DIMSTR, Any(bdyc.min));
    insertOption(_bcmax   , "bcmax"   , intVect , "Condition at upper boundary, dim=" DIMSTR, Any(bdyc.max));
    insertOption(_rmin    , "rmin"    , realVect, "Value of lower boundary, dim=" DIMSTR    , Any(domain.min));
    insertOption(_rmax    , "rmax"    , realVect, "Value of upper boundary, dim=" DIMSTR    , Any(domain.max));
    insertOption(_alfamin , "alfamin" , realVect, "Alfa value for mixed lower boundary, dim=" DIMSTR    , Any(alfa.min));
    insertOption(_alfamax , "alfamax" , realVect, "Alfa value for mixed upper boundary, dim=" DIMSTR    , Any(alfa.min));
    insertOption(_gbdymin , "gbdymin" , realVect, "Gbdy value for mixed lower boundary, dim=" DIMSTR    , Any(gbdy.min));
    insertOption(_gbdymax , "gbdymax" , realVect, "Gbdy value for mixed upper boundary, dim=" DIMSTR    , Any(gbdy.min));

    insertOption(_tolmax  , "tolmax"  , real    , "Maximum relative error"                  , Any(tolmax));
    insertOption(_delta   , "delta"   , real    , "Adaptive control parameter"              , Any(delta));
    insertOption(_maxcy   , "maxcy"   , integer , "Maximum number of cycles"                , Any(maxcy));
    insertOption(_iguess  , "iguess"  , boolean , "Enable initial guess start"              , Any(iguess));
    insertOption(_kcycle  , "kcycle"  , integer , "Number of cycles"                        , Any(kcycle));
    insertOption(_ipre    , "ipre"    , integer , "Number of pre-smoothing"                 , Any(ipre));
    insertOption(_ipost   , "ipost"   , integer , "Number of post-smoothing"                , Any(ipost));
    insertOption(_intpol  , "intpol"  , integer , "Interpolation type"                      , Any(intpol));
  }

  void MultiGridSolver::paramParsing()
  {
    parseOption(_gridsize, gridSize);

    parseOption(_bcmin   , bdyc.min);
    parseOption(_bcmax   , bdyc.max);

    parseOption(_rmin    , domain.min);
    parseOption(_rmax    , domain.max);

    parseOption(_alfamin , alfa.min);
    parseOption(_alfamax , alfa.max);
    parseOption(_gbdymin , gbdy.min);
    parseOption(_gbdymax , gbdy.max);

    parseOption(_tolmax  , tolmax);
    parseOption(_delta   , delta);
    parseOption(_maxcy   , maxcy);
    parseOption(_iguess  , iguess);
    parseOption(_kcycle  , kcycle);
    parseOption(_ipre    , ipre);
    parseOption(_ipost   , ipost);
    parseOption(_intpol  , intpol);
  }

  void MultiGridSolver::checkParam() const
  {
    checkMap(_bcmin , bcMap    , bdyc.min);
    checkMap(_bcmax , bcMap    , bdyc.max);
    checkMap(_kcycle, kcycleMap, kcycle);
    checkMap(_intpol, intpolMap, intpol);
  }

#undef ID

}
