/**************************************************************************
 *
 * $Id: test-solver.h,v 1.20 2012/11/09 13:51:55 patrick Exp $
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
 * solve the nonlienar 2D/3D problems
 *
 * nabla^2(phi) - exp(phi/Te) = -rho
 *
 * largely inspired from the linear multigrid method used in the MUDPACK
 * fortran package (http://www.scd.ucar.edu/css/software/mudpack/)
 *
 * Code uses the Blitz C++ library (http://oonumerics.org/blitz)
 *
 */

#if defined(HAVE_HDF4)
#include <mfhdf.h>
#if !defined(MAX_VAR_DIMS)
#define MAX_VAR_DIMS 32
#endif
#if !defined(MAX_NC_NAME)
#define MAX_NC_NAME 256
#endif
#endif

#if defined(POISSON_BOLTZMANN)
#include <poisson-boltzmann.h>
#define SOLVERCLASS PoissonBoltzmannSolver
#elif defined(POISSON)
#include <poisson.h>
#define SOLVERCLASS PoissonSolver
#else
#error macros POISSON_BOLTZMANN or POISSON must be defined
#endif

namespace mudfas {

#define ID "$Id: test-solver.h,v 1.20 2012/11/09 13:51:55 patrick Exp $"

  class TestSolver : public SOLVERCLASS {
  public:

    friend ostream& operator<<(ostream& os, const TestSolver &s);

    TestSolver(int nargs, char* args[]);
    ~TestSolver();

    virtual void initialise();
    void savePhiAndRhs();

  protected:

    enum func_enum { _f1=1, _f2, _f3, 
#if defined(HAVE_HDF4)
		_f4, 
#endif
		_f5, _last };
    enum parser_enum { _fun=SOLVERCLASS::next };

  private:

    int funId;
    LUT funMap;

    void f1(Field &phi, Field &rho);
    void f2(Field &phi, Field &rho);
    void f3(Field &phi, Field &rho);
#if defined(HAVE_HDF4)
    void f4(Field &phi, Field &rho);
#endif
    void f5(Field &phi, Field &rho);
    void (TestSolver::*f[_last])(Field &phi, Field &rho);

    void initParsing(int nargs, char *args[]);
    void paramParsing();
    void checkParam() const;
  };

  using blitz::cos;
  using blitz::pow2;

  using blitz::firstDim;
  using blitz::tensor::i;
  using blitz::tensor::j;

#if (DIM==3)

  using blitz::thirdDim;
  using blitz::tensor::k;
#endif

  using parser::header;
  using parser::map_elt;

  using std::endl;
  using std::ostream;


  ostream& operator<<(ostream& os, const TestSolver &s)
  {
    s.printOn(os);

    return os
           << header("Test solver setup")
           << "fun                 = " << map_elt(s.funMap,s.funId) << endl;
  }

  TestSolver::TestSolver(int nargs, char* args[])
    : SOLVERCLASS(nargs, args), funId(_f1)
  {
    f[_f1] = &TestSolver::f1;
    f[_f2] = &TestSolver::f2;
    f[_f3] = &TestSolver::f3;
#if defined(HAVE_HDF4)		
    f[_f4] = &TestSolver::f4;
#endif
    f[_f5] = &TestSolver::f5;
#if (DIM==2)

    funMap.insert(LUTPair(_f1, "f(x,y)            = 1+a*cos(kx*x+ky*y)"));
    funMap.insert(LUTPair(_f2, "f(x,y)            = 1, "
                          "f(1/3:1,1/3:2/3) *= 1.2"));
    funMap.insert(LUTPair(_f3, "f(x,y)            = 1 then 2"));
#if defined(HAVE_HDF4)
    funMap.insert(LUTPair(_f4, "f(x,y)            = load ./beamBS2d.hdf"));
#endif
    funMap.insert(LUTPair(_f5, "f(x,y)            = dn*exp(-x^2)"));
#elif (DIM==3)

    funMap.insert(LUTPair(_f1, "f(x,y,z)                  = 1+a*cos(kx*x+ky*y)"));
    funMap.insert(LUTPair(_f2, "f(x,y,z)                  = 1, "
                          "f(1/3:2/3,1/3:2/3,1/3:1) *= 1.2"));
    funMap.insert(LUTPair(_f3, "f(x,y,z)                  = 1 then 2"));
#if defined(HAVE_HDF4)
    funMap.insert(LUTPair(_f4, "f(x,y,z)                  = load ./beamBS3d.hdf"));
#endif
    funMap.insert(LUTPair(_f5, "f(x,y)                    = dn*exp(-x^2-z^2)"));
#endif

    initParsing(nargs, args);
    paramParsing();
    checkParam();
  }

  TestSolver::~TestSolver()
  {}

  void TestSolver::initParsing(int nargs, char *args[])
  {
    registerClass("TestSolver");
    registerPackage(PACKAGE, VERSION "\n" ID "\n", MUDFAS_COPYRIGHT);

    insertOption(_fun, "fun", parser::types::integer, "Function index", Any(funId));
  }

  void TestSolver::paramParsing()
  {
    parseOption(_fun, funId);
  }

  void TestSolver::checkParam() const
  {
    checkMap(_fun , funMap    , funId);
  }

  void TestSolver::savePhiAndRhs()
  {
    Range I(1, gridSize(0));
    Range J(1, gridSize(1));
#if (DIM==2)

    Field phi(mgPhi(numGrid-1)(I,J));
#elif (DIM==3)

    Range K(1, gridSize(2));
    Field phi(mgPhi(numGrid-1)(I,J,K));
#endif

    ofstream id_phi("phi.out");
    id_phi << "phi = " << phi << endl;

    Field rhs(mgRhs(numGrid-1));
    ofstream id_rhs("rhs.out");
    id_rhs << "rhs = " << rhs << endl;
  }

  void TestSolver::initialise()
  {

    Field phi(gridSize);
    Field rho(gridSize);

    (this->*f[funId])(phi, rho);

    SOLVERCLASS::initialise();

#if defined(POISSON)

    rho = -(rho-1.0);
#endif

    setPhiAndRhs(phi, rho);
  }

  void TestSolver::f1(Field &phi, Field &rho)
  {
    phi = 0.0;
    rho = 1.0;

    real a = 0.05;
    real kx = 10.0*M_PI;
    real ky = 8.0*M_PI;

    int nx = gridSize(0);
    int ny = gridSize(1);

    rho += a*cos(kx*i/(nx-1)+ky*j/(ny-1));
  }

  void TestSolver::f2(Field &phi, Field &rho)
  {
    phi = 0.0;
    rho = 1.0;

#if (DIM==2)

    int imn = gridSize(0)/3-1, imx = gridSize(0)-1;
    int jmn = gridSize(1)/3-1, jmx = 2*gridSize(1)/3-1;
    Range I(imn,imx), J(jmn,jmx);
    rho(I,J) *= 1.2;
#elif (DIM==3)

    int imn = gridSize(0)/3-1, imx = 2*gridSize(0)/3-1;
    int jmn = gridSize(1)/3-1, jmx = 2*gridSize(1)/3-1;
    int kmn = gridSize(2)/3-1, kmx = gridSize(2)-1;
    Range I(imn,imx), J(jmn,jmx), K(kmn,kmx);
    rho(I,J,K) *= 1.2;
#endif

  }

  void TestSolver::f3(Field &phi, Field &rho)
  {

    phi = 0.0;
    rho = 1.0;

    real a = 0.05;
    real kx = 10.0*M_PI;
    real ky = 8.0*M_PI;
    int nx = gridSize(0);
    int ny = gridSize(1);

    rho += a*cos(kx*i/(nx-1)+ky*j/(ny-1));

#if (DIM==2)

    int imn = gridSize(0)/3-1, imx = gridSize(0)-1;
    int jmn = gridSize(1)/3-1, jmx = 2*gridSize(1)/3-1;
    Range I(imn,imx), J(jmn,jmx);
    rho(I,J) *= 1.2;
#elif (DIM==3)

    int imn = gridSize(0)/3-1, imx = 2*gridSize(0)/3-1;
    int jmn = gridSize(1)/3-1, jmx = 2*gridSize(1)/3-1;
    int kmn = gridSize(2)/3-1, kmx = gridSize(2)-1;
    Range I(imn,imx), J(jmn,jmx), K(kmn,kmx);
    rho(I,J,K) *= 1.2;
#endif
  }

#if defined(HAVE_HDF4)

  void TestSolver::f4(Field &phi, Field &rho)
  {
    int32 sd_id, sds_id;
    intn status;
    int32 dim_sizes[MAX_VAR_DIMS];
    int32 start[MAX_VAR_DIMS], edges[MAX_VAR_DIMS];
    int32 rank, data_type, n_attrs;
    char name[MAX_NC_NAME];
    float *PtrHdfRho, *HdfRho;
    double lx, ly;
#if (DIM==2)

    static const char HdfFileName[] = "beamBS2d.hdf";
#elif (DIM==3)

    double lz;
    static const char HdfFileName[] = "beamBS3d.hdf";
#endif

    // Load rho from hdf file
    sd_id = SDstart(HdfFileName,DFACC_RDONLY);
    if (sd_id == FAIL) {
      HEpush((hdf_err_code_t)HEvalue(0),"fun4",__FILE__,(intn)__LINE__);
      HEprint(stderr, 0);
      throw(EXIT_FAILURE);
    }

    sds_id = SDfindattr(sd_id,"lx");
    SDreadattr(sd_id, sds_id, &lx);
    sds_id = SDfindattr(sd_id,"ly");
    SDreadattr(sd_id, sds_id, &ly);
#if (DIM==2)
#if 0

    cout<<"lx= " << lx << " ly= "<< ly << endl;
#endif
#elif (DIM==3)

    sds_id = SDfindattr(sd_id,"lz");
    SDreadattr(sd_id, sds_id, &lz);
#if 0

    cout<<"lx= " << lx << " ly= "<< ly << " lz= "<< lz << endl;
#endif
#endif

    sds_id = SDnametoindex(sd_id,"rho");
    sds_id = SDselect(sd_id,sds_id);
    status = SDgetinfo(sds_id,name,&rank,dim_sizes,&data_type,&n_attrs);
    cout<<"rank= "<<rank<<" dim_sizes= ";
    for (int i=0; i<rank; i++)
      cout<<dim_sizes[i]<<" ";
    cout<<endl;
    for (int i=0; i<rank; i++) {
      start[i] = 0;
      edges[i] = dim_sizes[i];
    }
#if (DIM==2)
    PtrHdfRho=HdfRho = new float[dim_sizes[0]*dim_sizes[1]];
#elif (DIM==3)

    PtrHdfRho=HdfRho = new float[dim_sizes[0]*dim_sizes[1]*dim_sizes[2]];
#endif

    status=SDreaddata(sds_id, start, 0, edges, static_cast<VOIDP>(HdfRho));

    status = SDendaccess (sds_id);
    status = SDend (sd_id);

    for (int i=0; i<rank; i++)
      gridSize(i) = dim_sizes[i];

    domain.min = 0.0;
#if (DIM==2)

    domain.max =  lx,  ly;
#elif (DIM==3)

    domain.max =  lx,  ly,  lz;
#endif

    phi.resize(gridSize);
    rho.resize(gridSize);

    phi = 0.0;

    for (int i=0; i<gridSize(0); i++)
      for (int j=0; j<gridSize(1); j++)
#if (DIM==2)

        rho(i,j) = static_cast<real>(*PtrHdfRho++);
#elif (DIM==3)

        for (int k=0; k<gridSize(2); k++)
          rho(i,j,k) = static_cast<real>(*PtrHdfRho++);
#endif

#if 0

    cout << "rho(Range(0,10),0)= " << rho(Range(0,10),0) << endl;
    cout << "rho(0,Range(0,10))= " << rho(0,Range(0,10)) << endl;
#endif

    delete[] HdfRho;
  }

#endif // defined(HAVE_HDF4)

  void TestSolver::f5(Field &phi, Field &rho)
  {

    phi = 0.0;
    rho = 1.0;

    real mn = 1.0;
    real mx = 2.0;
    real s  = 0.25;
    real s2 = pow2(s);

    int nx1 = rho.rows()-1;
    int ist = rho.lbound(firstDim);

#if (DIM==2)

    rho = (mn-mx)*exp(-pow2(2.0*(i-ist)/nx1-1.0)/(2.0*s2))+mn;
#elif (DIM==3)

    int nz1 = rho.depth()-1;
    int kst = rho.lbound(thirdDim);
    rho = (mn-mx)*exp(-(pow2(2.0*(i-ist)/nx1-1.0) +
                        pow2(2.0*(k-kst)/nz1-1.0))/(2.0*s2))+mn;
#endif
  }




#undef ID

}
