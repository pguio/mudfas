
EXTRA_DIST = test-solver.h

mudfas_srcdir = $(top_srcdir)/src
cpptools_srcdir = $(top_srcdir)/src/c++-tools

mudfas_libdir = $(top_srcdir)/lib

AM_CPPFLAGS = -I$(mudfas_srcdir) -I$(cpptools_srcdir)

check_PROGRAMS = Blitz Parser Integrate-float Integrate-double
bin_PROGRAMS = PoissonBoltzmann2d PoissonBoltzmann3d Poisson2d Poisson3d

if COMPILER_IS_COMPAQ_CXX
templates_mudfas2d_dir = -ptr $(mudfas_libdir)/mudfas2d
templates_mudfas3d_dir = -ptr $(mudfas_libdir)/mudfas3d
templates_mudfas_utils_dir = -ptr $(mudfas_libdir)/mudfas-utils
templates_local_dir = -ptr cxx_repository
endif

Blitz_SOURCES = Blitz.cpp

Parser_SOURCES = Parser.cpp 
Parser_LDFLAGS = $(templates_local_dir) $(templates_mudfas_utils_dir)
Parser_LDADD = $(mudfas_libdir)/libmudfas-utils.la

Integrate_float_SOURCES = Integrate.cpp
Integrate_float_CPPFLAGS = -DReal=float $(AM_CPPFLAGS)
Integrate_float_LDFLAGS = $(templates_local_dir) $(templates_mudfas_utils_dir)
Integrate_float_LDADD = $(mudfas_libdir)/libmudfas-utils.la

Integrate_double_SOURCES = Integrate.cpp
Integrate_double_CPPFLAGS = -DReal=double $(AM_CPPFLAGS)
Integrate_double_LDFLAGS = $(templates_local_dir) $(templates_mudfas_utils_dir)
Integrate_double_LDADD = $(mudfas_libdir)/libmudfas-utils.la

PoissonBoltzmann2d_SOURCES = TestSolver.cpp
PoissonBoltzmann2d_CPPFLAGS = -DPOISSON_BOLTZMANN -DDIM=2 $(AM_CPPFLAGS)
PoissonBoltzmann2d_LDFLAGS = $(templates_local_dir) $(templates_mudfas2d_dir) \
$(templates_mudfas_utils_dir)
PoissonBoltzmann2d_LDADD = $(mudfas_libdir)/libmudfas2d.la \
$(mudfas_libdir)/libmudfas-utils.la

PoissonBoltzmann3d_SOURCES = TestSolver.cpp
PoissonBoltzmann3d_CPPFLAGS = -DPOISSON_BOLTZMANN -DDIM=3 $(AM_CPPFLAGS)
PoissonBoltzmann3d_LDFLAGS = $(templates_local_dir) $(templates_mudfas3d_dir) \
$(templates_mudfas_utils_dir)
PoissonBoltzmann3d_LDADD = $(mudfas_libdir)/libmudfas3d.la \
$(mudfas_libdir)/libmudfas-utils.la

Poisson2d_SOURCES = TestSolver.cpp
Poisson2d_CPPFLAGS = -DPOISSON -DDIM=2 $(AM_CPPFLAGS)
Poisson2d_LDFLAGS = $(templates_local_dir) $(templates_mudfas2d_dir) \
$(templates_mudfas_utils_dir)
Poisson2d_LDADD = $(mudfas_libdir)/libmudfas2d.la \
$(mudfas_libdir)/libmudfas-utils.la

Poisson3d_SOURCES = TestSolver.cpp
Poisson3d_CPPFLAGS = -DPOISSON -DDIM=3 $(AM_CPPFLAGS)
Poisson3d_LDFLAGS = $(templates_local_dir) $(templates_mudfas3d_dir) \
$(templates_mudfas_utils_dir)
Poisson3d_LDADD = $(mudfas_libdir)/libmudfas3d.la \
$(mudfas_libdir)/libmudfas-utils.la

## Clean up template repositories
clean-local:
	-rm -rf cxx_repository ti_files ii_files

check:
	@echo Running test suite...
	@for prog in $(check_PROGRAMS); do \
		echo $$prog && (./$$prog > /dev/null 2>&1) || exit 1; \
	done
	@echo All tests passed.

