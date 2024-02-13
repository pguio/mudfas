
mudfas_srcdir = $(top_srcdir)/src
cpptools_srcdir = $(top_srcdir)/src/c++-tools
mudfas_libdir = $(top_srcdir)/lib

AM_CPPFLAGS = -I$(mudfas_srcdir) -I$(cpptools_srcdir)

utils_src = \
$(cpptools_srcdir)/parser.cpp \
$(cpptools_srcdir)/integrate.cpp \
$(cpptools_srcdir)/hdf-interface.cpp

mudfas_src = \
$(mudfas_srcdir)/init.cpp \
$(mudfas_srcdir)/linear-solver.cpp \
$(mudfas_srcdir)/mudfas.cpp \
$(mudfas_srcdir)/nonlinear-solver.cpp \
$(mudfas_srcdir)/poisson-boltzmann.cpp \
$(mudfas_srcdir)/poisson.cpp \
$(mudfas_srcdir)/prolon.cpp \
$(mudfas_srcdir)/restrict.cpp

lib_LTLIBRARIES = libmudfas2d.la libmudfas3d.la libmudfas-utils.la

if COMPILER_IS_COMPAQ_CXX
templates_mudfas2d_dir = -ptr $(mudfas_libdir)/mudfas2d
templates_mudfas3d_dir = -ptr $(mudfas_libdir)/mudfas3d
templates_mudfas_utils_dir = -ptr $(mudfas_libdir)/mudfas-utils
endif

## 2-dimensions library
libmudfas2d_la_SOURCES = $(mudfas_src)
libmudfas2d_la_CXXFLAGS = -DDIM=2 $(templates_mudfas2d_dir)
libmudfas2d_la_LDFLAGS = -version-info 0:0:0

## 3-dimensions library
libmudfas3d_la_SOURCES = $(mudfas_src)
libmudfas3d_la_CXXFLAGS = -DDIM=3 $(templates_mudfas3d_dir)
libmudfas3d_la_LDFLAGS = -version-info 0:0:0

## utility library
libmudfas_utils_la_SOURCES = $(utils_src)
libmudfas_utils_la_CXXFLAGS = $(templates_mudfas_utils_dir)
libmudfas_utils_la_LDFLAGS = -version-info 0:0:0

## Clean up template repositories
clean-local:
	-rm -rf mudfas2d  mudfas3d mudfas-utils ti_files ii_files

