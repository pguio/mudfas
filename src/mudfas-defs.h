/**************************************************************************
 *
 * $Id: mudfas-defs.h,v 1.47 2017/11/26 18:07:25 patrick Exp $
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

#ifndef MUDFAS_DEFS_H
#define MUDFAS_DEFS_H

#if defined(HAVE_CONFIG_H)
#include <mudfas-config.h>
#endif

#include <iostream>
#include <blitz/array.h>

#if defined(HAVE_MPI)
#include <mpi.h>
#endif

#include <parser.h>

namespace mudfas {

#define MUDFAS_COPYRIGHT \
"Copyright (c) 2000-2011 Patrick Guio <patrick.guio@gmail.com>\n\n"\
"This is free software; see the source for copying conditions.  There is NO\n"\
"warranty; not even for MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.\n"

#if defined(FLOAT_FIELD)
  typedef float real;
#elif defined(DOUBLE_FIELD)

  typedef double real;
#else
#error macro FLOAT_FIELD or DOUBLE_FIELD must be defined
#endif

  typedef blitz::Array<real,1> Array1dr;
  typedef blitz::Array<real,2> Array2dr;
  typedef blitz::Array<real,3> Array3dr;

  typedef blitz::TinyVector<real,2> Vector2dr;
  typedef blitz::TinyVector<real,3> Vector3dr;

  typedef blitz::TinyVector<int,2> Vector2di;
  typedef blitz::TinyVector<int,3> Vector3di;

  // Define DIM fpr Mudfas solver
#if defined(DIM) || defined(DIMR)
#if !defined(DIM)
#define DIM DIMR
#endif
#else
#error macro DIM or DIMR must be defined
#endif

  typedef blitz::Array<real,DIM> Field;
  typedef blitz::Array<int,DIM> iField;

  typedef blitz::TinyVector<Field,DIM> VectorField;

  typedef blitz::TinyVector<real,DIM> FieldVecr;
  typedef blitz::TinyVector<int,  DIM> FieldVeci;

  template<class T>
  std::ostream& operator<<(std::ostream &os, const blitz::TinyVector<T,DIM> &v)
  {
    os << v(0);
    for (int d=1; d<DIM; ++d)
      os << ',' << v(d);
    return os;
  }

#if (DIM == 2)
  const FieldVeci DEFAULT_GRIDSIZE(101,101);
  const FieldVecr DEFAULT_MIN(0.0,0.0);
  const FieldVecr DEFAULT_MAX(100.0,100.0);
#elif (DIM == 3)

  const FieldVeci DEFAULT_GRIDSIZE(37,37,37);
  const FieldVecr DEFAULT_MIN(0.0,0.0,0.0);
  const FieldVecr DEFAULT_MAX(36.0,36.0,36.0);
#endif

  // Define string DIMSTR
#if (DIM==2)
#define DIMSTR "2"
#elif (DIM==3)
#define DIMSTR "3"
#endif


#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#if defined(HAVE_MPI)
#if defined(FLOAT_FIELD)
#define MPI_FIELD MPI_FLOAT
#elif defined(DOUBLE_FIELD)
#define MPI_FIELD MPI_DOUBLE
#endif
#endif // HAVE_MPI

template<class T>
void saveMatlab(const std::string filename, T mode,
                const std::string varname, const real var)
{
  std::ofstream os(filename.c_str(), mode);
  os.setf(std::ios::scientific);
  os.precision(6);
  os << varname << " = " << var << ";\n";
}

template<class T>
void saveMatlab(const std::string filename, T mode,
                const std::string varname, const Array1dr var)
{
  std::ofstream os(filename.c_str(), mode);
  os.setf(std::ios::scientific);
  os.precision(6);
  os << varname << " = [ ";
  for (int i=0; i<var.size(); ++i) {
    os << var(i) << "; ";
    if (i != 0 && (i%10) == 0)  os << "...\n";
  }
  os << "];\n";
}

template<class T>
void saveMatlab(const std::string filename, T mode,
                const std::string varname, const Array2dr var)
{
  using blitz::firstDim;
  using blitz::secondDim;
  std::ofstream os(filename.c_str(), mode);
  os.setf(std::ios::scientific);
  os.precision(6);
  os << varname << " = [ ";
  for (int j=0; j<var.length(secondDim); ++j) {
    for (int i=0; i<var.length(firstDim); ++i) {
      os << var(i,j) << "; ";
      if (i != 0 && (i%10) == 0) 
			  os << "...\n";
    }
  }
  os << "];\n";
	os << varname << " = reshape(" << varname << "," 
	   << var.length(firstDim) << ","
	   << var.length(secondDim) << ");\n";
}

template<class T>
void saveMatlab(const std::string filename, T mode,
                const std::string varname, const Array3dr var)
{
  using blitz::firstDim;
  using blitz::secondDim;
  using blitz::thirdDim;
  std::ofstream os(filename.c_str(), mode);
  os.setf(std::ios::scientific);
  os.precision(6);
  for (int k=0; k<var.length(thirdDim); ++k) {
    for (int j=0; j<var.length(secondDim); ++j) {
      for (int i=0; i<var.length(firstDim); ++i) {
        os << var(i,j,k) << "; ";
        if (i != 0 && (i%10) == 0)  os << "...\n";
      }
		}
  }
  os << "];\n";
	os << varname << " = reshape(" << varname << "," 
	   << var.length(firstDim) << ","
	   << var.length(secondDim) << ","
	   << var.length(thirdDim) << ");\n";
}

}

#endif
