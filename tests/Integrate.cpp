/**************************************************************************
 *
 * $Id: Integrate.cpp,v 1.23 2011/03/26 12:56:40 patrick Exp $
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

#include <iostream>
#include <typeinfo>
using namespace std;

#include <parser.h>
using namespace parser;

#include <integrate.h>

#ifndef DIM
#define DIM 2 // dummy setting for including mudfas_defs.h
#endif
#include <mudfas-defs.h>

using blitz::tensor::i;
using blitz::tensor::j;
using blitz::tensor::k;
using blitz::pow2;
using blitz::pow3;
using blitz::Range;

template <class T_numtype>
void display_result(const char *title, T_numtype exact, T_numtype I,
                    T_numtype dI, T_numtype I0, T_numtype I1)
{
  cout << endl << title << " " << typeid(exact).name() <<  endl << endl
       << "Exact=" << exact << endl
       << "I    =" << I << " abs(I -Exact)=" << fabs(I -exact)
       << " dI=" << fabs(dI) << endl
       << "I0   =" << I0 << " abs(I0-Exact)=" << fabs(I0-exact) << endl
       << "I1   =" << I1 << " abs(I1-Exact)=" << fabs(I1-exact) << endl;
}

#define CMP_INT1D(type)                                  \
{                                                        \
blitz::Array<type,1> f(R0);                              \
blitz::TinyVector<type,1> dr;                            \
f = exp(-pow2(d*i/n1+a)/pow2(s));                        \
dr = d/n1;                                               \
type dI;                                                 \
type I  = integrate(f, dr, dI);                          \
type I0 = sum(f(R1))*product(dr);                        \
type I1 = sum(f(R2))*product(dr);                        \
type exact = int1d_value;                                \
display_result("1D integral", exact, I, dI, I0, I1);     \
}                                                        \
 
int main(int nargs, char *args[])
{
  double I, dI, I0, I1;
  double s(0.5), a(-1.0), b(1.0);
  int n(100);

#if defined(HAVE_MPI)
  MPI_Init(&nargs, &args);
#endif


  enum parser_enum { nb=1, std, lower, upper };
  Parser parser(nargs, args);
  parser.registerProgram(args[0]);
  parser.registerPackage(PACKAGE, VERSION, MUDFAS_COPYRIGHT);

  using parser::types::integer;
  using parser::types::real;

  parser.insertOption(nb   , "nb" , integer, "Number of points"  , Any(n));
  parser.insertOption(std  , "std", real   , "Standard deviation", Any(s));
  parser.insertOption(lower, "a"  , real   , "Lower limit"       , Any(a));
  parser.insertOption(upper, "b"  ,real   , "Upper limit"       , Any(b));

  if (parser.parseHelp() || parser.parseVersion())
    return 0;

  parser.parseOption(nb, n);
  parser.parseOption(std, s);
  parser.parseOption(lower, a);
  parser.parseOption(upper, b);

  if (a >= b) {
    cerr << "a must be strictly smaller than b" << endl;
    return EXIT_FAILURE;
  }

  const int n1 = n-1;
  const int n2 = n-2;

  const Range R0(Range(0,n1));
  const Range R1(Range(0,n2));
  const Range R2(Range(1,n1));

  const double d = b-a;

  const double int1d_value(1.0/M_2_SQRTPI*s*(erf(b/s)-erf(a/s)));

  cout.precision(8);
  cout.setf(ios::scientific);

#if 0
  {
    blitz::Array<double,1> f(R0);
    blitz::TinyVector<double,1> dr;
    f = exp(-pow2(d*i/n1+a)/pow2(s));
    dr = d/n1;
    I  = integrate(f, dr, dI);
    I0 = sum(f(R1))*product(dr);
    I1 = sum(f(R2))*product(dr);
    double exact = int1d_value;
    display_result("1D integral", exact, I, dI, I0, I1);
  }
#endif
  CMP_INT1D(double)
  CMP_INT1D(float)


  {
    blitz::Array<double,2> f(R0,R0);
    blitz::TinyVector<double,2> dr;
    f = exp(-(pow2(d*i/n1+a)+pow2(d*j/n1+a))/pow2(s));
    dr = d/n1;
    I  = integrate(f, dr, dI);
    I0 = sum(f(R1,R1))*product(dr);
    I1 = sum(f(R2,R2))*product(dr);
    double exact(pow2(int1d_value));
    display_result("2D integral",exact, I, dI, I0, I1);
  }

  {
    blitz::Array<double,3> f(R0,R0,R0);
    blitz::TinyVector<double,3> dr;
    f = exp(-(pow2(d*i/n1+a)+pow2(d*j/n1+a)+pow2(d*k/n1+a))/pow2(s));
    dr = d/n1;
    I = integrate(f, dr, dI);
    I0 = sum(f(R1,R1,R1))*product(dr);
    I1 = sum(f(R2,R2,R2))*product(dr);
    double exact(pow3(int1d_value));
    display_result("3D integral",exact, I, dI, I0, I1);
  }

#if defined(HAVE_MPI)
  MPI_Finalize();
#endif

  return 0;
}
