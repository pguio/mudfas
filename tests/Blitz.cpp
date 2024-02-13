/**************************************************************************
 *
 * $Id: Blitz.cpp,v 1.12 2011/03/26 12:56:40 patrick Exp $
 *
 * Copyright (c) 2000-2011 Patrick Guio <patrick.guio@gmail.com>
 *
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
using namespace std;

#include <blitz/array.h>
using namespace blitz;
using namespace blitz::tensor;

template <class Ti, class Tj>
void foo(Array<double,2> &fine,
         Ti &IM1, Ti &I, Ti &IP1, Tj &JM1, Tj &J, Tj &JP1,
         Array<double,2> &coarse, Ti &IC, Tj &JC)
{
  coarse(IC,JC) =
    (fine(IM1,JM1) + fine(IP1,JM1) + fine(IM1,JP1) + fine(IP1,JP1) +
     2.0*(fine(IM1,J) + fine(IP1,J) + fine(I,JM1) + fine(I,JP1)) +
     4.0*fine(I,J))*0.0625;
}

#define ADJFN(st,fn,stride) ( (fn) - (((fn)-(st))%(stride)) )
inline
int adjfn(int st, int fn, int stride)
{
  return fn-(fn-st)%stride;
}

int main(int nargs, char *args[])
{

  const int n=5;

  Array<double,2> A(n,n), B(n,n);
  Range IM1, I, IP1, IC;
  Range JM1, J, JP1, JC;

  A=0;
  B= 2.0 + tensor::i + tensor::j;

  cout << "B= " << B << endl;

  I.setRange(1,n-2,1);
  IM1=I-1;
  IP1=I+1;

  IC=I;

  J.setRange(1,n-2,1);
  JM1=J-1;
  JP1=J+1;

  JC=J;

  foo(B,IM1,I,IP1,JM1,J,JP1,A,IC,JC);
  cout << "A= " << A << endl;

  I.setRange(0,n-1,n-1);
  IM1.setRange(n-2,n-2,1);
  IP1.setRange(1,1,1);

  IC=I;

#if 0
  foo(B,IM1,I,IP1,JM1,J,JP1,A,IC,JC);
  cout << "A= " << A << endl;
#endif

  int i, im1, ip1, ic;
  i=0;
  im1=n-2;
  ip1=1;

  ic=i;
  foo(B,im1,i,ip1,JM1,J,JP1,A,ic,JC);

  i=n-1;

  ic=i;
  foo(B,im1,i,ip1,JM1,J,JP1,A,ic,JC);

  cout << "A= " << A << endl;

  cout << "Test adjfn...." << endl;

  I.setRange(1,ADJFN(1,8,2),2);
  cout << "With ADJFN = " << I << endl;
  I.setRange(1,adjfn(1,8,2),2);
  cout << "With adjfn = " << I << endl;

  return 0;
}
