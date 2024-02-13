/**************************************************************************
 *
 * $Id: boundary.h,v 1.14 2011/11/07 18:39:29 patrick Exp $
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

#ifndef BOUNDARY_H
#define BOUNDARY_H

#include <iostream>
#if defined(BLITZ2)
#include <blitz/tinyvec2.h>
#else
#include <blitz/tinyvec-et.h>
#endif

template <class T_numtype, int N_dim>
class Boundary {
  friend std::ostream& operator<<(std::ostream& os, const Boundary &boundary) {
    os << "min (";
    for (int d=0; d<N_dim-1; ++d) os << boundary.min(d) << ", ";
    os << boundary.min(N_dim-1) << ")\t";
    os << "max (";
    for (int d=0; d<N_dim-1; ++d) os << boundary.max(d) << ", ";
    os << boundary.max(N_dim-1) << ')';
    return os;
  }

public:
  Boundary(T_numtype _min=0, T_numtype _max=0) : min(_min), max(_max)
  {}
  template <class T_min, class T_max>
  Boundary(T_min _min, T_max _max) : min(_min), max(_max)
  {}
  template <class T_min, class T_max>
  Boundary(T_min *_min, T_max *_max) {
    for (unsigned i=0; i<min.length(); ++i)
      min(i) = _min[i];
    for (unsigned i=0; i<max.length(); ++i)
      max(i) = _max[i];

  }
  ~Boundary()
  {}

  template <class T_min, class T_max>
  void copy(T_min *_min, T_max *_max) {
    for (unsigned i=0; i<min.length(); ++i)
      _min[i] = min(i);
    for (unsigned i=0; i<max.length(); ++i)
      _max[i] = max(i);
  }

  Boundary& operator=(const T_numtype val) {
    this->min = val;
    this->max = val;
    return *this;
  }
  template <class P_numtype>
  Boundary& operator=(const P_numtype val) {
    this->min = val;
    this->max = val;
    return *this;
  }
  Boundary& operator+=(const Boundary &val) {
    this->min += val.min;
    this->max += val.max;
    return *this;
  }
  template <class P_numtype, int P_dim>
  Boundary& operator+=(const Boundary<P_numtype, P_dim> &val) {
    this->min += val.min;
    this->max += val.max;
    return *this;
  }

  blitz::TinyVector<T_numtype, N_dim> min;
  blitz::TinyVector<T_numtype, N_dim> max;
};

#endif // BOUNDARY_H
