/* ----------------------------------------------------------------------
    MS-STEM-FEM - Multislice Scanning & Fluctuation Transmission Electron 
    Microscopy Simulator
    Copyright (C) 2017 Nicholas Huebner Julian <njulian@ucla.edu>

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along
    with this program; if not, write to the Free Software Foundation, Inc.,
    51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

    See the README file in the top-level MS-STEM-FEM directory.
---------------------------------------------------------------------- */
// File: indexed_vector_magnitude_sqr.hpp
// Purpose:
//    - define a class to be used to obtain a vector of indices sorted by
//       the magnitude of the vectors associated with the indices,
//       e.g. sort indices (i,j) by the magnitudes of some vector 
//        [kx_{i,j}, ky_{i,j}]

#ifndef INDEXED_VECTOR_MAGNITUDE_SQR_HPP
#define INDEXED_VECTOR_MAGNITUDE_SQR_HPP

#include <cstdlib>

#include <algorithm> // sort()

#include <fftw3-mpi.h>
//#include "../include/fftw3-mpi.h"   // ptrdiff_t

using std::cerr;

namespace TEM_NS
{

class indexed_vector_magnitude_sqr
{
   public:
      ptrdiff_t i;
      ptrdiff_t j;
      double v_mag_sqr;
      double phi;

   // Constructor should take i, j, and the square of the magnitude of a two
   //  element array.
   indexed_vector_magnitude_sqr(  const ptrdiff_t& ii, const ptrdiff_t& jj, 
         const double& kk_sqr, const double& theta
         )
   {
      i = ii;
      j = jj;
      //v_mag_sqr = (v[0] * v[0] ) + (v[1] * v[1]);
      v_mag_sqr = kk_sqr;
      phi = theta;
   }

   // copy constructor
   //indexed_vector_magnitude( const indexed_vector_magnitude& in)
   //{
   //   i = in.i;
   //   j = in.j;
   //   v_mag_sqr = in.v_mag_sqr;
   //}

};// indexed_vector_magnitude

inline bool indexed_vector_magnitude_sqr_lt(
      const indexed_vector_magnitude_sqr& lhs,
      const indexed_vector_magnitude_sqr& rhs
      )
{
   return lhs.v_mag_sqr < rhs.v_mag_sqr;
}

} // TEM_NS
#endif
