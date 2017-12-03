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
// File: bw_limit.cpp
// Purpose:

#ifndef BW_LIMIT_CPP
#define BW_LIMIT_CPP

#include "bw_limit.hpp"

int TEM_NS::bw_limit( fftw_complex* in, 
      const ptrdiff_t& Nx_local, const double* const kx_local, 
      const ptrdiff_t& Ny, const double* const ky, 
      const double& cutoff) //, const int& mynode )
{
   const double cutoffsq = pow( cutoff, 2);

   for ( ptrdiff_t i = 0; i < Nx_local; i++ ) 
      for ( ptrdiff_t j = 0; j < Ny; j++ ) 
         if ( pow( kx_local[i], 2) + pow( ky[j], 2) >= cutoffsq ) 
         //if ( pow( kx_local[i], 2) + pow( ky[j], 2) > cutoffsq ) 
         {
            in[j + i*Ny][0] = 0.0;
            in[j + i*Ny][1] = 0.0;
         }

   return EXIT_SUCCESS;
}

int TEM_NS::bw_limit_renormalize( fftw_complex* in, 
      const ptrdiff_t& Nx, const double* const kx, 
      const ptrdiff_t& Ny, const double* const ky, 
      const double& norm,
      const double& cutoff) //, const int& mynode )
{
   const double cutoffsq = pow( cutoff, 2);

   for ( ptrdiff_t i = 0; i < Nx; i++ ) 
      for ( ptrdiff_t j = 0; j < Ny; j++ ) 
         //if ( pow( kx[i], 2) + pow( ky[j], 2) > cutoffsq ) 
         if ( pow( kx[i], 2) + pow( ky[j], 2) >= cutoffsq ) 
         {
            in[j + i*Ny][0] = 0.0;
            in[j + i*Ny][1] = 0.0;
         }
         else
         {
            in[j + i*Ny][0] = in[j + i*Ny][0] / norm;
            in[j + i*Ny][1] = in[j + i*Ny][1] / norm;
         }


   return EXIT_SUCCESS;
}

#endif
