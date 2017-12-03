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
// File: bw_limit.hpp
// Purpose:

#ifndef BW_LIMIT_HPP
#define BW_LIMIT_HPP

//#include <iostream>        // provides cout
//#include <fstream>
#include <cstdlib>         // provides EXIT_SUCCESS
#include <cmath>
//#include <cstdio>
//#include <string>
//#include <iomanip>         // provides setw()

#include <complex>

//#include <mpi.h>

#include <fftw3-mpi.h>
//#include "../include/fftw3-mpi.h"

namespace TEM_NS
{
   int bw_limit( fftw_complex* in, 
         const ptrdiff_t& Nx, const double* const kx,
         const ptrdiff_t& Ny, const double* const ky, 
         const double& cutoff ); //, const int& mynode );

   int bw_limit_renormalize( fftw_complex* in, 
         const ptrdiff_t& Nx, const double* const kx,
         const ptrdiff_t& Ny, const double* const ky, 
         const double& norm,
         const double& cutoff ); //, const int& mynode );
} // TEM_NS

#endif
