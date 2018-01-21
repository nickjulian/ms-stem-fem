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
// File: io_txt.hpp
// Purpose: to define tools for writing data to text files

#ifndef IO_TXT_HPP
#define IO_TXT_HPP

#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>   // EXIT_SUCCESS & EXIT_FAILURE 
//#include <cstdio>   
//#include <limits>   // std::numeric_limits::digits to find bits per sample
//#include <string>
//#include <string.h>  // required for memcpy
#include <iomanip>   // setprecision()
//#include <list>
#include <vector>

#include <mpi.h>

//#include <fftw3-mpi.h>
//#include "../include/fftw3-mpi.h"

//#include "radial_discretization.hpp"

using std::cout;
using std::endl;
using std::string;
using std::cerr;
//using std::setprecision;

namespace TEM_NS
{

int output_variance_to_txt( 
      const double* const outDataRaw,  // array of range elements
      const std::vector<double>& binning_boundaries, // domain elements
      const string& outFilePrefix
      );

//int output_psi_realspace_to_txt( 
//      const fftw_complex* const psi,  // array of range elements
//      const ptrdiff_t& local_alloc_size_fftw,
//      const ptrdiff_t Nx_local,
//      const double* const kx_joined, // domain elements
//      const ptrdiff_t Nx,
//      const double* const ky, // domain elements
//      const ptrdiff_t Ny,
//      const string& outFilePrefix,
//      const int& mynode,
//      const int& rootnode,
//      MPI_Comm comm
//      );

//int output_psi_reciprocalspace_to_txt( 
//      const fftw_complex* const psi,  // array of range elements
//      const ptrdiff_t& local_alloc_size_fftw,
//      const ptrdiff_t Nx_local,
//      const double* const kx_joined, // domain elements
//      const ptrdiff_t Nx,
//      const double* const ky, // domain elements
//      const ptrdiff_t Ny,
//      const string& outFilePrefix,
//      const int& mynode,
//      const int& rootnode,
//      MPI_Comm comm
//      );

//int output_psi_mag_reciprocalspace_to_txt( 
//      const double* const psi,  // array of range elements
//      const ptrdiff_t& local_alloc_size_fftw,
//      const ptrdiff_t Nx_local,
//      const double* const kx_joined, // domain elements
//      const ptrdiff_t Nx,
//      const double* const ky, // domain elements
//      const ptrdiff_t Ny,
//      const string& outFilePrefix,
//      const int& mynode,
//      const int& rootnode,
//      MPI_Comm comm
//      );

} // TEM_NS
#endif
