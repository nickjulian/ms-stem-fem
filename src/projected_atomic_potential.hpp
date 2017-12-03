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
// File: projected_atomic_potential.hpp
// Purpose:
//

// For the test, use a fake projected_atomic_potential function to retrieve
//  LUT parameters.

#ifndef PROJECTED_ATOMIC_POTENTIAL_HPP
#define PROJECTED_ATOMIC_POTENTIAL_HPP

#include <iostream>     // provides cout, cin, and cerr
#include <cstdlib>         // provides EXIT_SUCCESS, EXIT_FAILURE 
#include <vector>
#include <limits> 
#include <cmath>

#include <fftw3-mpi.h>
//#include "../include/fftw3-mpi.h"   // ptrdiff_t, fftw_complex types

//#include "scatterer.hpp" 
//#include "scatterer_param_LUT.hpp"


#ifndef PI
#define PI 3.14159265369
#endif

using std::cout;
using std::endl;

namespace TEM_NS
{

// sigmaVz() calls sigmaVz_single_scatterer(), which calls 
//  f_ej_single_scatterer

//int sigmaVz(
//      const std::vector<const scatterer*>& myScatterers,
//      const double& lambda, const double& gamma, const double& ab_inv,
//      const double& cutoff,
//      const double* const kxdomain, 
//      const ptrdiff_t& Nx,
//      const double* const kydomain, 
//      const ptrdiff_t& Ny,
//      fftw_complex* pap
//      );

int sigmaVz_single_scatterer(
      const double& lambda, const double& gamma, const double& ab_inv,
      const double& kx, const double& ky, // kz = 0, Fourier projection thm
      const double& x, const double& y,   // z integrated out
      //unsigned int& Z,  // atomic species
      const double* const a, const double* const b, 
      const double* const c, const double* const d,
      double& sigmaVz_re,
      double& sigmaVz_im
      //fftw_complex& sigmaVz
      );

double f_ej_single_scatterer(
      const double* const a, const double* const b, 
      const double* const c, const double* const d,
      const double& qsq
      );

//double projected_atomic_potential_single_scatterer(
//      const double* const a, const double* const b, 
//      const double* const c, const double* const d,
//      const double& rsq
//      );

} // TEM_NS

#endif
