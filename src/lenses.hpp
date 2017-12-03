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
// File: lenses.hpp
// Purpose:
//

#ifndef LENSES_HPP
#define LENSES_HPP

#include <iostream>     // provides cout, cin, and cerr
#include <cstdlib>         // provides EXIT_SUCCESS, EXIT_FAILURE 
#include <vector>
#include <limits> 

//#include "../include/fftw3-mpi.h"   // fftw_complex type

#include "scatterer.hpp" 

#ifndef PI
#define PI 3.14159265369
#endif

using std::cout;
using std::endl;

namespace TEM_NS
{

void objective_lens_transfer_function_coher(
          const double& ksqr, 
          const double& Cs3, 
          const double& defocus,
          const double& alpha_max_sqr,
          const double& lambda,
          const double& lambdasqr,
          double& H_0_re,
          double& H_0_im
          );

void objective_lens_transfer_function_wp(
          const double& ksqr, 
          const double& Cs3, 
          const double& defocus,
          const double& alpha_max_sqr,
          const double& lambda,
          const double& lambdasqr,
          double& H_0_re,
          double& H_0_im
          );

void objective_lens_transfer_function_wp_correctedtoCs5(
          const double& ksqr, 
          const double& Cs3, 
          const double& Cs5, 
          const double& defocus,
          const double& alpha_max_sqr,
          const double& lambda,
          const double& lambdasqr,
          double& H_0_re,
          double& H_0_im
          );

void objective_lens_transfer_function_wp_correctedtoCs5_spread(
          const double& ksqr,
          const double& Cs3,
          const double& Cs5,
          const double& defocus,
          const double& defocus_spread,   // \Delta_{0}
          const double& condenser_illumination_angle, // \lambda k_{s}
          const double& alpha_max_sqr,
          const double& lambda,
          const double& lambdasqr,
          double& H_0_re,
          double& H_0_im
          );

double aberration_function_uncorrected(
      const double& ksqr,
      const double& Cs3,
      const double& defocus,
      const double& lambda,
      const double& lambdasqr
      );

double aberration_function_correctedtoCs5(
      const double& ksqr,
      const double& Cs3,
      const double& Cs5,
      const double& defocus,
      const double& lambda,
      const double& lambdasqr
      );

} // TEM_NS

#endif
