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
// File: probe_wavefunction.hpp
// Purpose:
//
// Simulate the conventional TEM of a thick sample as viewed perpendicular
//  to the third axis (q[2]).

#ifndef PROBE_WAVEFUNCTION_HPP
#define PROBE_WAVEFUNCTION_HPP

#include <iostream>     // provides cout, cin, and cerr
#include <cstdlib>         // provides EXIT_SUCCESS, EXIT_FAILURE 
#include <vector>
#include <list>
#include <limits> 

#include <mpi.h>
#include <fftw3-mpi.h>
//#include "../include/fftw3-mpi.h"   // fftw_complex type

#include "scatterer.hpp"            // atom struct with fitting parameters
#include "projected_atomic_potential.hpp" // \sigma V_{z}
#include "func.hpp"                       // domain implementations
#include "write_tif.hpp"                  // tif image output functions
#include "bw_limit.hpp"                   // bandwidth limit functions
#include "lenses.hpp"                     // lense transfer functions
#include "slice.hpp"            

#ifndef PI
#define PI 3.14159265369
#endif

using std::cout;
using std::cerr;
using std::endl;
using std::string;

namespace TEM_NS
{

   // The wave function in realspace calculated via an inverse Fourier 
   //  transform:
   // Kirkland (2010) eq 5.47
   // \psi_{p}(\vec{x}, \vec{x}_{p}) = A_{p} FT^{-1}[
   //       \exp[ -i \Chi( \vec{k}) + 2 \pi i \vec{k} \cdot \vec{x}_{p}]
   //       A(\vec{k})
   //    ]
   //
   // where A(\vec{k}) is the aperture function (eq 5.28, =1 inside, else 0)
   // and A_{p} is a normalization constant chosen to yield
   // \int |\psi_{p}(\vec{x}, \vec{x}_{p})|^{2} d^{2}\vec{x} = 1

   // Kirkland (2016) eq 16 is identical to Kirkland (2010) eq 5.45, from 
   //  which Kirkland (2010) eq 5.47 is derived.

int probe_wavefunction_uncorrected_unnormalized(
      // Postcondition:
      //    - psi contains the values of eq 5.47 in reciprocal space, 
      //       lacking the factor A_{p}, as evaluated over kx and ky domain.
      const double& x_p,   // x component of \vec{x}_{p}
      const double& y_p,   // y component of \vec{x}_{p}
      const double* const kx,
      const ptrdiff_t& Nx,
      const double* const ky,
      const ptrdiff_t& Ny,
      const double& Cs3,      // spherical aberration
      const double& defocus,  // focal point distance from sample bottom
      const double& alpha_max_sqr, // objective aperture limiting angle
      const double& lambda,
      const double& lambda_sqr,
      fftw_complex* psi
      );

int probe_wavefunction_correctedtoCs5_unnormalized(
      // Postcondition:
      //    - psi contains the values of eq 5.47 in reciprocal space, 
      //       lacking the factor A_{p}, as evaluated over kx and ky domain.
      const double& x_p,   // x component of \vec{x}_{p}
      const double& y_p,   // y component of \vec{x}_{p}
      const double* const kx,
      const ptrdiff_t& Nx,
      const double* const ky,
      const ptrdiff_t& Ny,
      const double& Cs3,      // spherical aberration
      const double& Cs5,
      const double& defocus,  // focal point distance from sample bottom
      //const double& defocus_spread,  // \Delta_{0}
      const double& alpha_max_sqr, // objective aperture limiting angle
      const double& lambda,
      const double& lambda_sqr,
      fftw_complex* psi
      );

int probe_wavefunction_norm(
      // Precondition: 
      //    - psi contains the values of eq 5.47 in reciprocal space, 
      //       lacking the factor A_{p}, as evaluated over kx and ky domain.
      // Postcondition: 
      //    - Ap contains the normalizing factor A_{p} from eq 5.47
      double& Ap,
      const ptrdiff_t& Nx,
      const ptrdiff_t& Ny,
      fftw_complex* psi,
      //const int& mynode,
      //const int& rootnode,
      //const int& totalnodes,
      MPI_Comm comm
      );

int probe_wavefunction_norm(
      // Precondition: 
      //    - psi contains the values of eq 5.47 in reciprocal space, 
      //       lacking the factor A_{p}, as evaluated over kx and ky domain.
      // Postcondition: 
      //    - Ap contains the normalizing factor A_{p} from eq 5.47
      double& Ap,
      const ptrdiff_t& Nx,
      const ptrdiff_t& Ny,
      const double* const psi_re,
      const double* const psi_im,
      //const int& mynode,
      //const int& rootnode,
      //const int& totalnodes,
      MPI_Comm comm
      );

//int parallel_probe_wavefunction_uncorrected_unnormalized(
//      // Postcondition:
//      //    - psi contains a disc of values 1.0, with 0.0 elsewhere,
//      //       requiring normalization by A_{p}, as evaluated over the kx 
//      //       and ky domains.
//      const double& x_p,   // x component of \vec{x}_{p}
//      const double& y_p,   // y component of \vec{x}_{p}
//      const double* const kx,
//      const ptrdiff_t& Nx,
//      const double* const ky,
//      const ptrdiff_t& Ny,
//      const double& alpha_max_sqr, // objective aperture limiting angle
//      const double& lambda_sqr,
//      fftw_complex* psi
//      );

} // TEM_NS

#endif
