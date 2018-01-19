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
// File: adfstem.hpp
// Purpose:

#ifndef ADFSTEM_HPP
#define ADFSTEM_HPP

#include <iostream>     // provides cout, cin, and cerr
#include <cstdlib>         // provides EXIT_SUCCESS, EXIT_FAILURE 
#include <vector>
#include <list>
#include <limits> 
#include <cmath>  // floor

#include <mpi.h>
#include <fftw3-mpi.h>
//#include "../include/fftw3-mpi.h"   // fftw_complex type

#include "to_string.hpp"
#include "scatterer.hpp"            // atom struct with fitting parameters
#include "projected_atomic_potential.hpp" // \sigma V_{z}
#include "func.hpp"                       // domain implementations
#include "scherzer.hpp"                   // scherzer conditions
#include "write_tif.hpp"                  // tif image output functions
#include "bw_limit.hpp"                   // bandwidth limit functions
#include "lenses.hpp"                     // lense transfer functions
#include "slice.hpp"            
#include "probe_wavefunction.hpp"         // convergent probe wavefunction
#include "variance.hpp"
#include "io_netcdf.hpp"

#ifndef PI
#define PI 3.14159265369
#endif

using std::cout;
using std::cerr;
using std::endl;
using std::string;

namespace TEM_NS
{

int adfstem(
      const unsigned int& input_flag_fem,
      const unsigned int& input_flag_aberration_correction,
      const unsigned int& input_flag_complex_realspace_sum,
      const unsigned int& input_flag_image_output,
      const unsigned int& input_flag_netcdf_output,
      const unsigned int& input_flag_debug,
      // parameters taken from simulation of system evolution:
      //const std::list< slice* >& slicesOfScatterers,
      std::list< slice* >& slicesOfScatterers,
      //const size_t& total_population,
      const double& bwcutoff_pap,
      const double& bwcutoff_t,
      const double& azimuthal_binning_size_factor,
      //const double& xperiod, const double& yperiod, const double& zperiod,
      const double& xmin,      // lowest value of x domain
      const double& ymin,      // lowest value of y domain
      //const double& zmin,      // lowest value of z domain
      //const double& raster_spacing, // \AA units, x & y raster spacing
      const std::list<double>& x_p,
      const std::list<double>& y_p,
      const double& xperiod,
      const double& yperiod,
      const double& xperiod_duped,
      const double& yperiod_duped,
      const double* const xx_local,
      const double* const xx_joined,
      const double* const yy,
      const double* const kx_local,
      const double* const kx_joined,
      const ptrdiff_t& Nx,   // samples for each dimensional discretization
      const ptrdiff_t& Nx_local, 
      const double* const ky,
      const ptrdiff_t& Ny, 
      // parameters of the tem microscope:
      //const double& VV,       // accelerating voltage
      const double& Cs3,       // third order spherical aberration
      const double& Cs5,       // fifth order spherical aberration
      const double& defocus,  // focal point distance from sample bottom
      const double& defocus_spread,
      const double& alpha_max, // objective aperture limiting angle
      const double& detector_inner_angle, // detector dimension [\AA^{-1}]
      const double& detector_outer_angle, // detector dimension [\AA^{-1}]
      const double& lambda,
      const string& outFileName_prefix,
      const ptrdiff_t& local_alloc_size_fftw,
      const ptrdiff_t& local_0_start_fftw,
      //const ptrdiff_t& probe_idx_local_start_x,
      //const int& threads_ok,
      //const int& nthreads,
      const int& mynode,      // MPI rank of local node
      const int& rootnode,    // MPI rank of root node
      const int& totalnodes,
      MPI_Comm comm
      );

} // TEM_NS

#endif
