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
// File: variance.hpp
// Purpose:

#ifndef VARIANCE_HPP
#define VARIANCE_HPP

#include <cstdlib>
#include <cmath>
#include <string> 

#ifndef PI
#define PI 3.14159265369
#endif

// required for debug_output_complex_fftw_operand()
#include <mpi.h>
#include <iomanip> // provides setw
#include <vector>
//#include <list>

#include <fftw3-mpi.h>
//#include "../include/fftw3-mpi.h"

#include "indexed_vector_magnitude_sqr.hpp"
#include "func.hpp"
#include "radial_discretization.hpp"
#include "io_netcdf.hpp"
#include "io_txt.hpp"


namespace TEM_NS
{

   int variance_1D_STEM( 
         const double* const diffracted_wave_radial_intensity_sum,
         const double* const diffracted_wave_radial_intensity_sqr_sum,
         const std::vector<double>& binning_boundaries,
         const size_t& number_of_raster_points,
         const unsigned int& input_flag_netcdf_variance,
         const string& outFilePrefix,
         double* variance
         );

   int variance_2D_STEM(
         // average of intensity over beam position raster
         const double* const data2D_avgs_local,      // denominator
         // average of intensity^{2} over beam position raster
         const double* const data2D_sqr_avgs_local,  // numerator
         const unsigned int number_of_raster_points,
         //const std::vector<double>& binning_boundaries,
         const double* const kx_local, const size_t& Nx_local,
         const double* const kx_joined, // only non-null on rootnode
         const size_t& Nx_joined,
         const double* const ky, const size_t& Ny,
         const size_t& resolutionUnit_recip,
         const double& xResolution_recip, const double& yResolution_recip,
         const double& delta_r,
         const double& bwcutoff_t, //data known to be 0 beyond this |k|
         const string& outFilePrefix,
         const int& mynode,
         const int& rootnode,
         MPI_Comm comm
         );

   int variance_2D_BFCTEM(
         // average of intensity over beam position raster
         double* data2D_avgs_local, 
         // average of intensity^{2} over beam position raster
         double* data2D_sqr_avgs_local,
         const radial_discretization& sample_rotations,
         const string& outFilePrefix,
         const int& mynode,
         const int& rootnode,
         MPI_Comm comm
         );

   int integrate_out_theta_fftw(
         const fftw_complex* const psi,   // 2-D data input
         const double* const kx_local, const size_t& Nx_local,
         const double* const ky, const size_t& Ny,
         const std::vector<double>& binning_boundaries,
         int* bin_counts,
         double* f_of_k_magnitudes//,  // 1-D data output
         );

   int integrate_out_theta_double(
         const double* const psi,   // 2-D data input
         const double* const kx_local, const size_t& Nx_local,
         const double* const ky, const size_t& Ny,
         const std::vector<double>& binning_boundaries,
         double* f_of_k_magnitudes//,  // 1-D data output
         );

}
#endif
