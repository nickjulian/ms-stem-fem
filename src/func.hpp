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
// File: func.hpp
// Purpose:
//    - define functions used to generate data for the test.
//    - (x,y) == (row, column) == (vertical, horizontal) == (height, width)

#ifndef FUNC_HPP
#define FUNC_HPP

#include <cstdlib>
#include <cmath>

#ifndef PI
#define PI 3.14159265369
#endif

// required for debug_output_complex_fftw_operand()
#include <mpi.h>
#include <iomanip> // provides setw
#include <vector>
#include <list>

#include <fftw3-mpi.h>
//#include "../include/fftw3-mpi.h"

#include "to_string.hpp"
#include "write_tif.hpp"
#include "indexed_vector_magnitude_sqr.hpp"


namespace TEM_NS
{
   void domain_1D(
         // Precondition: 
         //    - x[N] has been instantiated (memory allocated).
         // Postcondition: 
         //    - x[N] has been initialized with values.
         const int N,
         const double P,
         double *x
         );

   void domain_2D(
         // Precondition: 
         //    - qx[Nx] and qy[Ny] have been instantiated (memory allocated)
         // Postcondition: 
         //    - qx[Nx] and qy[Ny] have been initialized with values.
         const int& Nx, const int& Ny,
         const double& qmin_x, const double& qmax_x,
         const double& qmin_y, const double& qmax_y,
         double *qx, double *qy
         );

   void domain_2D_periodic(
         // Precondition: 
         //    - qx[Nx] and qy[Ny] have been instantiated (memory allocated)
         // Postcondition: 
         //    - qx[Nx] and qy[Ny] have been initialized with values.
         //       qmin + qperiod_x == qmin 
         //       Nx * dx == qperiod_x \implies dx == qperiod_x/Nx
         const int& Nx, const int& Ny,
         const double& qperiod_x, const double& qperiod_y,
         const double& qmin_x, const double& qmin_y, 
         const double& dx, const double& dy,
         double *qx, double *qy
         );

   void assign_stem_raster_points(
         // Precondition: 
         //    - qx and qy have been instantiated as std::list<double>
         // Postcondition: 
         //    - qx and qy have been initialized with values.
         //       qmin + qperiod_x == qmin 
         //       Nx * dx == qperiod_x \implies dx == qperiod_x/Nx
         const int& Nx, const int& Ny,
         const double& qperiod_x,
         const double& qmin_x, 
         const double& qperiod_y,
         const double& qmin_y, 
         //double* qx, double* qy
         std::list<double>& qx, std::list<double>& qy
         );


   void domain_2D_periodic_recip(
         // Precondition: 
         //    - kx[Nx] and ky[Ny] have been instantiated (memory allocated)
         // Postcondition: 
         //    - kx[Nx] and ky[Ny] have been initialized with the values
         //       {0, dkx, ..., kmax_x,-kmax_x,-dkx] and  
         //       {0, dky, ..., kmax_y,-kmax_y,-dky] respectively
         //       where dkx = kmax_x*2/Nx, and dky = kmax_y*2/Ny
         const int& Nx, const int& Ny,
         const double& xperiod, const double& yperiod,
         const double& kxperiod, const double& kyperiod,
         const double& dkx, const double& dky,
         double *kx, double *ky
         );

   void domain_2D_periodic_recip_pos(
         // Precondition: 
         //    - kx[Nx] and ky[Ny] have been instantiated (memory allocated)
         // Postcondition: 
         //    - kx[Nx] and ky[Ny] have been initialized with the values
         //       {0, dkx, ..., kmax_x,-kmax_x,-dkx] and  
         //       {0, dky, ..., kmax_y,-kmax_y,-dky] respectively
         //       where dkx = kmax_x*2/Nx, and dky = kmax_y*2/Ny
         const int& Nx, const int& Ny,
         const double& xperiod,  // x period in real space
         const double& yperiod,
         double *kx, double *ky
         );

   void domain_3D_periodic(
         // Precondition: 
         //    - qx[Nx] and qy[Ny] have been instantiated (memory allocated)
         // Postcondition: 
         //    - qx[Nx] and qy[Ny] have been initialized with values.
         //       qmin + qperiod_x == qmin 
         //       Nx * dx == qperiod_x \implies dx == qperiod_x/Nx
         const int& Nx, const int& Ny, const int& Nz,
         const double& qperiod_x, const double& qmin_x, 
         const double& qperiod_y, const double& qmin_y, 
         const double& qperiod_z, const double& qmin_z, 
         double *qx, double *qy, double *qz
         );

   void domain_3D_periodic_recip(
         // Precondition: 
         //    - kx[Nx] and ky[Ny] have been instantiated (memory allocated)
         // Postcondition: 
         //    - kx[Nx] and ky[Ny] have been initialized with the values
         //       {0, dkx, ..., kmax_x,-kmax_x,-dkx] and  
         //       {0, dky, ..., kmax_y,-kmax_y,-dky] respectively
         //       where dkx = kmax_x*2/Nx, and dky = kmax_y*2/Ny
         const int& Nx, const int& Ny, const int& Nz,
         const double& xperiod,  // x period in real space
         const double& yperiod,
         const double& zperiod,
         double* kx, double* ky, double* kz
         );

   void step_2D(
         // Precondition: 
         //    - x[N] has been initialized with values.
         //    - f[N] has been instantiated (memory allocated).
         // Postcondition: 
         //    - f[N] has been initialized with values of 0 or 1 a function 
         //       of x[N]
         //       
         const int& Nx, const int& Ny,
         const double& q0x, const double& q0y,
         const double* const qx, const double* const qy,
         double *f
         );

   void cos_1D(
         // Precondition: 
         //    - W is an integer (to ensure continuity)
         //    - x[N] has been initialized with values.
         //    - f[N] has been initialized with values.
         // Postcondition: 
         //    - f[i] values have been incremented by A*cos( W * (x[i]-x0) )
         const int N,
         const double A,
         const double W,
         const double x0,
         double *x,
         double *f
         );

   void gaussian_1D(
         // Precondition: 
         //    - x[N] has been initialized with values.
         //    - f[N] has been initialized with values.
         // Postcondition: 
         //    - f[N] values have been incremented by a gaussian curve
         //       centered at mu, with stddev sigma and magnitude A.
         const int N,
         const double A,
         const double mu,
         const double sigma,
         double *x,
         double *f
         );

   void gaussian_2D(
         // Precondition: 
         //    - q[] & f[] are row-major arrays of 2-D data
         //    - the domain q[Nx*Ny] has been initialized with values.
         //    - the range f[Nx*Ny] has been initialized with values.
         // Postcondition: 
         //    - f[Nx*Ny] values have been incremented by a gaussian curve
         //       centered at (mux,muy), with stddev sigma and magnitude A.
         const int& Nx,
         const int& Ny,
         const double& A,
         const double& mux, const double& muy,
         const double& sigma,
         const double* const qx,
         const double* const qy,
         double *f
         );

   void gaussian_2D_lattice_square(
         // Precondition: 
         //    - q[] & f[] are row-major arrays of 2-D data
         //    - the domain q[Nx*Ny] has been initialized with values.
         //    - the range f[Nx*Ny] has been initialized with values.
         // Postcondition: 
         //    - f[Nx*Ny] values have been incremented by a gaussian curve
         //       centered at (mux,muy), with stddev sigma and magnitude A.
         const int& Nx,
         const int& Ny,
         const double& A,
         const int& Nkx, const int& Nky,
         const double& sigma,
         const double* const qx,
         const double* const qy,
         double *f
         );

   void positions_2D_lattice_square(
         const size_t& Nx,
         const size_t& Ny,
         const size_t& Nz,
         const double* const x,
         const double& x_size,
         const double* const y,
         const double& y_size,
         double* const q
         );

   void positions_3D_lattice_diamond_001(
         const size_t& x_cells,  // number of unit cells, x direction
         const size_t& y_cells,  // number of unit cells, y direction
         const size_t& z_cells,  // number of unit cells, z direction
         const double* const x,  // x domain
         const double& x_pixels, // number of points in discretized x domain
         const double* const y,  // y domain
         const double& y_pixels, // number of points in discretized y domain
         const double* const z,  // x domain
         const double& z_pixels, // number of points in discretized z domain
         double* const q         // lattice point positions (x,y,z), output
         );

   void positions_3D_lattice_diamond_001(
         const size_t& x_cells,  // number of unit cells, x direction
         const size_t& y_cells,  // number of unit cells, y direction
         const size_t& z_cells,  // number of unit cells, z direction
         const double* const x,  // x domain
         const double& x_pixels, // number of points in discretized x domain
         const double* const y,  // y domain
         const double& y_pixels, // number of points in discretized y domain
         const double* const z,  // x domain
         const double& z_pixels, // number of points in discretized z domain
         std::list<double> q     // lattice point positions (x,y,z), output
         );

   void positions_3D_lattice_diamond_011(
         const size_t& x_cells,  // number of unit cells, x direction
         const size_t& y_cells,  // number of unit cells, y direction
         const size_t& z_cells,  // number of unit cells, z direction
         const double* const x,  // x domain
         const double& x_pixels, // number of points in discretized x domain
         const double* const y,  // y domain
         const double& y_pixels, // number of points in discretized y domain
         const double* const z,  // x domain
         const double& z_pixels, // number of points in discretized z domain
         double* const q         // lattice point positions (x,y,z), output
         );

   void positions_3D_lattice_diamond_111_bad(
         const size_t& x_cells,  // number of unit cells, x direction
         const size_t& y_cells,  // number of unit cells, y direction
         const size_t& z_cells,  // number of unit cells, z direction
         const double* const x,  // x domain
         const double& x_pixels, // number of points in discretized x domain
         const double* const y,  // y domain
         const double& y_pixels, // number of points in discretized y domain
         const double* const z,  // x domain
         const double& z_pixels, // number of points in discretized z domain
         double* const q         // lattice point positions (x,y,z), output
         );

//   void positions_3D_lattice_diamond_111(
//         const size_t& x_cells,  // number of unit cells, x direction
//         const size_t& y_cells,  // number of unit cells, y direction
//         const size_t& z_cells,  // number of unit cells, z direction
//         const double* const x,  // x domain
//         const double& x_pixels, // number of points in discretized x domain
//         const double* const y,  // y domain
//         const double& y_pixels, // number of points in discretized y domain
//         const double* const z,  // x domain
//         const double& z_pixels, // number of points in discretized z domain
//         std::list<double> q     // lattice point positions (x,y,z), output
//         );

   void positions_3D_lattice_diamond_111(
         const size_t& x_cells,  // number of unit cells, x direction
         const size_t& y_cells,  // number of unit cells, y direction
         const size_t& z_cells,  // number of unit cells, z direction
         const double& lattice_constant,
         //const double* const x,  // x domain
         //const double& x_pixels, // number of points in discretized x domain
         //const double* const y,  // y domain
         //const double& y_pixels, // number of points in discretized y domain
         //const double* const z,  // x domain
         //const double& z_pixels, // number of points in discretized z domain
         std::list<double>& q     // lattice point positions (x,y,z), output
         );

   void positions_3D_lattice_zincblende_111(
         const size_t& x_cells,  // number of unit cells, x direction
         const size_t& y_cells,  // number of unit cells, y direction
         const size_t& z_cells,  // number of unit cells, z direction
         const double& lattice_constant,
         //const double* const x,  // x domain
         //const double& x_pixels, // number of points in discretized x domain
         //const double* const y,  // y domain
         //const double& y_pixels, // number of points in discretized y domain
         //const double* const z,  // x domain
         //const double& z_pixels, // number of points in discretized z domain
         const size_t& Z1, const size_t& Z2,// atom types to duplicate
         std::list<size_t>& Z,    // list of atom types, corresponding to 
                                 //  every third element of q
         std::list<double>& q     // lattice point positions (x,y,z), output
         );

   void positions_3D_lattice_zincblende_011(
         const size_t& x_cells,  // number of unit cells, x direction
         const size_t& y_cells,  // number of unit cells, y direction
         const size_t& z_cells,  // number of unit cells, z direction
         const double& lattice_constant,
         //const double* const x,  // x domain
         //const double& x_pixels, // number of points in discretized x domain
         //const double* const y,  // y domain
         //const double& y_pixels, // number of points in discretized y domain
         //const double* const z,  // x domain
         //const double& z_pixels, // number of points in discretized z domain
         const size_t& Z1, const size_t& Z2,// atom types to duplicate
         std::list<size_t>& Z,    // list of atom types, corresponding to 
                                 //  every third element of q
         std::list<double>& q     // lattice point positions (x,y,z), output
         );

   void positions_3D_lattice_zincblende_001(
         const size_t& x_cells,  // number of unit cells, x direction
         const size_t& y_cells,  // number of unit cells, y direction
         const size_t& z_cells,  // number of unit cells, z direction
         const double& lattice_constant,
         //const double* const x,  // x domain
         //const double& x_pixels, // number of points in discretized x domain
         //const double* const y,  // y domain
         //const double& y_pixels, // number of points in discretized y domain
         //const double* const z,  // x domain
         //const double& z_pixels, // number of points in discretized z domain
         const size_t& Z1, const size_t& Z2,// atom types to duplicate
         std::list<size_t>& Z,    // list of atom types, corresponding to 
                                 //  every third element of q
         std::list<double>& q     // lattice point positions (x,y,z), output
         );

   void debug_output_complex_fftw_operand(
         const fftw_complex* const complex_fftw_operand,
         const int& step,
         const ptrdiff_t& local_alloc_size_fftw,
         const ptrdiff_t& Nx_local,
         const ptrdiff_t& Nx,
         const ptrdiff_t& Ny,
         const size_t& resolutionUnit,
         const double& xResolution,
         const double& yResolution,
         const string& outFileName_prefix,
         const int* const psi_mag_strides,
         const int* const psi_mag_displacements,
         const int& mynode,
         const int& rootnode,
         MPI_Comm comm
         );

//   void output_stem_image(
//         const double* const stem_image,
//         const ptrdiff_t& Nx,
//         const ptrdiff_t& Ny,
//         const string& outFileName_prefix
//         );

//   double scherzer_defocus_bfctem_uncorrected(
//         const double& lambda,
//         const double& Cs3
//         );
//
//   double scherzer_alphamax_bfctem_uncorrected(
//         const double& lambda,
//         const double& Cs3
//         );
//
//   double scherzer_approxresolution_bfctem_uncorrected(
//         const double& lambda,
//         const double& Cs3
//         );
//
//   double scherzer_defocus_adfstem_uncorrected(
//         const double& lambda,
//         const double& Cs3
//         );
//
//   double scherzer_alphamax_adfstem_uncorrected(
//         const double& lambda,
//         const double& Cs3
//         );
//
//   double scherzer_approxresolution_adfstem_uncorrected(
//         const double& lambda,
//         const double& Cs3
//         );
//
//   double scherzer_Cs3_bfctem_correctedtoCs5(
//         const double& lambda,
//         const double& Cs5
//         );
//
//   double scherzer_defocus_bfctem_correctedtoCs5(
//         const double& lambda,
//         const double& Cs5
//         );
//
//   double scherzer_alphamax_bfctem_correctedtoCs5(
//         const double& lambda,
//         const double& Cs5
//         );
//
//   double scherzer_approxresolution_bfctem_correctedtoCs5(
//         const double& lambda,
//         const double& Cs5
//         );
//
//   double scherzer_Cs3_adfstem_correctedtoCs5(
//         const double& lambda,
//         const double& Cs5
//         );
//
//   double scherzer_defocus_adfstem_correctedtoCs5(
//         const double& lambda,
//         const double& Cs5
//         );
//
//   double scherzer_alphamax_adfstem_correctedtoCs5(
//         const double& lambda,
//         const double& Cs5
//         );
//
//   double scherzer_approxresolution_adfstem_correctedtoCs5(
//         const double& lambda,
//         const double& Cs5
//         );

   int avg_stddev_of_change_in_array_members(
         const double* const myArray,
         const size_t& number_of_elements,
         double& avg,
         double& stddev
         );

   double sum_detected_intensity(
         const fftw_complex* const psi,
         const double* const kx, const double* const ky,
         const double& k1, const double& k2, 
         const size_t& Nx, const size_t& Ny
         );

   // to_string has been moved to its own header file
   //template < typename T > string to_string( const T& n )
   //{
   //   std::ostringstream ss;
   //   ss << n;

   //   return ss.str();
   //}

}
#endif
