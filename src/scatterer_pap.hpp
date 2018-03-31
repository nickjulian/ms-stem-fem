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
// File: scatterer_pap.hpp
// Purpose:

#ifndef SCATTERER_PAP_HPP
#define SCATTERER_PAP_HPP

#include <cstdlib>
#include <fstream>
#include <cassert>

#include <fftw3-mpi.h>
//#include "../include/fftw3-mpi.h"   // contains ptrdiff_t, fftw_complex

#include "projected_atomic_potential.hpp" // sigmaVz_single_scatterer()
#include "scatterer_param_LUT.hpp"

using std::ifstream;
using std::cerr;
using std::endl;

namespace TEM_NS
{

class scatterer_pap
{
   // This is a struct because scatterer_pap_LUT will instantiate an array
   //  of these, ..., can one instantiate an array : double* params[]
   public:
      unsigned int Z;
      double* projected_atomic_potential_local_joined_re;
      double* projected_atomic_potential_local_joined_im;
   // LUT to be instantiated once and have all scatterers point to the 
   //    locations of their parameter values in the LUT

      // CONSTRUCTOR
      scatterer_pap( const unsigned int& Z_in,
                     scatterer_param_LUT& myScattererParamLUT,
                     const double& lambda, const double& gamma, 
                     const double& ab_inv,
                     const double& cutoff,
                     const double* const kx_local, 
                     const ptrdiff_t& Nx_local,
                     //const double* const kx_joined, 
                     const ptrdiff_t& Nx_joined,
                     const double& xmin,
                     const double* const ky, const ptrdiff_t& Ny,
                     const double& ymin,
                     const ptrdiff_t& local_alloc_size_fftw,
                     const int* const pap_strides,
                     const int* const pap_displacements,
                     const int& mynode, 
                     const int& rootnode, 
                     MPI_Comm comm )
      {
         Z = Z_in;

         const double cutoffsq = cutoff * cutoff;
         // fftw requires its operand to be an fftw_complex*
         fftw_complex* projected_atomic_potential_local_split;
         projected_atomic_potential_local_split
            = fftw_alloc_complex( local_alloc_size_fftw );  

         // Create an fftw plan on projected_atomic_potential_local_split
         fftw_plan pb_c2c_pap;
         pb_c2c_pap = fftw_mpi_plan_dft_2d( 
               Nx_joined, Ny,
               projected_atomic_potential_local_split,
               projected_atomic_potential_local_split,
               //comm, FFTW_BACKWARD, FFTW_ESTIMATE
               //comm, FFTW_BACKWARD, FFTW_PATIENT
               //comm, FFTW_BACKWARD, FFTW_EXHAUSTIVE
               comm, FFTW_BACKWARD, FFTW_MEASURE
               );

         // retrieve the scatterer parameters of Z_in from the LUT
         scatterer_param* paramsZ;
         paramsZ = myScattererParamLUT.get_paramPtr( Z_in );


         // Use sigmaVz_single_scatterer() to assign values to 
         //  projected_atomic_potential_local_split over the split domain.
         for ( ptrdiff_t i=0; i < Nx_local; ++i)
            for ( ptrdiff_t j=0; j < Ny; ++j)
            {
               projected_atomic_potential_local_split[j + i*Ny][0] = 0.0;
               projected_atomic_potential_local_split[j + i*Ny][1] = 0.0;

               if ( pow( kx_local[i], 2) + pow( ky[j], 2) < cutoffsq ) 
               {
                  // NOTE: using the cylindrical cutoff condition in 
                  //       reciprocal space is producing real space values
                  //       which appear to have 4-fold rotational symmetry
                  //       rather than cylindrical symmetry...
                  //       they look like snowflakes rather than 2-D 
                  //       gaussian distributions
                  sigmaVz_single_scatterer(  // adds to pap with +=
                        lambda, gamma, ab_inv,
                        kx_local[i], ky[j],  
                        xmin, ymin,   // realspace position of pap center
                        // 0.0, 0.0, // Put it at (0,0)? No. Where is (0,0)?
                        paramsZ->a,
                        paramsZ->b,
                        paramsZ->c,
                        paramsZ->d,
                        projected_atomic_potential_local_split[j + i*Ny][0],
                        projected_atomic_potential_local_split[j + i*Ny][1]
                        );

                  if ( (1.0/0.0) == // debug
                        projected_atomic_potential_local_split[j + i*Ny][0] // debug
                        ||  // debug
                        (1.0/0.0) == // debug
                        projected_atomic_potential_local_split[j + i*Ny][1]
                     ) // debug
                     cout << "Error, pap constructor, projected_atomic_potential_local_split[" << j + i*Ny << "][1] : " // debug
                     << projected_atomic_potential_local_split[j + i*Ny][1] // debug
                     << ", projected_atomic_potential_local_split[" << j + i*Ny << "][0] : " // debug
                     << projected_atomic_potential_local_split[j + i*Ny][0] // debug
                     << endl; // debug

               }
               //else
               //{
               //   projected_atomic_potential_local_split[j + i*Ny][0] = 0.0;
               //   projected_atomic_potential_local_split[j + i*Ny][1] = 0.0;
               //}
            }

         // Inverse fourier transform projected_atomic_potential_local_split
         fftw_execute( pb_c2c_pap );
         fftw_destroy_plan( pb_c2c_pap );

         // Scale it by 1/sqrt(Nx_joined * Ny) while splitting complex
         //  numbers into two doubles as required for the Allgather step.
         double sqrtNxNy = sqrt(Nx_joined * Ny);

         double* projected_atomic_potential_local_split_re;
         projected_atomic_potential_local_split_re 
            = new double[ local_alloc_size_fftw ];

         double* projected_atomic_potential_local_split_im;
         projected_atomic_potential_local_split_im 
            = new double[ local_alloc_size_fftw ];


         for (ptrdiff_t i=0; i < local_alloc_size_fftw ; ++i)
         {
               projected_atomic_potential_local_split_re[i]
                 = projected_atomic_potential_local_split[i][0] / sqrtNxNy;

               projected_atomic_potential_local_split_im[i]
                 = projected_atomic_potential_local_split[i][1] / sqrtNxNy;
         }


         fftw_free( projected_atomic_potential_local_split );

         // join the split projected atomic potential
         projected_atomic_potential_local_joined_re 
            = new double[Nx_joined * Ny];
         projected_atomic_potential_local_joined_im 
            = new double[Nx_joined * Ny];

         //for (size_t ii=0; ii < 3; ++ii)
         //{
         //   cout << "node " << mynode
         //      << " Nx_local: " << Nx_local << endl;
         //   cout << "node " << mynode
         //      << " pap_strides[" << ii << "]: "
         //      << pap_strides[ii ] << endl;
         //   cout << "node " << mynode
         //      << " pap_displacements[" << ii << "]: "
         //      << pap_displacements[ii ] << endl;
         //}

         //unsigned int Nx_local_uint = (unsigned int) Nx_local;
         //MPI_Allgather(
         //      &Nx_local_uint,
         //      1,
         //      MPI_UNSIGNED,
         //      Nx_local_all,
         //      totalnodes,
         //      MPI_UNSIGNED,
         //      comm
         //      );
         MPI_Allgatherv(
               projected_atomic_potential_local_split_re, // sendbuf
               Nx_local * Ny, //local_alloc_size_fftw, // sendcount
               MPI_DOUBLE, // sendtype
               projected_atomic_potential_local_joined_re,//*recvbuf
               pap_strides, // recvcounts[]
               pap_displacements, // displs[]
               MPI_DOUBLE, // recvtype
               comm);

         MPI_Allgatherv(
               projected_atomic_potential_local_split_im, 
               Nx_local * Ny, //local_alloc_size_fftw, 
               MPI_DOUBLE,
               projected_atomic_potential_local_joined_im,
               pap_strides, 
               pap_displacements, 
               MPI_DOUBLE,
               comm);

         
         delete[] projected_atomic_potential_local_split_im;
         delete[] projected_atomic_potential_local_split_re;
      }

      // DESTRUCTOR
      ~scatterer_pap()
      {
         delete[] projected_atomic_potential_local_joined_re;
         delete[] projected_atomic_potential_local_joined_im;
      }
};

}  // TEM_NS
#endif
