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
                     const double* const kx_split, 
                     const ptrdiff_t& Nx_split,
                     //const double* const kx_joined, 
                     const ptrdiff_t& Nx_joined,
                     const double& xmin,
                     const double* const ky, const ptrdiff_t& Ny,
                     const double& ymin,
                     const int& mynode, 
                     const int& rootnode, 
                     MPI_Comm comm )
      {
         Z = Z_in;
         //cout << "node " << mynode // debug
         //   << ", scatterer_pap constructor, Z_in : "  // debug
         //   << Z_in << endl; // debug

         const double cutoffsq = cutoff * cutoff;
         ptrdiff_t Nx_split_Ny = Nx_split * Ny;
         // fftw requires its operand to be an fftw_complex*
         fftw_complex* projected_atomic_potential_local_split;
         projected_atomic_potential_local_split
            = fftw_alloc_complex( Nx_split_Ny );  

         // Create an fftw plan on projected_atomic_potential_local_split
         fftw_plan pb_c2c_pap;
         pb_c2c_pap = fftw_mpi_plan_dft_2d( 
               Nx_joined, Ny,
               projected_atomic_potential_local_split,
               projected_atomic_potential_local_split,
               comm, FFTW_BACKWARD, FFTW_MEASURE
               );

         // retrieve the scatterer parameters of Z_in from the LUT
         scatterer_param* paramsZ;
         paramsZ = myScattererParamLUT.get_paramPtr( Z_in );

         //cout << "node " << mynode << ", pap constructor, *(paramLUT->a,b,c,d) :" // debug
         //   << *(paramsZ->a) << ", " // debug
         //   << *(paramsZ->b) << ", " // debug
         //   << *(paramsZ->c) << ", " // debug
         //   << *(paramsZ->d) << ", " // debug
         //   << endl; // debug
         //cout << "node " << mynode << " lambda, gamma, ab_inv, kx_split[0], ky[0], xmin, ymin : " // debug
         //              << lambda << ", " <<  gamma << ", " <<  ab_inv << ", " <<  kx_split[0] << ", " <<  ky[0] << ", " <<  xmin << ", " <<  ymin << endl; // debug

         // Use sigmaVz_single_scatterer() to assign values to 
         //  projected_atomic_potential_local_split over the split domain.
         for ( ptrdiff_t i=0; i < Nx_split; ++i)
            for ( ptrdiff_t j=0; j < Ny; ++j)
            {
               projected_atomic_potential_local_split[j + i*Ny][0] = 0.0;
               projected_atomic_potential_local_split[j + i*Ny][1] = 0.0;

               if ( pow( kx_split[i], 2) + pow( ky[j], 2) < cutoffsq ) 
               {
                  // NOTE: using the cylindrical cutoff condition in 
                  //       reciprocal space is producing real space values
                  //       which appear to have 4-fold rotational symmetry
                  //       rather than cylindrical symmetry...
                  //       they look like snowflakes rather than 2-D 
                  //       gaussian distributions
                  sigmaVz_single_scatterer(  // adds to pap with +=
                        lambda, gamma, ab_inv,
                        kx_split[i], ky[j],  
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
            = new double[Nx_split_Ny];

         double* projected_atomic_potential_local_split_im;
         projected_atomic_potential_local_split_im 
            = new double[Nx_split_Ny];

         //cout << "pap constructor, Nx_joined, Ny, sqrtNxNy : " // debug
         //   << Nx_joined << ", " << Ny << ", " << sqrtNxNy << endl;//debug

         for (ptrdiff_t i=0; i < Nx_split_Ny; ++i)
         {
               projected_atomic_potential_local_split_re[i]
                 = projected_atomic_potential_local_split[i][0] / sqrtNxNy;

               projected_atomic_potential_local_split_im[i]
                 = projected_atomic_potential_local_split[i][1] / sqrtNxNy;
         }

         //cout << "pap constructor, projected_atomic_potential_local_split_re[0] : " // debug
         //   << projected_atomic_potential_local_split_re[0] << endl; // debug
         //cout << "pap constructor, projected_atomic_potential_local_split_im[0] : " // debug
         //   << projected_atomic_potential_local_split_im[0] // debug
         //   << endl; // debug


         fftw_free( projected_atomic_potential_local_split );

         // join the split projected atomic potential
         projected_atomic_potential_local_joined_re 
            = new double[Nx_joined * Ny];
         projected_atomic_potential_local_joined_im 
            = new double[Nx_joined * Ny];

         MPI_Allgather(
               projected_atomic_potential_local_split_re, 
               Nx_split_Ny, MPI_DOUBLE,
               projected_atomic_potential_local_joined_re,
               Nx_split_Ny, MPI_DOUBLE,
               comm);

         MPI_Allgather(
               projected_atomic_potential_local_split_im, 
               Nx_split_Ny, MPI_DOUBLE,
               projected_atomic_potential_local_joined_im,
               Nx_split_Ny, MPI_DOUBLE,
               comm);

         // cleanup
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
