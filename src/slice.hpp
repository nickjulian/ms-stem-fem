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
// File: slice.hpp
// Purpose:
//    - define functions used to slice atomic systems and perform
//       multislice TEM simulations

#ifndef SLICE_HPP
#define SLICE_HPP

#include <cstdlib>
#include <cmath>
#include <vector>
#include <list>
#include <algorithm>    // sort()
#include <string>


// required for debug_output_complex_fftw_operand()
#include <mpi.h>
#include "scatterer.hpp"
#include "to_string.hpp"
#include "func.hpp"  // debug_output_complex_fftw_operand()
#include "bw_limit.hpp"
//#include "projected_atomic_potential.hpp"
#include <fftw3-mpi.h>
//#include "../include/fftw3-mpi.h"

using std::vector;
using std::list;
using std::string;

namespace TEM_NS
{

class slice
{
   // NOTE: The pointer elements of a slice do not have their destinations 
   //        memory managed within the slice constructor or destructor.
   //       *scatterers, *propagator, and *pap should have their memory 
   //        managed elsewhere.
   public:
      // - pointer to vector of scatterers
      vector< const scatterer* >* scatterers;
      // NOTE: each slice should have a unique vector< scatterer* > to which
      //    scatterers points to, so that the destructor does not find 
      //    itself trying to delete something twice.
      //
      // - lower bound
      double lower_bound;
      // - thickness
      double thickness;
      // - pointer to propagator function data
      //    Each unique slice thickness should induce a propagator function
      //     to be evaluated.
      //    Evaluation of the propagator function must be parallelized, 
      //     since the domain k is split among nodes.
      //    P(k, \Delta z) = \exp( -i \pi \lambda k^{2} \Delta z )
      //fftw_complex* propagator;
      double* propagator_x_re;
      double* propagator_x_im;
      double* propagator_y_re;
      double* propagator_y_im;
      // - projected atomic potential
      //fftw_complex* pap; // NOTE: pap is only used to evaluate exp_i_sigma_v
      // - transmission function 
      fftw_complex* exp_i_sigma_v;
      // - flag to determine if calculated members need to be updated
      bool update_propagator_flag;// = true;
      bool update_transmission_function_flag;// = true;

      int update_propagator(
         const double& lambda,
         const double* const kx,
         const ptrdiff_t& Nx,
         const double* const kydomain,
         const ptrdiff_t& Ny,
         const double& bwcutoff_t //,
         //const int& mynode,
         //const int& rootnode,
         //MPI_Comm comm
         );

      int propagate( 
            // This function requires separate evaluation of the 
            //  propagator.
            const double& lambda,
            const double* const kx,
            const ptrdiff_t& Nx,
            const double* const ky,
            const ptrdiff_t& Ny,
            //const double& sqrtNxNy,
            const double& bwcutoff_t,
            fftw_complex* ft_t_psi  // FT[ t_{n}(x,y) \psi_{n}(x,y) ]
            );

      int update_transmission_function(
            const double& lambda,
            const double& gamma,
            const double& ab_inv,
            const double& bwcutoff_t,
            //const double& bwcutoff_pap,
            //const scatterer_pap_LUT& myScattererPapLUT,
            const double* const kx_local,
            const ptrdiff_t& Nx_local,
            const double* const kx_joined,// reciprocal space x-domain
            const double* const xx_joined, // real space x-domain
            const ptrdiff_t& Nx_joined,
            const double* const ky, // reciprocal space y-domain
            const double* const yy,  // real space y-domain
            const ptrdiff_t& Ny,
            const double& NxNy,
            //const double& sqrtNxNy,
            const unsigned int& input_flag_pap_tif,
            const string& outFileName,//debug
            const size_t& sliceNumber,//debug
            const ptrdiff_t& local_alloc_size_fftw,//debug
            const ptrdiff_t& local_0_start_fftw,
            const int* const psi_mag_strides,
            const int* const psi_mag_displacements,
            const int& mynode,//debug
            const int& rootnode,//debug
            MPI_Comm comm
            );

   // CONSTRUCTORS
      slice( )
      {
         //scatterers = new vector< scatterer* >;
         scatterers = NULL;
         exp_i_sigma_v = NULL;
         //thickness = delta_z;
         lower_bound = 0.0;
         thickness = 0.0;
         propagator_x_re = NULL;
         propagator_x_im = NULL;
         propagator_y_re = NULL;
         propagator_y_im = NULL;
         update_propagator_flag = true;
         update_transmission_function_flag = true;
      }

      // DESTRUCTOR
      //~slice() 
      //{
      //   if ( scatterers != NULL )
      //      delete scatterers; // delete the vector< scatterer* > 
      //                         // pointed to by scatterers
      //   //vector< scatterer* >::iterator scatterersPtr_itr;
      //   //for ( scatterersPtr_itr = scatterers->begin();
      //   //   scatterersPtr_itr != scatterers->end();
      //   //   ++scatterersPtr_itr)
      //   //{
      //   //   delete *scatterersPtr_itr; // delete the vector< scatterer*>
      //   //                              //pointed to by scatterersPtr_itr
      //   //}//This is probably inappropriate, since the vector<scatterer*>
      //      //is not a member of slice, but only the pointer to it is.
      //}
};

// assign atoms to the slices
int assign_atoms_to_slices_z(
      const vector< scatterer >& allScatterers,
      const vector< double >& slice_locations_z,
      const double& xmin, const double& ymin, const double& zmin,
      const double& xperiod, const double& yperiod,
      const double& zperiod,
      list< slice* >& sliceList,
      const int& mynode,
      const int& rootnode,
      MPI_Comm comm
      );

int assign_atoms_to_slices_z_auto(
      const vector< scatterer >& allScatterers,
      const double& xmin, const double& ymin, const double& zmin,
      const double& xperiod, const double& yperiod,
      const double& zperiod,
      const double& min_thickness,
      list< slice* >& sliceList,
      const unsigned int& input_flag_debug,
      const int& mynode,
      const int& rootnode,
      MPI_Comm comm
      );

//int delete_members_of_slice_ctem( slice* mySlice );
int delete_slices_from_list( list< slice* > sliceList );
//int delete_slices_from_list_ctem( list< slice* > sliceList );
// Precondtion : 
//    -  No memory location is pointed to by more than one member of the 
//        slices contained in sliceList; each slice has unique pap, 
//        propagator, and scatterers.
//int delete_members_of_slice_stem( slice* mySlice );
//int delete_slices_from_list_stem( list< slice* > sliceList );
// Precondtion : 
//    -  No memory location is pointed to by more than one member of the 
//        slices contained in sliceList; each slice has unique pap, 
//        exp_i_sigma_v (transmission function), 
//        propagator_x_re, propagator_x_im, propagator_y_re, 
//        propagator_y_im, 
//        and scatterers.

//int evaluate_propagator_ms(
//      // Evaluate the propagator function in reciprocal space
//      // P_{n}(k_{x},k_{y},\Delta z_{n})
//      //    = \exp[-i \pi \lambda (k^{2}_{x} + k^{2}_{y}) \Delta z_{n}]
//      //
//      // Calculating the propagator separately and saving it as a slice 
//      //  member might only be useful for saving cpu time when slice 
//      //  thicknesses aren't unique.
//      //
//      //  Propagator may be factored into x & y components to save memory:
//      // P_{n} = P_{nx} P_{ny}
//      // P_{nx}(k_{x},k_{y},\Delta z_{n})
//      //    = \exp[-i \pi \lambda k^{2}_{x} \Delta z_{n}]
//      // P_{ny}(k_{x},k_{y},\Delta z_{n})
//      //    = \exp[-i \pi \lambda k^{2}_{y} \Delta z_{n}]
//      //
//      const double& lambda,
//      const double* const kxdomain,
//      const ptrdiff_t& Nx,
//      const double* const kydomain,
//      const ptrdiff_t& Ny,
//      const double& delta_z,
//      fftw_complex* propagator_x,
//      fftw_complex* propagator_y
//      );

//int propagate_ms( 
//      // This function does not require separate evaluation of the 
//      //  propagator, but instead incorporates the propagator calculation
//      //  into the propagation function.
//      // This might reduce memory usage, but increase cpu usage if any of
//      //  the propagators are duplicates (identical \Delta z).
//      // Calculating the propagator separately and saving it as a slice 
//      //  member might only be useful for saving cpu time when slice 
//      //  thicknesses aren't unique.
//      const double& lambda,
//      const double* const kxdomain,
//      const ptrdiff_t& Nx,
//      const double* const kydomain,
//      const ptrdiff_t& Ny,
//      //const double& sqrtNxNy,
//      const double& bwcutoff_t,
//      const double& delta_z,
//      fftw_complex* ft_t_psi  // FT[ t_{n}(x,y) \psi_{n}(x,y) ]
//      );

//int propagate_ms( 
//      // This function does not require separate evaluation of the 
//      //  propagator, but instead incorporates the propagator calculation
//      //  into the propagation function.
//      // This might reduce memory usage, but increase cpu usage if any of 
//      //  the propagators are duplicates (identical \Delta z).
//      // Calculating the propagator separately and saving it as a slice 
//      //  member might only be useful for saving cpu time when slice 
//      //  thicknesses aren't unique.
//      const double& lambda,
//      const double* const kxdomain,
//      const ptrdiff_t& Nx,
//      const double* const kydomain,
//      const ptrdiff_t& Ny,
//      //const double& sqrtNxNy,
//      const double& bwcutoff_t,
//      const double& delta_z,
//      double* propagator_x_re,
//      double* propagator_x_im,
//      double* propagator_y_re,
//      double* propagator_y_im,
//      fftw_complex* ft_t_psi  // FT[ t_{n}(x,y) \psi_{n}(x,y) ]
//      );

//int transmission_function_ms(
//      const ptrdiff_t& Nx,
//      const ptrdiff_t& Ny,
//      const double& sqrtNxNy,
//      const fftw_complex* const pap,
//      fftw_complex* exp_i_sigma_v
//      );

}
#endif
