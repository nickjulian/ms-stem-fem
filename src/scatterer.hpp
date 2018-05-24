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
// File: scatterer.hpp
// Purpose:
//
// Use: 
//    1. instantiate a scatterer_param_LUT
//    2. instantiate your scatterers using the scatterer_param_LUT from 1, 
//       its position q[3], and its species (Z) 
// or:
//    1. instantiate a scatterer_param_LUT
//    2. For each given scattering position q[3] and corresponding species 
//       (Z), increment the scattering factor by fetching parameter 
//       addresses from the scatterer_param_LUT.get_paramPtr( Z )->a ,b,c,d 

#ifndef SCATTERER_HPP
#define SCATTERER_HPP

#include <cstdlib>
#include <fstream>
#include <cassert>
#include <complex>
#include <fftw3-mpi.h>
//#include "fftw3-mpi.h"   // ptrdiff_t
//#include "scatterer_param_LUT.hpp"
#include "scatterer_pap_LUT.hpp"

using std::ifstream;
using std::cerr;
using std::cout;
using std::endl;
using std::complex;

namespace TEM_NS
{
class scatterer
{
   public:
      double* q;
      unsigned int Z;

      double* pap_to_translate_re;  // this projected atomic potential
      double* pap_to_translate_im;  // is centered at (0,0) and is shared 
                                    // among all atoms of type Z 

      // CONSTRUCTORS
      scatterer( double* const q_in,  const unsigned int& Z_in,
            const scatterer_pap_LUT& pap_lut
            )
      { // use scatterer( q, Z ) : pap( <address of pap from LUT> )
         q = q_in;
         Z = Z_in;
         
         // Select a pap from the scatterer_pap_LUT using Z and assign 
         //  pap_to_translate_{re,im} to point to it.
         std::list< scatterer_pap* >::const_iterator pap_lut_itr;
         for ( pap_lut_itr = pap_lut.pap_list.begin();
               pap_lut_itr != pap_lut.pap_list.end();
               ++pap_lut_itr)
         {
            if ( Z == (*pap_lut_itr)->Z )
            {
               pap_to_translate_re 
                  = (*pap_lut_itr)
                        ->projected_atomic_potential_local_joined_re;
               pap_to_translate_im
                  = (*pap_lut_itr)
                        ->projected_atomic_potential_local_joined_im;

               //cout << "scatterer constructor; Zs matched : Z == " // debug
               //   << Z << ", (*pap_lut_itr)-Z == " << (*pap_lut_itr)->Z // debug
               //   << endl; // debug
               break;
            }
            //else // debug
            //   cout << "scatterer constructor; Zs didn't match : Z == " // debug
            //      << Z << ", (*pap_lut_itr)-Z == " << (*pap_lut_itr)->Z // debug
            //      << endl; // debug
         }

         if ( pap_lut_itr == pap_lut.pap_list.end() ) 
         {
            cout << "Error, scatterer constructor could not match its Z"
               << " to any of those in the scatterer_pap_LUT." << endl;
            cout << "Z : " << Z << ", Z_in : " << Z_in << endl;
         }
         //std::cout <<  "scatterer constructor: pap_to_translate_re + i*pap_to_translate_im[0] : " // debug
         //     << pap_to_translate_re[0] << " + i*"  // debug
         //     << pap_to_translate_im[0] << endl; // debug
         //std::cout <<  "scatterer constructor: pap_to_translate_re : " << pap_to_translate_re  // debug
         //     << ", pap_to_translate_im : " << pap_to_translate_im // debug
         //     << "pap_lut_re[0] : " 
         //     << (*pap_lut_itr)->projected_atomic_potential_local_joined_re[0]
         //     << "pap_lut_im[0] : " 
         //     << (*pap_lut_itr)->projected_atomic_potential_local_joined_im[0]
         //     << endl; // debug
      }

      scatterer( double* const q_in,  const unsigned int& Z_in,
            const scatterer_pap_LUT* pap_lut
            )
      { // use scatterer( q, Z ) : pap( <address of pap from LUT> )
         q = q_in;
         Z = Z_in;
         
         // Select a pap from the scatterer_pap_LUT using Z and assign 
         //  pap_to_translate_{re,im} to point to it.
         std::list< scatterer_pap* >::const_iterator pap_lut_itr;
         for ( pap_lut_itr = pap_lut->pap_list.begin();
               pap_lut_itr != pap_lut->pap_list.end();
               ++pap_lut_itr)
         {
            if ( Z == (*pap_lut_itr)->Z )
            {
               pap_to_translate_re 
                  = (*pap_lut_itr)
                        ->projected_atomic_potential_local_joined_re;
               pap_to_translate_im
                  = (*pap_lut_itr)
                        ->projected_atomic_potential_local_joined_im;

               break;
            }
         }

         if ( pap_lut_itr == pap_lut->pap_list.end() ) 
         {
            cout << "Error, scatterer constructor could not match its Z"
               << " to any of those in the scatterer_pap_LUT." << endl;
            cout << "Z : " << Z << ", Z_in : " << Z_in << endl;
         }
      }
};

   // TODO: Idea: Overload scatterer_sorting_lt and sorting_metric() to 
   //        allow sorting along directions.
   //       Problem: sorting along directions requires special care for the
   //        periodicities of the domain.

   inline double scatterer_sorting_metric_x( const scatterer& aa) 
   { 
      return aa.q[0];
   }
   
   inline double scatterer_sorting_metric_y( const scatterer& aa) 
   { 
      return aa.q[1];
   }
   
   inline double scatterer_sorting_metric_z( const scatterer* const aa) 
   { 
      return aa->q[2];
   }

   inline bool scatterer_sorting_lt_x( 
         const scatterer& lhs, 
         const scatterer& rhs
         )
   {
      return 
         scatterer_sorting_metric_x(lhs) < scatterer_sorting_metric_x(rhs);
   }

   inline bool scatterer_sorting_lt_y( 
         const scatterer& lhs, 
         const scatterer& rhs
         )
   {
      return 
         scatterer_sorting_metric_y(lhs) < scatterer_sorting_metric_y(rhs);
   }
   
   inline bool scatterer_sorting_lt_z( 
         const scatterer* const lhs, 
         const scatterer* const rhs
         )
   {
      return 
         scatterer_sorting_metric_z(lhs) < scatterer_sorting_metric_z(rhs);
   }

// TODO : implement a function which rotates myScatterers and transmits
//       atoms across periodic boundaries

}  // TEM_NS

#endif
