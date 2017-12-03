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
// File: scatterer_pap_LUT.hpp
// Purpose:
//
// Use: 
//    1. instantiate a scatterer_param_LUT
//    2. instantiate your scatterers using the scatterer_param_LUT from 1, 
//       its position q[3], and its species (Z) 

#ifndef SCATTERER_PAP_LUT_HPP
#define SCATTERER_PAP_LUT_HPP

#include <cstdlib>
#include <fstream>
#include <cassert>
#include <list>

#include "scatterer_pap.hpp"

using std::ifstream;
using std::cerr;
using std::endl;

namespace TEM_NS
{

class scatterer_pap_LUT
{
   // LUT to be instantiated once and have all scatterers point to the 
   //    locations of their parameter values in the LUT
   //
   // I would rather this be a derived class from std::list<>, but when
   //  attempted I wasn't allowed to instantiate a 
   //  scatterer_pap_LUT::iterator .
   // The only purpose to having this as its own class is to specify the
   //  custom constructor, which could be implemented as a stand-alone 
   //  function if this were not a class.
   public:
      std::list< scatterer_pap* > pap_list;
      // CONSTRUCTOR
      scatterer_pap_LUT( const std::vector<unsigned int>& Z_list,
                     const double& lambda, const double& gamma, 
                     const double& ab_inv, 
                     const double& cutoff,
                     const double* const kx_split,
                     const ptrdiff_t& Nx_split,
                     //const double* kx_joined,
                     const ptrdiff_t& Nx_joined,
                     const double& xmin,
                     const double* const ky, const ptrdiff_t& Ny,
                     const double& ymin,
                     const int& mynode,
                     const int& rootnode,
                     MPI_Comm comm )
      {
         // The constructor will call the scatterer_pap constructor for 
         //  each member of Z_list.
         
         // Z_list shouldn't contain duplicates

         // Instantiate a scatterer_param_LUT to be used by pap 
         //  constructor.
         scatterer_param_LUT myScattererParamLUT;

         // Iterate over the reduced list of atomic species.
         for ( std::vector<unsigned int>::const_iterator 
                  Z_list_itr = Z_list.begin(); 
                  Z_list_itr != Z_list.end(); 
                  ++Z_list_itr )
         {
            // call the constructor of a pap using myScattererParamLUT,
            //  pushing the scatterer_pap onto the list<scatterer_pap*>

            pap_list.push_back(
                     new scatterer_pap( 
                              *Z_list_itr,
                              myScattererParamLUT,
                              lambda, gamma, ab_inv, cutoff,
                              kx_split, Nx_split,
                              //kx_joined,
                              Nx_joined,
                              xmin,
                              ky, Ny,
                              ymin,
                              mynode, rootnode, comm )
                  );
            //cout << "LUT constructor: pap_lut_re[0] : " // debug
            //   << pap_list.back()->projected_atomic_potential_local_joined_re[0] // debug
            //<< "LUT constructor: pap_lut_im[0] : " // debug
            //   << pap_list.back()->projected_atomic_potential_local_joined_im[0] << endl; // debug
         }
      }

      // DESTRUCTOR
      ~scatterer_pap_LUT()
      {
         for ( std::list<scatterer_pap*>::iterator 
               pap_list_itr = pap_list.begin(); 
               pap_list_itr != pap_list.end(); 
               ++pap_list_itr
               )
         {
            delete *pap_list_itr;
         }

         while ( ! pap_list.empty() )
         {
            pap_list.pop_back();
         }
      }
};

}  // TEM_NS

#endif
