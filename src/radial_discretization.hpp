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
// File: radial_discretization.hpp
// Purpose:
//    - define the radial_coordinate class having members : 
//       double r, vector<double[2]> cos_sin

#ifndef RADIAL_DISCRETIZATION_HPP
#define RADIAL_DISCRETIZATION_HPP

#ifndef PI
#define PI 3.14159265369
#endif

#include <cstdlib>
#include <cmath>
#include <vector>
#include <list>
#include <iostream>


// required for debug_output_complex_fftw_operand()
//#include <mpi.h>

using std::cerr;
using std::cout;
using std::endl;
using std::vector;
using std::list;

namespace TEM_NS
{

class cylindrical_coordinates
{
   public:
      double radius;
      vector< double > vector_of_angles;

   // CONSTRUCTORS
      cylindrical_coordinates( double rr, double delta_s )
      {
         radius = rr;

         if ( delta_s > 2 * PI * radius )
         {
            cerr << "Error: cylindrical_coordinates( rr: "
               << rr << ", delta_s: "
               << delta_s 
               << " ) "
               << "instantiation failed; delta_s > 2*PI*rr" << endl;
            return;
         }

         double delta_ss = 1.0;
         double number_of_delta_ss;
         double delta_phi;// = delta_ss/radius;
         number_of_delta_ss = floor( 2 * PI * rr / delta_s );
         if ( radius == 0.0 )
         {
            //delta_ss = 1.0;
            delta_phi = 2 * PI ;
         }
         else
         {
            delta_ss = 2 * PI * rr / number_of_delta_ss;
            delta_phi = delta_ss/radius;
         }


         for( double phi = 0.0; //delta_phi; 
               phi < 2 * PI ; 
               phi += delta_phi )
         {
            //cerr << "phi, 2*PI : " << phi << ", " << 2*PI << endl; // debug
            vector_of_angles.push_back( phi );
         }
      }

      // DESTRUCTOR
};

class radial_discretization : public list< cylindrical_coordinates > 
{
   public:

   // CONSTRUCTORS
      radial_discretization(  )
      {
         return;
      }

      radial_discretization( double r0, double r1, double delta_r, 
                           double delta_s )
      {
         // [r0, r1] : interval of radii to discretize
         // delta_r : approximate increment between subsequent radii
         // delta_s : arclength nearest to that which the circle will be
         //             split by 
         double delta_rr = floor( (r1 - r0)/delta_r ); 

         if ( delta_rr < 1.0 ) 
         {
            cerr << "Error : radial_discretization( r0: " 
               << r0 << ", r1: "
               << r1 << ", delta_r: "
               << delta_r << ", delta_s: "
               << delta_s << " ) "
               << "instantiation failed; delta_r > r1 - r0 " << endl;
            return;
         }

         //   double delta_ss;
         //   double number_of_delta_ss;
         //   number_of_delta_ss = floor( 2 * PI * rr / delta_s );
         //   delta_ss = 2 * PI * rr / number_of_delta_ss;
         //double phi;
         //if( delta_ss > 2 * PI * (r0 + delta_rr) )
         //{
         //   cerr << "Error : radial_discretization( r0: " 
         //      << r0 << ", r1: "
         //      << r1 << ", delta_rr: "
         //      << delta_rr << ", delta_ss: "
         //      << delta_ss << " ) "
         //      << "instantiation failed; delta_ss > 2 * PI * (r0 + delta_rr) " << endl;
         //   return;
         //}

         for( double rr = r0; rr <= r1; rr += delta_r )
         {
            if ( rr == 0.0 )
            {
               push_back( cylindrical_coordinates( 0.0, 0.0 ) );
               continue;
            }

            push_back( cylindrical_coordinates( rr, delta_s) );

            //for( double ss = 0.0; ss < 2 * PI * rr; ss += delta_ss )
            //{
            //   phi = ss / rr ;
            //   push_back( radial_coordinate( rr, phi ) );
            //}
         }
      }
};

} // TEM_NS
#endif
