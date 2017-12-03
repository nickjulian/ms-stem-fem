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
// File: func-lattice.cpp
// Purpose:
//    - implement function domain and atom lattice positions

#ifndef FUNC_LATTICE_CPP
#define FUNC_LATTICE_CPP

#include "func-lattice.hpp"

using namespace std;
using std::sqrt;
using std::cerr;
using std::endl;

void TEM_NS::domain_3D(
         const int& x_size, const int& y_size, const int& z_size,
         const double& qmin_x, const double& qmax_x, 
         const double& qmin_y, const double& qmax_y,
         const double& qmin_z, const double& qmax_z,
         double *qx, double *qy, double *qz
         )
{
   const double dx = (qmax_x - qmin_x)/(x_size-1);
   const double dy = (qmax_y - qmin_y)/(y_size-1);
   const double dz = (qmax_z - qmin_z)/(z_size-1);
   //cout << "dx, dy : " << dx << ", " << dy << endl; // debug

   for(int i=0; i<x_size; i++) qx[i] = qmin_x + ( i * dx) ;
   for(int i=0; i<y_size; i++) qy[i] = qmin_y + ( i * dy) ;
   for(int i=0; i<z_size; i++) qz[i] = qmin_z + ( i * dz) ;

   return ;
}

void TEM_NS::positions_3D_lattice_cube(
   const size_t& Nx,    // number of atomic lattice positions, x direction
   const size_t& Ny,    // number of atomic lattice positions, y direction
   const size_t& Nz,    // number of atomic lattice positions, z direction
   const double* const x,  // x domain
   const double& x_size,   // number of points in the discretized x domain
   const double* const y,  // y domain
   const double& y_size,   // number of points in the discretized y domain
   const double* const z,  // y domain
   const double& z_size,   // number of points in the discretized z domain
   double* const q         // lattice point positions (x,y,z), output
   )
{
   for( size_t i = 0; i < Nx; i++)
      for( size_t j = 0; j < Ny; j++)
         for( size_t k = 0; k < Nz; k++)
         {
            q[2*((i * Ny + j)*Nz + k) + 0 ] = x[ size_t( i * x_size / Nx )];
            q[2*((i * Ny + j)*Nz + k) + 1 ] = y[ size_t( j * y_size / Ny )];
            q[2*((i * Ny + j)*Nz + k) + 2 ] = z[ size_t( k * z_size / Nz )];
         }

   return;
}

#endif
