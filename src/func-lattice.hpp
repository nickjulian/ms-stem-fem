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
// File: func-lattice.hpp
// Purpose:
//    - prototype functions of domain and atom lattice positions

#ifndef FUNC_LATTICE_HPP
#define FUNC_LATTICE_HPP

#include <cstdlib>
#include <iostream>

using namespace std;
using std::sqrt;
using std::cerr;
using std::endl;

namespace TEM_NS
{
   void domain_3D(
            const int& x_size, const int& y_size, const int& z_size,
            const double& qmin_x, const double& qmax_x, 
            const double& qmin_y, const double& qmax_y,
            const double& qmin_z, const double& qmax_z,
            double *qx, double *qy, double *qz
            );
   // Precondition:
   //    - x_size, y_size, z_size contain values representing the number of 
   //    points in the discretized x, y, and z domains
   //    - qmin_x and qmax_x contain values for the boundary of the x 
   //       domain, same for qmin_y and qmax_y, ...
   //    - qx, qy, qz are pointers to allocated memory which will hold the 
   //       discretized domain points ( allocated with 'new double[Nx]' ...)
   // Postcondition:
   //    - qx, qy, qz have been initialized with values representing the
   //       3-D real-space domain
   
   void positions_3D_lattice_cube(
      const size_t& Nx,   // number of atomic lattice positions, x direction
      const size_t& Ny,   // number of atomic lattice positions, y direction
      const size_t& Nz,   // number of atomic lattice positions, z direction
      const double* const x, // x domain
      const double& x_size,  // number of points in the discretized x domain
      const double* const y, // y domain
      const double& y_size,  // number of points in the discretized y domain
      const double* const z, // y domain
      const double& z_size,  // number of points in the discretized z domain
      double* const q        // lattice point positions (x,y,z), output
      );
   // Precondition
   //    - x, y, z have been initialized with values representing the
   //       3-D real-space domain in row-major array format 
   //       ( x, y, z are the output of domain_3D() )
   //    - x_size, y_size, z_size contain values representing the number of 
   //       points in the discretized x, y, and z domains
   //    - Nx, Ny, Nz have been initialized with the number of atomic 
   //       lattice positions in each direction
   //    - q points to memory which has been allocated 
   //          ('new double[Nx*Ny*Nz]')
   // Postcondition
   //    - q points to a row-major array containing values representing the 
   //       locations of atoms in a cubic lattice

}

#endif
