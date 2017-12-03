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
// File: tem_io.hpp
// Purpose:
//

#ifndef TEM_IO_HPP
#define TEM_IO_HPP

#include <iostream>  // cout, cin, cerr
#include <iomanip>   // setw, setprecision
#include <fstream>   // ifstream, ofstream
#include <cstdlib>   // EXIT_SUCCESS, EXIT_FAILURE
#include <string>    
#include <vector>      
#include <list>      
#include <limits>    // numeric_limits<>
#include <algorithm> // transform()
//#include <locale>    // tolower() BINARY OPERATOR --> cannot transform()
#include <cctype>    // tolower() UNARY OPERATOR

#include <iomanip>   // setw(), setprecision()

#include "scatterer.hpp"
#include "func-lattice.hpp"
#include "func.hpp"

using std::ifstream;
using std::string;
using std::cout;
using std::cerr;
using std::endl;
using std::transform;

namespace TEM_NS
{
   int read_position_lammps_file( 
   // Precondition: 
   // - input is expected to be the same format as a lammps position
   //  input file, with an Atoms section and an enumerated list of 
   //  atoms types and coordinates
         const string& filename, 
         //scatterer_param_LUT& lut,
         //const ptrdiff_t& Nx, const ptrdiff_t& Ny, const ptrdiff_t& Nz,
         //std::vector<scatterer>& myScatterers,
         double*& qq_contig,
         unsigned int*& Z_contig,
         unsigned int& total_population,
         double& xlo, double& ylo, double& zlo, 
         double& xperiod, double& yperiod, double& zperiod,
         const unsigned int& input_flag_debug,
         const int& mynode,
         const int& rootnode,
         MPI_Comm comm
         );

   int read_position_xyz_file( 
   // Precondition: 
   // - input is expected to be the same format as a lammps position
   //  input file, with an Atoms section and an enumerated list of 
   //  atoms types and coordinates
         const string& filename, 
         //scatterer_param_LUT& lut,
         //const ptrdiff_t& Nx, const ptrdiff_t& Ny, const ptrdiff_t& Nz,
         //std::vector<scatterer>& myScatterers,
         double*& qq_contig,
         unsigned int*& Z_contig,
         unsigned int& total_population,
         double& xlo, double& ylo, double& zlo, 
         double& xperiod, double& yperiod, double& zperiod,
         const unsigned int& input_flag_debug,
         const int& mynode,
         const int& rootnode,
         MPI_Comm comm
         );

   //   TODO: implement function to read microscope parameters from a file
//   int read_microscope_file(
//         const string& filename,
//         int& Nx_int, int& Ny_int, int& Nz_int,
//         double& VV,
//         double& Cs3,
//         double& Cs5,
//         double& condenser_illumination_angle,
//         double& alphamax
//         );

   size_t atom_element_abbrev_to_Z( const string& element_name);

   int create_position_lammps_file_011diamond(
         const string& filename,
         const double& xlo, const double& ylo, const double& zlo, 
         const size_t& unit_cells_x,
         const size_t& unit_cells_y,
         const size_t& unit_cells_z,
         const int& Nx, const int& Ny, const int& Nz,
         const size_t& Z,
         const double& lattice_constant
         );
         
   int create_position_lammps_file_001diamond(
         const string& filename,
         const double& xlo, const double& ylo, const double& zlo, 
         const size_t& unit_cells_x,
         const size_t& unit_cells_y,
         const size_t& unit_cells_z,
         const int& Nx, const int& Ny, const int& Nz,
         const size_t& Z,
         const double& lattice_constant
         );
         
   int create_position_lammps_file_111diamond(
         const string& filename,
         const double& xlo, const double& ylo, const double& zlo, 
         const size_t& unit_cells_x,
         const size_t& unit_cells_y,
         const size_t& unit_cells_z,
         const int& Nx, const int& Ny, const int& Nz,
         const size_t& Z,
         const double& lattice_constant
         );
         
   int create_position_lammps_file_111zincblende(
         const string& filename,
         const double& xlo, const double& ylo, const double& zlo, 
         const size_t& unit_cells_x,
         const size_t& unit_cells_y,
         const size_t& unit_cells_z,
         const int& Nx, const int& Ny, const int& Nz,
         const size_t& Z1,
         const size_t& Z2,
         const double& lattice_constant
         );
         
   int create_position_lammps_file_011zincblende(
         const string& filename,
         const double& xlo, const double& ylo, const double& zlo, 
         const size_t& unit_cells_x,
         const size_t& unit_cells_y,
         const size_t& unit_cells_z,
         const int& Nx, const int& Ny, const int& Nz,
         const size_t& Z1,
         const size_t& Z2,
         const double& lattice_constant
         );
         
   int create_position_lammps_file_001zincblende(
         const string& filename,
         const double& xlo, const double& ylo, const double& zlo, 
         const size_t& unit_cells_x,
         const size_t& unit_cells_y,
         const size_t& unit_cells_z,
         const int& Nx, const int& Ny, const int& Nz,
         const size_t& Z1,
         const size_t& Z2,
         const double& lattice_constant
         );
         
}  // TEM_NS

#endif
