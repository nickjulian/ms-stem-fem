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
   struct input_flags
   {
      unsigned int m = 0;// microscope parameters
      //unsigned int pf = 0;// simulation parameter file given
      unsigned int o = 0;// output file name prefix
      unsigned int a = 0;// atom position and species input file
      unsigned int defocus = 0;   // defocus
      unsigned int spread = 0;  // defocus spread 
      unsigned int adfstem_corrected = 0;        // 
      unsigned int adfstem_uncorrected = 0;      // 
      unsigned int bfctem_corrected = 0;         // 
      unsigned int bfctem_uncorrected = 0;       // 
      unsigned int fem = 0;
      unsigned int gt17 = 0;
      unsigned int d1 = 0;
      unsigned int d2 = 0;
      unsigned int d3 = 0;
      unsigned int d4 = 0;
      unsigned int scherzer_defocus = 0;
      unsigned int scherzer_alphamax = 0;
      unsigned int scherzer_cs3 = 0;
      unsigned int cs3 = 0;
      unsigned int cs5 = 0;
      unsigned int alpha_max = 0;
      unsigned int aberration_correction = 0;
      unsigned int raster_spacing = 0;
      unsigned int pap_tif = 0;
      unsigned int dupe = 0;
      unsigned int image_output = 0;
      unsigned int netcdf_images = 0;
      unsigned int netcdf_variance = 0;
      unsigned int nx = 0;
      unsigned int ny = 0;
      unsigned int microscope_voltage = 0;
      unsigned int debug = 0;
   };

   int read_parameter_file(
      // sufficiency and conflicts of the input should be checked by caller
         const string& parameter_filename, 
         string& model_file,
         input_flags& flags,
         string& output_prefix,
         ptrdiff_t& Nx,
         ptrdiff_t& Ny,
         double& VV,
         double& defocus,
         double& alpha_max,
         double& defocus_spread,
         double& condenser_illumination_angle,
         double& Cs3,
         double& Cs5,
         double& raster_spacing,
         double& azimuthal_binning_size_factor,
         double& minSliceThickness,
         const int& mynode,
         const int& rootnode,
         MPI_Comm comm
      );

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
