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
// File: tem.cpp
// Purpose: 
//    To control execution of a multislice TEM simulation.
//    setup outline: 
//       x read positions and species
//       x broadcast positions and species to all nodes
//       - reduce species list to remove duplicates
//       x for each unique species, evaluate the pap centered at (0,0)
//       x Fourier transform each pap into realspace
//       - for each scatterer, add a translated copy of the pap of the 
//          corresponding species to the total pap
//          - original pap evaluated how?
//          - translation evaluated how?
//       - evaluate the transmission function (a function of the pap)
//
//    copies of the paps are made from a LUT in realspace


#include <iostream>     // provides cout, cin, and cerr
#include <fstream>
#include <string>
#include <cstdlib>         // provides EXIT_SUCCESS, EXIT_FAILURE 
#include <vector>
//#include <complex>
#include <list>
#include <algorithm>    // std::unique, std::distance
#include <iomanip> // provides setwidth, setprecision

#include <mpi.h>
#include <fftw3-mpi.h>
//#include "../include/fftw3-mpi.h"   // local fftw version

#include "tem_io.hpp"   
   // read_position_lammps_file()
   // read_position_xyz_file()
#include "scherzer.hpp"
   // functions which evaluate scherzer conditions

#include "adfstem.hpp"
//#include "bfctem.hpp"    // not yet implemented
#include "scatterer.hpp"   // atom struct with fitting parameters
#include "scatterer_pap_LUT.hpp"
#include "slice.hpp"            

using std::cout;
using std::endl;
using std::string;
using std::istringstream;
using std::setw;
using std::setprecision;
using namespace TEM_NS;

//int threads_ok; // global by request of fftw.org/doc/ ... section 6.11

// TODO: separate scherzer conditionals for defocus, alphamax, and Cs3
#define PRINT_USAGE cout << "Usage: " << argv[0] << " <options>" << endl << "OPTIONS : " << endl << "   -m <Nx> <Ny> <VV>" << endl << "   [--scherzer_defocus] calculate and use scherzer focus conditions" << endl << "   [--scherzer_alphamax] calculate and use scherzer focus conditions" << endl << "   [--scherzer_cs3] calculate and use scherzer Cs3 conditions" << endl << "   [--defocus <defocus>] (if not using --scherzer)" << endl << "   [--alphamax <alpha max>] (if not using --scherzer)" << endl << "   [--spread <defocus spread> <condenser_illumination_angle>] (only if using aberration correction)" << endl << "   [--cs3 <third order spherical aberration>] (ignored if using --scherzer with aberration correction)" << endl << "   [--cs5 <fifth order spherical aberration>] (only if using aberration correction)" << endl << "   [--rasterspacing <raster_spacing>] (Angstroms; default is 1.5 for STEM or 10 for STEM-FEM)" << endl << "   [--dupe <dupe_x> <dupe_y> <dupe_z>] (Periodically instantiate the given scatterers dupe_x, dupe_y, and dupe_z times in respective directions; fem only)" << endl << "   -a <lammps or xyz style position input file>" << endl << "   -o <output file prefix>" << endl << "   --paptif (output images of projected atomic potentials)"<< endl << "   [--adfstemcorrfem] simulate fluctuation microscopy using aberration corrected adfstem mode" << endl << "   [--adfstemuncorrfem] simulate fluctuation microscopy using adfstem mode without aberration correction" << endl << "   [--adfstemcorr] simulate aberration corrected adfstem" << endl << "   [--adfstemuncorr] simulate adfstem mode without aberration correction" << endl << "   [--bfctemcorr] simulate bright field TEM with aberration correction" << endl << "   [--bfctemuncorr] simulate bright field TEM without aberration correction" << endl << "   [--dr <azimuthal_binning_size_factor>] prefactor of sqrt(dx^2+dy^2) in bin size" << endl << "   [--minslice <minSliceThickness> minimum slice thickness in Angstroms, default is 1] " << endl << "   [--images] generate and save images" << endl << "   [--netcdfimages] save images as netcdf files" << endl << "   [--netcdfvariance] save 1-D variance as netcdf files" << endl << "   [--D1]" << endl << "   [--D2]" << endl << "   [--D3]" << endl << "   [--D4]" << endl << "   [--GT17]" << endl << "   [--debug] enable verbose debug output to files and stdout" << endl ;


int main( int argc, char* argv[])
{

   //////////////////////////////////////////////////////////////////
   // initializing MPI and fftw_mpi capabilities
   //////////////////////////////////////////////////////////////////
   
   //int provided; // used for fftw mpi + threads

   MPI_Init( &argc, &argv);
   //MPI_Init_thread( &argc, &argv, MPI_THREAD_FUNNELED, &provided );
   //threads_ok = ( provided >= MPI_THREAD_FUNNELED );

   int mynode, totalnodes;
   MPI_Comm_size( MPI_COMM_WORLD, &totalnodes );
   MPI_Comm_rank( MPI_COMM_WORLD, &mynode );
   int rootnode = 0;

   //if ( threads_ok )  threads_ok = fftw_init_threads();
   fftw_mpi_init();
   //const int nthreads = 1; // number of threads for fftw to use per process
   //if ( threads_ok ) fftw_plan_with_nthreads( nthreads );

   //////////////////////////////////////////////////////////////////
   // Microscope parameters
   //////////////////////////////////////////////////////////////////
   
   // TODO: Nx, Ny should be samples per Angstrom, based on xperiod ...
   //        but somehow the number of samples per dimension should be 
   //        restricted to a power of 2, for efficient fftw calculations.
   int Nx_int;
   int Ny_int;

   ptrdiff_t Nx;
   ptrdiff_t Ny;

   double NxNy;
   double NxNy_sqr;
   double sqrtNxNy; 

   double VV; 
   double Cs3; 
   double Cs3_corr; // Cs3 for use in abberation corrected routines
   double Cs5; 
   double defocus; 
   double defocus_spread = 1.0;

   unsigned int dupe_x;
   dupe_x = 1;
   unsigned int dupe_y;
   dupe_y = 1;
   unsigned int dupe_z;
   dupe_z = 1;
   
   // TODO : implement adfstem with defocus_spread

   double condenser_illumination_angle;
   // NOTE: Scherzer defocus is given by sqrt( C_s * lambda ),
   //       and the extended Scherzer defocus is 1.2 sqrt( C_s * lambda )
   double alpha_max_sqr; // = pow( strtod( argv[9], NULL), 2);
   double alpha_max;
   string outFileName_prefix; // = argv[2];
   //string input_position_lammps_file;

   // STEM and FEM raster point separation in both x and y directions
   double raster_spacing;

   //////////////////////////////////////////////////////////////////
   // Specimen domain
   //////////////////////////////////////////////////////////////////
   double xmin;
   double ymin;
   double zmin;
   double xperiod;
   double yperiod;
   double zperiod;

   double* q;
   double* q_duped;
   unsigned int* Z_array; // couldn't use size_t because it must have an 
                          //  MPI type
   unsigned int initial_population;
   unsigned int duped_population;

   //std::vector< double > slice_locations;

   //////////////////////////////////////////////////////////////////
   // read input parameters
   //////////////////////////////////////////////////////////////////
   std::vector<string> args( argv, argv + argc );
   unsigned int failflag = 0;
   unsigned int input_flag_m = 0;// microscope parameters
   unsigned int input_flag_mf = 0;// file containing microscope parameters
   unsigned int input_flag_o = 0;// output file name prefix
   unsigned int input_flag_a = 0;// atom position and species input file
   unsigned int input_flag_defocus = 0;   // defocus
   unsigned int input_flag_spread = 0;  // defocus spread 
   unsigned int input_flag_adfstem_corrected = 0;        // 
   unsigned int input_flag_adfstem_uncorrected = 0;      // 
   unsigned int input_flag_bfctem_corrected = 0;         // 
   unsigned int input_flag_bfctem_uncorrected = 0;       // 
   unsigned int input_flag_fem = 0;
   unsigned int input_flag_gt17 = 0;
   unsigned int input_flag_d1 = 0;
   unsigned int input_flag_d2 = 0;
   unsigned int input_flag_d3 = 0;
   unsigned int input_flag_d4 = 0;
   unsigned int input_flag_scherzer_defocus = 0;
   unsigned int input_flag_scherzer_alphamax = 0;
   unsigned int input_flag_scherzer_cs3 = 0;
   unsigned int input_flag_cs3 = 0;
   unsigned int input_flag_cs5 = 0;
   unsigned int input_flag_alpha_max = 0;
   unsigned int input_flag_aberration_correction = 0;
   unsigned int input_flag_raster_spacing = 0;
   unsigned int input_flag_pap_tif = 0;
   unsigned int input_flag_dupe = 0;
   unsigned int input_flag_image_output = 0;
   unsigned int input_flag_netcdf_images = 0;
   unsigned int input_flag_netcdf_variance = 0;
   unsigned int input_flag_debug = 0;

   //double azimuthal_binning_size_factor = 1.0;

   double azimuthal_binning_size_factor = 0.70710678118654746;// 1/sqrt(2)

   double minSliceThickness = 1.0;

   // Example input:
   // -m <Nx> <Ny> <VV> <condenser_illumination_angle> 
   // -o <output file prefix>
   // -a <lammps or xyz style position input file>
   // --defocus <defocus> 
   // --alphamax <alphamax>
   // --spread <defocus spread>
   // --cs3 <objective lens third order spherical aberation>
   // --adfstemuncorrfem
   // --scherzer
   // --dupe <dupe_x> <dupe_y> <dupe_z>
   // --images
   // --netcdf
   // --debug
   // --paptif
   
   //double tmp_double;// for accepting input from operator>>() and 
   //                      pushing onto a vector<double> 

   for ( size_t idx=1; idx < args.size(); idx++)
   {
      if ( args[idx] == "-o" )
      {
         if (idx + 1 < args.size()) 
            outFileName_prefix = string(args[idx + 1]);//=string(argv[2]);
         input_flag_o = 1;
         idx += 1;
      }
      else if ( args[idx] == "-m" )
      {
         if (idx + 1 < args.size()) 
            istringstream(args[idx + 1]) >> Nx_int;
         if (idx + 2 < args.size()) 
            istringstream(args[idx + 2]) >> Ny_int;
         //if (idx + 3 < args.size()) 
         // istringstream(args[idx + 3]) >> Nz_int;
         if (idx + 3 < args.size()) istringstream(args[idx + 3]) >> VV;

         Nx = (ptrdiff_t) Nx_int;
         Ny = (ptrdiff_t) Ny_int;
         // Nz = (ptrdiff_t) Nz_int;

         NxNy = Nx * Ny;
         NxNy_sqr = NxNy * NxNy;
         if ( Nx == Ny ) sqrtNxNy =  Nx ;   
         else sqrtNxNy = sqrt( NxNy );

         input_flag_m = 1;
         idx += 3;
      }
      // TODO: specify an input file structure
      else if (args[idx] == "--mf")   // microscope parameter file name
      {
         string microscope_file_name;
         if (idx + 1 < args.size()) 
            microscope_file_name = string(args[idx + 1]);
            //istringstream( args[idx + 1] ) >> microscope_file_name;

         if ( mynode == rootnode )
         {
            cerr << "Error, file could not be read : " 
               << microscope_file_name << endl;
            cerr << "read_microscope_file() has not yet been implemented" 
               << endl;
         }
         failflag = 1 ;
         //microscope_file_name = string(argv[2]);

         //if (
         //   read_microscope_file(
         //      microscope_file_name, 
         //      Nx_int, Ny_int, Nz_int,
         //      VV,
         //      Cs3,
         //      Cs5,
         //      condenser_illumination_angle,
         //      alpha_max
         //      )
         //   )
         //{
         //   cerr << "Error, file could not be read : " << args[idx + 1] 
         //      << endl;
         //   failflag = 1 ;
         //}

         //alpha_max_sqr = pow( alpha_max, 2);

         //Nx = (ptrdiff_t) Nx_int;
         //Ny = (ptrdiff_t) Ny_int;
         //Nz = (ptrdiff_t) Nz_int;

         input_flag_mf = 1;
         idx += 1;
      }
      else if ( args[idx] == "-a" )
      {
    // TODO: would it be more efficient to read only on rootnode and Bcast?
         //if ( mynode == rootnode )
         //{
         string input_filename;
         //istringstream( args[idx + 1] ) >> input_filename;
         if (idx + 1 < args.size()) 
            input_filename = string( args[ idx + 1] );
         
         if ( 
               input_filename.substr( 
                  input_filename.find_last_of(".") + 1 
                  )
               == "xyz" 
            )
         {
            if ( mynode == rootnode )
            {
               cout << "Reading xyz format file: "
                 << input_filename  << endl; // debug
            }
            if(  
                  read_position_xyz_file(
                     input_filename,
                     q, // will be allocated within this call
                     Z_array, // will be allocated within this call
                     initial_population,
                     xmin, ymin, zmin,
                     xperiod, yperiod, zperiod,
                     input_flag_debug,
                     mynode, rootnode, MPI_COMM_WORLD
                     )
              )
            {
               if ( mynode == rootnode )
               {
                  cerr << "Error, file could not be read : " 
                     << args[idx +1] << endl;
               }
               failflag = 1 ;
            }
         }
         else
         {
            if ( mynode == rootnode )
            {
               cout << "Reading lammps format file: "
                 << input_filename  << endl; // debug
            }
            if(  
                  read_position_lammps_file(
                     input_filename,
                     q,       // will be allocated within this call
                     Z_array, // will be allocated within this call
                     initial_population,
                     xmin, ymin, zmin,
                     xperiod, yperiod, zperiod,
                     input_flag_debug,
                     mynode, rootnode, MPI_COMM_WORLD
                     )
              )
            {
               if ( mynode == rootnode )
               {
                  cerr << "Error, file could not be read : " 
                     << args[idx +1] << endl;
               }
               failflag = 1 ;
            }
         }

         double periods_and_mins[6];
         if ( mynode == rootnode )
         {
            periods_and_mins[0] = xperiod;
            periods_and_mins[1] = yperiod;
            periods_and_mins[2] = zperiod;
            periods_and_mins[3] = xmin;
            periods_and_mins[4] = ymin;
            periods_and_mins[5] = zmin;
         }

         MPI_Bcast( periods_and_mins, 6, MPI_DOUBLE, 
               rootnode, MPI_COMM_WORLD);

         //// multiple mpi calls are probably slow
         //MPI_Bcast( &xperiod, 1, MPI_DOUBLE, rootnode, MPI_COMM_WORLD);
         //MPI_Bcast( &yperiod, 1, MPI_DOUBLE, rootnode, MPI_COMM_WORLD);
         //MPI_Bcast( &zperiod, 1, MPI_DOUBLE, rootnode, MPI_COMM_WORLD);
         //MPI_Bcast( &xmin, 1, MPI_DOUBLE, rootnode, MPI_COMM_WORLD);
         //MPI_Bcast( &ymin, 1, MPI_DOUBLE, rootnode, MPI_COMM_WORLD);
         //MPI_Bcast( &zmin, 1, MPI_DOUBLE, rootnode, MPI_COMM_WORLD);
         if ( mynode != rootnode )
         {
            xperiod = periods_and_mins[0];
            yperiod = periods_and_mins[1];
            zperiod = periods_and_mins[2];
            xmin = periods_and_mins[3];
            ymin = periods_and_mins[4];
            zmin = periods_and_mins[5];
         }


         input_flag_a = 1;
         idx += 1;
      }
      else if ( args[idx] == "--scherzer_defocus")
      {
         input_flag_scherzer_defocus = 1;
      }
      else if ( args[idx] == "--scherzer_alphamax")
      {
         input_flag_scherzer_alphamax = 1;
      }
      else if ( args[idx] == "--scherzer_cs3")
      {
         input_flag_scherzer_cs3 = 1;
      }
      else if ( args[idx] == "--defocus" )
      {
         if (idx + 1 < args.size()) 
            istringstream( args[idx + 1] ) >> defocus;
         input_flag_defocus= 1;
         idx += 1;
      }
      else if ( args[idx] == "--alphamax" )
      {
         if (idx + 1 < args.size()) 
            istringstream( args[idx + 1] ) >> alpha_max;
         input_flag_alpha_max = 1;
         idx += 1;
      }
      else if ( args[idx] == "--spread" )
      {
         if (idx + 1 < args.size()) 
            istringstream( args[idx + 1] ) >> defocus_spread;
         if (idx + 2 < args.size()) 
            istringstream( args[idx + 2] ) >> condenser_illumination_angle;
         input_flag_spread = 1;
         idx += 2;
      }
      else if ( args[idx] == "--cs3" )
      {
         if (idx + 1 < args.size()) 
            istringstream( args[idx + 1] ) >> Cs3;
         input_flag_cs3 = 1;
         idx += 1;
      }
      else if ( args[idx] == "--cs5" )
      {
         if (idx + 1 < args.size()) 
            istringstream( args[idx + 1] ) >> Cs5;
         input_flag_cs5 = 1;
         idx += 1;
      }
      else if ( args[idx] == "--rasterspacing" )
      {
         if ( idx + 1 < args.size())
            istringstream( args[idx + 1] ) >> raster_spacing;
         input_flag_raster_spacing = 1;
         idx += 1;
      }
      else if ( args[idx] == "--adfstemuncorrfem" )
      {
         input_flag_adfstem_uncorrected = 1;
         input_flag_fem = 1;
      }
      else if ( args[idx] == "--adfstemcorrfem" )
      {
         input_flag_adfstem_corrected = 1;
         input_flag_fem = 1;
      }
      else if ( args[idx] == "--GT17" )
      {
         input_flag_gt17 = 1;
      }
      else if ( args[idx] == "--D1" )
      {
         input_flag_d1 = 1;
      }
      else if ( args[idx] == "--D2" )
      {
         input_flag_d2 = 1;
      }
      else if ( args[idx] == "--D3" )
      {
         input_flag_d3 = 1;
      }
      else if ( args[idx] == "--D4" )
      {
         input_flag_d4 = 1;
      }
      else if ( args[idx] == "--adfstemcorr" )
      {
         input_flag_adfstem_corrected = 1;
      }
      else if ( args[idx] == "--adfstemuncorr" )
      {
         input_flag_adfstem_uncorrected = 1;
      }
      else if ( args[idx] == "--bfctemcorr" )
      {
         input_flag_bfctem_corrected = 1;
      }
      else if ( args[idx] == "--bfctemuncorr" )
      {
         input_flag_bfctem_uncorrected = 1;
      }
      else if ( args[idx] == "--paptif" )
      {
          input_flag_pap_tif = 1;
      }
      else if ( args[idx] == "--dupe" )
      {
         if (idx + 3 < args.size()) 
            istringstream( args[idx + 1] ) >> dupe_x;
            istringstream( args[idx + 2] ) >> dupe_y;
            istringstream( args[idx + 3] ) >> dupe_z;
         input_flag_dupe = 1;
         idx += 3;
      }
      else if ( args[idx] == "--dr" )
      {
         if ( idx + 1 < args.size())
            istringstream( args[idx + 1] ) 
               >> azimuthal_binning_size_factor;
         idx += 1;
      }
      else if ( args[idx] == "--minslice" )
      {
         if ( idx + 1 < args.size())
            istringstream( args[idx + 1] ) 
               >> minSliceThickness;
         idx += 1;
      }
      else if ( args[idx] == "--images" )
      {
         input_flag_image_output = 1;
      }
      else if ( args[idx] == "--netcdfimages" )
      {
         input_flag_netcdf_images = 1;
      }
      else if ( args[idx] == "--netcdfvariance" )
      {
         input_flag_netcdf_variance = 1;
      }
      else if ( args[idx] == "--debug" )
      {
         input_flag_debug  = 1;
      }
      else
      {
         if ( mynode == rootnode )
         {
            cerr << "Error, unexpected argument : " << args[idx] << endl;
         }
         failflag = 1;
      }
   }
      if ( (! input_flag_gt17 ) 
            && (! input_flag_d1 ) && (! input_flag_d2 ) 
            && (! input_flag_d3 ) && (! input_flag_d4 ) )
      {
         input_flag_d1 = 1; // default variance calculation mode
      }

      // debug
      if( mynode == rootnode  && input_flag_debug )
      {
         cout << 
         "input_flag_m " << 
         input_flag_m << endl
         << "input_flag_mf " << 
         input_flag_mf << endl <<
         "input_flag_o " << 
         input_flag_o << endl <<
         "input_flag_a " << 
         input_flag_a << endl <<
         "input_flag_defocus " << 
         input_flag_defocus << endl <<
         "input_flag_spread " << 
         input_flag_spread << endl <<
         "input_flag_dupe " << 
         input_flag_dupe << endl <<
         "input_flag_adfstem_corrected " << 
         input_flag_adfstem_corrected << endl <<
         "input_flag_adfstem_uncorrected  " << 
         input_flag_adfstem_uncorrected  << endl <<
         "input_flag_bfctem_corrected  " << 
         input_flag_bfctem_corrected  << endl <<
         "input_flag_bfctem_uncorrected  " << 
         input_flag_bfctem_uncorrected  << endl <<
         "input_flag_fem  " << 
         input_flag_fem  << endl <<
         "input_flag_gt17 " << 
         input_flag_gt17  << endl <<
         "input_flag_d1 " << 
         input_flag_d1  << endl <<
         "input_flag_d2 " << 
         input_flag_d2  << endl <<
         "input_flag_d3 " << 
         input_flag_d3  << endl <<
         "input_flag_d4 " << 
         input_flag_d4  << endl <<
         "input_flag_scherzer_defocus " << 
         input_flag_scherzer_defocus << endl <<
         "input_flag_scherzer_alphamax " << 
         input_flag_scherzer_alphamax << endl <<
         "input_flag_scherzer_cs3 " << 
         input_flag_scherzer_cs3 << endl <<
         "input_flag_cs3 " << 
         input_flag_cs3 << endl <<
         "input_flag_cs5 " << 
         input_flag_cs5 << endl <<
         "input_flag_alpha_max " << 
         input_flag_alpha_max << endl <<
         "input_flag_aberration_correction " << 
         input_flag_aberration_correction << endl <<
         "input_flag_raster_spacing " << 
         input_flag_raster_spacing << endl <<
         "input_flag_debug " << 
         input_flag_debug 
         << endl;
      }
      // end debug

   if (
         !(
            input_flag_o   // output file name prefix
            &&
            input_flag_a   // atom position and species input file
            &&
            ( // command line xor file input of microscope parameters
               // exclusive or:
               (! input_flag_m) != (! input_flag_mf) 
            )
            &&
            // if aberration correction is used with bfctem, require 
            //    input_flag_spread, input_flag_cs3, input_flag_cs5
            !(
               (
                  input_flag_adfstem_corrected 
                  || 
                  input_flag_bfctem_corrected
               )
               &&
               (!
                  (
                     //input_flag_spread
                     //&&
                                                // defocus_spread is 
                                                // currently only used
                                                // in bfctem
                     input_flag_cs5
                     &&
                     (
                        input_flag_scherzer_cs3
                        ||
                        input_flag_cs3
                     )
                  )
               )
            )
            &&
            !(
               (// require specifying Cs3 if using uncorrected TEM
                  input_flag_adfstem_uncorrected
                  ||
                  input_flag_bfctem_uncorrected
               )
               && 
               (! input_flag_cs3)
            )
            &&
            (
               input_flag_scherzer_defocus 
               || 
               input_flag_defocus
            )
            &&
            (
               input_flag_scherzer_alphamax
               || 
               input_flag_alpha_max
            )
         )
      )
   {
      if ( mynode == rootnode )
      {
         cout << argv[0] << " : lacking required parameters" << endl;
         // TODO: the following two if(){} statements duplicate tests above
         //       Perhaps implement these error messages inside the above
         //       tests or vice versa.
         if (
              (input_flag_adfstem_corrected || input_flag_bfctem_corrected)
               && 
               ! (
                  //input_flag_spread
                  //&&
                                          // defocus_spread is 
                                          // currently only used 
                                          // in bfctem
                  (input_flag_cs3 || input_flag_scherzer_cs3)
                  && 
                  input_flag_cs5
               )
                
            )
         {
            cout << "Use of aberraction correction requires"
               << " specifying --defocus_spread and --cs5. If not using"
               << " --scherzer_cs3, then --cs3 is also required." << endl;
         }
         if (
               (// require specifying Cs3 if using uncorrected TEM
                  input_flag_adfstem_uncorrected
                  ||
                  input_flag_bfctem_uncorrected
               )
               && 
               (! input_flag_cs3)
            )
         {
            cout << "Use of uncorrected aberration"
             << " requires specifying --cs3." << endl;
         }

         if (
               !(
                  input_flag_scherzer_defocus 
                  || 
                  input_flag_defocus
               )
            )
         {
            cout << "You must specify either --defocus <value [angstrom]>"
               << " or --scherzer_defocus" << endl;
         }
         if (
               !(
                  input_flag_scherzer_alphamax
                  || 
                  input_flag_alpha_max
               )
            )
         {
            cout << "You must specify either --alphamax <value [radians]>"
               << " or --scherzer_alphamax" << endl;
         }
      }
      failflag = 1;
   }

   if (
         failflag == 0
         &&
         ! input_flag_adfstem_uncorrected
         &&
         ! input_flag_adfstem_corrected
         &&
         ! input_flag_bfctem_corrected
         &&
         ! input_flag_bfctem_uncorrected
      )
   {
      if ( mynode == rootnode )
      {
         cout << argv[0] << " : you must specify at least one of the following" << endl
            << "   [--adfstemcorrfem] simulate fluctuation microscopy using aberration corrected adfstem mode" 
            << endl 
            << "   [--adfstemuncorrfem] simulate fluctuation microscopy using adfstem mode without aberration correction" 
            << endl 
            << "   [--adfstemcorr] simulate aberration corrected adfstem" 
            << endl 
            << "   [--adfstemuncorr] simulate adfstem mode without aberration correction" 
            << endl 
            << "   [--bfctemcorr] simulate bright field TEM with aberration correction" 
            << endl 
            << "   [--bfctemuncorr] simulate bright field TEM without aberration correction" 
            << endl;
      }
      failflag = 1;
   }
   if ( input_flag_dupe && ! input_flag_fem && mynode == rootnode )
   {
      cout << "Note: --dupe flag only duplicates the sample for FTEM."
         << " The sample will not be duplicated this time." 
         << endl;
   }

   MPI_Bcast( &failflag, 1, MPI_UNSIGNED, rootnode, MPI_COMM_WORLD );
   if ( failflag ) 
   {
      fftw_mpi_cleanup();
      MPI_Finalize();
      if ( mynode == rootnode)
      {
         cout << "Exiting tem() call due to failure." << endl;
         PRINT_USAGE // macro defined above
      }
      return EXIT_FAILURE;
   }

   // At this point rootnode should have arrays of all scatterer 
   //  positions and species Z_array.
   
   // Each node has a scatterer_pap_LUT, and should instantiate its 
   //  own vector of scatterers upon recieving q[][3] and Z_array[] 
   //  from rootnode.

   //if ( mynode != rootnode ) 
   //{
   //   q = new double[ 3 * initial_population ];
   //   Z_array = new unsigned int[ initial_population ];
   //}
   //MPI_Bcast( q, 3 * initial_population, MPI_DOUBLE,    // 3-D
   //           rootnode, MPI_COMM_WORLD);
   //MPI_Bcast( Z_array, initial_population, MPI_UNSIGNED, 
   //           rootnode, MPI_COMM_WORLD);
   
   // TODO: ensure that all required parameters have been assigned at
   //        every node

   //////////////////////////////////////////////////////////////////
   // Physical constants
   //////////////////////////////////////////////////////////////////
   const double m_e0 = 9.109e-31; // electron rest mass [kg]
   const double cc = 299792458; // speed of light in vacuum [\frac{m}{s}]
   const double ee = 1.60217662e-19; // charge of an electron [C]
   const double pi = 3.14159265369;
   const double hh = 6.62607004e-34; // [\frac{m^{2} kg}{s}]

   // Electron wavelength
   const double lambda = hh * cc * pow(10, 10)  // [\AA]
                / sqrt(ee * VV * (2.0 * m_e0 * cc * cc + ee * VV));

   // Lorentz factor for reciprocal space V_z calculation
   const double gamma
      = (hh * pow(10,10)/ lambda) * (hh * pow(10,10)/ lambda)
         * (m_e0 * cc * cc + ee * VV) 
         / ((2 * m_e0 * cc * cc + ee * VV) * ee * VV * m_e0);

   if ( mynode == rootnode && input_flag_debug )
   {
      // sigma: interaction parameter (for comparisons to book value)
      const double sigma
         = (2 * pi /(lambda * VV))
         * (m_e0 * pow(cc,2) + ee * VV)
         / (2 * m_e0 * pow(cc,2) + ee * VV); // [1/(V \AA)]
      cout << "Physical constants --------------------" << endl;
      cout << "VV : " << VV << " [V]" << endl;
      cout << "m_e0 : " << m_e0 << " [kg]" << endl;
      cout << "cc : " << cc << " [m/s]" << endl;
      cout << "ee : " << ee << " [C]" << endl;
      cout << "hh : " << hh << " [m^{2} kg/s]" << endl;
      cout << "lambda : " << lambda << " [AA]" << endl;
      cout << "sigma : " << sigma << " [1/(V AA)]" << endl;
      cout << "gamma : " << gamma << endl;
      cout << "---------------------------------------" << endl;
      cout << endl;
   }


   //////////////////////////////////////////////////////////////////
   // Allocate local variables for fftw_mpi
   //////////////////////////////////////////////////////////////////
   ptrdiff_t Nx_local;
   ptrdiff_t local_alloc_size_fftw;
   ptrdiff_t local_0_start_fftw;

   local_alloc_size_fftw
      = fftw_mpi_local_size_2d(     // fftw will be 2-D here
            Nx, Ny,
            MPI_COMM_WORLD,
            &Nx_local,
            &local_0_start_fftw );

   if ( input_flag_debug )
      cout << "node " << mynode 
         << ", (Nx_local, local_0_start_fftw) : (" 
         << Nx_local << ", " << local_0_start_fftw << ")"<< endl; // debug

   if ( local_alloc_size_fftw != Nx_local * Ny ) 
   {
      cerr << "local_alloc_size_fftw != Nx_local * Ny " << endl;
      return( EXIT_FAILURE );
   }

   //////////////////////////////////////////////////////////////////
   // Retrieve fftw wisdom 
   //////////////////////////////////////////////////////////////////
   // TODO: ensure wisdom file is customized to local machine type, 
   //        put machine type in the wisdom file name

   string wisdom_filename
      = "wisdom_file_c2c2c_cifb_" + TEM_NS::to_string(Nx_local)
         + "x" + TEM_NS::to_string(Ny) + "_MPI" 
         + TEM_NS::to_string(totalnodes);
   
   int flag_wisdomFile = 1;
   if ( mynode == rootnode )
   {
      flag_wisdomFile = 
         fftw_import_wisdom_from_filename( wisdom_filename.c_str() );
         // returns 0 (== false) on failure
      if ( ! flag_wisdomFile )
      {
         cout << "Could not open wisdom file: " << wisdom_filename
            << " , will create it after execution ..." << endl;
      }
   }

   //////////////////////////////////////////////////////////////////
   // Modify parameters sensitive to sample duplication
   //////////////////////////////////////////////////////////////////
   // NOTE: change Nx, Ny, Nz ?  No? 
   //       change xperiod, ..., ? For some cases, yes.
   //       The period of the cell containing duplicates is needed
   //       for the inverse area factor, but not needed for assigning
   //       raster points.
   
   double xperiod_duped, yperiod_duped, zperiod_duped;
   double duped_periods[3]; // Array of periods so that only one Bcast 
                            //  is used.
   if ( mynode == rootnode )
   {
      xperiod_duped = xperiod * dupe_x;
      yperiod_duped = yperiod * dupe_y;
      zperiod_duped = zperiod * dupe_z;

      duped_periods[0] = xperiod_duped;
      duped_periods[1] = yperiod_duped;
      duped_periods[2] = zperiod_duped;
   }

   //MPI_Bcast(&xperiod_duped, 1, MPI_DOUBLE, rootnode, MPI_COMM_WORLD);
   //MPI_Bcast(&yperiod_duped, 1, MPI_DOUBLE, rootnode, MPI_COMM_WORLD);
   //// !!!!!!!!! THIS Bcast isn't working for zperiod_duped ...
   //MPI_Bcast(&zperiod_duped, 1, MPI_DOUBLE, rootnode, MPI_COMM_WORLD);
   
   MPI_Bcast(duped_periods, 3, MPI_DOUBLE, rootnode, MPI_COMM_WORLD);

   if ( mynode != rootnode )
   {
      xperiod_duped = duped_periods[0];
      yperiod_duped = duped_periods[1];
      zperiod_duped = duped_periods[2];
      cout << " node " << mynode 
         << " (duped_periods[2], zperiod_duped) : (" 
         << duped_periods[2] << ", " << zperiod_duped 
         << ")" <<  endl;
   }

   if ( input_flag_debug )
      cout << "Check MPI_Bcast() quality:"
         << "  node " << mynode 
         << ", (xperiod, yperiod, zperiod) : ("
         << xperiod << ", " << yperiod << ", " 
         << zperiod << ")" <<  endl
         << "  node " << mynode 
         << ", (xperiod_duped, yperiod_duped, zperiod_duped) : ("
         << xperiod_duped << ", " << yperiod_duped << ", " 
         << zperiod_duped << ")" <<  endl;

   //////////////////////////////////////////////////////////////////
   // Establish bandwidth limiting cutoff values
   //////////////////////////////////////////////////////////////////
   // Kirkland page 91: 
   // "This bandwidth limit should be applied to both the projected atomic
   //  potential and the transmission function ..." 
   //   derived from it, 
   //  ... 
   //  and the wave function with which it will be multiplied by in 
   //  realspace.
   //
   //  page 147:
   // "When the electron wave function is multiplied by the transmission 
   //  function (at each step or slice of the multislice method) its 
   //  bandwidth doubles and is then reduced to the required maximum by 
   //  convolution with the propagator function (multiplication in 
   //  reciprocal space)."

   double bwcutoff_pap;
   if ( Nx/(2 * xperiod_duped) < Ny/(2 * yperiod_duped) ) 
      bwcutoff_pap = Nx/(2 * xperiod_duped);
      //bwcutoff_pap = Nx/(3 * xperiod_duped);
   else 
      bwcutoff_pap = Ny/(2 * yperiod_duped);
      //bwcutoff_pap = Ny/(3 * yperiod_duped);

   double bwcutoff_t; 
   bwcutoff_t = (2.0/3.0) * bwcutoff_pap; // transmission function
   //bwcutoff_t = bwcutoff_pap; // transmission function


   // TODO: delete the following two lines after running cmd_046.txt
   //bwcutoff_t = 0.25 * bwcutoff_t;
   //bwcutoff_pap = 0.25 * bwcutoff_pap;

   double bwcutoff_t_sqr; 
   bwcutoff_t_sqr = bwcutoff_t * bwcutoff_t;

   if ( mynode == rootnode && input_flag_debug )
   {
      cout << "Bandwidth limits " << endl // debug
         << "  projected atomic potential bandwidth : " // debug
         << bwcutoff_pap << endl // debug
         << "  transmission function bandwidth : " // debug
         << bwcutoff_t << endl; // debug
   }

   //////////////////////////////////////////////////////////////////
   // Create domains
   //////////////////////////////////////////////////////////////////

   // Only splitting the kx domain amongst nodes, not x, y, or ky

   double* kx_joined; kx_joined = new double[ Nx ];
   double* xx_joined; xx_joined = new double[ Nx ];
   double* kx_local; kx_local = new double[ Nx_local ];
   double* xx_local; xx_local = new double[ Nx_local ];// debug
   double* ky; ky = new double[ Ny ];
   double* yy; yy = new double[ Ny ];
   double delta_x, delta_y, delta_kx, delta_ky;
   delta_x = xperiod_duped / Nx; 
   delta_y = yperiod_duped / Ny; 
   delta_kx = 1.0/xperiod_duped;
   delta_ky = 1.0/yperiod_duped;

   const double kxperiod = Nx / xperiod_duped; 
   const double kyperiod = Ny / yperiod_duped;

   if ( mynode == rootnode )
   {
      // The reciprocal space domain my be restricted to 2-D, since the
      //  Fourier projection theorem is being used.

      domain_2D_periodic_recip( Nx, Ny, 
            xperiod_duped, yperiod_duped, 
            kxperiod, kyperiod, 
            delta_kx, delta_ky,
            kx_joined, ky
            );

      domain_2D_periodic( Nx, Ny, 
            xperiod_duped, yperiod_duped, 
            xmin, ymin, 
            delta_x, delta_y,
            xx_joined, yy
            );

      if ( input_flag_debug )
      {
         int precision = 7;
         int width = precision + 2;
         cout << "Periodicity from input file: " << endl 
            << " xperiod_duped = "
            << setw(width) << setprecision(precision)
            << xperiod_duped << " [A]"<< endl
            << " yperiod_duped = " 
            << setw(width) << setprecision(precision)
            << yperiod_duped << " [A]" <<  endl
            << " kxperiod = "
            << setw(width) << setprecision(precision)
            << kxperiod << " [A^{-1}]" << endl
            << " kyperiod = "
            << setw(width) << setprecision(precision)
            << kyperiod << " [A^{-1}]" << endl;
         cout << "Domains :" << endl;
         cout << " xx_joined: [" 
            << setw(width) << setprecision(precision)
            << xx_joined[0] << ": " 
            << setw(width) << setprecision(precision)
            << delta_x << ": "
            << setw(width) << setprecision(precision)
            <<  xx_joined[Nx-1] << "]" << endl;//debug
         cout << " yy: [" 
            << setw(width) << setprecision(precision)
            << yy[0] << ": " 
            << setw(width) << setprecision(precision)
            << delta_y << ": "
            << setw(width) << setprecision(precision)
            <<  yy[Ny-1] << "]" << endl;//debug

         cout << " kx_joined: [" 
            << setw(width) << setprecision(precision)
            << kx_joined[(Nx/2)-1] << ": " 
            << setw(width) << setprecision(precision)
            << delta_kx << ": "
            << setw(width) << setprecision(precision)
            <<  kx_joined[Nx/2] << "]" << endl;//debug
         cout << " ky: [" 
            << setw(width) << setprecision(precision)
            << ky[(Ny/2)-1] << ": " 
            << setw(width) << setprecision(precision)
            << delta_ky << ": "
            << setw(width) << setprecision(precision)
            <<  ky[Ny/2] << "]" << endl;//debug

      }
   }
   // TODO: note that the following is just a test and should be deleted
   MPI_Scatter( xx_joined, Nx_local, MPI_DOUBLE, // 2-D
               xx_local, Nx_local, MPI_DOUBLE,
               rootnode, MPI_COMM_WORLD);
   // TODO: note that the above is just a test and should be deleted

   MPI_Scatter( kx_joined, Nx_local, MPI_DOUBLE, // 2-D
               kx_local, Nx_local, MPI_DOUBLE,
               rootnode, MPI_COMM_WORLD);
   MPI_Bcast( ky, Ny, MPI_DOUBLE, rootnode, MPI_COMM_WORLD);
   MPI_Bcast( yy, Ny, MPI_DOUBLE, rootnode, MPI_COMM_WORLD);
   MPI_Bcast( kx_joined, Nx, MPI_DOUBLE, rootnode, MPI_COMM_WORLD);
   MPI_Bcast( xx_joined, Nx, MPI_DOUBLE, rootnode, MPI_COMM_WORLD);
   
   //////////////////////////////////////////////////////////////////
   // Strategy for efficiently assigning the projected atomic potential:
   //  - At rootnode: reduce Z_array[] to a list of unique atomic 
   //     numbers in a std::vector<unisgned int> Z_vector
   //  - broadcast Z_array elements to each node
   //  - instantiate a scatterer_pap_LUT:vector<scatterer*>
   //    - for each unique species, evaluate the pap centered at 
   //       (xmin, ymin)
   //    - transform the LUT of paps into realspace
   //    - MPI_Allgather each pap into locations which will be 
   //       pointed to by a pointer member of scatterer objects
   //////////////////////////////////////////////////////////////////
   
   // debug
   //  Determine the value of idx_local_start_x such that
   //   kx_joined[i + idx_local_start_x] == kx_local[i]
   //   and thus
   //   pap_joined_x[i + idx_local_start_x] == pap_split_x[i]
   size_t idx_local_start_x; 
   size_t idx_tmp; idx_tmp = 0;
   while ( idx_tmp < Nx)
   {
      // NOTE: comparing reciprocal space domains to determine real space 
      //  splitting location ...
      if ( kx_joined[idx_tmp] == kx_local[0] )
      {
         idx_local_start_x = idx_tmp;
         break;
      }
      else 
         ++idx_tmp;
   }
   if ( idx_tmp == Nx)
   {
      cout << "Error : update_transmission_function() failed;" 
         << " could not identify appropriate idx_local_start" << endl;
      return EXIT_FAILURE;
   }
   if ( idx_local_start_x != local_0_start_fftw )
   {
      cout << "Error : idx_local_start_x != local_0_start_fftw"
         << endl << " idx_local_start_x: " << idx_local_start_x << endl
         << endl << " local_0_start_fftw: " << local_0_start_fftw << endl
         << " Shifting of cached projected atomic potentials and probes"
         << " will be erroneous." << endl;
      return( EXIT_FAILURE );
   }
   // end debug

   //////////////////////////////////////////////////////////////////
   // Create a vector containing a list of the unique species types Z
   //////////////////////////////////////////////////////////////////
   
   std::vector<unsigned int> Z_vector;
   unsigned int Z_vector_size;

   // TODO: parallelize this process of elliminating duplicate Zs
   if ( mynode == rootnode )
   {
      // push Z_array elements onto Z_vector, then sort and remove 
      //  duplicates from Z_vector
      for (unsigned int i = 0; i < initial_population; i++)
      {
         Z_vector.push_back( Z_array[i] );
      }

      std::sort( Z_vector.begin(), Z_vector.end() );

      std::vector<unsigned int>::iterator itr;
      itr = std::unique( Z_vector.begin(), Z_vector.end() );
      Z_vector.resize( std::distance( Z_vector.begin(), itr) );

      Z_vector_size = Z_vector.size();
   }

   MPI_Bcast( &Z_vector_size, 1, 
         MPI_UNSIGNED, rootnode, MPI_COMM_WORLD);
   if ( input_flag_debug )
      cout << "node " << mynode // debug
         << ", Z_vector_size : " << Z_vector_size << endl; // debug

   if ( mynode != rootnode ) Z_vector.resize( Z_vector_size );
   
   //MPI_Bcast( &Z_vector[0], Z_vector.size(),
   MPI_Bcast( &Z_vector[0], Z_vector_size,
               MPI_UNSIGNED, rootnode, MPI_COMM_WORLD);

   if ( input_flag_debug )
      cout << "node " << mynode  // debug
         << ", Z_vector.size (resized) : " 
         << Z_vector.size() << endl;// debug

   //////////////////////////////////////////////////////////////////
   // Initialize a look up table containing Zs paired with projected 
   //  atomic potentials centered at (0,0) (I thought is was (xmin,ymin)?)
   //////////////////////////////////////////////////////////////////

   scatterer_pap_LUT myScattererPapLUT(
         Z_vector,
         lambda, gamma, 1/(xperiod_duped * yperiod_duped), 
         bwcutoff_pap,
         kx_local, Nx_local,
         //kx_joined, // kx_joined is not necessary
         Nx,
         xmin,
         ky, Ny, ymin,
         mynode, rootnode, MPI_COMM_WORLD
         );

   if ( input_flag_debug )
   {
      cout << "List of unique Zs after scatter_pap_LUT instantiation,"
        << " node " << mynode << " : ";
      for ( std::vector<unsigned int>::iterator 
               Z_vector_itr = Z_vector.begin();
             Z_vector_itr != Z_vector.end(); 
             ++Z_vector_itr )
      {
         cout << *Z_vector_itr << ", ";
      }
      cout << endl;
      //cout << "Value of myScattererPapLUT members"
      //   for ( std::list<scatterer_pap*>::iterator 
      //         pap_list_itr = pap_list.begin(); 
      //         pap_list_itr != pap_list.end(); 
      //         ++pap_list_itr
      //         )
      //   {
      //   }
   }

   // TODO: Check that the pap_LUT is identical on each node
   if ( input_flag_debug )
   {
      for ( std::list< scatterer_pap* >::const_iterator
            pap_list_itr = myScattererPapLUT.pap_list.begin();
            pap_list_itr != myScattererPapLUT.pap_list.end();
            pap_list_itr++ 
          )
      {
         cout << "node " << mynode << ", myScattererPapLUT.pap_list"
           << " contents : Z = "
            << (*pap_list_itr)->Z
            << " , pap...joined[0] = " << 
            (*pap_list_itr)->projected_atomic_potential_local_joined_re[0]
            << " + i " << 
            (*pap_list_itr)->projected_atomic_potential_local_joined_im[0]
            << endl;
      }
   }
   // end debug

   //////////////////////////////////////////////////////////////////
   // Create the scatterer objects containing position, species, and
   //  corresponding pointers to template projected atomic potentials
   //  centered at (0,0) which will be translated to the scatterer
   //  position by slice.update_transmission_function() .
   //////////////////////////////////////////////////////////////////
   std::vector<scatterer> myScatterers;

   // Create duplicate scatterers to increase FTEM V(|k|) resolution
   if ( input_flag_dupe && input_flag_fem 
        &&
        (dupe_x > 1 || dupe_y > 1 || dupe_z > 1)
        &&
        (dupe_x > 0 && dupe_y > 0 && dupe_z > 0)
      )
   {
      if ( mynode == rootnode && input_flag_debug ) // debug
      {
         cout << "Instantiating the given scatterer positions " << endl
            << "   " << dupe_x << " times in the x direction, " << endl
            << "   " << dupe_y << " times in the y direction, " << endl
            << "   " << dupe_z << " times in the z direction" << endl;
      } // debug

      duped_population 
         = initial_population * dupe_x * dupe_y * dupe_z;
      q_duped = new double[ 3 * duped_population ];

      for ( unsigned int ii = 0; ii < dupe_x; ++ii)
         for ( unsigned int jj = 0; jj < dupe_y; ++jj)
            for ( unsigned int kk = 0; kk < dupe_z; ++kk)
               for (size_t i = 0; i < initial_population; i++)
               {
                  q_duped[
                        3 * (
                          i + initial_population * (
                             kk + dupe_z * (jj + dupe_y * ii)
                          )
                     )
                    ]
                     = q[3*i]     + ii * xperiod;

                  q_duped[
                    1 + 3 * (
                          i + initial_population * (
                             kk + dupe_z * (jj + dupe_y * ii)
                          )
                       )
                     ]
                     = q[3*i + 1] + jj * yperiod;

                  q_duped[
                    2 + 3 * (
                          i + initial_population * (
                             kk + dupe_z * (jj + dupe_y * ii)
                             )
                       )
                    ]
                     = q[3*i + 2] + kk * zperiod;

                  myScatterers.push_back(
                        scatterer(
                           &(q_duped[
                              3 * ( i + initial_population * (
                                     kk + dupe_z * (jj + dupe_y * ii)
                                     )
                                 )
                           ]),
                           Z_array[i],
                           myScattererPapLUT
                           )
                        );
               }
   }
   else
   {
      q_duped = q;
      duped_population = initial_population;

      for (size_t i = 0; i < initial_population; i++)
      {
         //cout << "node " << mynode 
         //   << ", pushing scatterers onto myScatterers : " // debug
         //   << "Z_vector[" << i << "] : " << Z_vector[i] // debug
         //   << ", q[0,1,2] : " // debug
         //   << q[3*i]  // debug
         //   << ", " << q[3*i + 1] // debug
         //   << ", " << q[3*i + 2] // debug
         //   << endl;  // debug

         myScatterers.push_back(
               scatterer( 
                  &(q[3*i]),
                  Z_array[i],
                  myScattererPapLUT
                  )
            );
      }
   }

   delete[] Z_array; // NOTE: be wary of the consequences of freeing this

   // debug
   //cout << "List of unique atomic species after myScatterers instantiation, node " 
   //   << mynode << " : ";
   //for ( std::vector<unsigned int>::iterator 
   //         Z_vector_itr = Z_vector.begin();
   //       Z_vector_itr != Z_vector.end(); 
   //       ++Z_vector_itr )
   //{
   //   cout << *Z_vector_itr << ", ";
   //}
   //cout << endl;
   // end debug

   //////////////////////////////////////////////////////////////////
   // Split the sample along the transmission axis into slices
   //  to be used by all tem simulations
   //////////////////////////////////////////////////////////////////
   // Notes: 
   // "If the slice thickness in the multislice calculation matches 
   // the natural periodicity in the specimen (i.e., an integer 
   // number of slices in the repeat length of the specimen) then the
   // multislice simulation will reproduce the higher order Laue 
   // zones. If the slice thickness does not match the specimen 
   // periodicity then beating can occur between the slice thickness
   // and the specimen periodicity to produce artifacts in the image.
   // The slice thickness can produce an artificial periodicity in 
   // the specimen if it does not match the natural periodicity of 
   // the specimen. If the multislice slice thickness is much larger
   // than the natural periodicity of the speciment then false HOLZ 
   // lines can be created at: 
   // k = \sqrt(\frac{2}{\Delta z \lambda}) eqn 6.95, 
   // where \Delta z is the multislice slice thickness, \lambda 
   // is the electron wavelength and \alpha = \lambda k is the 
   // electron scattering angle. The slice thickness effectively 
   // takes the place of the normal crystal lattice spacing along z.
   // An overly large slice thickness can produce serious artifacts 
   // in a simulated image and should be avoided." 
   // ...
   // "The standard multislice method [(6.92) and (6.93)] is only 
   // accurate to order \Delta z." ... "As long as the slice 
   // thickness is large compared to the effective range of the 
   // atomic potential this approximation is valid. However, if the 
   // slice thickness \Delta z is made less than the range of the 
   // atomic potential (about 1 \AA) then this approximation is no 
   // longer valid. Therefore, if the total projected atomic 
   // potential is used (as is typical) the minimum slice thickness 
   // is about 1 \AA, which sets a limit on the maximum achievable 
   // accuracy of the multislice calculation. There is some 
   // question if 1 \AA is thin enough for heavy atoms, such as gold,
   // at 100 keV, however higher voltage or lower atomic number 
   // should be all right (Watanabe [372]). The electron wavelength 
   // is typically of order 0.03 \AA, which is much smaller than 
   // the minimum slice thickness. If the full wave function were 
   // used then the slice thickness would have to be much smaller 
   // than the electron wavelength and the multislice method would 
   // not work at all. However, only the slowly varying part of the 
   // wave function is calculated so a slice thickness of several 
   // Angstroms is acceptable."
   // - Kirkland (2010)

   std::list< slice* > slicesOfScatterers;

   // TODO: parallelize assign_atoms_to_slices_z_auto() or compute on
   // a single node and broadcast the results.
   assign_atoms_to_slices_z_auto(
         myScatterers,
         xmin,ymin,zmin,
         xperiod_duped, yperiod_duped, zperiod_duped,
         minSliceThickness,  // minimum slice thickness [\AA]
         //1.0,  // minimum slice thickness [\AA]
         slicesOfScatterers,
         input_flag_debug,
         mynode,
         rootnode,
         MPI_COMM_WORLD
         );

   // Projected atomic potentials of each slice are not yet 
   // calculated. The next step should use the 
   // scatterer.pap_to_translate_{re,im} members to calculate the 
   // total pap of each slice 
   //    (slicesOfScatterers[i].slice.exp_i_sigma_v[0,1]).
   // This is done within 
   //   (*sliceList_itr)->update_transmission_function()

   // TODO: if ( <atoms in a slice have been altered> )
   //       <re-evaluate relevant slices, else reuse previous slices>
   
   // NOTE: Scatterers contained in myScatterers are pointed to from
   //       within members of slicesOfScatterers, so DO NOT release
   //       memory of the scatterers pointed to by myScatterers 
   //       members until all use of slicesOfScatterers is concluded.

   //////////////////////////////////////////////////////////////////
   // Initialize slice members
   //////////////////////////////////////////////////////////////////

   if ( mynode == rootnode && input_flag_debug )
      cout << "Initializing slice members ..." << endl;

   size_t sliceNumber = 0;//debug
   for (std::list< slice* >::iterator sliceList_itr 
         = slicesOfScatterers.begin();
         sliceList_itr != slicesOfScatterers.end();
         ++sliceList_itr)
   {
      ++sliceNumber;//debug

      (*sliceList_itr)->propagator_x_re = new double[Nx_local];
      (*sliceList_itr)->propagator_x_im = new double[Nx_local];
      (*sliceList_itr)->propagator_y_re = new double[Ny];
      (*sliceList_itr)->propagator_y_im = new double[Ny];
      // NOTE: this allocation is deleted within 
      //    delete_slices_from_list()
      // Evaluate the propagator member of each slice
      //if ( mynode == rootnode )  // debug
      //   cout << "Initializing propagator." << endl; // debug
      //if ( (*sliceList_itr)->update_propagator_flag )
      //{
      (*sliceList_itr)->update_propagator(
            lambda, 
            kx_local, Nx_local, ky, Ny,
            bwcutoff_t
            );
      //}
   
      // Evaluate the transmission function of each slice
      //if ( mynode == rootnode )  // debug
      //   cout << "Initializing transmission function." // debug
      //      << endl; // debug

      (*sliceList_itr)->exp_i_sigma_v 
         = fftw_alloc_complex( local_alloc_size_fftw );

      //if ( (*sliceList_itr)->update_transmission_function_flag )
      //{
      if (
            (*sliceList_itr)->update_transmission_function(
               lambda, gamma, 
               //1/(xperiod * yperiod), // (ab)^{-1}
               1/(xperiod_duped * yperiod_duped), // (ab)^{-1}
               bwcutoff_t,
               kx_local, Nx_local, 
               kx_joined, xx_joined, Nx, 
               ky, yy, Ny, 
               NxNy,
               //sqrtNxNy,
               input_flag_pap_tif,
               outFileName_prefix,
               sliceNumber,
               local_alloc_size_fftw,  
               local_0_start_fftw,
               mynode,
               rootnode,
               MPI_COMM_WORLD
               )
            == EXIT_FAILURE
         )
      {
         cout << "Exiting program." << endl;
         return EXIT_FAILURE;
      }
      //}
   }
   if ( mynode == rootnode && input_flag_debug )
      cout << "Evaluated slice members." << endl;
   // debug
   if ( mynode == rootnode && input_flag_debug )
   {
      sliceNumber = 0;
      for (std::list< slice* >::iterator sliceList_itr 
            = slicesOfScatterers.begin();
            sliceList_itr != slicesOfScatterers.end();
            ++sliceList_itr)
      {
         ++sliceNumber;//debug
         cout << "slice : " << sliceNumber << endl;
         cout << "exp_i_sigma_v[50][0,1] : "
            << (*sliceList_itr)->exp_i_sigma_v[50][0] 
            << "+ i" << (*sliceList_itr)->exp_i_sigma_v[50][1]
            << endl;
         cout << "propagator_x_re[50] : "
            << (*sliceList_itr)->propagator_x_re[50]
            << endl;
      }
   }
   // end debug

   //////////////////////////////////////////////////////////////////
   // Run the requested microscope simulations
   //////////////////////////////////////////////////////////////////


   // TODO: use the following flags to decide which simulations to 
   //        run :
   // unsigned int input_flag_adfstem_corrected = 0;        // 
   // unsigned int input_flag_adfstem_uncorrected = 0;      // 
   // unsigned int input_flag_bfctem_corrected = 0;         // 
   // unsigned int input_flag_bfctem_uncorrected = 0;       // 
   // unsigned int input_flag_fem = 0;
   
   // run adfstem()
   if ( 
         input_flag_adfstem_corrected 
         || input_flag_adfstem_uncorrected 
      )
   {
      double detector_inner_angle;
      double detector_outer_angle;
      // TODO: allow input to specify adf detector inner and outer 
      //       angles
      if ( VV >= 200000 )
      {
         // 40 mrad inner angle for 200 keV & 400 keV
         detector_inner_angle = 0.040;
         // TODO: make detector dimensions an input variable

         // 200 mrad outer angle for 200 keV & 400 keV
         detector_outer_angle = 0.200;
      }
      else
      {
         detector_inner_angle = 0.045; // 45 mrad for 100 keV
         detector_outer_angle = 0.200; // 200 mrad for 100 keV
      }

      // Create default raster positions if they weren't given by
      //  the parent call.
      if ( ! input_flag_raster_spacing )
      {
         if ( input_flag_fem )
            raster_spacing  = 10; // \AA // TODO: are these appropriate?
         else
            raster_spacing = 1.5; // \AA
      }

      // Create points at which the rastered stem beam will evaluated
      std::list<double> x_p;
      std::list<double> y_p;
      if (input_flag_adfstem_uncorrected || input_flag_adfstem_corrected)
      {

         // NOTE: do not increase number of raster points when using
         //        --dupe
         size_t number_of_raster_points_x
            = floor( xperiod/raster_spacing );  // * dupe_x ??? No.

         size_t number_of_raster_points_y
            = floor( yperiod/raster_spacing );  // * dupe_y ??? No.

         if ( mynode == rootnode && input_flag_debug )
            cout << "Determining STEM raster points with parameters" 
               << endl
               << " (xmin, ymin) : (" 
               << xmin << ", " << ymin << ")" << endl
               << " (xperiod, yperiod) : (" 
               << xperiod << ", " << yperiod << ")" << endl;

         assign_stem_raster_points(
               number_of_raster_points_x,
               number_of_raster_points_y,
               xperiod, xmin,
               yperiod, ymin,
               x_p, y_p
               );
         if ( mynode == rootnode && input_flag_debug )
            cout << "Determined STEM raster points." << endl
               << " (x_p.front(), y_p.front()) : ("
               << x_p.front() << ", " << y_p.front() << ")" << endl;
      }

      //   // TODO: insert units [radians] [angstrom]
      //   cout << "Running adfstem with parameters : " << endl;//debug
      //   cout << "   fem : " << input_flag_fem << endl
      //      << "   aberration correction : " << 0 << endl
      //      //<< "   initial_population : " 
      //      //<< initial_population << endl
      //      //<< "   xperiod, yperiod, zperiod : " 
      //      //<< xperiod << ", " << yperiod << ", " << zperiod << endl
      //      //<< "   xmin, ymin, zmin : " 
      //      //<< xmin << ", " << ymin << ", " << zmin << endl
      //      << "   Nx, Ny : " 
      //      << Nx << ", " << Ny << ", " 
      //      //<< Nz 
      //      << endl
      //      << "   raster position spacing x, y : "
      //      << xperiod / number_of_raster_points_x << ", "
      //      << yperiod / number_of_raster_points_y << ", " << endl
      //      << "   Accelerating voltage : " << VV << endl
      //      << "   Cs3 : " << Cs3 << endl
      //      << "   defocus : " 
      //      << defocus << endl
      //      << "   alpha max : " << alpha_max << endl
      //      //<< "   alpha max squared : " << alpha_max_sqr << endl
      //      << "   detector inner angle : " << detector_inner_angle 
      //      << endl
      //      << "   detector outer angle : " << detector_outer_angle 
      //      << endl
      //      << "   output file name prefix : " 
      //      << outFileName_prefix + "_adfstem_uncorr" << endl;

      //   adfstem(
      //         input_flag_fem,
      //         0,                         // aberration correction?
      //         0,
      //         input_flag_image_output,
      //         input_flag_debug,
      //         slicesOfScatterers,        // atom properties
      //         bwcutoff_pap,
      //         bwcutoff_t,
      //         x_p, y_p,
      //         xperiod, yperiod,
      //         kx_local, kx_joined, // sample domain, x direction
      //         Nx,  Nx_local, // number of samples
      //         ky, Ny, // sample domain, y direction
      //         Cs3, Cs5,               // spherical aberration
      //         defocus,                
      //         defocus_spread,   // TODO : implement adfstem with defocus_spread
      //         alpha_max_sqr,
      //         detector_inner_angle, detector_outer_angle,
      //         lambda, 
      //         outFileName_prefix + "_adfstem_uncorr",
      //         local_alloc_size_fftw,
      //         mynode, rootnode, totalnodes,
      //         MPI_COMM_WORLD
      //       );
      //   if ( mynode == rootnode ) 
      //      cout << "root node finished running adfstem" << endl; // debug
      //   //cout << "node " << mynode // debug
      //   //   << " Finished running adfstem " // debug
      //   //   << endl; // debug
      //}

      if ( mynode == rootnode && input_flag_debug )
      {
         cout << "Optimal adfstem parameters without aberration"
            << " compensation with fixed Cs3 (Kirkland 2016):" << endl;
         cout << "defocus : " << scherzer_defocus_adfstem_uncorrected(
               lambda, Cs3) << endl;
         cout << "alpha_max : " 
            << scherzer_alphamax_adfstem_uncorrected( lambda, Cs3) 
            << endl;
         if ( input_flag_cs5 )
         {
            cout << "Optimal adfstem parameters with aberration"
               << " compensation with fixed Cs5 :" << endl;
            cout << "defocus : " 
               << scherzer_defocus_adfstem_correctedtoCs5( lambda, Cs5) 
               << endl;
            cout << "alpha_max : " 
               << scherzer_alphamax_adfstem_correctedtoCs5( lambda, Cs5) 
               << endl;
         }

         if ( input_flag_adfstem_uncorrected )
            cout << "approx. min. resolution [AA], adfstem uncorrected : " 
               << scherzer_approxresolution_adfstem_uncorrected(
                  lambda, Cs3) << endl;

         if ( input_flag_adfstem_corrected )
            cout << "approx. min. resolution [AA], adfstem corrected : " 
               << scherzer_approxresolution_adfstem_correctedtoCs5(
                  lambda, Cs5) << endl;
      }


      // NOTE: input_flag_adfstem_corrected and 
      //       input_flag_adfstem_uncorrected are not necessarily 
      //       opposites.
      //       We should allow for both corrected and uncorrected 
      //       simulations within the same run.


      // TODO: modify the defocus by the sample thickness to ensure
      //       that sample defocus is w.r.t. the bottom of the sample.
      if ( input_flag_adfstem_uncorrected )
      {
         if ( input_flag_scherzer_defocus  )
            defocus = scherzer_defocus_adfstem_uncorrected( 
                  lambda, Cs3
                  );
         if ( input_flag_scherzer_alphamax  )
            alpha_max = scherzer_alphamax_adfstem_uncorrected(
                  lambda, Cs3
                  );

         //alpha_max_sqr = alpha_max * alpha_max;

         // TODO: insert units [radians] [angstrom]
         if ( mynode == rootnode && input_flag_debug )
         {
            cout << "Running adfstem with parameters : " << endl;//debug
            cout << "   fem : " << input_flag_fem << endl
               << "   aberration correction : " 
               << 0 << endl
               << "   Nx, Ny : " 
               << Nx << ", " << Ny << ", " 
               << endl
               << "   Accelerating voltage : " << VV << endl
               << "   Cs3 : " << Cs3 << ", " << endl
               << "   defocus : " << defocus << endl
               << "   alpha max : " << alpha_max << " [radians]" << endl
               << "   detector inner angle : " 
               << detector_inner_angle << endl
               << "   detector outer angle : " 
               << detector_outer_angle << endl
               << "   output file name prefix : " 
               << outFileName_prefix + "_adfstem_uncorr" << endl;
         }

         adfstem(
               input_flag_fem,
               input_flag_gt17,
               input_flag_d1, input_flag_d2, input_flag_d3, input_flag_d4,
               0,                         // aberration correction?
               0, // calculate using complex average over real space?
               input_flag_image_output,
               input_flag_netcdf_images,
               input_flag_netcdf_variance,
               input_flag_debug,
               slicesOfScatterers,        // atom properties
               bwcutoff_pap,
               bwcutoff_t,
               azimuthal_binning_size_factor,
               xmin, ymin,
               x_p, y_p,
               xperiod, yperiod,
               xperiod_duped, yperiod_duped,
               xx_local, xx_joined, yy,
               kx_local, kx_joined, // sample domain, x direction
               Nx,  Nx_local, // number of samples
               ky, Ny, // sample domain, y direction
               Cs3, Cs5,               // spherical aberration
               defocus,                
               defocus_spread,   // TODO : implement adfstem with defocus_spread
               alpha_max,
               detector_inner_angle, detector_outer_angle,
               lambda, 
               outFileName_prefix + "_adfstem_uncorr",
               local_alloc_size_fftw,
               local_0_start_fftw,
               mynode, rootnode, totalnodes,
               MPI_COMM_WORLD
             );
         
         if ( mynode == rootnode && input_flag_debug )
            cout << "root node finished running adfstem"
               << " without aberration correction" << endl; // debug
      }

      if ( input_flag_adfstem_corrected )
      {
         if ( input_flag_scherzer_defocus )
            defocus 
               = scherzer_defocus_adfstem_correctedtoCs5(lambda,Cs5);

         if ( input_flag_scherzer_alphamax )
            alpha_max 
               = scherzer_alphamax_adfstem_correctedtoCs5(lambda,Cs5);

         if ( input_flag_scherzer_cs3 )
            Cs3_corr = scherzer_Cs3_adfstem_correctedtoCs5(lambda, Cs5);
         else
            Cs3_corr = Cs3;

         //alpha_max_sqr = alpha_max * alpha_max;

         // TODO: insert units [radians] [angstrom]
         if ( mynode == rootnode && input_flag_debug )
         {
            cout << "Running adfstem with parameters : " << endl;//debug
            cout << "   fem : " << input_flag_fem << endl
               << "   aberration correction : " 
               << 1 << endl
               << "   Nx, Ny : " 
               << Nx << ", " << Ny << ", " << endl
               << "   Accelerating voltage : " << VV << endl
               << "   Cs3 : " << Cs3_corr << ", " << endl
               << "   Cs5 : " << Cs5 << endl 
               << "   defocus : " << defocus << endl
               << "   defocus spread : " << defocus_spread << endl
               << "   alpha max : " << alpha_max << " [radians]" << endl
               << "   detector inner angle : " << detector_inner_angle 
               << endl
               << "   detector outer angle : " << detector_outer_angle 
               << endl
               << "   output file name prefix : " 
               << outFileName_prefix + "_adfstem_corr" << endl;
         }

         adfstem(
               input_flag_fem,
               input_flag_gt17,
               input_flag_d1, input_flag_d2, input_flag_d3, input_flag_d4,
               1, // aberration correction?
               0, // calculate using complex average over real space?
               input_flag_image_output,
               input_flag_netcdf_images,
               input_flag_netcdf_variance,
               input_flag_debug,
               slicesOfScatterers,        // atom properties
               bwcutoff_pap,
               bwcutoff_t,
               azimuthal_binning_size_factor,
               xmin, ymin,
               x_p, y_p,
               xperiod, yperiod,
               xperiod_duped, yperiod_duped,
               xx_local, xx_joined, yy,
               kx_local, kx_joined, // sample domain, x direction
               Nx,  Nx_local, // number of samples
               ky, Ny, // sample domain, y direction
               Cs3_corr, Cs5,               // spherical aberration
               defocus,                
               defocus_spread,   // TODO : implement adfstem with defocus_spread
               alpha_max,
               detector_inner_angle, detector_outer_angle,
               lambda, 
               outFileName_prefix + "_adfstem_corr",
               local_alloc_size_fftw,
               local_0_start_fftw,
               mynode, rootnode, totalnodes,
               MPI_COMM_WORLD
             );
         
         if ( mynode == rootnode && input_flag_debug )
            cout << "root node finished running adfstem"
               << " with aberration correction" << endl; // debug
      }

   }
   if ( 
         input_flag_bfctem_corrected 
         || input_flag_bfctem_uncorrected 
      )
   {
      cout << "Error, bfctem() has not yet been integrated into this"
        << " version of the program." << endl;
      // run bfctem()
   }
   //////////////////////////////////////////////////////////////////
   // CTEM
   //bfctem_ms_uncorrected(
   //      myScatterers,                 // atom properties
   //      initial_population,             // total number of atoms
   //      xperiod, yperiod, zperiod,    // periodic boundary widths
   //      xmin, ymin, zmin,             // lowest values of x y z domains
   //      Nx, Ny, Nz,                   // resolution
   //      slice_locations,  
   //      VV, Cs3, defocus,
   //      alpha_max_sqr,  // microscope parameters
   //      outFileName_prefix,              
   //      threads_ok, nthreads, mynode, rootnode, totalnodes,
   //      MPI_COMM_WORLD
   //    );
   //bfctem_ms_uncorrected_autoslice(
   //      myScatterers,                 // atom properties
   //      initial_population,             // total number of atoms
   //      xperiod, yperiod, zperiod,    // periodic boundary widths
   //      xmin, ymin, zmin,             // lowest values of x y z domains
   //      Nx, Ny, Nz,                   // resolution
   //      VV, Cs3, defocus,
   //      alpha_max_sqr,  // microscope parameters
   //      outFileName_prefix,              
   //      threads_ok, nthreads, mynode, rootnode, totalnodes,
   //      MPI_COMM_WORLD
   //    );



   //////////////////////////////////////////////////////////////////
   // Clean up
   //////////////////////////////////////////////////////////////////

   delete_slices_from_list( slicesOfScatterers );
     // ^^ calls fftw_free() on paps, exp_i_sigma_v, and deletes 
     //  propagator_x_re, propagator_x_im, propagator_y_re,
     //  and propagator_y_im
   
   //for (std::list< slice* >::iterator sliceList_itr 
   //      = slicesOfScatterers.begin();
   //      sliceList_itr != slicesOfScatterers.end();
   //      ++sliceList_itr)
   //{
   //   delete[] (*sliceList_itr)->propagator_x_re;
   //   delete[] (*sliceList_itr)->propagator_x_im;
   //   delete[] (*sliceList_itr)->propagator_y_re;
   //   delete[] (*sliceList_itr)->propagator_y_im;

   //   fftw_free( (*sliceList_itr)->exp_i_sigma_v )
   //}

   delete[] q;
   //delete[] Z_array;
   if ( input_flag_dupe && input_flag_fem 
        &&
        (dupe_x > 1 || dupe_y > 1 || dupe_z > 1)
        &&
        (dupe_x > 0 && dupe_y > 0 && dupe_z > 0)
      )
   {
      delete[] q_duped;
      //delete[] Z_array_duped;
   }

   delete[] kx_local;
   delete[] kx_joined;
   delete[] xx_joined;

   delete[] ky;
   delete[] yy;

   for (size_t i = 0; i < myScatterers.size() ; i++) 
   {
   // NOTE: Scatterers contained in myScatterers are pointed to from
   //       within members of slicesOfScatterers, so DO NOT release
   //       memory of the scatterers pointed to by myScatterers 
   //       members until all use of slicesOfScatterers is concluded.
      myScatterers.pop_back();
   }

   // Save wisdom to a file if it wasn't successfully opened above
   if ( mynode == rootnode )
      if ( ! flag_wisdomFile )
      {
         if ( input_flag_debug )
            cout << "Exporting wisdom to file: " 
               << wisdom_filename << endl;
         if ( 
               ! fftw_export_wisdom_to_filename( 
                  wisdom_filename.c_str()
                  )
            )
         {
            cerr << "Error - could not export wisdom to file: "
              << wisdom_filename << endl;
         }
      }

   fftw_mpi_cleanup();
   MPI_Finalize();

   return EXIT_SUCCESS;
}


