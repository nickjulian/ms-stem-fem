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
// File: ftem_lammps.cpp
// Purpose: 
//    To control execution of a multislice FTEM simulation.
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

//#include <algorithm> // transform()
//#include <cctype> // tolower() UNARY OPERATOR

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
#include "to_string.hpp"

// lammps library headers
#include "lammps.h"
//#include "lmptype.h"
#include "input.h"
#include "atom.h"
#include "comm.h"
#include "error.h"
#include "memory.h"
#include "library.h"

typedef int tagint;

using std::cout;
using std::endl;
using std::string;
using std::istringstream;
using std::setw;
using std::setprecision;
using namespace TEM_NS;

//int threads_ok; // global by request of fftw.org/doc/ ... section 6.11
#define PRINT_USAGE cout << "Usage: " << argv[0] << " <options>" << endl << "OPTIONS : " << endl << "   --parameter_file <file name>" << endl << "      Contents of this file may specify any of the other arguments, but will" << endl << "      be superceded by those on the command line." << endl << "   -m <samples_x> <samples_y> <VV>" << endl << "      number of samples in the 2-D discretization and microscope voltage" << endl << "      WARNING: if samples_x != samples_y, or if they are odd, bad data" << endl << "               may creep from the edges of reciprocal space boundaries" << endl << "   [--scherzer_defocus]" << endl << "      calculate and use scherzer focus conditions" << endl << "   [--scherzer_alphamax]" << endl << "      calculate and use scherzer focus conditions" << endl << "   [--scherzer_cs3]" << endl << "      calculate and use scherzer Cs3 conditions" << endl << "   [--defocus <defocus>]" << endl << "      required if not using --scherzer_defocus" << endl << "   [--alphamax <alpha max>]" << endl << "      required if not using --scherzer_alphamax" << endl << "   [--spread <defocus spread> <condenser_illumination_angle>]" << endl << "      applicable only if using aberration correction" << endl << "   [--cs3 <third order spherical aberration>]" << endl << "      ignored if using --scherzer with aberration correction" << endl << "   [--cs5 <fifth order spherical aberration>]" << endl << "     applicable only if using aberration correction" << endl << "   [--detectorangles <inner_angle> <outer_angle>]" << endl << "     specify in radians; applicable only if using adfstem;" << endl << "     default is 0.04 0.2 for V >= 200kV, else 0.045 and 0.2" << endl << "        Note: bandwidth limit will reduce detected signal if" << endl << "              the outer angle is greater than it, being" << endl << "              (1/3)*(number of x or y samples)/(x or y realspace period)" << endl << "   [--rasterspacing <raster_spacing>]" << endl << "      units: angstroms; default is 1.5 for STEM or 10 for STEM-FEM" << endl << "   [--dupe <dupe_x> <dupe_y> <dupe_z>]" << endl << "      Periodically instantiate the given scatterers dupe_x, dupe_y, and " << endl << "      dupe_z times in respective directions; fem only" << endl << "   -a <lammps or xyz style position input file>" << endl << "   --output_prefix <output file prefix>" << endl << "   [--paptif]" << endl << "      output images of projected atomic potentials" << endl << "   [--adfstemcorrfem]" << endl << "      simulate fluctuation microscopy using aberration corrected" << endl << "      adfstem mode" << endl << "   [--adfstemuncorrfem]" << endl << "      simulate fluctuation microscopy using" << endl << "      adfstem mode without aberration correction" << endl << "   [--adfstemcorr]" << endl << "      simulate aberration corrected adfstem" << endl << "   [--adfstemuncorr]" << endl << "      simulate adfstem mode without aberration correction" << endl << "   [--dr <azimuthal_binning_size_factor>]" << endl << "      prefactor of sqrt(dx^2+dy^2) in azimuthal integration bin size" << endl << "   [--minslice <minSliceThickness>" << endl << "      minimum slice thickness in Angstroms, default is 1 " << endl << "   [--images]" << endl << "      generate and save images" << endl << "   [--V_omega]" << endl << "      FEM: normalized variance of the annular mean" << endl << "   [--Vbar_r]" << endl << "      FEM: mean of normalized variances of rings" << endl << "   [--V_re]" << endl << "      FEM: normalized variance of ring ensemble" << endl << "   [--Omega_vimage]" << endl << "      FEM: annular mean of variance image" << endl << "   [--V_gt17]" << endl << "      FEM: ratio of annular means of mean diffraction squared and the square of" << endl << "      mean diffraction" << endl << "   [--V_rbar]" << endl << "      FEM: angular normalized variance of the average diffraction" << endl << "   [--correlograph]" << endl << "      FEM: 2-D average of autocorrelation along azimuthal phi axis" << endl << "   [--correlograph_everyimage]" << endl << "      FEM: 2-D tif image output of autocorrelation along azimuthal phi axis for" << endl << "      every STEM raster point" << endl << "   [--correlograph_everytxt]" << endl << "      FEM: 2-D txt output of autocorrelation along azimuthal phi axis for every" << endl << "      STEM raster point" << endl << "   [--correlograph_variance]" << endl << "      FEM: 2-D variance of autocorrelation along azimuthal phi axis" << endl << "   [--debug]" << endl << "      enable verbose debug output to files and stdout" << endl ;
//"   [--netcdfimages]" << endl << "      save images as netcdf files" << endl << "   [--netcdfvariance]" << endl << "      save 1-D variance as netcdf files" << endl << 
//<< endl << "   [--bfctemcorr]" << endl << "      simulate bright field TEM with aberration correction" << endl << "   [--bfctemuncorr] simulate bright field TEM without aberration correction" 


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
   double detector_inner_angle;
   double detector_outer_angle;

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
   double alpha_max;
   string output_prefix; // = argv[2];
   //string input_position_lammps_file;

   // STEM and FEM raster point separation in both x and y directions
   double raster_spacing;

   //////////////////////////////////////////////////////////////////
   // Specimen domain
   //////////////////////////////////////////////////////////////////
   double xmin, xmax;
   double ymin, ymax;
   double zmin, zmax;
   double xperiod;
   double yperiod;
   double zperiod;

   double* q;
   double* q_duped;
   unsigned int* Z_array; // couldn't use size_t because it must have an 
                          //  MPI type
   unsigned int* uniqueZs; // same content as Z_vector
   uniqueZs = NULL; // to ensure conditional deletion doesn't break
   // TODO: replace Z_vector with uniqueZs regardless of model file fmt
   unsigned int numberOfSpecies;
   unsigned int initial_population;
   unsigned int duped_population;

   //std::vector< double > slice_locations;

   //////////////////////////////////////////////////////////////////
   // read input parameters
   //////////////////////////////////////////////////////////////////
   std::vector<string> args( argv, argv + argc );
   string model_file_name;
   //unsigned int failflag = 0;

   std::vector< double > mtf_1D;
   std::vector< double > mtf_domain_sqr;
   double mtf_resolution;

   std::list<double> x_p;
   std::list<double> y_p;

   //double azimuthal_binning_size_factor = 1.0;

   double azimuthal_binning_size_factor = 0.70710678118654746;// 1/sqrt(2)

   double minSliceThickness = 1.0;

   input_flags flags;

   string lammps_preTEM_file_name;
   unsigned int lammps_TEM_steps;
   unsigned int lammps_TEM_samples;

   //double tmp_double;// for accepting input from operator>>() and 
   //                      pushing onto a vector<double> 

   if ( read_cmdline_options( 
         args,
         model_file_name,
         flags,
         output_prefix,
         Nx,
         Ny,
         VV,
         defocus,
         alpha_max,
         defocus_spread,
         condenser_illumination_angle,
         Cs3,
         Cs5,
         detector_inner_angle,
         detector_outer_angle,
         mtf_1D,
         mtf_domain_sqr, // not squared until later
         mtf_resolution,
         raster_spacing,
         azimuthal_binning_size_factor,
         minSliceThickness,
         dupe_x,
         dupe_y,
         dupe_z,
         lammps_preTEM_file_name,
         lammps_TEM_steps,
         lammps_TEM_samples,
         mynode,
         rootnode,
         MPI_COMM_WORLD
         ) == EXIT_FAILURE 
   )
   {
      if ( flags.debug && mynode == rootnode) 
         cerr << "Failed to read command line arguments" << endl;
      //failflag = 1 ;
      flags.fail = 1;
   }

   //MPI_Bcast( &flags.fail, 1, MPI_UNSIGNED, rootnode, MPI_COMM_WORLD );
   //if ( flags.fail) 
   //{
   //   if ( mynode == rootnode)
   //   {
   //      cout << "Exiting ms-stem-fem"
   //         << endl;
   //      PRINT_USAGE // macro defined above
   //   }
   //   fftw_mpi_cleanup();
   //   MPI_Finalize();
   //   return EXIT_FAILURE;
   //}

   NxNy = Nx * Ny;
   NxNy_sqr = NxNy * NxNy;
   if ( Nx == Ny ) sqrtNxNy =  Nx ;   
   else sqrtNxNy = sqrt( NxNy );

   size_t mtf_domain_size; // = mtf_domain.size();
   if ( flags.mtf_file && flags.mtf_resolution )
   {
      mtf_domain_size = mtf_domain_sqr.size();
      for ( size_t mtf_domain_idx = 0; 
            mtf_domain_idx < mtf_domain_size;
            ++mtf_domain_idx)
      { // mtf_file was in units of Nyquist freq. = 0.5 * sampling freq.
         mtf_domain_sqr[mtf_domain_idx]
            = pow(
               (0.5 * mtf_resolution * mtf_domain_sqr[mtf_domain_idx]),
                  2);
      }
   }

   // Read model file
   unsigned int read_xyz_flag = 0;
   //if ( flags.a )
   //{
   if ( model_file_name.substr( 
            model_file_name.find_last_of(".") + 1 
            ) == "xyz" )
   {
      if ( mynode == rootnode && flags.debug)
      {
         cout << "Reading xyz format file: "
           << model_file_name  << endl; // debug
      }
      if ( 
            read_position_xyz_file(
               model_file_name,   // only valid on root node
               q, // will be allocated within this call
               Z_array, // will be allocated within this call
               initial_population,
               xmin, ymin, zmin,
               xperiod, yperiod, zperiod,
               flags.debug,
               mynode, rootnode, MPI_COMM_WORLD
               )
         )
      {
         if ( mynode == rootnode )
         {
            cerr << "Error, file could not be read : " 
               << model_file_name << endl;
         }
         flags.fail = 1;
      }
   }
   else
   {
      if ( mynode == rootnode && flags.debug )
      {
         cout << "Reading lammps data format file: "
           << model_file_name  << endl; // debug
      }
      if(  
            read_position_lammps_data_file(
               model_file_name,   // only valid on root node
               q,       // will be allocated within this call
               Z_array, // will be allocated within this call
               numberOfSpecies,
               uniqueZs, // will be allocated within this call
               initial_population,
               xmin, ymin, zmin,
               xperiod, yperiod, zperiod,
               flags.debug,
               mynode, rootnode, MPI_COMM_WORLD
               ) != EXIT_SUCCESS
        )
      {
         if ( mynode == rootnode )
         {
            cerr << "Error, file could not be read : " 
               << model_file_name << endl;
         }
         flags.fail = 1;
      }
   }
   cout << "node: " << mynode << ",2, flags.fail: " << flags.fail << endl;

   if ( 
         flags.bfctem_corrected 
         || flags.bfctem_uncorrected 
      )
   {
      cout << "Error, bfctem() has not yet been integrated into"
            << " this version of the program." << endl;
      flags.fail = 1;
   }


   if ( check_runtime_flags(
         flags,
         args[0],
         mynode,
         rootnode
         ) == EXIT_FAILURE)
   {
      flags.fail = 1;
   }

   MPI_Bcast( &flags.fail, 1, MPI_UNSIGNED, rootnode, MPI_COMM_WORLD );
   if ( flags.fail != 0)
   {
      if ( mynode == rootnode)
      {
         cout << "Exiting ms-stem-fem due to failure." << endl;
         PRINT_USAGE // macro defined above
      }
      fftw_mpi_cleanup();
      MPI_Finalize();
      return EXIT_FAILURE;
   }

   //////////////////////////////////////////////////////////////////
   // read lammps_preTEM_file 
   //////////////////////////////////////////////////////////////////
   // clear positions and element types from ms-stem-fem parameters
   if ( flags.a ) 
   {
      delete[] q;
      delete[] Z_array;
   }

   std::list< string > lammps_preTEM_commands;
   int lammps_preTEM_size;
   int lammps_preTEM_line_length;
   char lammps_preTEM_line_c_str[1024];
   if ( mynode == rootnode )
   {
      ifstream lammps_preTEM_file( lammps_preTEM_file_name.c_str() );
      string lammps_preTEM_line;

      while (  getline( lammps_preTEM_file, lammps_preTEM_line )
                     && lammps_preTEM_file.good()) 
      {
         lammps_preTEM_commands.push_back( lammps_preTEM_line );
      }

      lammps_preTEM_size = lammps_preTEM_commands.size();
      MPI_Bcast( &lammps_preTEM_size, 1, 
                  MPI_INT, rootnode, MPI_COMM_WORLD);
      std::list< string >::iterator
         lammps_preTEM_cmd_itr = lammps_preTEM_commands.begin();
      for ( size_t ii=0; ii < lammps_preTEM_size; ++ii )
      {
         strcpy( lammps_preTEM_line_c_str, 
                  (*lammps_preTEM_cmd_itr).c_str() );
         lammps_preTEM_line_length = strlen(lammps_preTEM_line_c_str)+1;
         MPI_Bcast( &lammps_preTEM_line_length, 1, 
                     MPI_INT, rootnode, MPI_COMM_WORLD);
         MPI_Bcast( lammps_preTEM_line_c_str, lammps_preTEM_line_length, 
                     MPI_CHAR, rootnode, MPI_COMM_WORLD);
         ++lammps_preTEM_cmd_itr;
      }
      lammps_preTEM_file.close();
   } // rootnode
   else
   {
      MPI_Bcast( &lammps_preTEM_size, 1, 
                  MPI_INT, rootnode, MPI_COMM_WORLD);
      for ( size_t ii=0; ii < lammps_preTEM_size; ++ii )
      {
         MPI_Bcast( &lammps_preTEM_line_length, 1, 
                     MPI_INT, rootnode, MPI_COMM_WORLD);
         MPI_Bcast( lammps_preTEM_line_c_str, lammps_preTEM_line_length, 
                     MPI_CHAR, rootnode, MPI_COMM_WORLD);

         lammps_preTEM_commands.push_back( 
               string( lammps_preTEM_line_c_str ) );
      }
   }
   // debug
   //cout << "node " << mynode << " lammps_preTEM_line contents: " << endl;
   //for ( std::list< string >::iterator 
   //      lammps_preTEM_cmd_itr = lammps_preTEM_commands.begin(); 
   //      lammps_preTEM_cmd_itr != lammps_preTEM_commands.end(); 
   //      ++lammps_preTEM_cmd_itr )
   //{
   //   cout << *lammps_preTEM_cmd_itr << endl;
   //}
   //cout << "node " << mynode << " end of lammps_preTEM_line contents." << endl;
   // end debug

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

   if ( mynode == rootnode && flags.debug )
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
   // Running lammps preTEM file commands
   //////////////////////////////////////////////////////////////////

   //bool lammps_node;
   //if ( mynode < lammps_node_count ) lammps_node = true; 
   //// check truth of lammps_node before each lammps command

   LAMMPS_NS::LAMMPS *lmp = NULL;
   lmp = new LAMMPS_NS::LAMMPS(0, NULL, MPI_COMM_WORLD );
   
   std::vector< string > element_name_list;
   size_t found_pair_coeff;
   size_t found_pair_style;
   size_t found_pair_style_sw;
   size_t found_pair_style_tersoff;
   size_t found_pair_style_eam;

   bool tersoff_flag = false;
   bool sw_flag = false;
   bool eam_flag = false;
   for ( std::list< string >::iterator 
         lammps_preTEM_cmd_itr = lammps_preTEM_commands.begin(); 
         lammps_preTEM_cmd_itr != lammps_preTEM_commands.end(); 
         ++lammps_preTEM_cmd_itr )
   {
      found_pair_style
         = (*lammps_preTEM_cmd_itr).find("pair_style");
      if ( found_pair_style != string::npos ) 
      {
         found_pair_style_sw
            = (*lammps_preTEM_cmd_itr).find("sw");

         if ( found_pair_style_sw != string::npos ) 
            sw_flag = true;

         found_pair_style_tersoff
            = (*lammps_preTEM_cmd_itr).find("tersoff");

         if ( found_pair_style_tersoff != string::npos ) 
         {
            tersoff_flag = true;
         }

         found_pair_style_eam
            = (*lammps_preTEM_cmd_itr).find("eam");

         if ( found_pair_style_eam != string::npos ) 
         {
            eam_flag = true;
         }
      }
      
      found_pair_coeff = (*lammps_preTEM_cmd_itr).find("pair_coeff");
      if ( found_pair_coeff != string::npos )
      {
         // read element names from pair_coeff line 
         //  and associate their Z to the lammps type number
         // NOTE: the following is tailored to tersoff and 
         //  stillinger-weber pair_coeff commands and need to be 
         //  altered if using a different type
         istringstream pair_coeff_line(
               (*lammps_preTEM_cmd_itr).substr(found_pair_coeff)
               );
         string tmp_str;
         if ( tersoff_flag || sw_flag || eam_flag )
         {
            pair_coeff_line >> tmp_str >> tmp_str >> tmp_str; 
            pair_coeff_line >> tmp_str;   // file name
            while ( pair_coeff_line.good()  )
            {
               pair_coeff_line >> tmp_str;
               element_name_list.push_back( tmp_str );
            }
         }
         else
         {
            if (lmp->comm->me == rootnode)
              lmp->error->warning(FLERR,
                    "Error in reading lammps_preTEM_file_name : only tersoff, stillinger-weber, and eam pair_coeff commands currently allowed");
         }
         //if ( flags.debug )
         //   for ( size_t ii=0; ii < element_name_list.size(); ++ii)
         //   {
         //      cout << "element_name_list[" << ii << "] : " 
         //         << element_name_list[ii] << endl;
         //   }
      }
      strcpy( lammps_preTEM_line_c_str, (*lammps_preTEM_cmd_itr).c_str());
      lammps_command(lmp, lammps_preTEM_line_c_str);
   }

   // create the list of unique Zs to match with lammps type numbers
   // NOTE: this is only needed when the model is an xyz format file,
   //    otherwise it duplicates uniqueZs and numberOfSpecies
   std::vector<unsigned int> Z_vector;
   unsigned int Z_vector_size;
   string element_name;
   for ( std::vector<string>::iterator 
         itr = element_name_list.begin();
         itr != element_name_list.end();
         ++itr)
   {
      //transform( (*itr).begin(), (*itr).end(),
      //      (*itr).begin(), (int(*)(int))tolower );
      element_name = *itr;
      //cout << "element_name : " << element_name << endl;
      transform( element_name.begin(), element_name.end(), 
                 element_name.begin(), (int(*)(int))tolower );
      //cout << "element_name (lower case) : " << element_name << endl;

      Z_vector.push_back(
            atom_element_abbrev_to_Z( element_name )
            );
   }
   if ( flags.debug )
   {
      cout << "Unique list of Zs: ";
      for (size_t ii=0; ii < Z_vector.size(); ++ii)
         cout << Z_vector[ii] << ", ";
      cout << endl;
   }

   /////////////////////////////////////////////////////////////////////
   // copy domain info from lammps to ms-stem-fem variables
   /////////////////////////////////////////////////////////////////////
   unsigned int lammps_population; 
   lammps_population = static_cast<unsigned int>( lmp->atom->natoms );
   initial_population = lammps_population ;
   duped_population = dupe_x * dupe_y * dupe_z * lammps_population;
   
   if ( flags.adfstem_uncorrected )
   {
      if ( flags.scherzer_defocus  )
         defocus = scherzer_defocus_adfstem_uncorrected( 
               lambda, Cs3
               );
      if ( flags.scherzer_alphamax  )
         alpha_max = scherzer_alphamax_adfstem_uncorrected(
               lambda, Cs3
               );
      // TODO: insert units [radians] [angstrom]
      if ( mynode == rootnode && flags.debug )
      {
         cout << "Running adfstem with parameters : " << endl;//debug
         cout << "   fem : " << flags.fem << endl
            << "   aberration correction : " 
            << "no" << endl
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
            << output_prefix + "_adfstem_uncorr" << endl;
      }

      flags.aberration_correction = 0;
      flags.complex_realspace_sum = 0;
   } else if ( flags.adfstem_corrected )
   {
      if ( flags.scherzer_defocus )
         defocus 
            = scherzer_defocus_adfstem_correctedtoCs5(lambda,Cs5);

      if ( flags.scherzer_alphamax )
         alpha_max 
            = scherzer_alphamax_adfstem_correctedtoCs5(lambda,Cs5);

      if ( flags.scherzer_cs3 )
         Cs3_corr = scherzer_Cs3_adfstem_correctedtoCs5(lambda, Cs5);
      else
         Cs3_corr = Cs3;
      if ( mynode == rootnode && flags.debug )
      {
         cout << "Running adfstem with parameters : " << endl;//debug
         cout << "   fem : " << flags.fem << endl
            << "   aberration correction : " 
            << "yes" << endl
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
            << output_prefix + "_adfstem_corr" << endl;
      }

      flags.aberration_correction = 1;
      flags.complex_realspace_sum = 0;
   }

   if ( mynode == rootnode && flags.debug )
   {
      cout << endl
         << "  minSliceThickness : " << minSliceThickness << endl
         << "  azimuthal_binning_size_factor : " 
         << azimuthal_binning_size_factor << endl
         << "  azimuthal_binning_size_factor : " 
         << azimuthal_binning_size_factor << endl
         << "  dupe : " << dupe_x << ", " << dupe_y << ", " << dupe_z 
         << endl
         << "  raster_spacing : " << raster_spacing << endl
         << "  Cs3 : " << Cs3 << endl
         << "  Cs5 : " << Cs5 << endl
         << "  defocus_spread : " << defocus_spread << endl
         << "  condenser_illumination_angle : " 
         << condenser_illumination_angle << endl
         << "  alpha_max : " << alpha_max << endl
         << "  defocus : " << defocus << endl
         << "  model_file_name : " << model_file_name << endl
         << "  Nx : " << Nx << endl
         << "  Ny : " << Ny << endl
         << "  output_prefix : " << output_prefix << endl;
   }
   if ( mynode == rootnode && flags.debug 
         && (  flags.adfstem_corrected || flags.adfstem_uncorrected )
         )
   {
      cout << "Optimal adfstem parameters without aberration"
         << " compensation with fixed Cs3 (Kirkland 2016):" << endl;
      cout << "defocus : " << scherzer_defocus_adfstem_uncorrected(
            lambda, Cs3) << endl;
      cout << "alpha_max : " 
         << scherzer_alphamax_adfstem_uncorrected( lambda, Cs3) 
         << endl;
      if ( flags.cs5 )
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

      if ( flags.adfstem_uncorrected )
         cout << "approx. min. resolution [AA], adfstem uncorrected : " 
            << scherzer_approxresolution_adfstem_uncorrected(
               lambda, Cs3) << endl;

      if ( flags.adfstem_corrected )
         cout << "approx. min. resolution [AA], adfstem corrected : " 
            << scherzer_approxresolution_adfstem_correctedtoCs5(
               lambda, Cs5) << endl;
   }


   //////////////////////////////////////////////////////////////////
   // Instantiate local variables for fftw_mpi
   //////////////////////////////////////////////////////////////////
   ptrdiff_t Nx_local;
   ptrdiff_t local_alloc_size_fftw;
   ptrdiff_t local_0_start_fftw;

   MPI_Comm comm;
   comm = MPI_COMM_WORLD;

   local_alloc_size_fftw
      = fftw_mpi_local_size_2d(     // fftw will be 2-D here
            Nx, Ny,
            MPI_COMM_WORLD,
            &Nx_local,
            &local_0_start_fftw );

   //if ( flags.debug )
   //{
   //   cout << "node " << mynode 
   //      << ", (Nx_local, local_0_start_fftw, local_alloc_size_fftw)"
   //      << " : (" 
   //      << Nx_local << ", " << local_0_start_fftw << ", " 
   //      << local_alloc_size_fftw << ")"<< endl; // debug
   //}

   if ( local_alloc_size_fftw != Nx_local * Ny ) 
   {
      // TODO: fix
      //if ( flags.debug )
      cout << "WARNING: local_alloc_size_fftw != Nx_local * Ny " 
         << ", this may yield bad data at the edges of images" << endl;
      //return( EXIT_FAILURE );
   }

   //////////////////////////////////////////////////////////////////
   // Gather the values of Nx_local and local_alloc_size_fftw into 
   //    arrays on root node
   //////////////////////////////////////////////////////////////////
   int* Nx_strides;
   int* Nx_displacements;
   int* psi_mag_strides;
   int* psi_mag_displacements;
   Nx_strides = new int[ totalnodes ];
   Nx_displacements = new int[ totalnodes ];
   // psi_mag_strides and displacements are required for MPI_Allgatherv
   psi_mag_strides = new int[ totalnodes ];
   psi_mag_displacements = new int[ totalnodes ];

   int Nx_local_int;
   Nx_local_int = (int) Nx_local;   // problem here if domain is gigantic
   MPI_Allgather(
         &Nx_local_int,   // const void* sendbuf,
         1,                // int sendcount,
         MPI_INT,          // MPI_Datatype sendtype,
         Nx_strides,       // void *recvbuf,
         1,//totalnodes,       // int recvcount
         MPI_INT,     // MPI_Datatype recvtype
         //rootnode,         // int root
         MPI_COMM_WORLD    // MPI_Comm comm
         );
   if ( flags.debug )
   {
      if ( mynode == rootnode )
      {
         cout << "rootnode Nx_strides[]: ";
         for (size_t ii=0; ii < totalnodes; ++ii)
            cout << Nx_strides[ii] << ", ";
         cout << endl;
      }
   }

   //int local_alloc_size_fftw_int;
   //local_alloc_size_fftw_int = (int) local_alloc_size_fftw;
   int Nx_local_Ny_int = (int) (Nx_local * Ny);
   MPI_Allgather(
         &Nx_local_Ny_int,//&local_alloc_size_fftw_int, // *sendbuf,
         1,                // int sendcount,
         MPI_INT,          // MPI_Datatype sendtype,
         psi_mag_strides,  // void *recvbuf,
         1,//totalnodes,       // int recvcount
         MPI_INT,          // MPI_Datatype recvtype
         //rootnode,         // int root
         MPI_COMM_WORLD    // MPI_Comm comm
         );

   int Nx_offset = 0;
   int psi_mag_offset = 0;
   for ( size_t ii=0; ii < totalnodes; ++ii)
   {
      Nx_displacements[ii] = Nx_offset;
      Nx_offset += Nx_strides[ii];

      psi_mag_displacements[ii] = psi_mag_offset;
      psi_mag_offset += psi_mag_strides[ii];

      if ( (mynode == rootnode) && flags.debug )
            cout << "node " << ii << ", Nx_strides[" << ii << "]: " 
               << Nx_strides[ii] 
               << endl
               << "node " << ii << ", local_alloc_size_fftw: " 
               << psi_mag_strides[ii] << endl;
   }

   //////////////////////////////////////////////////////////////////
   // Retrieve fftw wisdom 
   //////////////////////////////////////////////////////////////////
   // TODO: ensure wisdom file is customized to local machine type, 
   //        put machine type in the wisdom file name

   string wisdom_file_name
      = "wisdom_file_c2c2c_cifb_" + TEM_NS::to_string(Nx)
         + "x" + TEM_NS::to_string(Ny)
         + "_"
         + TEM_NS::to_string(totalnodes)
         + "nodes";
   
   flags.wisdomFile = 1;
   if ( mynode == rootnode )
   {
      flags.wisdomFile =
         fftw_import_wisdom_from_filename( wisdom_file_name.c_str() );
         // returns 0 (== false) on failure
      if ( ! flags.wisdomFile && flags.debug )
      {
         cout << "Could not open wisdom file: " << wisdom_file_name
            << " , will create it after execution ..." << endl;
      }
   }

   MPI_Bcast( &flags.wisdomFile, 1, MPI_UNSIGNED, rootnode, 
         MPI_COMM_WORLD );

   //////////////////////////////////////////////////////////////////
   // allocate domains and domain related variables
   //////////////////////////////////////////////////////////////////
   double* kx_joined; kx_joined = new double[ Nx ];
   double* xx_joined; xx_joined = new double[ Nx ];
   double* kx_local; kx_local = new double[ Nx_local ];
   double* xx_local; xx_local = new double[ Nx_local ];// debug
   double* ky; ky = new double[ Ny ];
   double* yy; yy = new double[ Ny ];
   double periods_and_mins[6];
   double bwcutoff_t_sqr; 
   double bwcutoff_t; // transmission function bandwidth cutoff
   double bwcutoff_pap;
   double delta_x, delta_y, delta_kx, delta_ky;
   double xperiod_duped, yperiod_duped, zperiod_duped;
   double duped_periods[3]; // Array of periods so that only one Bcast 
                            //  is used.
   size_t number_of_raster_points_x, number_of_raster_points_y;
   double kxperiod;
   double kyperiod;

   scatterer_pap_LUT* myScattererPapLUT;
   bool papLUTexistsFlag; papLUTexistsFlag = false;

   std::vector< scatterer > myScatterers;
   std::list< slice* > slicesOfScatterers;

   fftw_complex* init_probe;
   init_probe = fftw_alloc_complex( local_alloc_size_fftw );
   fftw_plan pb_c2c_probe_split;
   pb_c2c_probe_split = fftw_mpi_plan_dft_2d( 
                           Nx, Ny,
                           init_probe, init_probe,
                           //comm, FFTW_BACKWARD, FFTW_ESTIMATE );
                           //comm, FFTW_BACKWARD, FFTW_PATIENT );
                           //comm, FFTW_BACKWARD, FFTW_EXHAUSTIVE );
                           comm, FFTW_BACKWARD, FFTW_MEASURE );

   double* probe_split_re;
   probe_split_re = new double[ local_alloc_size_fftw ];
   double* probe_split_im;
   probe_split_im = new double[ local_alloc_size_fftw ];

   double* probe_joined_re;
   probe_joined_re = new double[Nx * Ny];
   double* probe_joined_im;
   probe_joined_im = new double[Nx * Ny];

   fftw_complex* psi;
   psi = fftw_alloc_complex( local_alloc_size_fftw );

   size_t mtf_1D_size = mtf_1D.size();
   double rsqr, domain1, domain2, sincx, sincy, sincarg, xxx, yyy;
   fftw_complex* mtf_2D_split; 
   fftw_plan pf_c2c_mtf, pb_c2c_mtf;

   unsigned int number_of_k_bins;  // azimuthal integration subintervals
   unsigned int number_of_phi_bins = 0;

   double* diffracted_wave_mag_sum;
   double* diffracted_wave_mag_sum_sqr;
   double* diffracted_wave_mag_sqr_sum;
   double* diffracted_wave_re_sum;
   double* diffracted_wave_im_sum;

   //double* diffracted_wave_mag_avg_sqr; // used by gt17, d3, and d4 fem
   //double* diffracted_wave_mag_sqr_avg; // used by gt17, d3, and d4 fem

   double* diffracted_wave_mag_sqr_sum_radial_intensity_local;// gt17 & d3
   double* diffracted_wave_mag_sqr_sum_radial_intensity_total;// gt17 & d3

   double* diffracted_wave_mag_sum_sqr_radial_intensity_local;//gt17 & rva
   double* diffracted_wave_mag_sum_sqr_radial_intensity_total;//gt17 & rva

   double* diffracted_wave_mag_sum_radial_intensity_local; // d3 & rva
   double* diffracted_wave_mag_sum_radial_intensity_total; // d3 & rva

   double* diffracted_wave_mag_sqr;  // d1 & d2
   double* diffracted_wave_mag_radial_intensity_local; // d1 & d2 & corr.
   double* diffracted_wave_mag_radial_intensity_total; // d1 & d2 & corr.

   double* diffracted_wave_mag_radial_intensity_total_sum; // d1
   //double* diffracted_wave_mag_radial_intensity_sqr_local; // d1
   double* diffracted_wave_mag_radial_intensity_total_sqr_sum; // d1

   double* diffracted_wave_mag_sqr_radial_intensity_local; // d2
   double* diffracted_wave_mag_sqr_radial_intensity_total; // d2
   double* diffracted_wave_variance_sum;

   double* diffracted_wave_mag_variance; // d4
   double* diffracted_wave_mag_variance_radial_intensity_local; // d4
   double* diffracted_wave_mag_variance_radial_intensity_total; // d4

   double* diffracted_wave_mag_in_radial_coords_local; // correlograph
   double* diffracted_wave_mag_in_radial_coords_total; // correlograph
   double* diffracted_wave_mag_in_radial_coords_sum; // correlograph

   double* correlograph;
   double* correlograph_sum; 
   double* correlograph_sqr_sum; // correlograph_variance
   //double* correlograph_variance;

   double* ftem_gt17;
   double* ftem_d1;
   //double* ftem_d2;
   double* ftem_d3;
   double* ftem_d4_local;
   double* ftem_d4;
   double* ftem_rva;

   int* k_bin_counts_local;
   int* k_bin_counts_aggregated;
   int* k_bin_counts_sqr_sum_aggregated;
   int* k_bin_counts_sum_sqr_aggregated;
   
   int* phi_bin_counts_local;
   int* phi_bin_counts_aggregated;

   fftw_complex* diffracted_wave_fftw_complex_sum;// only used for 1 image
   double* diffracted_wave_complex_sum_mag;

   double delta_k;   // length of aziumthal integration subintervals
   std::vector<double> k_binning_boundaries;
   std::vector<double> phi_binning_boundaries;
   double delta_phi; // TODO: make this an input parameter
   vector<indexed_vector_magnitude_sqr> indexed_magnitudes;

   size_t probe_idx_shift_x; probe_idx_shift_x = 0;
   size_t probe_idx_shift_y; probe_idx_shift_y = 0;
   size_t probe_idx_x; probe_idx_x = 0;
   size_t probe_idx_y; probe_idx_y = 0;

   /////////////////////////////////////////////////////////////
   // Create forward and reverse fftw plans for psi
   /////////////////////////////////////////////////////////////
   fftw_plan pf_c2c_psi, pb_c2c_psi;//, pf_c2c_t, pb_c2c_t;//,pb_c2c_pap;

   // c2c in-place forward fft
   pf_c2c_psi = fftw_mpi_plan_dft_2d( // c2c in-place fft,
                           Nx, Ny, 
                           psi, psi, 
                           //comm, FFTW_BACKWARD, FFTW_ESTIMATE );
                           //comm, FFTW_BACKWARD, FFTW_PATIENT );
                           //comm, FFTW_BACKWARD, FFTW_EXHAUSTIVE );
                           comm, FFTW_FORWARD, FFTW_MEASURE );

   // c2c in-place reverse fft
   pb_c2c_psi = fftw_mpi_plan_dft_2d( 
                           Nx, Ny,
                           psi, psi, 
                           //comm, FFTW_BACKWARD, FFTW_ESTIMATE );
                           //comm, FFTW_BACKWARD, FFTW_PATIENT );
                           //comm, FFTW_BACKWARD, FFTW_EXHAUSTIVE );
                           comm, FFTW_BACKWARD, FFTW_MEASURE );


   bool first_probe_flag; first_probe_flag  = true;
   double detected_intensity, detected_intensity_reduced;
   size_t pixel_number_x = 0;
   size_t pixel_number_y = 0;
   size_t number_of_raster_points; //= x_p.size() * y_p.size();
   double half_delta_x;// = 0.5 * (xx_joined[1] - xx_joined[0]);
   double half_delta_y;// = 0.5 * (yy[1] - yy[0]);

   if ( flags.mtf_file && flags.mtf_resolution )
   {
      mtf_1D_size = mtf_1D.size();
      mtf_2D_split =  fftw_alloc_complex( local_alloc_size_fftw );
      pf_c2c_mtf = fftw_mpi_plan_dft_2d( // c2c in-place fft,
                                 Nx, Ny, 
                                 mtf_2D_split, mtf_2D_split,
                                 //comm, FFTW_FORWARD, FFTW_ESTIMATE );
                                 //comm, FFTW_FORWARD, FFTW_PATIENT );
                                 //comm, FFTW_FORWARD, FFTW_EXHAUSTIVE );
                                 comm, FFTW_FORWARD, FFTW_MEASURE );

      pb_c2c_mtf = fftw_mpi_plan_dft_2d( // c2c in-place fft,
                                 Nx, Ny, 
                                 mtf_2D_split, mtf_2D_split,
                                 //comm, FFTW_BACKWARD, FFTW_ESTIMATE );
                                 //comm, FFTW_BACKWARD, FFTW_PATIENT );
                                 //comm, FFTW_BACKWARD, FFTW_EXHAUSTIVE );
                                 comm, FFTW_BACKWARD, FFTW_MEASURE );
      
      for ( size_t i=0; i < local_alloc_size_fftw ; ++i)
      {
         mtf_2D_split[i][0] = 0; // default value
         mtf_2D_split[i][1] = 0; // default value
      }
   }

   //// broadcast domains from ms-stem-fem input
   //// TODO: remove if positions and domains are input from lammps script
   //if ( mynode == rootnode )
   //{
   //   periods_and_mins[0] = xperiod;
   //   periods_and_mins[1] = yperiod;
   //   periods_and_mins[2] = zperiod;
   //   periods_and_mins[3] = xmin;
   //   periods_and_mins[4] = ymin;
   //   periods_and_mins[5] = zmin;
   //}

   //MPI_Bcast( periods_and_mins, 6, MPI_DOUBLE, 
   //      rootnode, MPI_COMM_WORLD);

   //if ( mynode != rootnode )
   //{
   //   xperiod = periods_and_mins[0];
   //   yperiod = periods_and_mins[1];
   //   zperiod = periods_and_mins[2];
   //   xmin = periods_and_mins[3];
   //   ymin = periods_and_mins[4];
   //   zmin = periods_and_mins[5];
   //} // TODO: remove if positions and domains are input from lammps script
   //////////////////////////////////////////////////////////////////
   // gather lammps model into ms-stem-fem variables
   ////////////////////////////////////////////////////////////////////
   // Notes: 
   // - Gathering all scatterers without sorting them simultaneously is
   //    inefficient. They will be sorted by 
   //    TEM_NS::assign_atoms_to_slices_z_auto(), which is also not 
   //    parallel.
   //    This serial sorting and slicing operation will be very 
   //     inefficient.

   if ( *((int *) lammps_extract_global(lmp, (char *) "triclinic")) == 0)
   {
      xmin = *((double *) lammps_extract_global(lmp, (char *) "boxxlo"));
      ymin = *((double *) lammps_extract_global(lmp, (char *) "boxylo"));
      zmin = *((double *) lammps_extract_global(lmp, (char *) "boxzlo"));

      xmax = *((double *) lammps_extract_global(lmp, (char *) "boxxhi"));
      ymax = *((double *) lammps_extract_global(lmp, (char *) "boxyhi"));
      zmax = *((double *) lammps_extract_global(lmp, (char *) "boxzhi"));

      xperiod = xmax - xmin;
      yperiod = ymax - ymin;
      zperiod = zmax - zmin;

      //double mins[3];
      //double maxs[3];
      //int periodicity[3];
      //double xy, xz, yz;
      //int box_change;
      //lammps_extract_box( lmp, mins, maxs, 
      //                     &xy, &yz, &xz, 
      //                     periodicity, &box_change);
      //xmin = mins[0];
      //ymin = mins[1];
      //zmin = mins[2];
      //xperiod = maxs[0] - xmin;
      //yperiod = maxs[1] - ymin;
      //zperiod = maxs[2] - zmin;

   } else {
      if ( mynode == rootnode )
         cout << "Error in ms-stem-fem-md: triclinic not supported";
      return EXIT_FAILURE;
   }


   int offset;
   ////double** q_lammps;
   ////q_lammps = new double[ 3* lammps_population ];
   q = new double[ 3* lammps_population ];
   Z_array = new unsigned int[ lammps_population ];
   if ( flags.dupe && flags.fem 
        &&
        (dupe_x > 1 || dupe_y > 1 || dupe_z > 1)
        &&
        (dupe_x > 0 && dupe_y > 0 && dupe_z > 0)
      )
   {
      q_duped = new double[ 3 * duped_population ];
   }

   ////int* Z_array_int;
   ////Z_array_int = new int[ lammps_population ];

   //////lammps_gather_atoms( lmp, (char *) "x", 1, 3, &q_lammps[0][0]);
   ////lammps_gather_atoms( lmp, (char *) "x", 1, 3, q);
   ////cout << "node " << mynode << " q[0], q[1], q[4] : " << q[0] << ", " << q[1] << ", " << q[4] << endl;

   ////lammps_gather_atoms( lmp, (char *) "type", 1, 1, Z_array_int);
   ////for (size_t ii=0; ii < lammps_population; ++ii)
   ////   Z_array[ii] = static_cast<unsigned int>(Z_array_int[ii]);
   ////cout << "node " << mynode << " Z_array[0] : " << Z_array[0] << endl;
   ////delete[] Z_array_int;

   //// error if tags are not defined or not consecutive
   //// NOTE: test that name = image or ids is not a 64-bit int in code?

   int errorFlag = 0;
   //if (lmp->atom->tag_enable == 0 || lmp->atom->tag_consecutive() == 0)
   //  errorFlag = 1;
   //if (lmp->atom->natoms > MAXSMALLINT) errorFlag = 1;
   //if (errorFlag) {
   //  if (lmp->comm->me == rootnode)
   //    lmp->error->warning(FLERR,"error in ms-stem-fem/gather");
   //  // TODO: clean up memory before this return
   //  return EXIT_FAILURE;
   //}

   ////  Notes: 
   ////   What is atom type stored as?
   ////    atom_vec_atomic.cpp: avec->data_atom() converts the 
   ////    second string token of an atoms section of an input file
   ////    to an integer. The atom type is described by an integer.
   ////    I think it's stored in int atom->avec->type[].
   ////    Positions are stored in atom->avec->x[][].

   //// TODO: make persistent and adapt to changing value of natoms
   ////        using memory->grow()
   //// types = Natom length vector of per-atom types
   int* types_sendbuf;   // needed when not using MPI_IN_PLACE
   //types_sendbuf  = lmp->memory->create(types_sendbuf, lmp->atom->natoms,
   //                "ms-stem-fem:types_sendbuf");//  storage
   //for (size_t ii=0; ii < lmp->atom->natoms; ++ii) types_sendbuf[ii] = 0; 
   //                                          // make types[i]=0 
   //                                          // if using MPI_IN_PLACE

   char name_type[] = "type";

   void *types_vptr;// = lmp->atom->extract(name_type);

   //if (types_vptr  == NULL) {
   //   lmp->error->warning(FLERR,
   //     "ms-stem-fem: unable to extract property, atom->extract('type')");
   //  // TODO: clean up memory before this return
   //   return EXIT_FAILURE;
   //}

   //// MPI_Allreduce with MPI_SUM to merge into data, ordered by atom ID
   //// TODO: is ordering by atom ID necessary? May each node differ in order?

   int *types_vector = NULL;

   //types_vector = (int *) types_vptr ;  // because types_vptr is a void*

   //// use atom ID to insert each atom's values into types[]
   tagint *tag;// = lmp->atom->tag;
   int nlocal;// = lmp->atom->nlocal;

   //for (size_t ii = 0; ii < nlocal; ii++) 
   //   types_sendbuf[tag[ii]-1] = types_vector[ii];

   //MPI_Allreduce(
   //      //MPI_IN_PLACE,         // sendbuf
   //      types_sendbuf,      // sendbuf
   //      Z_array,         // recvbuf
   //      lammps_population,       // count
   //      MPI_INT,      // datatype
   //      MPI_SUM,      // op
   //      lmp->world);  // comm

   //lmp->memory->destroy(types_sendbuf);

   //// send atom positions
   char name_x[] = "x";
   void *vptr_positions;// = lmp->atom->extract(name_x);
   double* positions_sendbuf;   // needed when not using MPI_IN_PLACE
   //positions_sendbuf = lmp->memory->create(positions_sendbuf, 
   //                     3 * lmp->atom->natoms,
   //                "ms-stem-fem:positions_sendbuf");// storage

   //offset = 0;
   //for (size_t ii=0; ii < 3* lammps_population; ++ii) 
   //   positions_sendbuf[offset++] = 0.0;

   double **positions_array = NULL;
   //positions_array = (double **) vptr_positions;

   //for (size_t ii = 0; ii < nlocal; ii++) {
   //  offset = 3*(tag[ii]-1);
   //  for (size_t jj = 0; jj < 3; jj++)
   //    positions_sendbuf[offset++] = positions_array[ii][jj];
   //}

   //MPI_Allreduce(
   //    //MPI_IN_PLACE, only available when all nodes are in a single group
   //      positions_sendbuf,
   //      q,
   //      3* lammps_population,
   //      MPI_DOUBLE,
   //      MPI_SUM,
   //      lmp->world);

   //lmp->memory->destroy(positions_sendbuf);

   //if ( flags.debug && mynode == rootnode )
   //{
   //   cout << "Finished gathering initial positions and types"
   //      << " from lammps"
   //      << endl;
   //}

   //// change Z_array numbers from lammps type numbers to Z values
   //for ( size_t ii=0; ii < lammps_population; ++ii)
   //   Z_array[ii] = Z_vector[ Z_array[ii] -1 ];

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

   ////////////////////////////////////////////////////////////////////
   //// Create a vector containing a list of the unique species types Z
   ////////////////////////////////////////////////////////////////////
   //
   //// TODO: parallelize this process of eliminating duplicate Zs
   //if ( mynode == rootnode )
   //{
   //   // push Z_array elements onto Z_vector, then sort and remove 
   //   //  duplicates from Z_vector
   //   for (unsigned int i = 0; i < lammps_population; i++)
   //   {
   //      Z_vector.push_back( Z_array[i] );
   //   }

   //   std::sort( Z_vector.begin(), Z_vector.end() );

   //   std::vector<unsigned int>::iterator itr;
   //   itr = std::unique( Z_vector.begin(), Z_vector.end() );
   //   Z_vector.resize( std::distance( Z_vector.begin(), itr) );

   //   Z_vector_size = Z_vector.size();
   //}

   //MPI_Bcast( &Z_vector_size, 1, 
   //      MPI_UNSIGNED, rootnode, MPI_COMM_WORLD);
   //if ( flags.debug )
   //   cout << "node " << mynode // debug
   //      << ", Z_vector_size : " << Z_vector_size << endl; // debug

   //if ( mynode != rootnode ) Z_vector.resize( Z_vector_size );
   //
   ////MPI_Bcast( &Z_vector[0], Z_vector.size(),
   //MPI_Bcast( &Z_vector[0], Z_vector_size,
   //            MPI_UNSIGNED, rootnode, MPI_COMM_WORLD);

   //if ( flags.debug )
   //   cout << "node " << mynode  // debug
   //      << ", Z_vector.size (resized) : " 
   //      << Z_vector.size() << endl;// debug
   // TODO: The above is probably too reliant on communication, 
   //       the list of Zs is probably small enough that communication is
   //       slower than having duplicate calls to sort on each node.

   //////////////////////////////////////////////////////////////////
   // identify the STEM probe locations
   //////////////////////////////////////////////////////////////////
   if ( 
         flags.adfstem_corrected 
         || flags.adfstem_uncorrected 
      )
   {
      if ( ! flags.adfstem_detector_angles )
      {  // set default detector angles
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
      }

      // Create default raster positions if they weren't given by
      //  the parent call.
      if ( ! flags.raster_spacing )
      {
         if ( flags.fem )
            raster_spacing  = 10; // \AA // TODO: are these appropriate?
         else
            raster_spacing = 1.5; // \AA
      }

      // Create points at which the rastered stem beam will evaluated
      if (flags.adfstem_uncorrected || flags.adfstem_corrected)
      {

         // NOTE: do not increase number of raster points when using
         //        --dupe
         number_of_raster_points_x
            = floor( xperiod/raster_spacing );

         number_of_raster_points_y
            = floor( yperiod/raster_spacing );

         if ( mynode == rootnode && flags.debug )
            cout << "Determining STEM raster points with parameters" 
               << endl
               << " (xmin, ymin) : (" 
               << xmin << ", " << ymin << ")" << endl
               << " (xperiod, yperiod) : (" 
               << xperiod << ", " << yperiod << ")" << endl
               << "(number_of_raster_points_x, number_of_raster_points_y)"
               << " : (" << number_of_raster_points_x << ", "
               << number_of_raster_points_y << ")" << endl;

         assign_stem_raster_points(
               number_of_raster_points_x,
               number_of_raster_points_y,
               xperiod, xmin,
               yperiod, ymin,
               x_p, y_p
               );
         if ( mynode == rootnode && flags.debug )
            cout << "Determined STEM raster points." << endl
               << " (x_p.front(), y_p.front()) : ("
               << x_p.front() << ", " << y_p.front() << ")" << endl;

            number_of_raster_points = x_p.size() * y_p.size();
      }
   }

   const double lambda_sqr = lambda * lambda;
   const double alpha_max_sqr = alpha_max * alpha_max;
   
   size_t resolutionUnit, resolutionUnit_recip;
   double xResolution, yResolution,             // for real space images
          xResolution_stem, yResolution_stem,   // for real space images
          xResolution_recip, yResolution_recip, // for diffraction images
          kResolution_correlograph, phiResolution_correlograph;
   resolutionUnit = 3;  // 1: none, 2: inches, 3: centimeters
   resolutionUnit_recip = 1;  // 1: none, 2: inches, 3: centimeters
   double bcast_resolutions[6];
   // XRESOLUTION : "The number of pixels per RESOLUTIONUNIT in the 
   //                IMAGEWIDTH direction."
   if ( mynode == rootnode ) 
   {
      xResolution_stem =              // pixels per centimeter
         (x_p.size() / xperiod)  // angstrom^{-1}
               * 1.0e8;          // cm^{-1}

      yResolution_stem =              // pixels per centimeter
         (y_p.size() / yperiod)  // angstrom^{-1}
               * 1.0e8;          // cm^{-1}

      xResolution = (Nx / xperiod_duped) * 1.0e8;
      yResolution = (Ny / yperiod_duped) * 1.0e8;
   
      // xResolution_recip  = Nx / kxperiod
      // kxperiod = Nx / xperiod
      // \implies
      // xResolution_recip = Nx / (Nx / xperiod) = xperiod
      xResolution_recip =        // pixels per centimeter^{-1}  
         xperiod_duped * 1.0e-8;
      yResolution_recip =        // pixels per centimeter^{-1}  
         yperiod_duped * 1.0e-8;

      bcast_resolutions[0] = xResolution;
      bcast_resolutions[1] = yResolution;
      bcast_resolutions[2] = xResolution_stem;
      bcast_resolutions[3] = yResolution_stem;
      bcast_resolutions[4] = xResolution_recip;
      bcast_resolutions[5] = yResolution_recip;
   }
   MPI_Bcast(bcast_resolutions, 6, MPI_DOUBLE, rootnode, comm);

   double diffraction_scale_factor;


   ///////////////////////////////////////////////////////////////////
   // Instantiate the STEM image
   ///////////////////////////////////////////////////////////////////
   double* stem_image;
   if ( mynode == rootnode )
   {
      size_t raster_point_count;
      raster_point_count = x_p.size() * y_p.size();

      stem_image = new double[ raster_point_count ];
      for ( size_t ii=0; ii < raster_point_count; ++ii )
         stem_image[ii] = 0;
   }

   // variable to which diffraction will accumulate over lammps runs
   fftw_complex* psi_accumulation;
   psi_accumulation = fftw_alloc_complex( local_alloc_size_fftw );

   // flags for controlling domain dependent quantity updates
   bool update_domain_x; update_domain_x = true;
   bool update_domain_y; update_domain_y = true;

   bool fresh_raster_point; fresh_raster_point = true;

   //////////////////////////////////////////////////////////////////
   // loop over raster points
   //////////////////////////////////////////////////////////////////

   // create dumps of each configuration sampled by TEM
   //string lammps_dump_cmd;
   //lammps_dump_cmd = "dump 2 all atom " 
   //                     + to_string(lammps_TEM_samples)
   //                     + " " + to_string(output_prefix) 
   //                     + ".dump.gz";
   //lmp->input->one( lammps_dump_cmd.c_str() );

   for( std::list<double>::const_iterator 
         x_itr = x_p.begin();
         x_itr != x_p.end(); 
         ++x_itr )
   {
      pixel_number_y = 0;
      for( std::list<double>::const_iterator 
            y_itr = y_p.begin();
            y_itr != y_p.end(); 
            ++y_itr )
      { // loop over raster points
         fresh_raster_point = true;
         for ( size_t ii=0; ii < local_alloc_size_fftw; ++ii)
         {
            psi_accumulation[ii][0] = 0;
            psi_accumulation[ii][1] = 0;
         }

         if ( flags.debug ) 
            cout << "Running lammps at probe position : ("
               << *x_itr << ", " << *y_itr << ")" << endl;
         string lammps_run_cmd;
         lammps_run_cmd = "run " + to_string(lammps_TEM_steps);
         //cout << "lammps_run_cmd : " << lammps_run_cmd << endl;//debug
         for ( size_t lammps_step_idx=0; 
               lammps_step_idx < lammps_TEM_samples; 
               ++lammps_step_idx)
         {
            // loop over lammps runs  (closes near line 2878)

            if ( lammps_TEM_steps > 0 )
               lmp->input->one( lammps_run_cmd.c_str() );

            ////////////////////////////////////////////////////////////
            // update population and reallocate positions and Zs
            ////////////////////////////////////////////////////////////
            // lammps_population = 
            // initial_population = lammps_population
            lammps_population 
                  = static_cast<unsigned int>( lmp->atom->natoms );
            
            if ( lammps_population != initial_population )
            {
               initial_population = lammps_population;
               duped_population = dupe_x * dupe_y * dupe_z 
                                    * lammps_population;
               delete[] q;
               q = new double[ 3* lammps_population ];
               delete[] Z_array;
               Z_array = new unsigned int[ lammps_population ];
               
               if ( flags.dupe && flags.fem 
                    &&
                    (dupe_x > 1 || dupe_y > 1 || dupe_z > 1)
                    &&
                    (dupe_x > 0 && dupe_y > 0 && dupe_z > 0)
                  )
               {
                  delete[] q_duped;
                  q_duped = new double[ 3 * duped_population ];
               }

            }

            ////////////////////////////////////////////////////////////
            // update periods and domain from lammps
            ////////////////////////////////////////////////////////////
            if ( flags.debug && mynode == rootnode )
               cout << "updating periods and domain" << endl;//debug
            if ( *((int *) lammps_extract_global(lmp, 
                           (char *) "triclinic")) == 0)
            {
               double xmin_next, ymin_next, zmin_next;
               double xmax_next, ymax_next, zmax_next;
               double xperiod_next, yperiod_next, zperiod_next;

               xmin_next = *((double *) lammps_extract_global(lmp, 
                              (char *) "boxxlo"));
               ymin_next = *((double *) lammps_extract_global(lmp, 
                              (char *) "boxylo"));
               zmin_next = *((double *) lammps_extract_global(lmp, 
                              (char *) "boxzlo"));

               xmax_next = *((double *) lammps_extract_global(lmp, 
                              (char *) "boxxhi"));
               ymax_next = *((double *) lammps_extract_global(lmp, 
                              (char *) "boxyhi"));
               zmax_next = *((double *) lammps_extract_global(lmp, 
                              (char *) "boxzhi"));

               xperiod_next = xmax_next - xmin_next;
               yperiod_next = ymax_next - ymin_next;
               zperiod_next = zmax_next - zmin_next;

               zperiod = zperiod_next;

               if ( xperiod != xperiod_next ) 
               {
                  update_domain_x = true;
                  xperiod = xperiod_next;
               }
               //else 
               //{
               //   update_domain_x = false;
               //}
               if ( yperiod != yperiod_next ) 
               {
                  update_domain_y = true;
                  yperiod = yperiod_next;
               }
               //else 
               //{
               //   update_domain_y = false;
               //}

            } else {
               if ( mynode == rootnode )
                  cout << "Error in ms-stem-fem-md: triclinic not supported";
               return EXIT_FAILURE;
            }

            ////////////////////////////////////////////////////////////
            // update arrays of positions and Zs from lammps
            ////////////////////////////////////////////////////////////
            if ( flags.debug && mynode == rootnode )
               cout << "updating ms-stem-fem positions and Zs" 
                  << endl;//debug

            errorFlag = 0;
            if (lmp->atom->tag_enable == 0 
                  || lmp->atom->tag_consecutive() == 0)
              errorFlag = 1;
            if (lmp->atom->natoms > MAXSMALLINT) errorFlag = 1;
            if (errorFlag) {
              if (lmp->comm->me == rootnode)
                lmp->error->warning(FLERR,"error in ms-stem-fem/gather");
              // TODO: clean up memory before this return
              return EXIT_FAILURE;
            }

            // TODO: make persistent and adapt to changing value of natoms
            //        using memory->grow()
            // types = Natom length vector of per-atom types
            //int* types_sendbuf;   // needed when not using MPI_IN_PLACE
            types_sendbuf  = lmp->memory->create(types_sendbuf, 
                                 lmp->atom->natoms,
                                 "ms-stem-fem:types_sendbuf");//  storage
            for (size_t ii=0; ii < lmp->atom->natoms; ++ii) 
               types_sendbuf[ii] = 0; 
                                                 // make types[i]=0 
                                                 // if using MPI_IN_PLACE

            //char name_type[] = "type";
            types_vptr = lmp->atom->extract(name_type);
            if (types_vptr  == NULL) {
               lmp->error->warning(FLERR,
                 "ms-stem-fem: unable to extract property, atom->extract('type')");
              // TODO: clean up memory before this return
               return EXIT_FAILURE;
            }


            // MPI_Allreduce with MPI_SUM to merge into data, ordered by atom ID
            // TODO: is ordering by atom ID necessary? 
            //       May each node differ in order?

            types_vector = NULL;
            types_vector = (int *) types_vptr ; // because types_vptr 
                                                // is a void*

            // use atom ID to insert each atom's values into types[]
            tag = lmp->atom->tag;
            nlocal = lmp->atom->nlocal;

            for (size_t ii = 0; ii < nlocal; ii++) 
               types_sendbuf[tag[ii]-1] = types_vector[ii];

            MPI_Allreduce(
                  //MPI_IN_PLACE,         // sendbuf
                  types_sendbuf,      // sendbuf
                  Z_array,         // recvbuf
                  lammps_population,       // count
                  MPI_INT,      // datatype
                  MPI_SUM,      // op
                  lmp->world);  // comm

            lmp->memory->destroy(types_sendbuf);

            // send atom positions
            //char name_x[] = "x";
            vptr_positions = lmp->atom->extract(name_x);
            //double* positions_sendbuf;   // needed when not using MPI_IN_PLACE
            positions_sendbuf = lmp->memory->create(positions_sendbuf, 
                                 3 * lmp->atom->natoms,
                            "ms-stem-fem:positions_sendbuf");// storage
            offset = 0;
            for (size_t ii=0; ii < 3* lammps_population; ++ii) 
               positions_sendbuf[offset++] = 0.0;

            //double **positions_array = NULL;
            positions_array = (double **) vptr_positions;

            for (size_t ii = 0; ii < nlocal; ii++) {
              offset = 3*(tag[ii]-1);

              // enforce periodic boundary conditions
              while ( positions_array[ii][0] <= xmin )
                 positions_array[ii][0] += xperiod;
              while ( positions_array[ii][0] > xmax )
                 positions_array[ii][0] -= xperiod;
              while ( positions_array[ii][1] <= ymin )
                 positions_array[ii][1] += yperiod;
              while ( positions_array[ii][1] > ymax )
                 positions_array[ii][1] -= yperiod;

              for (size_t jj = 0; jj < 3; jj++)
                positions_sendbuf[offset++] = positions_array[ii][jj];
            }

            MPI_Allreduce(
                //MPI_IN_PLACE, only available when all nodes are in a single group
                  positions_sendbuf,
                  q,
                  3* lammps_population,
                  MPI_DOUBLE,
                  MPI_SUM,
                  lmp->world);

            lmp->memory->destroy(positions_sendbuf);

            if ( flags.debug && (mynode == rootnode) )
               cout << "updating Z_array values" << endl;//debug
            // change Z_array numbers from lammps type numbers to Z values
            for ( size_t ii=0; ii < lammps_population; ++ii)
               Z_array[ii] = Z_vector[ Z_array[ii] -1 ];

            if ( flags.debug && mynode == rootnode )
            {
               cout << "Finished gathering positions and types"
                  << " from lammps"
                  << endl;
            }

            ////////////////////////////////////////////////////////////
            // update domain dependent quantities
            ////////////////////////////////////////////////////////////

            //if ( xperiod or yperiod have changed ) // TODO
            //debug
            //update_domain_x = true;
            //update_domain_y = true;
            //end debug
            if ( update_domain_x || update_domain_y )
            {   
               // broadcast domain periods and minima
               // TODO: why bother with broadcasting from rootnode,
               //       when periods and mins were grabbed from lammps
               //       at all nodes?
               //if ( mynode == rootnode )
               //{
               //   periods_and_mins[0] = xperiod;
               //   periods_and_mins[1] = yperiod;
               //   periods_and_mins[2] = zperiod;
               //   periods_and_mins[3] = xmin;
               //   periods_and_mins[4] = ymin;
               //   periods_and_mins[5] = zmin;
               //}

               //MPI_Bcast( periods_and_mins, 6, MPI_DOUBLE, 
               //      rootnode, MPI_COMM_WORLD);

               //if ( mynode != rootnode )
               //{
               //   xperiod = periods_and_mins[0];
               //   yperiod = periods_and_mins[1];
               //   zperiod = periods_and_mins[2];
               //   xmin = periods_and_mins[3];
               //   ymin = periods_and_mins[4];
               //   zmin = periods_and_mins[5];
               //}

               /////////////////////////////////////////////////////////
               // Modify parameters sensitive to sample duplication
               /////////////////////////////////////////////////////////
               // NOTE: change Nx, Ny, Nz ?  No? 
               //       change xperiod, ..., ? For some cases, yes.
               //       The period of the cell containing duplicates is 
               //       needed for the inverse area factor, but not 
               //       needed for assigning raster points.
               
               //if ( mynode == rootnode )
               //{
               xperiod_duped = xperiod * dupe_x;
               yperiod_duped = yperiod * dupe_y;
               zperiod_duped = zperiod * dupe_z;

               //   duped_periods[0] = xperiod_duped;
               //   duped_periods[1] = yperiod_duped;
               //   duped_periods[2] = zperiod_duped;
               //}

               //MPI_Bcast(duped_periods, 3, MPI_DOUBLE, 
               //            rootnode, MPI_COMM_WORLD);

               //if ( mynode != rootnode )
               //{
               //   xperiod_duped = duped_periods[0];
               //   yperiod_duped = duped_periods[1];
               //   zperiod_duped = duped_periods[2];
               //}

               /////////////////////////////////////////////////////////
               // Create domains
               /////////////////////////////////////////////////////////
               // Only splitting the kx domain amongst nodes, not x, y, 
               //  or ky

               delta_x = xperiod_duped / Nx; 
               delta_y = yperiod_duped / Ny; 
               delta_kx = 1.0/xperiod_duped;
               delta_ky = 1.0/yperiod_duped;

               kxperiod = Nx / xperiod_duped; 
               kyperiod = Ny / yperiod_duped;

               //if ( mynode == rootnode )
               //{
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

               if ( flags.debug  && (mynode == rootnode) )
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
               //}

               if ( flags.debug && (mynode == rootnode))
               {
                  for ( size_t ii=0; ii< totalnodes; ++ii)
                  {
                     cout << "Nx_strides[" << ii 
                        << "], Nx_displacements[ "
                        << ii << "]: " << Nx_strides[ii] << ", " 
                        << Nx_displacements[ii] << endl;
                     cout << "Scatterving xx and kx ..." << endl;
                  }
               }

               MPI_Scatterv( 
                           xx_joined,        // const void *sendbuf
                           Nx_strides,       // const int sendcounts[]
                           Nx_displacements, // const int displacements[]
                           MPI_DOUBLE,       // MPI_Datatype sendtype
                           xx_local,         // void *recvbuf
                           Nx_local,         // int recvcount
                           MPI_DOUBLE,       // MPI_Datatype recvtype
                           rootnode,         // int root
                           MPI_COMM_WORLD);  // MPI_Comm comm

               MPI_Scatterv( 
                           kx_joined,        // const void *sendbuf
                           Nx_strides,       // const int sendcounts[]
                           Nx_displacements, // const int displacements[]
                           MPI_DOUBLE,       // MPI_Datatype sendtype
                           kx_local,         // void *recvbuf
                           Nx_local,         // int recvcount
                           MPI_DOUBLE,       // MPI_Datatype recvtype
                           rootnode,         // int root
                           MPI_COMM_WORLD);  // MPI_Comm comm

               if ( flags.debug && (mynode == rootnode))
                  cout << "Scattered xx and kx domains ..." << endl;

               // TODO: consider bundling the follow data into a 
               //  single Bcast
               MPI_Bcast( ky, Ny, MPI_DOUBLE, rootnode, MPI_COMM_WORLD);
               MPI_Bcast( yy, Ny, MPI_DOUBLE, rootnode, MPI_COMM_WORLD);
               MPI_Bcast( kx_joined, Nx, MPI_DOUBLE, 
                           rootnode, MPI_COMM_WORLD);
               MPI_Bcast( xx_joined, Nx, MPI_DOUBLE, 
                           rootnode, MPI_COMM_WORLD);

               /////////////////////////////////////////////////////////
               // Establish bandwidth limiting cutoff values
               /////////////////////////////////////////////////////////

               if ( Nx/(2 * xperiod_duped) < Ny/(2 * yperiod_duped) ) 
                  bwcutoff_pap = Nx/(2 * xperiod_duped);
               else 
                  bwcutoff_pap = Ny/(2 * yperiod_duped);

               bwcutoff_t = (2.0/3.0) * bwcutoff_pap;

               bwcutoff_t_sqr = bwcutoff_t * bwcutoff_t;

               if ( mynode == rootnode && flags.debug )
               {
                  cout << "Bandwidth limits " << endl // debug
                     << "  projected atomic potential bandwidth : "//debug
                     << bwcutoff_pap << endl // debug
                     << "  transmission function bandwidth : " // debug
                     << bwcutoff_t << endl; // debug
               }
               
               half_delta_x = 0.5 * (xx_joined[1] - xx_joined[0]);
               half_delta_y = 0.5 * (yy[1] - yy[0]);

               if ( (flags.adfstem_corrected || flags.adfstem_uncorrected)
                     && flags.adfstem_detector_angles )
               {
                  if (detector_inner_angle < 0.0) 
                  {
                     if (mynode == rootnode)
                        cout << "Error: detector_inner_angle < 0" << endl;
                     flags.fail = 1;
                  }
                  if (detector_outer_angle < 0.0) 
                  {
                     if (mynode == rootnode)
                        cout << "Error: detector_outer_angle < 0" << endl;
                     flags.fail = 1;
                  }
                  if (detector_inner_angle > bwcutoff_t * lambda)
                  {
                     if (mynode == rootnode)
                     {
                        cout << "Error: detector inner angle falls beyond"
                          << " bandwidth limits" << endl
                          << " detector angles: " << detector_inner_angle
                          << ", " << detector_outer_angle << endl
                          << " bandwidth limit: " << bwcutoff_t * lambda 
                          << endl;
                     }
                     flags.fail = 1;
                  }
               } // end adfstem

               // update mtf components which are dependent upon domain 
               //   changes
               if ( flags.mtf_file && flags.mtf_resolution )
               {
                  ///////////////////////////////////////////////////////
                  // Evaluate the two dimensional representations of 
                  //  the modulation transfer function
                  ///////////////////////////////////////////////////////
                  for ( size_t i=0; i < Nx_local; ++i)
                  {
                     for ( size_t j=0; j < Ny; ++j)
                     {
                        xxx = pow( xperiod_duped, 2) * kx_local[i] / Nx;
                        yyy = pow( yperiod_duped, 2) * ky[j] / Ny;
                        rsqr = pow( xxx, 2) + pow( yyy, 2);
                        //rsqr = pow( pow(xperiod_duped, 2) 
                        //          *kx_local[i]/Nx, 2) 
                        //     + pow( pow(yperiod_duped, 2)* ky[j]/Ny, 2);

                        mtf_2D_split[j + i * Ny][0] = 0; // default value
                        mtf_2D_split[j + i * Ny][1] = 0; // default value
                        for ( size_t mtf_idx = 0; 
                              mtf_idx < mtf_1D_size -1; 
                              ++mtf_idx)
                        {
                           if ( mtf_1D_size == 1 )
                           {
                              mtf_2D_split[j + i * Ny][0] 
                                 = mtf_1D[mtf_idx];
                           }
                           if (
                                 (rsqr >= mtf_domain_sqr[mtf_idx]) 
                                && (rsqr < mtf_domain_sqr[mtf_idx+1]) 
                              )
                           {
                              domain1 = sqrt(mtf_domain_sqr[mtf_idx]);
                              domain2 = sqrt(mtf_domain_sqr[mtf_idx+1]);
                              mtf_2D_split[j + i * Ny][0] = 
                                mtf_1D[mtf_idx] 
                                + (mtf_1D[mtf_idx+1] - mtf_1D[mtf_idx])
                                 * (sqrt(rsqr) - domain1)
                                  /(domain2 - domain1);

                              // Krause  2013
                              sincarg = xxx / mtf_resolution;
                              // Thust 2009
                              //sincarg = 0.5 * PI * xxx ;

                              if ( sincarg != 0 )
                                 sincx = sin( sincarg ) / sincarg;
                              else 
                                 sincx = 1;

                              // Krause 2013
                              sincarg = yyy / mtf_resolution;
                              // Thust 2009
                              //sincarg = 0.5 * PI * yyy ;

                              if ( sincarg != 0 )
                                 sincy = sin( sincarg ) / sincarg;
                              else
                                 sincy = 1;

                              mtf_2D_split[j + i * Ny][0] = 
                               mtf_2D_split[j + i * Ny][0] 
                                * sincx * sincy;
                           }
                        }
                     }
                  }
                  //bw_limit(
                  //      mtf_2D_split,
                  //      Nx_local, kx_local, Ny, ky,
                  //      //xperiod_duped/3
                  //      1.5*bwcutoff_t
                  //      //bwcutoff_t
                  //      );
                  fftw_execute( pb_c2c_mtf );
                  //if ( flags.debug && flags.image_output )
                  //{
                  //   cout << "saving 2-D modulation transfer function" 
                  //     <<    " to tiff" << endl;

                  //   debug_output_complex_fftw_operand(
                  //         mtf_2D_split,
                  //         0,
                  //         local_alloc_size_fftw,
                  //         Nx_local, Nx, Ny,
                  //         resolutionUnit,
                  //         xResolution, yResolution,
                  //         output_prefix 
                  //         + "_mtf_2D_reciprocalspace",
                  //         psi_mag_strides,
                  //         psi_mag_displacements,
                  //         mynode, rootnode, comm);
                  //}
                  bw_limit(
                        mtf_2D_split,
                        Nx_local, kx_local, Ny, ky,
                        bwcutoff_t
                        );
                  //if ( flags.debug && flags.image_output )
                  //{
                  //   cout << "saving 2-D modulation transfer function"
                  //      << " to tiff" << endl;

                  //   debug_output_complex_fftw_operand(
                  //         mtf_2D_split,
                  //         0,
                  //         local_alloc_size_fftw,
                  //         Nx_local, Nx, Ny,
                  //         resolutionUnit,
                  //         xResolution, yResolution,
                  //         output_prefix 
                  //         + "_mtf_2D_realspace_bwlimited",
                  //         psi_mag_strides,
                  //         psi_mag_displacements,
                  //         mynode, rootnode, comm);
                  //}

                  fftw_execute( pf_c2c_mtf );

                  // renormalize mtf_2D
                  for (ptrdiff_t i=0; i < local_alloc_size_fftw; ++i)
                  {
                     mtf_2D_split[i][0] = mtf_2D_split[i][0] / NxNy;
                     mtf_2D_split[i][1] = mtf_2D_split[i][1] / NxNy;
                  }

                  if ( flags.debug && flags.image_output )
                  {
                     cout << "saving 2-D modulation transfer function"
                       << " to tiff" << endl;

                     debug_output_complex_fftw_operand(
                           mtf_2D_split,
                           0,
                           local_alloc_size_fftw,
                           Nx_local, Nx, Ny,
                           resolutionUnit,
                           xResolution, yResolution,
                           output_prefix 
                           + "_mtf_2D",
                           psi_mag_strides,
                           psi_mag_displacements,
                           mynode, rootnode, comm);
                  }
               } // end update of mtf components w.r.t. domain changes

               //update_domain_x = false;
               //update_domain_y = false;
            } // end update domains

            if ( flags.fail == 1)
            {
               // TODO: update this
               // clean up all allocated memory and exit
               delete[] Nx_strides;
               delete[] Nx_displacements;
               delete[] psi_mag_strides;
               delete[] psi_mag_displacements;
               delete[] kx_joined;
               delete[] xx_joined;
               delete[] kx_local;
               delete[] xx_local;
               delete[] ky;
               delete[] yy;

               if ( flags.mtf_file && flags.mtf_resolution )
               {
                  fftw_destroy_plan( pf_c2c_mtf );
                  fftw_destroy_plan( pb_c2c_mtf );
                  fftw_free( mtf_2D_split );
               }

               if ( mynode == rootnode ) 
                  cerr << " Exiting in error after attempting to"
                    << " update domains" << endl;

               fftw_mpi_cleanup();
               MPI_Finalize();
               return EXIT_FAILURE;
            }

            //  Allocate FEM relevant variables
            if ( flags.fem && first_probe_flag )
                  // TODO:
                  //  If the domain changes during MD, how should the 
                  //   following be altered?
            { // Initalize k_binning_boundaries for azimuthal integration.

               delta_k =   
                  azimuthal_binning_size_factor * sqrt( 
                    (kx_local[1] - kx_local[0]) 
                     * (kx_local[1] - kx_local[0])
                      + (ky[1] - ky[0]) * (ky[1] - ky[0])
                     );

               unsigned int* number_of_bins;
               number_of_bins = new unsigned int[2];

               if ( mynode == rootnode ) 
               {
                  // NOTE: why not duplicate the following on all nodes?
                  //       I want to be certain that the number of bins 
                  //       are all identical.
                  unsigned int phi_bdy_count = 0;
                  number_of_k_bins = bwcutoff_t / delta_k;

                  for ( unsigned int ii=0; 
                        ii < number_of_k_bins + 1; 
                        ++ii)
                  {  // count the number of k bins between 0.2 and 1.0,
                     //  and use that as the number of phi bins
                     double k_bdy = ii * delta_k;
                     k_binning_boundaries.push_back( k_bdy );
                     if ( flags.correlograph ) // is this inefficient?
                        if (( k_bdy >= 0.2) && (k_bdy <= 1.0) ) 
                           ++phi_bdy_count;
                  }

                  if ( flags.correlograph ) // is this inefficient?
                  {
                     if (phi_bdy_count > 0)
                        number_of_phi_bins = phi_bdy_count - 1;
                     else
                        number_of_phi_bins = 0;
                  }

                  number_of_bins[0] = number_of_k_bins;
                  number_of_bins[1] = number_of_phi_bins;
                  if ( flags.debug && (mynode == rootnode) )
                  {
                     cout << "number_of_bins[] : (" << number_of_bins[0] 
                        << ", " << number_of_bins[1] << ")" << endl;
                  }
               }

               MPI_Bcast( number_of_bins, 2, MPI_UNSIGNED, 
                           rootnode, comm);

               if ( mynode != rootnode ) 
               {
                  number_of_k_bins = number_of_bins[0];
                  number_of_phi_bins = number_of_bins[1];
               }

               if ( flags.debug && (mynode != rootnode) )
               {
                  cout << "node : " << mynode << "number_of_bins[] : (" 
                     << number_of_bins[0] << ", " << number_of_bins[1] 
                     << ")" << endl;
               }

               delete[] number_of_bins;

               if ( flags.correlograph ) // is this inefficient?
                  delta_phi = 2*PI/number_of_phi_bins;

               if ( mynode != rootnode )
               {
                  for ( unsigned int ii=0; 
                        ii < number_of_k_bins + 1; 
                        ++ii)
                  {
                     k_binning_boundaries.push_back( ii * delta_k );
                  }
               }
               if ( flags.correlograph ) // is this inefficient?
                  for ( unsigned int ii=0; 
                        ii < number_of_phi_bins + 1; 
                        ++ii)
                     phi_binning_boundaries.push_back( ii * delta_phi );
               ////NOTE: it's probably faster to push_back on each node 
               //          rather than broadcast the values, even though 
               //          it's possible that they may slightly vary 
               //          between nodes.
               //}
               //else
               //{
               //   k_binning_boundaries.resize(number_of_k_bins + 1);
               //   phi_binning_boundaries.resize(number_of_phi_bins + 1);
               //}

               //MPI_Bcast( &k_binning_boundaries[0], number_of_k_bins +1,
               //            MPI_DOUBLE, rootnode, comm);
               //MPI_Bcast( &phi_binning_boundaries[0], 
               //             number_of_phi_bins + 1, 
               //            MPI_DOUBLE, rootnode, comm);

               // Initialize the vector of indexed_vector_magitude_sqr 
               //  for re-use
               for ( ptrdiff_t i=0; i<Nx_local; ++i)
                  for ( ptrdiff_t j=0; j<Ny; ++j)
                  {
                     indexed_magnitudes.push_back(
                           indexed_vector_magnitude_sqr( 
                                 i, j,
                                 kx_local[i] * kx_local[i] 
                                  + ky[j] * ky[j],
                                 PI + atan2(ky[j], kx_local[i]) 
                                  // TODO: improve?
                              )
                           );
                  }
               // sort the indexed_magnitudes by the magnitude of 
               //  their |k|^2 values
               std::sort(
                     indexed_magnitudes.begin(),
                     indexed_magnitudes.end(),
                     indexed_vector_magnitude_sqr_lt // pointer 
                                                     // to "<" function 
                     );
               //cout << "min |k|^2 of indexed_magnitudes : "  // debug
               //  << indexed_magnitudes.front().v_mag_sqr << endl;//debug
               //cout << "max |k|^2 of indexed_magnitudes : "  // debug
               //   << indexed_magnitudes.back().v_mag_sqr << endl;//debug

               // Allocate variables for the various FTEM modes
               if ( flags.d1 || flags.d2  || flags.d3 
                     || flags.correlograph
                     || flags.d4 || flags.gt17 || flags.rva )
               {  
                  k_bin_counts_local = new int[ number_of_k_bins ];
                  if ( flags.correlograph ) // is this inefficient?
                     phi_bin_counts_local = new int[ number_of_phi_bins ];
                  //k_bin_counts_sqr_local = new int[ number_of_k_bins ];
                  // Initiallization to 0 will be done by 
                  //    integrate_out_phi_...()
                  //for ( ptrdiff_t ii=0; ii < number_of_k_bins; ++ii )
                  //{
                  //   k_bin_counts_local[ii] = 0;
                  //   k_bin_counts_sqr_local[ii] = 0;
                  //}
                  if ( mynode == rootnode )
                  {
                     k_bin_counts_aggregated
                        = new int[ number_of_k_bins ];
                     k_bin_counts_sqr_sum_aggregated
                        = new int[ number_of_k_bins ];
                     k_bin_counts_sum_sqr_aggregated
                        = new int[ number_of_k_bins ];
                     for ( size_t ii=0; ii < number_of_k_bins; ++ii)
                     {
                        k_bin_counts_aggregated[ii] = 0;
                        k_bin_counts_sqr_sum_aggregated[ii] = 0;
                        k_bin_counts_sum_sqr_aggregated[ii] = 0;
                     }
                     if ( flags.correlograph )
                     {
                        phi_bin_counts_aggregated
                           = new int[ number_of_phi_bins ];
                        for ( size_t ii=0; ii < number_of_phi_bins; ++ii)
                        {
                           phi_bin_counts_aggregated[ii] = 0;
                        }
                     }
                  }
               }

               if ( flags.d1 || flags.d2 || flags.correlograph)
               {
                  diffracted_wave_mag_radial_intensity_local
                     = new double[number_of_k_bins];
                  //k_bin_counts_mag_local = new int[number_of_k_bins];
                  for ( ptrdiff_t ii=0; ii < number_of_k_bins; ++ii )
                  {
                     //k_bin_counts_mag_local[ii] = 0;
                     diffracted_wave_mag_radial_intensity_local[ii] = 0.0;
                  }

                  if ( flags.correlograph )
                  {
                     // total of this variable is held at the root node
                     //  but every node needs to accumulate its own under 
                     //  the same name
                     diffracted_wave_mag_in_radial_coords_local
                        = new double[
                              number_of_k_bins * number_of_phi_bins];
                     correlograph 
                        = new double[
                              number_of_phi_bins * number_of_k_bins];
                     correlograph_sum
                        = new double[
                              number_of_phi_bins * number_of_k_bins];
                     if ( flags.correlograph_variance )
                        correlograph_sqr_sum
                         = new double[
                              number_of_phi_bins * number_of_k_bins];
                     for ( ptrdiff_t ii=0; 
                           ii < number_of_k_bins * number_of_phi_bins;
                           ++ii)
                     {
                        correlograph_sum[ii] = 0.0;
                        if ( flags.correlograph_variance )
                           correlograph_sqr_sum[ii] = 0.0;
                     }
                  }

                  if ( flags.d2 )
                  {
                     diffracted_wave_mag_sqr
                        = new double[local_alloc_size_fftw];
                     diffracted_wave_mag_sqr_radial_intensity_local
                        = new double[ number_of_k_bins ];
                     //k_bin_counts_mag_sqr_local 
                     // = new int[number_of_k_bins];
                     for ( ptrdiff_t ii=0; ii < number_of_k_bins; ++ii )
                     {
                        //k_bin_counts_mag_sqr_local[ii] = 0;
                        diffracted_wave_mag_sqr_radial_intensity_local[ii]
                           = 0.0;
                     }
                 //for (ptrdiff_t ii=0; ii < local_alloc_size_fftw; ++ii)
                     //{
                     //   diffracted_wave_mag_sqr[ii] = 0.0;
                     //}
                  }

                  if ( mynode == rootnode )
                  {
                     diffracted_wave_mag_radial_intensity_total//temporary
                                                         // storage to
                        = new double[ number_of_k_bins ];// allow squaring
                     //k_bin_counts_mag_aggregated 
                     //   = new int[ number_of_k_bins ];
                     for ( ptrdiff_t ii=0; ii < number_of_k_bins; ++ii )
                     { // quantities common to both d1 & d2
                        diffracted_wave_mag_radial_intensity_total[ii] 
                           = 0.0;
                        //k_bin_counts_mag_aggregated[ii] = 0;
                     }

                     if ( flags.d1 )
                     {
                        ftem_d1 = new double[ number_of_k_bins ];
                        diffracted_wave_mag_radial_intensity_total_sum
                           = new double[ number_of_k_bins ];
                        diffracted_wave_mag_radial_intensity_total_sqr_sum
                           = new double[ number_of_k_bins ];
                        for ( ptrdiff_t ii=0; ii < number_of_k_bins; ++ii)
                        {
                       diffracted_wave_mag_radial_intensity_total_sum[ii]
                              = 0.0;
                    diffracted_wave_mag_radial_intensity_total_sqr_sum[ii]
                              = 0.0;
                        }
                     }
                     if ( flags.d2 )
                     {
                        diffracted_wave_mag_sqr_radial_intensity_total
                           = new double[ number_of_k_bins ];
                        diffracted_wave_variance_sum
                           = new double[ number_of_k_bins ];
                        for ( ptrdiff_t ii=0; ii < number_of_k_bins; ++ii)
                        {
                       diffracted_wave_mag_sqr_radial_intensity_total[ii]
                              = 0.0;
                           diffracted_wave_variance_sum[ii] 
                              = 0.0;
                        }
                     }
                     if ( flags.correlograph )
                     {
                        // average correlograph will be computed at root 
                        //  node
                        diffracted_wave_mag_in_radial_coords_sum
                         = new double[
                                 number_of_k_bins * number_of_phi_bins ];
                        diffracted_wave_mag_in_radial_coords_total
                         = new double[
                                 number_of_k_bins * number_of_phi_bins ];
                        for ( ptrdiff_t ii=0; 
                              ii < number_of_k_bins * number_of_phi_bins;
                              ++ii)
                        {
                       diffracted_wave_mag_in_radial_coords_sum[ii]= 0.0;
                     diffracted_wave_mag_in_radial_coords_total[ii]= 0.0;
                        }
                     }
                  }
               } // ( flags.d1 || flags.d2 || flags.correlograph )
 
               if ( flags.gt17 || flags.d3 || flags.d4 
                     || flags.rva ) 
               {  // accumulate I(\vec{k}) and I^{2}(\vec{k})

                  // used for 2-D variance image and GT17 & d3 & d4
                  diffracted_wave_mag_sum 
                     = new double[ local_alloc_size_fftw ];
                  diffracted_wave_mag_sum_sqr
                     = new double[ local_alloc_size_fftw ];
                  if ( flags.gt17 || flags.d3 || flags.d4) 
                  {
                     diffracted_wave_mag_sqr_sum 
                        = new double[ local_alloc_size_fftw ];
                  }

                  //diffracted_wave_mag_sqr_avg
                  //   = new double[ local_alloc_size_fftw ];

                  //if ( flags.gt17 || flags.d4 )
                  //   diffracted_wave_mag_avg_sqr
                  //      = new double[ local_alloc_size_fftw ];

                  //if ( flags.d3 )
                  //   diffracted_wave_mag_avg_sqr_d3
                  //      = new double[ local_alloc_size_fftw ];

                  if (flags.complex_realspace_sum) 
                  {
                     // used for 2-D variance image when adding phase in 
                     //  the image plane
                     diffracted_wave_re_sum 
                        = new double[ local_alloc_size_fftw ];
                     diffracted_wave_im_sum 
                        = new double[ local_alloc_size_fftw ];
                     //diffracted_wave_complex_variance_2D
                     //   = new double[ local_alloc_size_fftw ];
                  }
                  
                  for ( ptrdiff_t ii=0; ii < local_alloc_size_fftw; ++ii )
                  {// ensure that they are initialized to 0.0
                     // 2-D image
                     diffracted_wave_mag_sum[ii] = 0.0;
                     if ( flags.gt17 || flags.d3 || flags.d4) 
                     {
                        diffracted_wave_mag_sqr_sum[ii] = 0.0;
                     }

                     if (flags.complex_realspace_sum) 
                     {
                        diffracted_wave_re_sum[ii] = 0.0;
                        diffracted_wave_im_sum[ii] = 0.0;
                        //diffracted_wave_complex_variance_2D[ii] = 0.0;
                     }
                  }

                  if ( flags.gt17 || flags.d3 ) 
                  {
                     diffracted_wave_mag_sqr_sum_radial_intensity_local
                        = new double[ number_of_k_bins ];
                     //k_bin_counts_sqr_sum_local
                     //   = new double[ number_of_k_bins ];
                     for ( size_t ii=0; ii < number_of_k_bins; ++ii)
                     {
                    diffracted_wave_mag_sqr_sum_radial_intensity_local[ii]
                           = 0.0;
                        //k_bin_counts_sqr_sum_local = 0;
                     }
                     if (mynode == rootnode)
                     {
                        diffracted_wave_mag_sqr_sum_radial_intensity_total
                           = new double[ number_of_k_bins ];
                        //k_bin_counts_mag_sqr_sum_aggregated
                        //   = new double[ number_of_k_bins ];
                        for ( size_t ii=0; ii < number_of_k_bins; ++ii)
                        {
                   diffracted_wave_mag_sqr_sum_radial_intensity_total[ii]
                              = 0.0;
                           //k_bin_counts_sqr_sum_aggregated = 0;
                        }
                     }
                  } // ( flags.gt17 || flags.d3 ) 

                  if ( flags.gt17 )
                  {
                     diffracted_wave_mag_sum_sqr_radial_intensity_local
                        = new double[ number_of_k_bins ];
                     //k_bin_counts_sum_sqr_local   // TODO: rename
                     //   = new double[ number_of_k_bins ];
                     //for ( size_t ii=0; ii < number_of_k_bins; ++ii)
                     //{
                 //diffracted_wave_mag_sum_sqr_radial_intensity_local[ii]
                     //      = 0.0;
                     //   k_bin_counts_sum_sqr_local[ii] = 0;
                     //}
                     if ( mynode == rootnode )
                     {
                        ftem_gt17 = new double[ number_of_k_bins ];
                        diffracted_wave_mag_sum_sqr_radial_intensity_total
                           = new double[ number_of_k_bins ];
                        //k_bin_counts_sum_sqr_aggregated // TODO: rename
                        //   = new double[ number_of_k_bins ];
                        for ( size_t ii=0; ii < number_of_k_bins; ++ii)
                        {
                    diffracted_wave_mag_sum_sqr_radial_intensity_total[ii]
                              = 0.0;
                           //k_bin_counts_sum_sqr_aggregated[ii] = 0;
                        }
                     }
                  } // flags.gt17 

                  if ( flags.d3 || flags.rva )
                  {
                     diffracted_wave_mag_sum_radial_intensity_local
                        = new double[ number_of_k_bins ];
                     //k_bin_counts_sum_local // TODO: rename
                     //   = new double[ number_of_k_bins ];
                     for ( size_t ii=0; ii < number_of_k_bins; ++ii)
                     {
                      diffracted_wave_mag_sum_radial_intensity_local[ii] 
                           = 0.0;
                        //k_bin_counts_sum_local[ii] = 0;
                     }
                     if ( mynode == rootnode )
                     {
                        diffracted_wave_mag_sum_radial_intensity_total
                           = new double[ number_of_k_bins ];
                        //k_bin_counts_sum_aggregated // TODO: rename
                        //   = new double[ number_of_k_bins ];
                        for ( size_t ii=0; ii < number_of_k_bins; ++ii)
                        {
                       diffracted_wave_mag_sum_radial_intensity_total[ii]
                              = 0.0;
                           //k_bin_counts_sum_aggregated[ii] = 0;
                        }
                        if ( flags.d3 )
                        {
                           ftem_d3 = new double[ number_of_k_bins ];
                        }
                        if ( flags.rva )
                        {
                           ftem_rva = new double[ number_of_k_bins ];
                        }
                     }
                  } // flags.d3 || flags.rva

                  if ( flags.d4 )
                  {
                     diffracted_wave_mag_variance
                        = new double[ local_alloc_size_fftw ];
                     diffracted_wave_mag_variance_radial_intensity_local
                        = new double[ number_of_k_bins ];
                     //bin_count_variance_local
                     //   = new double[ number_of_k_bins ];
                     //for ( size_t ii=0; ii < number_of_k_bins; ++ii)
                     //{
                     //   diffracted_wave_mag_variance[ii]
                     //      = 0.0;
                     //   //bin_count_variance_local[ii] = 0;
                     //}

                     if ( mynode == rootnode )
                     {
                        ftem_d4 = new double[ number_of_k_bins ];
                       diffracted_wave_mag_variance_radial_intensity_total
                           = new double[ number_of_k_bins ];
                        //bin_count_variance_aggregated 
                        //   = new int[ number_of_k_bins ];
                        for ( size_t ii=0; ii < number_of_k_bins; ++ii)
                        {
                   diffracted_wave_mag_variance_radial_intensity_total[ii]
                              = 0.0;
                           //bin_count_variance_aggregated[ii] = 0;
                        }
                     }
                  } // flags.d4 
               } // ( flags.gt17 || flags.d3 || flags.d4) 
            } // flags.fem && first_probe_flag

            ////////////////////////////////////////////////////////////
            // update the look up table containing Zs paired with 
            //  projected atomic potentials centered at (0,0) 
            //   (I thought is was (xmin,ymin)?)
            ////////////////////////////////////////////////////////////

            // TODO: check for change in xperiod_duped and yperiod_duped
            //if ( xperiod_duped or yperiod_duped have changed )
            if ( update_domain_x || update_domain_y )
            {
               if ( papLUTexistsFlag == true ) delete myScattererPapLUT;

               myScattererPapLUT = new scatterer_pap_LUT(
                     Z_vector,
                     lambda, gamma, 1/(xperiod_duped * yperiod_duped), 
                     bwcutoff_pap,
                     kx_local, Nx_local,
                     Nx,
                     xmin,
                     ky, Ny, ymin,
                     local_alloc_size_fftw,
                     psi_mag_strides,
                     psi_mag_displacements,
                     mynode, rootnode, MPI_COMM_WORLD
                     );

               papLUTexistsFlag = true;
            }


            ///////////////////////////////////////////////////////////
            // update scatterer positions
            
            // clear any previous scatterers from the vector<scatterer>
            myScatterers.clear();
            // the above doesn't free q or the pap LUT member of 
            //  scatterers

            ////////////////////////////////////////////////////////////
            // Create the scatterer objects containing position, species, 
            //  and corresponding pointers to template projected atomic 
            //  potentials centered at (0,0) which will be translated to 
            //  the scatterer position by 
            //  slice.update_transmission_function() .
            ////////////////////////////////////////////////////////////
            // Create duplicate scatterers to increase FTEM V(|k|) 
            //  resolution
            if ( flags.dupe && flags.fem 
                 &&
                 (dupe_x > 1 || dupe_y > 1 || dupe_z > 1)
                 &&
                 (dupe_x > 0 && dupe_y > 0 && dupe_z > 0)
               )
            {
               if ( mynode == rootnode && flags.debug )
               {
                  cout << "Instantiating the given scatterer positions " 
                     << endl
                     << "   " << dupe_x << " times in the x direction, " 
                     << endl
                     << "   " << dupe_y << " times in the y direction, " 
                     << endl
                     << "   " << dupe_z << " times in the z direction" 
                     << endl;
               }

               duped_population 
                  = initial_population * dupe_x * dupe_y * dupe_z;

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
                                    *myScattererPapLUT
                                    )
                                 );
                        }
            }
            else
            {
               //q_duped = q;
               duped_population = initial_population;

               for (size_t i = 0; i < initial_population; i++)
               {
                  myScatterers.push_back(
                        scatterer( 
                           &(q[3*i]),
                           Z_array[i],
                           *myScattererPapLUT
                           )
                     );
               }
            }

            ///////////////////////////////////////////////////////////
            // update slices and slice members
            ///////////////////////////////////////////////////////////

            //delete_slices_from_list( slicesOfScatterers );
            //// the above frees propagator, transmission function,
            ////  and vector< scatterers* >* members
            ////  of all slices in slicesOfScatterers
            //slicesOfScatterers.clear();

            // Split the sample along the transmission axis into slices
            //  to be used by all tem simulations

            // TODO: parallelize assign_atoms_to_slices_z_auto() or 
            //       compute on a single node and broadcast the results.
            if ( (mynode == rootnode) && flags.debug)
               cout << "assigning atoms to slices ..." << endl;

            assign_atoms_to_slices_z_auto(
                  myScatterers,
                  xmin,ymin,zmin,
                  xperiod_duped, yperiod_duped, zperiod_duped,
                  minSliceThickness,  // [\AA]
                  slicesOfScatterers,
                  flags.debug,
                  mynode,
                  rootnode,
                  MPI_COMM_WORLD
                  );

            if ( (mynode == rootnode) && flags.debug)
               cout << "assigned atoms to slices" << endl;

            // Projected atomic potentials of each slice are not yet 
            // calculated. The next step should use the 
            // scatterer.pap_to_translate_{re,im} members to calculate the
            // total pap of each slice 
            //    (slicesOfScatterers[i].slice.exp_i_sigma_v[0,1]).
            // This is done within 
            //   (*sliceList_itr)->update_transmission_function()

            // TODO: if ( <atoms in a slice have been altered> )
            //   <re-evaluate relevant slices, else reuse previous slices>
            
            // NOTE: Scatterers contained in myScatterers are pointed to 
            //       from within members of slicesOfScatterers, so DO 
            //       NOT release memory of the scatterers pointed to by 
            //       myScatterers members until use of slicesOfScatterers
            //       is concluded.

            ////////////////////////////////////////////////////////////
            // Initialize slice members
            ////////////////////////////////////////////////////////////

            if ( mynode == rootnode && flags.debug )
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

               //if ((*sliceList_itr)->update_transmission_function_flag )
               //{
               //if ( flags.debug && (mynode == rootnode)) 
               //{
               //   cout << "updating transmission function, slice "
               //      << sliceNumber << endl;
               //   cout << " lambda " << lambda
               //      << endl << "gamma " << gamma 
               //      << endl <<  "1/(xperiod_duped * yperiod_duped) " << 1/(xperiod_duped * yperiod_duped) 
               //      << endl << "bwcutoff_t " << bwcutoff_t 
               //      << endl << "kx_local[0] " << kx_local[0]
               //      << endl << "Nx_local " << Nx_local
               //      << endl <<  "kx_joined[0] " << kx_joined[0]
               //         << endl <<  "xx_joined[0] "<<  xx_joined[0]
               //         << endl <<  "Nx " << Nx
               //         << endl <<  "ky[0] " << ky[0]
               //         << endl <<  "yy[0] " << yy[0]
               //         << endl <<  "Ny " << Ny 
               //         << endl <<  "NxNy " << NxNy
               //         << endl <<  "local_alloc_size_fftw " << local_alloc_size_fftw
               //         << endl <<  "local_0_start_fftw " << local_0_start_fftw
               //         << endl <<  "psi_mag_strides[0] " << psi_mag_strides[0]
               //         << endl <<  "psi_mag_displacements[0[ " << psi_mag_displacements[0]
               //         << endl;//debug
               //}
               if (
                     (*sliceList_itr)->update_transmission_function(
                        lambda, gamma, 
                        1/(xperiod_duped * yperiod_duped), // (ab)^{-1}
                        bwcutoff_t,
                        kx_local, Nx_local, 
                        kx_joined, xx_joined, Nx, 
                        ky, yy, Ny, 
                        NxNy,
                        //sqrtNxNy,
                        flags.pap_tif,
                        output_prefix,
                        sliceNumber,
                        local_alloc_size_fftw,  
                        local_0_start_fftw,
                        psi_mag_strides,
                        psi_mag_displacements,
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
            if ( mynode == rootnode && flags.debug )
               cout << "Evaluated slice members." << endl;

            //if ( xperiod_duped or yperiod_duped have changed or this is the first pass)
            if ( update_domain_x || update_domain_y || first_probe_flag )
            {
               //////////////////////////////////////////////////////
               // set up the STEM probe centered at (0,0)
               //////////////////////////////////////////////////////

               if ( flags.adfstem_corrected )
               {
                  if ( mynode == rootnode && flags.debug )
                     cout 
                        << "calling"
                      << " probe_wavefunction_correctedtoCs5_unnormalized"
                        << endl << "   Cs3 : " << Cs3
                        << endl << "   Cs5 : " << Cs5
                        << endl << "   defocus : " << defocus
                        << endl << "   alpha_max_sqr : "
                        << alpha_max_sqr 
                        << endl << "   lambda, lambda_sqr : "
                        << lambda << ", " << lambda_sqr
                        << endl << "   position : ("
                        << xx_joined[0] << ", " << yy[0] << ")" 
                        << endl << "   Nx_local, Ny : "
                        << Nx_local << ", " << Ny
                        << endl;

                  probe_wavefunction_correctedtoCs5_unnormalized(
                        xx_joined[0], yy[0],
                        kx_local, Nx_local, ky, Ny, 
                        Cs3, Cs5, defocus, alpha_max_sqr, 
                        lambda, lambda_sqr,
                        init_probe
                        );

                  if ( flags.debug )
                     cout << "Evaluated probe ..." << endl;// debug
               }
               else  // if ( flags.adfstem_uncorrected )
               {
                  if ( mynode == rootnode 
                        && flags.debug )
                     cout << "calling"
                        << " probe_wavefunction_uncorrected_unnormalized"
                        << endl << "   Cs3 : " << Cs3
                        << endl << "   defocus : " << defocus
                        << endl << "   alpha_max_sqr : "
                        << alpha_max_sqr 
                        << endl << "   lambda, lambda_sqr : "
                        << lambda << ", " << lambda_sqr
                        << endl << "   position : ("
                        << xx_joined[0] << ", " << yy[0] << ")" 
                        << endl << "   Nx_local, Ny: "
                        << Nx_local << ", " << Ny
                        << endl;

                  probe_wavefunction_uncorrected_unnormalized(
                        xx_joined[0], yy[0],
                        kx_local, Nx_local, ky, Ny, 
                        Cs3, defocus, alpha_max_sqr, 
                        lambda, lambda_sqr,
                        init_probe
                        );

                  if ( flags.debug )
                     cout << "Evaluated probe ..." << endl;// debug
               }

               // Transform the split probe into realspace
               if ( mynode == rootnode && flags.debug )
                  cout << "Transforming probe to realspace ..." << endl;

               fftw_execute( pb_c2c_probe_split );

               // Gather the probe realspace wavefunction onto each node

               for (ptrdiff_t i=0; i < local_alloc_size_fftw; ++i)
               {  // reassign probe values for MPI_Allgather while 
                  //  removing fftw scale factor
                  probe_split_re[i] = init_probe[i][0] / sqrtNxNy;
                  probe_split_im[i] = init_probe[i][1] / sqrtNxNy;
               }

               if ( mynode == rootnode && flags.debug )
                  cout << "Allgathering probe ..." << endl;

               MPI_Allgatherv(
                     probe_split_re,         // sendbuf
                     Nx_local * Ny, //local_alloc_size_fftw,  // sendcount
                     MPI_DOUBLE,             // sendtype
                     probe_joined_re,        // recvbuf
                     psi_mag_strides,        // recvcounts[]
                     psi_mag_displacements,  // displs[]
                     MPI_DOUBLE,             // recvtype
                     comm);                  // MPI_Comm comm

               MPI_Allgatherv(
                     probe_split_im,         // sendbuf
                     Nx_local * Ny, //local_alloc_size_fftw,  // sendcount
                     MPI_DOUBLE,             // sendtype
                     probe_joined_im,        // recvbuf
                     psi_mag_strides,        // recvcounts[]
                     psi_mag_displacements,  // displs[]
                     MPI_DOUBLE,             // recvtype
                     comm);                  // MPI_Comm comm

               // Normalize the probe in real space
               //  Why not normalize before MPI_AllGatherv ?
               //  Because the split versions may overlap.

               // Periodic boundary conditions ensure that the probe norm
               //  is invariant with probe position.
               if ( mynode == rootnode && flags.debug )
                  cout << "Calculating probe norm in real space ... " 
                     << endl;

               double Ap; Ap = 0.0; // probe normalization factor
               probe_wavefunction_norm(
                     Ap, 
                     Nx_local, Ny, 
                     probe_joined_re,
                     probe_joined_im,
                     comm
                     );
               //if (mynode == rootnode && flags.debug ) // debug
               //   cout << "Probe norm in real space : " << Ap << endl;

               for ( ptrdiff_t i=0; i< Nx * Ny; ++i)
               {
                  probe_joined_re[i] = Ap * probe_joined_re[i];
                  probe_joined_im[i] = Ap * probe_joined_im[i];
               }
            } // end of updating unshifted probe in response to domain 
            //    changes

            //if (update_domain_x || update_domain_y || fresh_raster_point)
            //{
            ////////////////////////////////////////////////////////
            // translate the initial STEM probe w.r.t. the updated 
            //  domain
            //  and assign it to psi
            ////////////////////////////////////////////////////////
            // TODO: the following is inefficient, 
            //        in the future save the probe_idx_shift_... in an
            //        array to avoid recalculation
            // Determine indices of xx_joined[] and yy[] which 
            //  represent  the point closest to (*x_itr, *y_itr).
            // The following assumes that 
            //  xx_joined[0] == min(xx_joined[]),
            //  and that the probe is centered somewhere on the sample
            //       *x_itr \in [xx_joined[0], xx_joined[Nx - 1]

            for (size_t i=0; i < Nx; ++i)
            {
               if ( xx_joined[i] - *x_itr >= half_delta_x )
               {
                  probe_idx_shift_x = i-1;
                  break;
               }
            }

            for (size_t j=0; j < Ny; ++j)
            {
               if ( yy[j] - *y_itr >= half_delta_y )
               {
                  probe_idx_shift_y = j-1;
                  break;
               }
            }
            // TODO: test that the probe center is within the domains.

            //while ((idx_shift_x + 1 < Nx) 
            //         && (xx_joined[ idx_shift_x + 1] < *x_itr))
            //   ++idx_shift_x;
            //// state: xx_joined[ idx_shift_x + 1] >= *x_itr or idx_shift_x + 1 >= Nx
            //if ( ( *x_itr - xx_joined[ idx_shift_x ]) 
            //      > 
            //      (xx_joined[ idx_shift_x + 1 ] - *x_itr) )
            //{
            //   ++idx_shift_x;
            //}
            //if ( idx_shift_x + 1 >= Nx ) idx_shift_x = Nx - 1;

            // xx_joined[ idx_shift_x ] is nearest to *x_itr
            // yy[ idx_shift_y ] is nearest to *y_itr

            // copy probe_joined_re,im to psi[][0,1] with 
            //  shifted x and y indices 

            // Those indices should be the same as the number of 
            //  positions which the cached probe array should be 
            //  shifted.

            // Iterate over the domains, assign probe locally to
            //  the variable which will be used in the fft.
            for ( size_t i=0; i < Nx_local; ++i)
            {
               // enforce periodic boundary condition, x direction
               if ( i + local_0_start_fftw < probe_idx_shift_x )
                  probe_idx_x = i + Nx + local_0_start_fftw 
                                 - probe_idx_shift_x;
               else
                  probe_idx_x = i + local_0_start_fftw 
                                 - probe_idx_shift_x;

               for ( size_t j=0; j < Ny; ++j)
               {
                  // enforce periodic boundary condition, y direction 
                  if ( j < probe_idx_shift_y ) 
                     probe_idx_y = j + Ny - probe_idx_shift_y;
                  else
                     probe_idx_y = j - probe_idx_shift_y;

                  // debug
                  if (
                        probe_idx_x < 0 || probe_idx_x > Nx
                        || probe_idx_y < 0 || probe_idx_y > Ny
                        || probe_idx_y + probe_idx_x * Ny > Nx * Ny
                     )
                     cout << "Error, (probe_idx_x, probe_idx_y,"
                        << " probe_idx_y + probe_idx_x*Ny,"
                        << " local_alloc_size_fftw): ("
                        << probe_idx_x << ", " << probe_idx_y 
                        << ", " << probe_idx_y + probe_idx_x * Ny 
                        << ", " << local_alloc_size_fftw
                        << ") exceeds bounds (Nx,Ny)"
                        << Nx << ", " << Ny << ")" << endl;
                  // end debug

                  psi[j + i * Ny][0] = 
                     probe_joined_re[ probe_idx_y 
                                             + probe_idx_x * Ny ];

                  psi[j + i * Ny][1] = 
                     probe_joined_im[ probe_idx_y 
                                             + probe_idx_x * Ny ];
               }
            }

            if ( flags.debug && flags.image_output )
            {
               if ( mynode == rootnode && flags.debug )  // debug
                  cout << "saving initial probe to tiff" << endl; 

               debug_output_complex_fftw_operand(//TODO: rename this function
                  psi,
                  0, 
                  local_alloc_size_fftw,
                  Nx_local, Nx, Ny,
                  resolutionUnit,
                  xResolution, yResolution,
                  output_prefix 
                  + "_initial_probe_realspace",
                  psi_mag_strides,
                  psi_mag_displacements,
                  mynode, rootnode, comm);
            }

            // transform translated beam into reciprocal space for 
            //  slice loop
            fftw_execute( pf_c2c_psi ); 

            for ( ptrdiff_t i=0; i<local_alloc_size_fftw; ++i)
            {
               psi[i][0] = psi[i][0] / sqrtNxNy;
               psi[i][1] = psi[i][1] / sqrtNxNy;
            }
            //}// end of updating shifted probe 

            ///////////////////////////////////////////////////////////
            // propagate
            ///////////////////////////////////////////////////////////
            sliceNumber = 0; //debug
            for(
                std::list< slice* >::iterator sliceList_itr
                  = slicesOfScatterers.begin();
                sliceList_itr != slicesOfScatterers.end(); 
                ++sliceList_itr )
            {
               ++sliceNumber;//debug
               fftw_execute( pb_c2c_psi );
               if ( mynode == rootnode && flags.debug )
               {
                  cout << "multiplying psi by "
                     << "transmission function, slice : "
                     << sliceNumber
                     << ", node : " << mynode << endl;
               }
               // limit the probe in realspace, to prevent aliasing 
               //  artifacts from appearing in reciprocal space
               for ( ptrdiff_t i=0; i < Nx_local; i++)
               {
                  for ( ptrdiff_t j=0; j < Ny; j++)
                  {
                     psi[(j + i * Ny)][0]
                        = (
                           (*sliceList_itr)->exp_i_sigma_v[j + i * Ny][0] 
                           * psi[(j + i * Ny)][0]
                           - 
                           (*sliceList_itr)->exp_i_sigma_v[j + i * Ny][1] 
                           * psi[(j + i * Ny)][1]
                          ) / sqrtNxNy;// NxNy; // / NxNy_sqr; 
                     // normalizing due to pf_c2c_t, pb_c2c_t, pf_c2c_psi,
                     //  pb_c2c_psi
                     // NOTE: transmission_function_ms() already has a 
                     //       factor of 1/sqrtNxNy

                     psi[(j + i * Ny)][1]
                        = (
                           (*sliceList_itr)->exp_i_sigma_v[j + i * Ny][0] 
                           * psi[(j + i * Ny)][1]
                           + 
                           (*sliceList_itr)->exp_i_sigma_v[j + i * Ny][1] 
                           * psi[(j + i * Ny)][0]
                          ) / sqrtNxNy;// NxNy; // / NxNy_sqr; 

                  }
               }
               fftw_execute( pf_c2c_psi );
               for ( ptrdiff_t i=0; i < local_alloc_size_fftw; i++)
               {
                  psi[i][0] = psi[i][0] /sqrtNxNy;
                  psi[i][1] = psi[i][1] /sqrtNxNy;
               }

               bw_limit(
                     psi,
                     Nx_local, kx_local, Ny, ky,
                     bwcutoff_t   
                     );

               (*sliceList_itr)->propagate(
                     lambda,
                     kx_local, Nx_local, ky, Ny, 
                     bwcutoff_t,
                     psi
                     );

            } // end of iteration over slices

            // Assume that the slices will need to be rebuilt after lammps
            //  and clear their memory pre-emptively
            delete_slices_from_list( slicesOfScatterers );
            // the above frees propagator, transmission function,
            //  and vector< scatterers* >* members
            //  of all slices in slicesOfScatterers
            slicesOfScatterers.clear();


            if ( mynode == rootnode && flags.debug )
               cout << "finished propagating through slices,"
                 << " probe centered at : (x,y) = ( " 
                 << *x_itr << ", " << *y_itr << ")" << endl;
            
            if ( flags.debug && first_probe_flag && flags.image_output )
            {
               diffraction_scale_factor = 1.0e-5;
               output_diffraction(
                     psi,
                     diffraction_scale_factor,
                     local_alloc_size_fftw,
                     Nx_local, Nx, Ny,
                     resolutionUnit_recip,
                     xResolution_recip, yResolution_recip,
                     output_prefix + "_first",
                     psi_mag_strides,
                     psi_mag_displacements,
                     mynode, rootnode, comm
                     );
               diffraction_scale_factor = 1.0e-1;
               output_diffraction(
                     psi,
                     diffraction_scale_factor,
                     local_alloc_size_fftw,
                     Nx_local, Nx, Ny,
                     resolutionUnit_recip,
                     xResolution_recip, yResolution_recip,
                     output_prefix + "_first",
                     psi_mag_strides,
                     psi_mag_displacements,
                     mynode, rootnode, comm
                     );
               diffraction_scale_factor = 1.0e0;
               output_diffraction(
                     psi,
                     diffraction_scale_factor,
                     local_alloc_size_fftw,
                     Nx_local, Nx, Ny,
                     resolutionUnit_recip,
                     xResolution_recip, yResolution_recip,
                     output_prefix + "_first",
                     psi_mag_strides,
                     psi_mag_displacements,
                     mynode, rootnode, comm
                     );
               diffraction_scale_factor = 1.0e10;
               output_diffraction(
                     psi,
                     diffraction_scale_factor,
                     local_alloc_size_fftw,
                     Nx_local, Nx, Ny,
                     resolutionUnit_recip,
                     xResolution_recip, yResolution_recip,
                     output_prefix + "_first",
                     psi_mag_strides,
                     psi_mag_displacements,
                     mynode, rootnode, comm
                     );
               diffraction_scale_factor = 1.0e+20;//1e-20;
               output_diffraction(
                     psi,
                     diffraction_scale_factor,
                     local_alloc_size_fftw,
                     Nx_local, Nx, Ny,
                     resolutionUnit_recip,
                     xResolution_recip, yResolution_recip,
                     output_prefix + "_first",
                     psi_mag_strides,
                     psi_mag_displacements,
                     mynode, rootnode, comm
                     );
            }
            first_probe_flag = false;

            //--------------------------------------------------------
            // End of beam transmission through individual slices
            //--------------------------------------------------------

            ///////////////////////////////////////////////////////////
            // accumulate detected wavefunction

            //////////////////////////////////////////////////////////
            // STEM: integrate the diffracted intensity impingent upon 
            //       the annular detector
            //////////////////////////////////////////////////////////

            if ( flags.image_output )
            {
               detected_intensity = 
                  sum_detected_intensity( 
                        psi, 
                        kx_local, ky,
                        detector_inner_angle / lambda,
                        detector_outer_angle / lambda,
                        Nx_local, Ny
                        );

               MPI_Reduce( 
                     &detected_intensity, 
                     &detected_intensity_reduced,
                     1, MPI_DOUBLE,
                     MPI_SUM, rootnode, comm );

               if( mynode == rootnode )
               {
                  stem_image[pixel_number_y + pixel_number_x * y_p.size()]
                     += detected_intensity_reduced / lammps_TEM_samples;
               }
            }

            ////////////////////////////////////////////////////////////
            // Apply modulation transfer function of the image detector
            ////////////////////////////////////////////////////////////
            if ( flags.mtf_file && flags.mtf_resolution )
            {
               // transform psi into its magnitude
               //  and apply beam blocker to direct beam
               double ksqr;
               for ( size_t i=0; i < Nx_local; ++i)
               {
                  for ( size_t j=0; j < Ny; ++j)
                  {
                     ksqr = pow(kx_local[i], 2) + pow(ky[j], 2);
                     if ( lambda_sqr * ksqr > alpha_max_sqr )
                     {
                        psi[j + i*Ny][0] = 
                           sqrt( (psi[j + i*Ny][0] *psi[j + i*Ny][0])
                                 + (psi[j + i*Ny][1] * psi[j + i*Ny][1]));
                     }
                     else
                     {
                        psi[j + i*Ny][0] =  0.0;
                     }
                     psi[j + i*Ny][1] = 0.0;
                  }
               }

               //bw_limit(
               //      psi,
               //      Nx_local, kx_local, Ny, ky,
               //      //0.5*bwcutoff_t
               //      bwcutoff_t
               //      );
               
               // transform psi into real-space
               fftw_execute( pb_c2c_psi );

               if ( flags.debug && flags.image_output)
               {
                  diffraction_scale_factor = 1.0e+10;
                  output_diffraction_append(
                        psi,
                        diffraction_scale_factor,
                        lambda_sqr,
                        alpha_max_sqr,
                        local_alloc_size_fftw,
                        Nx_local, 
                        kx_local,
                        Nx, 
                        Ny,
                        ky,
                        resolutionUnit_recip,
                        xResolution_recip, yResolution_recip,
                        output_prefix + "_psi_mtf2D_recipspace",
                        psi_mag_strides,
                        psi_mag_displacements,
                        mynode, rootnode, comm
                        );
               }

               
               for (size_t i=0; i < Nx_local; ++i)
               {
                  for ( size_t j=0; j < Ny; ++j)
                  {
                     // multiply psi by the mtf
                     psi[j + i * Ny][0]
                        = (
                           (psi[j + i * Ny][0])
                           *(mtf_2D_split[j + i * Ny][0])
                          - (psi[j + i * Ny][1])
                            *(mtf_2D_split[j + i * Ny][1])
                          )/sqrtNxNy;

                     psi[j + i * Ny][1]
                        = (
                           (psi[j + i * Ny][0])
                            *(mtf_2D_split[j + i * Ny][1])
                           + (psi[j + i * Ny][1])
                            *(mtf_2D_split[j + i * Ny][0])
                          )/sqrtNxNy;
                  }
               }
               //if ( flags.debug && flags.image_output )
               //{
               //   diffraction_scale_factor = 1.0e+10;
               //   output_diffraction(
               //         psi,
               //         diffraction_scale_factor,
               //         //lambda_sqr,
               //         //alpha_max_sqr,
               //         local_alloc_size_fftw,
               //         Nx_local, //kx_local, 
               //         Nx, 
               //         //kx_joined, 
               //         Ny, //ky,
               //         resolutionUnit_recip,
               //         xResolution_recip, yResolution_recip,
               //         output_prefix + "_psi_mtf2D_realspace",
               //         psi_mag_strides,
               //         psi_mag_displacements,
               //         mynode, rootnode, comm
               //         );
               //}

               //  transform psi back into reciprocal space
               fftw_execute( pf_c2c_psi );
               // renormalize
               for (size_t i=0; i < local_alloc_size_fftw; ++i)
               {
                     psi[i][0] = psi[i][0] / sqrtNxNy;
                     psi[i][1] = psi[i][1] / sqrtNxNy;
               }
               bw_limit(
                     psi,
                     Nx_local, kx_local, Ny, ky,
                     bwcutoff_t
                     );
            } // end of mtf application

            if ( flags.diffraction_output )
            {
               if ( flags.debug && mynode == rootnode )
                  cout << "appending diffraction to tiff" << endl;//debug

               diffraction_scale_factor = 1.0e+13;
               output_diffraction_append(
                     psi,
                     diffraction_scale_factor,
                     lambda_sqr,
                     alpha_max_sqr,
                     local_alloc_size_fftw,
                     Nx_local, 
                     kx_local,
                     Nx, 
                     Ny,
                     ky,
                     resolutionUnit_recip,
                     xResolution_recip, yResolution_recip,
                     output_prefix 
                        + "_raster_position_"
                        + to_string(pixel_number_x) 
                        + "_"
                        + to_string(pixel_number_y) ,
                     psi_mag_strides,
                     psi_mag_displacements,
                     mynode, rootnode, comm
                     );
                     
               diffraction_scale_factor = 1.0;
               output_diffraction_append(
                     psi,
                     diffraction_scale_factor,
                     lambda_sqr,
                     alpha_max_sqr,
                     local_alloc_size_fftw,
                     Nx_local, 
                     kx_local,
                     Nx, 
                     Ny,
                     ky,
                     resolutionUnit_recip,
                     xResolution_recip, yResolution_recip,
                     output_prefix 
                        + "_raster_position_"
                        + to_string(pixel_number_x) 
                        + "_"
                        + to_string(pixel_number_y) ,
                     psi_mag_strides,
                     psi_mag_displacements,
                     mynode, rootnode, comm
                     );
                     
               //if ( flags.netcdf_images )
               //{
               //   if ( mynode == rootnode && flags.debug)
               //      cout << "saving initial probe to netCDF" << endl;

               //   if( output_psi_realspace_to_netcdf(
               //         psi,
               //         local_alloc_size_fftw,
               //         Nx_local, kx_joined, Nx, ky, Ny,
               //         output_prefix 
               //            + "_diffraction_raster_position_"
               //            + to_string(pixel_number_x) 
               //            + "_"
               //            + to_string(pixel_number_y) ,
               //         psi_mag_strides,
               //         psi_mag_displacements,
               //         mynode, rootnode, comm
               //         ) == EXIT_FAILURE)
               //      cout << "output_psi_realspace_to_netcdf() failed" 
               //         << endl;
               //}
            }

            for (size_t i=0; i < Nx_local; ++i)
            {
               for ( size_t j=0; j < Ny; ++j)
               {
                  psi_accumulation[j + i * Ny][0] 
                     += psi[j + i * Ny][0] / lammps_TEM_samples;
                  psi_accumulation[j + i * Ny][1] 
                     += psi[j + i * Ny][1] / lammps_TEM_samples;
               }
            }

            update_domain_x = false;
            update_domain_y = false;
            fresh_raster_point = false;

         } // end of loop over lammps runs

         // save accumulated diffraction
         if ( flags.diffraction_output )
         {
            if ( flags.debug && mynode == rootnode )
               cout << "appending diffraction to tiff" << endl;//debug

            diffraction_scale_factor = 1.0e+13;
            output_diffraction_append(
                  psi,
                  diffraction_scale_factor,
                  lambda_sqr,
                  alpha_max_sqr,
                  local_alloc_size_fftw,
                  Nx_local, 
                  kx_local,
                  Nx, 
                  Ny,
                  ky,
                  resolutionUnit_recip,
                  xResolution_recip, yResolution_recip,
                  output_prefix + "_accumulated",
                  psi_mag_strides,
                  psi_mag_displacements,
                  mynode, rootnode, comm
                  );
                  
            diffraction_scale_factor = 1.0;
            output_diffraction_append(
                  psi,
                  diffraction_scale_factor,
                  lambda_sqr,
                  alpha_max_sqr,
                  local_alloc_size_fftw,
                  Nx_local, 
                  kx_local,
                  Nx, 
                  Ny,
                  ky,
                  resolutionUnit_recip,
                  xResolution_recip, yResolution_recip,
                  output_prefix + "_accumulated",
                  psi_mag_strides,
                  psi_mag_displacements,
                  mynode, rootnode, comm
                  );
                  
            //if ( flags.netcdf_images )
            //{
            //   if ( mynode == rootnode && flags.debug)
            //      cout << "saving initial probe to netCDF" << endl;

            //   if( output_psi_realspace_to_netcdf(
            //         psi,
            //         local_alloc_size_fftw,
            //         Nx_local, kx_joined, Nx, ky, Ny,
            //         output_prefix 
            //            + "_diffraction_raster_position_"
            //            + to_string(pixel_number_x) 
            //            + "_"
            //            + to_string(pixel_number_y) ,
            //         psi_mag_strides,
            //         psi_mag_displacements,
            //         mynode, rootnode, comm
            //         ) == EXIT_FAILURE)
            //      cout << "output_psi_realspace_to_netcdf() failed" 
            //         << endl;
            //}
         }

         ///////////////////////////////////////////////////////////////
         // accumulate FTEM quantities
         if ( flags.fem )
         {
            if ( (mynode == rootnode) && flags.debug )
               cout << "calculating and accumulating FTEM quantities ..." 
                  << endl;
            if ( flags.d1 || flags.d2 || flags.correlograph )
            {
               // Azimuthally integrate diffraction for 1-D variance.
               if (integrate_out_phi_fftw(
                     psi_accumulation,
                     flags,
                     kx_local, Nx_local,
                     ky, Ny,
                     indexed_magnitudes,
                     k_binning_boundaries, 
                     k_bin_counts_local,
                     phi_binning_boundaries, 
                     phi_bin_counts_local,
                     diffracted_wave_mag_radial_intensity_local,
                     diffracted_wave_mag_in_radial_coords_local
                     ) == EXIT_FAILURE)
               {
                  cerr << "integrate_out_phi_fftw failed" << endl;
                  return EXIT_FAILURE;
               }

               // zero direct beam region of azimuthally integrated diffraction
               size_t number_of_phi_bins = phi_binning_boundaries.size() -1;
               for ( size_t ii=0; ii < number_of_k_bins; ++ii)
               {
                  if ( lambda * k_binning_boundaries[ ii + 1] <= alpha_max)
                  {
                     diffracted_wave_mag_radial_intensity_local[ ii] = 0.0;
                     if ( flags.correlograph)
                     {
                        for ( size_t jj=0; jj < number_of_phi_bins; ++jj)
                        {
                           // zero correlograph of direct beam region
                           diffracted_wave_mag_in_radial_coords_local[
                              jj + ii*number_of_phi_bins ] = 0.0;
                        }
                     }
                  }
               }

               // TODO: is the following initialization necessary?
               // No.
               //if ( mynode == rootnode )
               //   for (unsigned int ii=0; ii < number_of_k_bins; ++ii)
               //      diffracted_wave_radial_intensity_sum_d1_tmp[ii]=0.0;
               if ( number_of_k_bins != k_binning_boundaries.size() - 1)//debug
               {//debug
                  cerr << "Error - number_of_k_bins != "
                     << "k_binning_boundaries.size() - 1" << endl;//debug
                  // TODO: make a function to release memory for a clean 
                  //        exit
               }//debug

               MPI_Reduce(
                     diffracted_wave_mag_radial_intensity_local,
                     diffracted_wave_mag_radial_intensity_total,
                     number_of_k_bins,
                     MPI_DOUBLE, MPI_SUM,
                     rootnode, comm
                     );

               MPI_Reduce(
                     k_bin_counts_local,
                     k_bin_counts_aggregated,
                     number_of_k_bins,
                     MPI_INT, MPI_SUM,
                     rootnode, comm
                     );

               if ( flags.correlograph )
               {
                  MPI_Reduce(
                        diffracted_wave_mag_in_radial_coords_local,
                        diffracted_wave_mag_in_radial_coords_total,
                        number_of_k_bins * number_of_phi_bins,
                        MPI_DOUBLE, MPI_SUM,
                        rootnode, comm
                        );

                  MPI_Reduce(
                        phi_bin_counts_local,
                        phi_bin_counts_aggregated,
                        number_of_phi_bins,
                        MPI_INT, MPI_SUM,
                        rootnode, comm
                        );
               }

               if ( mynode == rootnode )
               {
                  if ( flags.d1 )
                  {
                     for ( unsigned int ii=0; ii < number_of_k_bins; ++ii)
                     {
                        // division by k_bin_counts_aggregated implements 
                        //  azimuthal averaging 
                    diffracted_wave_mag_radial_intensity_total_sqr_sum[ii]
                       += (diffracted_wave_mag_radial_intensity_total[ii]
                        * diffracted_wave_mag_radial_intensity_total[ii])
                                 / (k_bin_counts_aggregated[ii] 
                                       * k_bin_counts_aggregated[ii]);
                        
                    diffracted_wave_mag_radial_intensity_total_sum[ii]
                       += diffracted_wave_mag_radial_intensity_total[ii] 
                              / k_bin_counts_aggregated[ii];
                     }
                  }

                  if ( flags.correlograph ) // at rootnode
                  {
                     // k_binning_boundaries,
                     // k_bin_counts_aggregated, // # of points in rings
                     // phi_binning_boundaries,
                     // phi_bin_counts_aggregated,// points of small areas
                     // diffracted_wave_mag_radial_intensity_local
                     // diffracted_wave_mag_in_radial_coords_local
                     // jj indexes phi, ii indexes k
                     size_t idx; // TODO: is this efficient?
                     size_t c_idx;
                     double tmpdbl;
                     for (size_t ii=0; ii < number_of_k_bins; ++ii)
                     {
                        for (size_t jj=0; jj < number_of_phi_bins; ++jj)
                        {
                           idx = jj + ii * number_of_phi_bins;
                           c_idx = ii + jj * number_of_k_bins;
                           tmpdbl 
                                 = (k_bin_counts_aggregated[ii] 
                        * diffracted_wave_mag_in_radial_coords_total[idx]
                        / diffracted_wave_mag_radial_intensity_total[ii]);

                           if ( (tmpdbl >= 1) 
                                 && (k_bin_counts_aggregated[ii] !=0) 
                      && (diffracted_wave_mag_in_radial_coords_total[idx] 
                                    !=0)
                      && (diffracted_wave_mag_radial_intensity_total[ii] 
                                    !=0)
                           )
                           {
                              correlograph[c_idx]
                                 = (k_bin_counts_aggregated[ii] 
                         * diffracted_wave_mag_in_radial_coords_total[idx]
                         / diffracted_wave_mag_radial_intensity_total[ii])
                                 - 1;
                           }
                           else
                              correlograph[c_idx] = 0.0;

                           // accumulate average correlograph
                           correlograph_sum[c_idx] += correlograph[c_idx];
                           if ( flags.correlograph_variance )
                              correlograph_sqr_sum[c_idx] 
                                 += (correlograph[c_idx]) 
                                    * (correlograph[c_idx]);
                        }
                     }
                     kResolution_correlograph //limiting domain to (0.2,1)
                        = ( number_of_k_bins / (1.0 - 0.2)) * 1.0e8;
                     phiResolution_correlograph 
                        = number_of_phi_bins / (360);// degrees for output

                     if ( (mynode == rootnode) 
                           && flags.correlograph_everyimage )
                     {
                        output_correlograph_image(
                              correlograph,
                              number_of_phi_bins,
                              number_of_k_bins,
                              resolutionUnit,
                              phiResolution_correlograph,
                              kResolution_correlograph,
                              output_prefix 
                              + "_" + to_string(pixel_number_x) 
                              + "_" + to_string(pixel_number_y),
                              0
                              );
                     }

                     if ( (mynode == rootnode) 
                           && flags.correlograph_everytxt )
                     {
                        output_correlograph_to_txt(
                              correlograph,
                              phi_binning_boundaries,
                              k_binning_boundaries,
                              output_prefix 
                                 + "_" + to_string(pixel_number_x) 
                                 + "_" + to_string(pixel_number_y)
                              );
                     }
                  }
               }

               if ( flags.d2 )
               {
                  for ( ptrdiff_t ii=0; ii < local_alloc_size_fftw; ++ii)
                  {
                        diffracted_wave_mag_sqr[ii] 
                           = (psi_accumulation[ii][0] 
                                 * psi_accumulation[ii][0]) 
                              + (psi_accumulation[ii][1] 
                                 * psi_accumulation[ii][1]);
                  }

                  integrate_out_phi_double(
                        diffracted_wave_mag_sqr,
                        kx_local, Nx_local,
                        ky, Ny,
                        indexed_magnitudes,
                        k_binning_boundaries, 
                        k_bin_counts_local,
                        diffracted_wave_mag_sqr_radial_intensity_local
                        );

                  // zero the region of the direct beam
                  for ( size_t ii=0; ii < number_of_k_bins; ++ii)
                  {
                     if ( lambda * k_binning_boundaries[ ii + 1] <= alpha_max)
                     {
                        diffracted_wave_mag_sqr_radial_intensity_local[ ii] = 0.0;
                     }
                  }

                  MPI_Reduce(
                        diffracted_wave_mag_sqr_radial_intensity_local,
                        diffracted_wave_mag_sqr_radial_intensity_total,
                        number_of_k_bins,
                        MPI_DOUBLE, MPI_SUM,
                        rootnode, comm
                        );

                  MPI_Reduce(
                        k_bin_counts_local,
                        k_bin_counts_sqr_sum_aggregated,
                        number_of_k_bins,
                        MPI_INT, MPI_SUM,
                        rootnode, comm
                        );
                  if ( mynode == rootnode )
                  {
                     for ( unsigned int ii=0; ii<number_of_k_bins; ++ii)
                     {
                        diffracted_wave_variance_sum[ii] += 
                   ((diffracted_wave_mag_sqr_radial_intensity_total[ii]
                           / k_bin_counts_sqr_sum_aggregated[ii])
                           /
                           (
                   (diffracted_wave_mag_radial_intensity_total[ii]
                     * diffracted_wave_mag_radial_intensity_total[ii])
                           / (k_bin_counts_aggregated[ii] 
                              * k_bin_counts_aggregated[ii])
                           )) - 1.0;

                        k_bin_counts_aggregated[ii] = 0;
                        k_bin_counts_sqr_sum_aggregated[ii] = 0;
                     }
                  }
               } // flags.d2
            } // flags.d1 || flags.d2 || flags.correlograph

            if ( flags.gt17 || flags.d3 || flags.d4
                  || flags.rva ) 
            {
               // Accumulate values for 2-D variance image.
               for ( ptrdiff_t ii=0; ii < local_alloc_size_fftw; ++ii )
               {
                  diffracted_wave_mag_sum[ii] 
                     += 
                     sqrt(psi_accumulation[ii][0] 
                          * psi_accumulation[ii][0] 
                            + psi_accumulation[ii][1] 
                              * psi_accumulation[ii][1]);

                  if ( flags.gt17 || flags.d3 || flags.d4) 
                  {
                     diffracted_wave_mag_sqr_sum[ii] 
                        += psi_accumulation[ii][0] 
                           * psi_accumulation[ii][0] 
                             + psi_accumulation[ii][1] 
                               * psi_accumulation[ii][1];
                  }
                     //(psi_accumulation[ii][0] 
                     //  * psi_accumulation[ii][0] 
                     //   + psi_accumulation[ii][1] 
                     //     * psi_accumulation[ii][1])
                     //   * (psi_accumulation[ii][0] 
                     //     * psi_accumulation[ii][0] 
                     //   + psi_accumulation[ii][1] 
                     //     * psi_accumulation[ii][1]);
                  
                  if (flags.complex_realspace_sum) 
                  {
                     diffracted_wave_re_sum[ii] 
                        += psi_accumulation[ii][0] ;
                     diffracted_wave_im_sum[ii] 
                        += psi_accumulation[ii][1] ;
                  }
               }
            } // flags.gt17 || flags.d3 || flags.d4
         } // end of flags.fem dependent block

         // TODO:

         ++pixel_number_y; 
      }
      ++pixel_number_x; 
   } // end of loop over raster points

   //=================================================================
   // End of STEM beam position rastering
   //=================================================================


   ////////////////////////////////////////////////////////////////
   // Calculation of FTEM statistical measures
   ////////////////////////////////////////////////////////////////
   if ( flags.fem )
   {
      if ((mynode == rootnode) && flags.debug)
         cout << "calculating post-raster FTEM quantities ..." 
            << endl;
      if ( flags.d1 )
      {
         // V(|\vec{k}|,K_{aperture})
         //  = \frac{ \langle I^{2}(|\vec{k}|,K_{aperture}) \rangle }
         //          { \langle I(|\vec{k}|,K_{aperture}) \rangle^{2} } -1
         // where I() and I^{2}() are 
         // diffracted_wave_mag_radial_intensity_total_sum and 
         // diffracted_wave_mag_radial_intensity_total_sqr_sum, 
         // respectively.

         if ( mynode == rootnode )
         {
            for ( size_t idx=0; idx < number_of_k_bins; ++idx)
            {
               if ( diffracted_wave_mag_radial_intensity_total_sum[idx] 
                     == 0.0 )
               {
                  ftem_d1[idx] = 0.0;
               }
               else
               {
                  ftem_d1[idx] 
                     = ( 
                          number_of_raster_points 
                          *  
                 diffracted_wave_mag_radial_intensity_total_sqr_sum[idx]
                     / ( 
                   diffracted_wave_mag_radial_intensity_total_sum[idx]
                 * diffracted_wave_mag_radial_intensity_total_sum[idx]
                    )
                       ) - 1.0;
               }

               if ( ftem_d1[idx] < 0.0 )
               {
                  cout << " WARNING, fem --V_omega, ftem_d1[" 
                     << idx << "] < 0.0 : "
                     << ftem_d1[idx] << endl;
               }
            }

            // write the V_omega variance to a file
            //if( flags.netcdf_variance )
            //{
            //   if ( 
            //         output_variance_to_netcdf(
            //            ftem_d1,
            //            k_binning_boundaries,
            //            output_prefix + "_fem_V_omega"
            //         ) != EXIT_SUCCESS)
            //   {
            //      cout << " failed to write variance data to a"
            //        << " netCDF file : " << output_prefix << endl;
            //   }
            //} else {
            if ( 
               output_variance_to_txt(
                     ftem_d1,
                     k_binning_boundaries,
                     output_prefix + "_fem_V_omega"
                  ) != EXIT_SUCCESS)
            {
               cout << " failed to write variance data to a"
                  << " txt file : " << output_prefix << endl;
            }
            //}
         } // rootnode
      } // flags.d1

      if ( flags.d2 )
      {
         if ( mynode == rootnode )
         {
            for (size_t ii=0; ii < number_of_k_bins; ++ii)
            {
               diffracted_wave_variance_sum[ii]
                  = diffracted_wave_variance_sum[ii] 
                     / number_of_raster_points;

               if ( diffracted_wave_variance_sum[ii] < 0.0 )
               {
                  cout << " WARNING, fem --Vbar_r,"
                     << " diffracted_wave_variance_sum[" 
                     << ii << "] < 0.0 : "
                     << diffracted_wave_variance_sum[ii] << endl;
               }
            }
            // write the d2 variance to a file
            //if( flags.netcdf_variance )
            //{
            //   if (
            //         output_variance_to_netcdf(
            //            diffracted_wave_variance_sum,
            //            k_binning_boundaries,
            //            output_prefix + "_fem_Vbar_r"
            //         ) != EXIT_SUCCESS)
            //   {
            //      cout << " failed to write variance data to a"
            //        << " netCDF file : " << output_prefix << endl;
            //   }
            //} else {
            if (
               output_variance_to_txt(
                     diffracted_wave_variance_sum,
                     k_binning_boundaries,
                     output_prefix + "_fem_Vbar_r"
                  ) != EXIT_SUCCESS)
            {
               cout << " failed to write --Vbar_r data to a"
                  << " txt file : " << output_prefix << endl;
            }
            //}
         } // rootnode
      } // flags.d2 

      if (flags.gt17 || flags.d3 || flags.d4 
            || flags.rva )
      {
         ////////////////////////////////////////////////////////////////
         // Calculate the FTEM 2-D variance V(\vec{k},K_{aperture})
         ////////////////////////////////////////////////////////////////
         // V(\vec{k},K_{aperture})
         //  = \frac{ \langle I^{2}(\vec{k},K_{aperture}) \rangle }
         //          { \langle I(\vec{k},K_{aperture}) \rangle^{2} } -1
         // where I() and I^{2}() are diffracted_wave_mag_sum and 
         // diffracted_wave_mag_sqr_sum, respectively.

         // * save diffracted_wave_mag_sum and 
         //    diffracted_wave_mag_sqr_sum as images
         // * calculate the averages of diffracted_wave_mag_sum and
         //    diffracted_wave_mag_sqr_sum as a function of |\vec{k}|
         // * calculate the variance from the averages
         //    - determine a proper \Delta k for the variance curve 
         //       resolution
         //    * instantiate a variance domain using \Delta k
         //    * instantiate a variance[] array
         //    * for each k_{i} in the domain of variance, integrate 
         //       diffracted_wave_mag_sum and diffracted_wave_mag_sqr_sum
         //       between k_{i} and k_{i+1}
         

         // Divide the sum by the number of points to obtain the average
         //  2-D diffracted intensity.
         if ( flags.gt17 || flags.d3 || flags.d4
               || flags.rva )
         {
            for ( ptrdiff_t ii=0; ii < local_alloc_size_fftw; ++ii )
            {
               diffracted_wave_mag_sum[ii] 
                  = diffracted_wave_mag_sum[ii] 
                     / number_of_raster_points; 

               diffracted_wave_mag_sum_sqr[ii] 
                  = diffracted_wave_mag_sum[ii] 
                     * diffracted_wave_mag_sum[ii];

               if ( flags.gt17 || flags.d3 || flags.d4 )
               {
                  diffracted_wave_mag_sqr_sum[ii] 
                     = diffracted_wave_mag_sqr_sum[ii] 
                        / number_of_raster_points;
               }
               if (flags.complex_realspace_sum)
               {
                  diffracted_wave_re_sum[ii] 
                     = diffracted_wave_re_sum[ii] 
                        / number_of_raster_points;

                  diffracted_wave_im_sum[ii] 
                     = diffracted_wave_im_sum[ii] 
                        / number_of_raster_points;
               }
            }

            //// begin debug
            //if (diffracted_wave_mag_sum_max < diffracted_wave_mag_sum[ii])
            //   diffracted_wave_mag_sum_max =  diffracted_wave_mag_sum[ii];

            //if (
            //      diffracted_wave_mag_sqr_sum_max 
            //      < diffracted_wave_mag_sqr_sum[ii] 
            //   )
            //   diffracted_wave_mag_sqr_sum_max 
            //      =  diffracted_wave_mag_sqr_sum[ii];
            //// end debug

            // Combine re_sum and im_sum to obtain the magnitudes
            //  for 2-D variance images

            // TODO: modify to obtain 1-D variance curve using the complex sum
            if (flags.complex_realspace_sum) 
            {
               diffracted_wave_complex_sum_mag 
                  = new double[ local_alloc_size_fftw ];

               diffracted_wave_fftw_complex_sum 
                  = fftw_alloc_complex( local_alloc_size_fftw );

               for ( size_t ii=0; ii < local_alloc_size_fftw; ++ii )
               {
                  diffracted_wave_complex_sum_mag[ii] 
                     = sqrt(
                           (diffracted_wave_re_sum[ii] 
                              * diffracted_wave_re_sum[ii])
                           + (diffracted_wave_im_sum[ii] 
                             * diffracted_wave_im_sum[ii])
                       );

                  //cout << "diffracted_wave_im_sum[" << ii << "] : " // debug
                  //   << diffracted_wave_im_sum[ii] << endl; // debug
                  //cout << "diffracted_wave_complex_sum[" << ii << "] : " // debug
                  //   << diffracted_wave_complex_sum[ii] << endl; // debug

                  // Store the sums as fftw_complex variables for 
                  //    input into output_psi_reciprocalspace_to_netcdf
                  // TODO : check to see if this is necessary, or even useful
                  // TODO : do these need to be divided by 
                  //          number_of_raster_points? they already have been
                  diffracted_wave_fftw_complex_sum[ii][0] 
                     = diffracted_wave_re_sum[ii];
                  diffracted_wave_fftw_complex_sum[ii][1] 
                     = diffracted_wave_im_sum[ii];

                  // TODO: remove this averaging when performed by 
                  //      variance_2D_STEM. ... why would it be done there?
                  diffracted_wave_complex_sum_mag[ii]
                     = diffracted_wave_complex_sum_mag[ii] 
                        / number_of_raster_points;
               }
            } // flags.complex_realspace_sum
         }  // flags.gt17 || flags.d3 || flags.d4 
            //  || flags.rva



         ////////////////////////////////////////////////////////////////
         // Save average of 2-D diffraction magnitude taken from all 
         //  raster points
         ////////////////////////////////////////////////////////////////
         if ( flags.image_output )
         {
            //diffraction_scale_factor = 1.0e15 ;//1e-53;
            diffraction_scale_factor = 1.0e13 ;//1e-53;
            output_diffraction(
                  diffracted_wave_mag_sum,
                  diffraction_scale_factor,
                  local_alloc_size_fftw,
                  Nx_local, Nx, Ny,
                  resolutionUnit_recip,
                  xResolution_recip, yResolution_recip,
                  output_prefix + "_diffracted_wave_avg",
                  psi_mag_strides,
                  psi_mag_displacements,
                  mynode, rootnode, comm
                  );
            if ( flags.gt17 || flags.d3 || flags.d4) 
            {
               diffraction_scale_factor = 1.0e13 ;//1e-53;
               output_diffraction(
                     diffracted_wave_mag_sqr_sum,
                     diffraction_scale_factor,
                     local_alloc_size_fftw,
                     Nx_local, Nx, Ny,
                     resolutionUnit_recip,
                     xResolution_recip, yResolution_recip,
                     output_prefix + "_diffracted_wave_sqr_avg",
                     psi_mag_strides,
                     psi_mag_displacements,
                     mynode, rootnode, comm
                     );
            }
         }

         if (flags.image_output && flags.complex_realspace_sum) 
         {
            diffraction_scale_factor = 1.0e13 ;//1e-53;
            output_diffraction(
                  diffracted_wave_re_sum,
                  diffraction_scale_factor,
                  local_alloc_size_fftw,
                  Nx_local, Nx, Ny,
                  resolutionUnit_recip,
                  xResolution_recip, yResolution_recip,
                  output_prefix + "_diffracted_wave_re_sum",
                  psi_mag_strides,
                  psi_mag_displacements,
                  mynode, rootnode, comm
                  );

            output_diffraction(
                  diffracted_wave_im_sum,
                  diffraction_scale_factor,
                  local_alloc_size_fftw,
                  Nx_local, Nx, Ny,
                  resolutionUnit_recip,
                  xResolution_recip, yResolution_recip,
                  output_prefix + "_diffracted_wave_im_sum",
                  psi_mag_strides,
                  psi_mag_displacements,
                  mynode, rootnode, comm
                  );

            output_diffraction(
                  diffracted_wave_complex_sum_mag,
                  diffraction_scale_factor,
                  local_alloc_size_fftw,
                  Nx_local, Nx, Ny,
                  resolutionUnit_recip,
                  xResolution_recip, yResolution_recip,
                  output_prefix 
                     + "_diffracted_wave_complex_sum_mag_avg",
                  psi_mag_strides,
                  psi_mag_displacements,
                  mynode, rootnode, comm
                  );

         } // flags.image_output && flags.complex_realspace_sum

         //if (flags.complex_realspace_sum && flags.netcdf_images)
         //{
         //  if(
         //     output_psi_reciprocalspace_to_netcdf(
         //            diffracted_wave_fftw_complex_sum,
         //            local_alloc_size_fftw,
         //            Nx_local,
         //            kx_joined, Nx,
         //            ky, Ny,
         //            output_prefix 
         //               + "_diffracted_wave_fftw_complex_sum",
         //            psi_mag_strides,
         //            psi_mag_displacements,
         //            mynode, rootnode, comm
         //            ) == EXIT_FAILURE
         //  )
         //  {
         //     cout << "output_psi_reciprocalspace_to_netcdf() failed" 
         //        << endl;
         //      fftw_free( diffracted_wave_fftw_complex_sum ); // debug

         //   }
         //}

         if ( flags.gt17 || flags.d3 )
         {
            // Azimuthally average and MPI_Reduce the 2-D averages 
            //    of intensity and intensity squared.

            integrate_out_phi_double(
                  diffracted_wave_mag_sqr_sum,
                  kx_local, Nx_local,
                  ky, Ny,
                  indexed_magnitudes,
                  k_binning_boundaries,
                  k_bin_counts_local,
                  diffracted_wave_mag_sqr_sum_radial_intensity_local
                  );

            // zero the region of the direct beam
            for ( size_t ii=0; ii < number_of_k_bins; ++ii)
            {
               if ( lambda * k_binning_boundaries[ ii + 1] <= alpha_max)
               {
                  diffracted_wave_mag_sqr_sum_radial_intensity_local[ ii] = 0.0;
               }
            }

            MPI_Reduce(
                  k_bin_counts_local,
                  k_bin_counts_sqr_sum_aggregated,
                  number_of_k_bins,
                  MPI_INT, MPI_SUM,
                  rootnode, comm
                  );

            MPI_Reduce(
                  diffracted_wave_mag_sqr_sum_radial_intensity_local,
                  diffracted_wave_mag_sqr_sum_radial_intensity_total,
                  number_of_k_bins,
                  MPI_DOUBLE, MPI_SUM,
                  rootnode, comm
                  );
         }

         if ( flags.gt17 || flags.rva )
         {
            integrate_out_phi_double(
                  diffracted_wave_mag_sum_sqr,
                  kx_local, Nx_local,
                  ky, Ny,
                  indexed_magnitudes,
                  k_binning_boundaries,
                  k_bin_counts_local,
                  diffracted_wave_mag_sum_sqr_radial_intensity_local
                  );

            // zero the region of the direct beam
            for ( size_t ii=0; ii < number_of_k_bins; ++ii)
            {
               if ( lambda * k_binning_boundaries[ ii + 1] <= alpha_max)
               {
                  diffracted_wave_mag_sum_sqr_radial_intensity_local[ ii] = 0.0;
               }
            }

            MPI_Reduce(
                  k_bin_counts_local,
                  k_bin_counts_sum_sqr_aggregated,
                  number_of_k_bins,
                  MPI_INT, MPI_SUM,
                  rootnode, comm
                  );

            MPI_Reduce(
                  diffracted_wave_mag_sum_sqr_radial_intensity_local,
                  diffracted_wave_mag_sum_sqr_radial_intensity_total,
                  number_of_k_bins,
                  MPI_DOUBLE, MPI_SUM,
                  rootnode, comm
                  );

            if ( mynode == rootnode && flags.gt17 )
            {
               for( size_t ii=0; ii < number_of_k_bins; ++ii)
               {
                  // Caculate the ratio  to obtain 
                  //  \frac{\langle \langle I^{2}(\vec{k})\rangle_{\vec{r}} \rangle_{\phi}}{\langle \langle I(\vec{k})\rangle_{\vec{r}}^{2} \rangle_{\phi}}
                  ftem_gt17[ii] 
            = (   diffracted_wave_mag_sqr_sum_radial_intensity_total[ii]
                           / k_bin_counts_sqr_sum_aggregated[ii])
                           / 
              (   diffracted_wave_mag_sum_sqr_radial_intensity_total[ii]
                               / k_bin_counts_sum_sqr_aggregated[ii]);

                  if ( ! flags.rva )
                  {
                     k_bin_counts_aggregated[ii] = 0;
                  }
                  if (! flags.d3  )
                  {
                     k_bin_counts_sqr_sum_aggregated[ii] = 0;
                     // NOTE: k_bin_counts_sqr_sum_aggregated and
                     //       k_bin_counts_aggregated 
                     //       should be equal, and might not be 
                     //       needed here.
                  }
               }

               // save output of gt17
               //if( flags.netcdf_variance )
               //{
               //   if ( 
               //         output_variance_to_netcdf(
               //            ftem_gt17,
               //            k_binning_boundaries,
               //            output_prefix 
               //               + "_fem_V_gt17"
               //         ) != EXIT_SUCCESS)
               //   {
               //      cout << " failed to write V_gt17 variance data to"
               //        << " netCDF file : " 
               //        << output_prefix 
               //            + "_fem_V_gt17.nc" 
               //        << endl;
               //   }
               //} else {
                  // debug
                  //output_variance_to_txt(
                  //   diffracted_wave_radial_intensity_sqr_avg_gt17,
                  //   k_binning_boundaries,
                  //   output_prefix 
                  //      + "_radial_intensity_sqr_avg_gt17"
                  //   );
                  //output_variance_to_txt(
                  //   diffracted_wave_radial_intensity_avg_sqr_gt17,
                  //   k_binning_boundaries,
                  //   output_prefix 
                  //      + "_radial_intensity_avg_sqr_gt17"
                  //   );
                  // end debug
                  if ( 
                     output_variance_to_txt(
                           ftem_gt17,
                           k_binning_boundaries,
                           output_prefix 
                              + "_fem_V_gt17"
                        ) != EXIT_SUCCESS)
                  {
                     cout << " failed to write V_gt17 variance data"
                        << " to a txt file : " 
                        << output_prefix 
                           + "_fem_V_gt17.txt"
                        << endl;
                  } // gt17 text output block
               //} 
            } // mynode == rootnode && flags.gt17 
         } // flags.gt17 || flags.rva

         if ( flags.d3 || flags.rva )
         {
            integrate_out_phi_double(
                  diffracted_wave_mag_sum,
                  kx_local, Nx_local,
                  ky, Ny,
                  indexed_magnitudes,
                  k_binning_boundaries,
                  k_bin_counts_local,
                  diffracted_wave_mag_sum_radial_intensity_local
                  );

            // zero the region of the direct beam
            for ( size_t ii=0; ii < number_of_k_bins; ++ii)
            {
               if ( lambda * k_binning_boundaries[ ii + 1] <= alpha_max)
               {
                  diffracted_wave_mag_sum_radial_intensity_local[ ii] = 0.0;
               }
            }

            MPI_Reduce(
                  k_bin_counts_local,
                  k_bin_counts_aggregated,
                  number_of_k_bins,
                  MPI_INT, MPI_SUM,
                  rootnode, comm
                  );

            MPI_Reduce(
                  diffracted_wave_mag_sum_radial_intensity_local,
                  diffracted_wave_mag_sum_radial_intensity_total,
                  number_of_k_bins,
                  MPI_DOUBLE, MPI_SUM,
                  rootnode, comm
                  );

            if ( mynode == rootnode )
            {
               for (size_t ii=0; ii < number_of_k_bins; ++ii)
               {
                  if ( flags.d3 ) 
                  {
                     ftem_d3[ii] =
                        ((
              diffracted_wave_mag_sqr_sum_radial_intensity_total[ii]
                     / k_bin_counts_sqr_sum_aggregated[ii]
                        ) / (
              diffracted_wave_mag_sum_radial_intensity_total[ii]
              * diffracted_wave_mag_sum_radial_intensity_total[ii]
              / (
                 k_bin_counts_aggregated[ii] * k_bin_counts_aggregated[ii]
                )
                        )) - 1.0;
                  }
                  if ( flags.rva ) 
                  {
                     ftem_rva[ii] 
                = ((diffracted_wave_mag_sum_sqr_radial_intensity_total[ii]
                     / k_bin_counts_sum_sqr_aggregated[ii])  
                  / (
                        diffracted_wave_mag_sum_radial_intensity_total[ii]
                        *
                        diffracted_wave_mag_sum_radial_intensity_total[ii]
                        / (k_bin_counts_aggregated[ii] 
                           * k_bin_counts_aggregated[ii])
                    )) - 1.0;
                  }
               }
               // save output of d3
               if ( flags.d3 ) 
               {
                  //if( flags.netcdf_variance )
                  //{
                  //   if ( 
                  //         output_variance_to_netcdf(
                  //            ftem_d3,
                  //            k_binning_boundaries,
                  //            output_prefix 
                  //               + "_fem_V_re"
                  //         ) != EXIT_SUCCESS)
                  //   {
                  //      cout << " failed to write --V_re variance data to"
                  //        << " netCDF file : " 
                  //        << output_prefix 
                  //            + "_fem_V_re.nc" 
                  //        << endl;
                  //   }
                  //} else {
                     if ( 
                        output_variance_to_txt(
                              ftem_d3,
                              k_binning_boundaries,
                              output_prefix 
                                 + "_fem_V_re"
                           ) != EXIT_SUCCESS)
                     {
                        cout << " failed to write --V_re variance data"
                           << " to a txt file : " 
                           << output_prefix 
                              + "_fem_V_re.txt"
                           << endl;
                     }
                  //} // d3 text output block
               } // flags.d3

               if ( flags.rva ) 
               {
                  //if( flags.netcdf_variance )
                  //{
                  //   if ( 
                  //         output_variance_to_netcdf(
                  //            ftem_rva,
                  //            k_binning_boundaries,
                  //            output_prefix 
                  //               + "_fem_V_rbar"
                  //         ) != EXIT_SUCCESS)
                  //   {
                  //      cout << " failed to write --V_rbar variance data"
                  //        << " to netCDF file : " 
                  //        << output_prefix 
                  //            + "_fem_V_rbar.nc" 
                  //        << endl;
                  //   }
                  //} else {
                     if ( 
                        output_variance_to_txt(
                              ftem_rva,
                              k_binning_boundaries,
                              output_prefix 
                                 + "_fem_V_rbar"
                           ) != EXIT_SUCCESS)
                     {
                        cout << " failed to write --V_rbar variance data"
                           << " to a txt file : " 
                           << output_prefix 
                              + "_fem_V_rbar.txt"
                           << endl;
                     }
                  //} // rva text output block
               } // flags.rva
            } // root node 
         } // flags.d3 || flags.rva

         if ( flags.d4 )
         {
            for ( size_t ii=0; ii < local_alloc_size_fftw; ++ii)
            {
               diffracted_wave_mag_variance[ii] = 
                  (
                     diffracted_wave_mag_sqr_sum[ii] 
                     /
                     (
                        diffracted_wave_mag_sum[ii]
                        * diffracted_wave_mag_sum[ii]
                     )
                  ) - 1.0;
            }
            // zero the direct-beam region of the 2-D variance image
            double ksqr;
            for ( size_t i=0; i < Nx_local; ++i)
            {
               for ( size_t j=0; j < Ny; ++j)
               {
                  ksqr = pow(kx_local[i], 2) + pow(ky[j], 2);
                  if ( lambda_sqr * ksqr <= alpha_max_sqr )
                  {
                     diffracted_wave_mag_variance[ j + Ny*i] = 0.0;
                  }
               }
            }
            if ( flags.image_output )
            {
               diffraction_scale_factor = 1.0e0;
               output_diffraction(
                     diffracted_wave_mag_variance,
                     diffraction_scale_factor,
                     local_alloc_size_fftw,
                     Nx_local, Nx, Ny,
                     resolutionUnit_recip,
                     xResolution_recip, yResolution_recip,
                     output_prefix + "_variance_2-D",
                     psi_mag_strides,
                     psi_mag_displacements,
                     mynode, rootnode, comm
                     );
               diffraction_scale_factor = 1.0e-10;
               output_diffraction(
                     diffracted_wave_mag_variance,
                     diffraction_scale_factor,
                     local_alloc_size_fftw,
                     Nx_local, Nx, Ny,
                     resolutionUnit_recip,
                     xResolution_recip, yResolution_recip,
                     output_prefix + "_variance_2-D",
                     psi_mag_strides,
                     psi_mag_displacements,
                     mynode, rootnode, comm
                     );
               //diffraction_scale_factor = 1.0e-25;
               //output_diffraction(
               //      psi,
               //      diffraction_scale_factor,
               //      local_alloc_size_fftw,
               //      Nx_local, Nx, Ny,
               //      resolutionUnit_recip,
               //      xResolution_recip, yResolution_recip,
               //      output_prefix + "_variance_2-D",
               //      mynode, rootnode, comm
               //      );
            }

            // azimuthally average 2-D variance to obtain d4
            integrate_out_phi_double(
                  diffracted_wave_mag_variance,
                  kx_local, Nx_local,
                  ky, Ny,
                  indexed_magnitudes,
                  k_binning_boundaries,
                  k_bin_counts_local,
                  diffracted_wave_mag_variance_radial_intensity_local
                  );

            // zero the region of the direct beam
            for ( size_t ii=0; ii < number_of_k_bins; ++ii)
            {
               if ( lambda * k_binning_boundaries[ ii + 1] <= alpha_max)
               {
                  diffracted_wave_mag_variance_radial_intensity_local[ ii] = 0.0;
               }
            }

            MPI_Reduce(
                  diffracted_wave_mag_variance_radial_intensity_local,
                  ftem_d4,
                  number_of_k_bins,
                  MPI_DOUBLE, MPI_SUM,
                  rootnode, comm
                  );

            MPI_Reduce(
                  k_bin_counts_local,
                  k_bin_counts_aggregated,
                  number_of_k_bins,
                  MPI_INT, MPI_SUM,
                  rootnode, comm
                  );

            if ( mynode == rootnode )
            {
               for (size_t ii=0; ii < number_of_k_bins; ++ii)
               {
                  ftem_d4[ii] = ftem_d4[ii] / k_bin_counts_aggregated[ii];
               }
               // save output of d4
               //if( flags.netcdf_variance )
               //{
               //   if ( 
               //         output_variance_to_netcdf(
               //            ftem_d4,
               //            k_binning_boundaries,
               //            output_prefix 
               //               + "_fem_Omega_vimage"
               //         ) != EXIT_SUCCESS)
               //   {
               //      cout << " failed to write --Omega_vimage variance data to"
               //        << " netCDF file : " 
               //        << output_prefix 
               //            + "_fem_Omega_vimage.nc" 
               //        << endl;
               //   }
               //} else {
                  if ( 
                     output_variance_to_txt(
                           ftem_d4,
                           k_binning_boundaries,
                           output_prefix 
                              + "_fem_Omega_vimage"
                        ) != EXIT_SUCCESS)
                  {
                     cout << " failed to write --Omega_vimage variance data"
                        << " to a txt file : " 
                        << output_prefix 
                           + "_fem_Omega_vimage.txt"
                        << endl;
                  } // d4 text output block
               //}
            }
         }

      } // flags.gt17 || flags.d3 || flags.d4 
         // flags.rva

      if ( flags.correlograph && (mynode == rootnode) )
      {
         size_t idx;
         size_t c_idx;
         size_t d_idx;
         for ( size_t ii=0; ii < number_of_k_bins; ++ii)
         {
            for ( size_t jj=0; jj < number_of_phi_bins; ++jj)
            {
               //idx = jj + ii * number_of_phi_bins;
               //c_idx = jj + ii * number_of_phi_bins;
               c_idx = ii + jj * number_of_k_bins;
               correlograph_sum[c_idx]
                  = correlograph_sum[c_idx] / number_of_raster_points;
            }
         }
         output_correlograph_to_txt(
               correlograph_sum,
               phi_binning_boundaries,
               k_binning_boundaries,
               output_prefix 
                  + "_correlograph_avg"
               );
         output_correlograph_image(
               correlograph_sum,
               number_of_phi_bins,
               number_of_k_bins,
               resolutionUnit,
               phiResolution_correlograph,
               kResolution_correlograph,
               output_prefix + "_correlograph_avg",
               0
               );
         // iterate through k and phi bins to determine dimensions of 
         //  the variance image having k \in [0.2, 1.0]
         if ( flags.correlograph_variance )
         {
            if ( (mynode == rootnode) && flags.debug)
               cout << "calculating variance of correlographs ..." 
                  << endl;
            size_t k_correlograph_start_idx, k_correlograph_end_idx;
            k_correlograph_start_idx=0;
            k_correlograph_end_idx = number_of_k_bins - 1;
            for ( size_t ii=0; ii < number_of_k_bins; ++ii)
            {
               //if ( flags.debug )
               //   cout << "k_binning_boundaries[" << ii << "]: "
               //      << k_binning_boundaries[ii] 
               //      << " k_binning_boundaries[" << ii + 1 << "]: "
               //      << k_binning_boundaries[ii + 1] << endl;

               if ( (k_binning_boundaries[ii] <= 0.2) 
                     && (k_binning_boundaries[ii+1] > 0.2))
                  k_correlograph_start_idx = ii;

               if ( (k_binning_boundaries[ii] <= 1.0) 
                     && (k_binning_boundaries[ii+1] > 1.0))
                  k_correlograph_end_idx = ii;
            }
            double* correlograph_variance;
            double tmpdbl;
            size_t number_of_reducedk_bins;
            if ( k_correlograph_end_idx >= k_correlograph_start_idx )
               number_of_reducedk_bins
                  = (k_correlograph_end_idx 
                     - k_correlograph_start_idx );
            else 
               number_of_reducedk_bins = 0;

            correlograph_variance 
               = new double[number_of_reducedk_bins * number_of_phi_bins];
            //for ( size_t ii = 0; ii < number_of_k_bins; ++ii)
            for ( size_t ii = k_correlograph_start_idx; 
                  ii < k_correlograph_end_idx ; 
                  ++ii)
            {
               for ( size_t jj=0; jj < number_of_phi_bins; ++jj)
               {
                  //idx = jj + ii * number_of_phi_bins;
                  c_idx = ii + jj * number_of_k_bins;
                  //d_idx = jj + (ii - k_correlograph_start_idx) 
                  //               * number_of_phi_bins;
                  d_idx = (ii - k_correlograph_start_idx)
                           + jj * number_of_reducedk_bins;
                  //correlograph_variance[idx]
                  tmpdbl = 
                         correlograph_sqr_sum[c_idx]
                           / (
                              number_of_raster_points 
                           * correlograph_sum[c_idx] 
                           * correlograph_sum[c_idx]);
                  if ( (correlograph_sqr_sum[c_idx] != 0)
                        && (number_of_raster_points != 0)
                        && (correlograph_sum[c_idx] != 0)
                        && (tmpdbl >= 1))
                  {
                     correlograph_variance[ d_idx ] = tmpdbl - 1;
                     //   jj + c_idx * number_of_phi_bins]
                     //correlograph_sqr_sum[idx]
                        //(
                        // correlograph_sqr_sum[c_idx]
                        //   / (
                        //      number_of_raster_points 
                        //   * correlograph_sum[idx] * correlograph_sum[idx]
                        //   )
                        //) - 1;
                  }
                  else
                     correlograph_variance[ d_idx ] = 0;
                     //correlograph_sqr_sum[idx] = 0;
               }
            }
            output_correlograph_to_txt(
                  correlograph_sqr_sum,   // correlograph_variance
                  phi_binning_boundaries,
                  k_binning_boundaries,
                  output_prefix 
                     + "_correlograph_variance"
                  );
            output_correlograph_image(
                  //correlograph_sqr_sum,   // correlograph_variance
                  correlograph_variance,   // correlograph_variance
                  //number_of_k_bins,
                  number_of_phi_bins,
                  number_of_reducedk_bins,
                  resolutionUnit,
                  phiResolution_correlograph,
                  kResolution_correlograph,
                  output_prefix + "_correlograph_variance",
                  1  // flag to indicate logscale
                  );
            output_correlograph_image(
                  correlograph_variance,   // correlograph_variance
                  number_of_phi_bins,
                  number_of_reducedk_bins,
                  resolutionUnit,
                  phiResolution_correlograph,
                  kResolution_correlograph,
                  output_prefix + "_correlograph_variance",
                  0  // flag to indicate logscale
                  );
            delete[] correlograph_variance;
         }
      }
   }// end fem specific commands

   ////////////////////////////////////////////////////////////////
   // Save stem_image as a tif
   ////////////////////////////////////////////////////////////////
   if ( mynode == rootnode && flags.image_output )
   {
      output_stem_image(
            stem_image,
            x_p.size(), y_p.size(),
            resolutionUnit,
            xResolution_stem, yResolution_stem,
            output_prefix 
            //mynode, rootnode, comm
            );
   }

   ////////////////////////////////////////////////////////////////
   // Save diffraction as tif
   ////////////////////////////////////////////////////////////////
   
   if ( flags.image_output )
   {
      //output_diffraction_with_renormalization(
      //diffraction_scale_factor = 1.0e13;//1e-20;
      //diffraction_scale_factor = 1.0e+20;//1e-20;
      diffraction_scale_factor = 1.0e+0;
      output_diffraction(
            psi_accumulation,
            diffraction_scale_factor,
            local_alloc_size_fftw,
            Nx_local, Nx, Ny,
            resolutionUnit_recip,
            xResolution_recip, yResolution_recip,
            output_prefix + "_final",
            psi_mag_strides,
            psi_mag_displacements,
            mynode, rootnode, comm
            );

      //debug_output_complex_fftw_operand(
      //      psi,
      //      5, 
      //      local_alloc_size_fftw,
      //      Nx_local, Nx, Ny,
      //      resolutionUnit,
      //      xResolution, yResolution,
      //      output_prefix + "_diffraction_" 
      //         + TEM_NS::to_string(sliceNumber),
      //      psi_mag_strides,
      //      psi_mag_displacements,
      //      mynode, rootnode, comm
      //      );
   }


   //////////////////////////////////////////////////////////////////
   // Clean up
   //////////////////////////////////////////////////////////////////

   delete lmp;

   if ( mynode == rootnode && flags.debug )
            cout << "cleaning up memory ..." << endl; // debug

   //delete_slices_from_list( slicesOfScatterers );
     // ^^ calls fftw_free() on paps, exp_i_sigma_v, and deletes 
     //  propagator_x_re, propagator_x_im, propagator_y_re,
     //  and propagator_y_im, and also deletes the vector<scatterer*>*
     //  member of each slice

   delete[] q;
   delete[] Z_array;
   if ( uniqueZs != NULL) delete[] uniqueZs;
   if ( flags.dupe && flags.fem 
        &&
        (dupe_x > 1 || dupe_y > 1 || dupe_z > 1)
        &&
        (dupe_x > 0 && dupe_y > 0 && dupe_z > 0)
      )
   {
      delete[] q_duped;
      //delete[] Z_array_duped;
   }

   if ( flags.mtf_file && flags.mtf_resolution )
   {
     //   mtf_1D;
     //   mtf_domain;
     //   // TODO: does the mtf need to be split among nodes into local sub-intervals? Not the 1-D MTF, but the 2-D MTF does.
      fftw_destroy_plan( pf_c2c_mtf );
      fftw_destroy_plan( pb_c2c_mtf );

      fftw_free( mtf_2D_split );
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

   delete[] Nx_strides;
   delete[] Nx_displacements;
   delete[] psi_mag_strides;
   delete[] psi_mag_displacements;

   // Save wisdom to a file if it wasn't successfully opened above
   if ( ! flags.wisdomFile )
   {
      fftw_mpi_gather_wisdom( MPI_COMM_WORLD); // gathers to rank 0 proc
      if ( mynode == rootnode)
      {
         if ( flags.debug)
            cout << "Exporting wisdom to file: "
               << wisdom_file_name << endl;

         if ( ! fftw_export_wisdom_to_filename(
                wisdom_file_name.c_str()))
         {
            cerr << "Error - could not export wisdom to file: "
              << wisdom_file_name << endl;
         }
      }

      delete[] stem_image;
   }

   fftw_destroy_plan( pb_c2c_probe_split );
   fftw_destroy_plan( pf_c2c_psi );
   fftw_destroy_plan( pb_c2c_psi );
   //fftw_destroy_plan( pf_c2c_t );
   //fftw_destroy_plan( pb_c2c_t );

   fftw_free( init_probe );
   fftw_free( psi );

   delete[] probe_split_re;
   delete[] probe_split_im;

   delete[] probe_joined_re;
   delete[] probe_joined_im;

   if ( flags.fem )
   {
      indexed_magnitudes.clear();
      if ( flags.d1 || flags.d2  || flags.d3 
            || flags.d4 || flags.gt17 || flags.rva )
      {
         delete[] k_bin_counts_local;
         //delete[] k_bin_counts_sqr_local;
         if ( mynode == rootnode )
         {
            delete[] k_bin_counts_aggregated;
            delete[] k_bin_counts_sqr_sum_aggregated;
            delete[] k_bin_counts_sum_sqr_aggregated;
         }
      }
      if ( flags.d1 || flags.d2 || flags.correlograph )
      {
         delete[] diffracted_wave_mag_radial_intensity_local;
         if ( flags.correlograph )
         {
            delete[] diffracted_wave_mag_in_radial_coords_local;
            delete[] correlograph;
            delete[] correlograph_sum;
            if ( flags.correlograph_variance )
               delete[] correlograph_sqr_sum;
         }
         if ( flags.d2 )
         {
            delete[] diffracted_wave_mag_sqr_radial_intensity_local;
         }
         if ( mynode == rootnode )
         {
            delete[] diffracted_wave_mag_radial_intensity_total;
            if ( flags.correlograph )
            {
               delete[] diffracted_wave_mag_in_radial_coords_sum;
               delete[] diffracted_wave_mag_in_radial_coords_total;
            }
            if ( flags.d1 )
            {
               delete[] ftem_d1;
               delete[] diffracted_wave_mag_radial_intensity_total_sum;
               delete[] 
                  diffracted_wave_mag_radial_intensity_total_sqr_sum;
            }
            if ( flags.d2 )
            {
               delete[] diffracted_wave_mag_sqr;
               delete[] diffracted_wave_mag_sqr_radial_intensity_total;
               delete[] diffracted_wave_variance_sum;
            }
         }
      }
      if ( flags.gt17 || flags.d3 || flags.d4
            || flags.rva ) 
      {
         delete[] diffracted_wave_mag_sum;
         delete[] diffracted_wave_mag_sum_sqr;
         if ( flags.gt17 || flags.d3 || flags.d4 ) 
         {
            delete[] diffracted_wave_mag_sqr_sum;
         }
         if (flags.complex_realspace_sum) 
         {
            delete[] diffracted_wave_re_sum;
            delete[] diffracted_wave_im_sum;
            //delete[] diffracted_wave_complex_variance_2D;
         }
         if ( flags.gt17 || flags.d3 || flags.rva ) 
         {
            if ( flags.gt17 || flags.d3 ) 
            {
               delete[] 
                  diffracted_wave_mag_sqr_sum_radial_intensity_local;
               if ( mynode == rootnode )
               {
                  delete[] 
                     diffracted_wave_mag_sqr_sum_radial_intensity_total;
               }
               if ( flags.gt17 )
               {
                  delete[] 
                     diffracted_wave_mag_sum_sqr_radial_intensity_local;
                  if ( mynode == rootnode )
                  {
                     delete[] ftem_gt17;
                     delete[] 
                       diffracted_wave_mag_sum_sqr_radial_intensity_total;
                  }
               }
               if ( flags.d3 || flags.rva )
               {
                  delete[] diffracted_wave_mag_sum_radial_intensity_local;
                  if ( mynode == rootnode )
                  {
                     if ( flags.d3 )
                     {
                        delete[] ftem_d3;
                     }
                     if ( flags.rva )
                     {
                        delete[] ftem_rva;
                     }
                     delete[] 
                        diffracted_wave_mag_sum_radial_intensity_total;
                  }
               }
            }
         }
         if ( flags.d4 )
         {
            delete[] diffracted_wave_mag_variance;
            delete[] diffracted_wave_mag_variance_radial_intensity_local;
            if ( mynode == rootnode )
            {
               delete[] ftem_d4;
               delete[] 
                  diffracted_wave_mag_variance_radial_intensity_total;
            }
         }
      }
   }

   fftw_mpi_cleanup();
   MPI_Barrier(MPI_COMM_WORLD);
   MPI_Finalize();

   return EXIT_SUCCESS;
}


