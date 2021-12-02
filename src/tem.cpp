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
#include "to_string.hpp"

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
   double alpha_max_sqr; // = pow( strtod( argv[9], NULL), 2);
   double alpha_max;
   string output_prefix; // = argv[2];
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

   std::vector< double > mtf;
   std::vector< double > mtf_domain;
   double mtf_resolution;

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
         mtf,
         mtf_domain,
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
      if ( flags.debug ) 
         cerr << "Failed to read command line arguments" << endl;
      //failflag = 1 ;
   }
   //MPI_Bcast( &flags.fail, 1, MPI_UNSIGNED, rootnode, MPI_COMM_WORLD );
   //if ( flags.fail) 
   //{
   //   fftw_mpi_cleanup();
   //   MPI_Finalize();
   //   if ( mynode == rootnode)
   //   {
   //      cout << "Exiting ms-stem-fem due to failure." << endl;
   //      PRINT_USAGE // macro defined above
   //   }
   //   return EXIT_FAILURE;
   //}


   // debug
   //if ( flags.debug )
   //{
   //   cout << "node " << mynode << "MTF size : " << mtf.size() << endl;
   //   if (mynode == rootnode)
   //   {
   //      cout << "MTF file contents:" << endl;
   //   }
   //   for ( std::vector<double>::iterator 
   //            itr = mtf.begin();
   //            itr != mtf.end();
   //            ++itr)
   //   {
   //      cout << "node " << mynode << ", mtf : " << *itr << endl;
   //   }
   //   for ( std::vector<double>::iterator 
   //            itr = mtf_domain.begin();
   //            itr != mtf_domain.end();
   //            ++itr)
   //   {
   //      cout << "node " << mynode << ", mtf_domain : " << *itr << endl;
   //   }
   //}
   // debug

   if ( flags.mtf_file && flags.mtf_resolution )
   {
      size_t mtf_domain_size = mtf_domain.size();
      for ( size_t mtf_domain_idx = 0; 
            mtf_domain_idx < mtf_domain_size;
            ++mtf_domain_idx)
      { // mtf_file was in units of Nyquist freq. = 0.5 * sampling freq.
         mtf_domain[mtf_domain_idx]
            = pow((0.5 * mtf_resolution * mtf_domain[mtf_domain_idx]), 2);
      }
   }

   // Read model file
   unsigned int read_xyz_flag = 0;
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
         flags.fail = 1 ;
      }
   }
   else
   {
      if ( mynode == rootnode && flags.debug )
      {
         cout << "Reading lammps data format file: "
           << model_file_name  << endl; // debug
      }
      if (
            read_position_lammps_data_file(
               model_file_name, // only valid on root node
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
         flags.fail = 1 ;
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
      fftw_mpi_cleanup();
      MPI_Finalize();
      if ( mynode == rootnode)
      {
         cout << "Exiting ms-stem-fem due to failure." << endl;
         PRINT_USAGE // macro defined above
      }
      return EXIT_FAILURE;
   }

   NxNy = Nx * Ny;
   NxNy_sqr = NxNy * NxNy;
   if ( Nx == Ny ) sqrtNxNy =  Nx ;   
   else sqrtNxNy = sqrt( NxNy );


   // At this point rootnode should have arrays of all scatterer 
   //  positions and species Z_array.
   
   // Each node has a scatterer_pap_LUT, and should instantiate its 
   //  own vector of scatterers upon recieving q[][3] and Z_array[] 
   //  from rootnode.

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


   //////////////////////////////////////////////////////////////////
   // Instantiate local variables for fftw_mpi
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

   //if ( flags.debug )
   //{
   //   cout << "node " << mynode 
   //      << ", (Nx_local, local_0_start_fftw, local_alloc_size_fftw) : (" 
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
   }

   MPI_Bcast( &(flags.wisdomFile), 1, MPI_INT, rootnode, MPI_COMM_WORLD);

   if ( flags.wisdomFile == 0)
   {
      if (( mynode == rootnode) && flags.debug )
         cout << "Could not open wisdom file: " << wisdom_file_name
            << " ," << endl
            << " will create it after fftw plan ..." << endl;
   }
   else
   {
      // send wisdom from file to all processors
      fftw_mpi_broadcast_wisdom( MPI_COMM_WORLD);
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

   MPI_Bcast(duped_periods, 3, MPI_DOUBLE, rootnode, MPI_COMM_WORLD);

   if ( mynode != rootnode )
   {
      xperiod_duped = duped_periods[0];
      yperiod_duped = duped_periods[1];
      zperiod_duped = duped_periods[2];
      //cout << " node " << mynode 
      //   << " (duped_periods[2], zperiod_duped) : (" 
      //   << duped_periods[2] << ", " << zperiod_duped 
      //   << ")" <<  endl;
   }

   if ( flags.debug )
      cout << "Checking MPI_Bcast() quality:"
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

   double bwcutoff_t; // transmission function bandwidth cutoff
   bwcutoff_t = (2.0/3.0) * bwcutoff_pap;

   double bwcutoff_t_sqr; 
   bwcutoff_t_sqr = bwcutoff_t * bwcutoff_t;

   if ( mynode == rootnode && flags.debug )
   {
      cout << "Bandwidth limits " << endl // debug
         << "  projected atomic potential bandwidth : " // debug
         << bwcutoff_pap << endl // debug
         << "  transmission function bandwidth : " // debug
         << bwcutoff_t << endl; // debug
   }
   
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
              << " bandwidth limit: " << bwcutoff_t * lambda << endl;
         }
         flags.fail = 1;
      }
   }

   if ( flags.fail == 1)
   {
      // clean up all allocated memory and exit
      delete[] Nx_strides;
      delete[] Nx_displacements;
      delete[] psi_mag_strides;
      delete[] psi_mag_displacements;
      
      if ( mynode == rootnode ) 
         cerr << " Exiting in error after evaluating input" << endl;

      MPI_Finalize();
      return EXIT_FAILURE;
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

      if ( flags.debug )
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
   if ( flags.debug && (mynode == rootnode))
   {
      for ( size_t ii=0; ii< totalnodes; ++ii)
      {
         cout << "Nx_strides[" << ii << "], Nx_displacements[ " 
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
   //MPI_Scatter( xx_joined, Nx_local, MPI_DOUBLE,
   //            xx_local, Nx_local, MPI_DOUBLE,
   //            rootnode, MPI_COMM_WORLD);

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
   //MPI_Scatter( kx_joined, Nx_local, MPI_DOUBLE,
   //            kx_local, Nx_local, MPI_DOUBLE,
   //            rootnode, MPI_COMM_WORLD);

   // TODO: consider bundling the follow data into a single Bcast
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
   // NOTE: this is only needed when the model is an xyz format file,
   //  otherwise it duplicates uniqueZs and numberOfSpecies.
   
   std::vector<unsigned int> Z_vector;
   unsigned int Z_vector_size;

   // TODO: parallelize this process of eliminating duplicate Zs
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
   if ( flags.debug )
      cout << "node " << mynode // debug
         << ", Z_vector_size : " << Z_vector_size << endl; // debug

   if ( mynode != rootnode ) Z_vector.resize( Z_vector_size );
   
   //MPI_Bcast( &Z_vector[0], Z_vector.size(),
   MPI_Bcast( &Z_vector[0], Z_vector_size,
               MPI_UNSIGNED, rootnode, MPI_COMM_WORLD);

   if ( flags.debug )
      cout << "node " << mynode  // debug
         << ", Z_vector.size (resized) : " 
         << Z_vector.size() << endl;// debug
   // TODO: The above is probably too reliant on communication, 
   //       the list of Zs is probably small enough that communication is
   //       slower than having duplicate calls to sort on each node.

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
         local_alloc_size_fftw,
         psi_mag_strides,
         psi_mag_displacements,
         mynode, rootnode, MPI_COMM_WORLD
         );

   if ( flags.debug )
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
   if ( flags.debug )
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
   if ( flags.dupe && flags.fem 
        &&
        (dupe_x > 1 || dupe_y > 1 || dupe_z > 1)
        &&
        (dupe_x > 0 && dupe_y > 0 && dupe_z > 0)
      )
   {
      if ( mynode == rootnode && flags.debug ) // debug
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
         //if ( (mynode == rootnode) && flags.debug)
         //   cout << "node " << mynode 
         //      << ", pushing scatterers onto myScatterers : " // debug
         //      << "Z_vector[" << i << "] : " << Z_vector[i] // debug
         //      << ", q[0,1,2] : " // debug
         //      << q[3*i]  // debug
         //      << ", " << q[3*i + 1] // debug
         //      << ", " << q[3*i + 2] // debug
         //      << endl;  // debug

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
   if ( uniqueZs != NULL) delete[] uniqueZs;

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

   std::list< slice* > slicesOfScatterers;

   // TODO: parallelize assign_atoms_to_slices_z_auto() or compute on
   // a single node and broadcast the results.
   if ( (mynode == rootnode) && flags.debug)
      cout << "assigning atoms to slices ..." << endl;
   assign_atoms_to_slices_z_auto(
         myScatterers,
         xmin,ymin,zmin,
         xperiod_duped, yperiod_duped, zperiod_duped,
         minSliceThickness,  // minimum slice thickness [\AA]
         //1.0,  // minimum slice thickness [\AA]
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
   //       <re-evaluate relevant slices, else reuse previous slices>
   
   // NOTE: Scatterers contained in myScatterers are pointed to from
   //       within members of slicesOfScatterers, so DO NOT release
   //       memory of the scatterers pointed to by myScatterers 
   //       members until all use of slicesOfScatterers is concluded.

   //////////////////////////////////////////////////////////////////
   // Initialize slice members
   //////////////////////////////////////////////////////////////////

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
   // debug
   if ( mynode == rootnode && flags.debug )
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
   // unsigned int flags.adfstem_corrected = 0;        // 
   // unsigned int flags.adfstem_uncorrected = 0;      // 
   // unsigned int flags.bfctem_corrected = 0;         // 
   // unsigned int flags.bfctem_uncorrected = 0;       // 
   // unsigned int flags.fem = 0;
   
   // run adfstem()
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
      std::list<double> x_p;
      std::list<double> y_p;
      if (flags.adfstem_uncorrected || flags.adfstem_corrected)
      {

         // NOTE: do not increase number of raster points when using
         //        --dupe
         size_t number_of_raster_points_x
            = floor( xperiod/raster_spacing );

         size_t number_of_raster_points_y
            = floor( yperiod/raster_spacing );

         if ( mynode == rootnode && flags.debug )
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
         if ( mynode == rootnode && flags.debug )
            cout << "Determined STEM raster points." << endl
               << " (x_p.front(), y_p.front()) : ("
               << x_p.front() << ", " << y_p.front() << ")" << endl;
      }

      //   // TODO: insert units [radians] [angstrom]
      //   cout << "Running adfstem with parameters : " << endl;//debug
      //   cout << "   fem : " << flags.fem << endl
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
      //      << output_prefix + "_adfstem_uncorr" << endl;

      //   adfstem(
      //         flags.fem,
      //         0,                         // aberration correction?
      //         0,
      //         flags.image_output,
      //         flags.debug,
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
      //         output_prefix + "_adfstem_uncorr",
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

      if ( mynode == rootnode && flags.debug )
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


      // NOTE: flags.adfstem_corrected and 
      //       flags.adfstem_uncorrected are not necessarily 
      //       opposites.
      //       We should allow for both corrected and uncorrected 
      //       simulations within the same run.


      // TODO: modify the defocus by the sample thickness to ensure
      //       that sample defocus is w.r.t. the bottom of the sample.
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

         //alpha_max_sqr = alpha_max * alpha_max;

         // TODO: insert units [radians] [angstrom]
         if ( mynode == rootnode && flags.debug )
         {
            cout << "Running adfstem with parameters : " << endl;//debug
            cout << "   fem : " << flags.fem << endl
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
               << output_prefix + "_adfstem_uncorr" << endl;
         }

         flags.aberration_correction = 0;
         flags.complex_realspace_sum = 0;
         adfstem(
               flags,
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
               mtf, // modulation transfer function
               mtf_domain, // has been squared (=mtf_domain_sqr)
               mtf_resolution,
               lambda, 
               output_prefix + "_adfstem_uncorr",
               local_alloc_size_fftw,
               local_0_start_fftw,
               Nx_strides,
               Nx_displacements,
               psi_mag_strides,
               psi_mag_displacements,
               mynode, rootnode, totalnodes,
               MPI_COMM_WORLD
             );
         
         if ( mynode == rootnode && flags.debug )
            cout << "root node finished running adfstem"
               << " without aberration correction" << endl; // debug
      }

      if ( flags.adfstem_corrected )
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

         //alpha_max_sqr = alpha_max * alpha_max;

         // TODO: insert units [radians] [angstrom]
         if ( mynode == rootnode && flags.debug )
         {
            cout << "Running adfstem with parameters : " << endl;//debug
            cout << "   fem : " << flags.fem << endl
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
               << output_prefix + "_adfstem_corr" << endl;
         }

         flags.aberration_correction = 1;
         flags.complex_realspace_sum = 0;
         adfstem(
               flags,
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
               mtf, // modulation transfer function
               mtf_domain, // mtf_domain_sqr
               mtf_resolution,
               lambda, 
               output_prefix + "_adfstem_corr",
               local_alloc_size_fftw,
               local_0_start_fftw,
               Nx_strides,
               Nx_displacements,
               psi_mag_strides,
               psi_mag_displacements,
               mynode, rootnode, totalnodes,
               MPI_COMM_WORLD
             );
         
         if ( mynode == rootnode && flags.debug )
            cout << "root node finished running adfstem"
               << " with aberration correction" << endl; // debug
      }

   }
   if ( 
         flags.bfctem_corrected 
         || flags.bfctem_uncorrected 
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
   //      output_prefix,              
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
   //      output_prefix,              
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
   //if ( flags.mtf_file )
   //{
   ////   mtf;
   ////   mtf_domain;
   ////   // TODO: does the mtf need to be split among nodes into local sub-intervals? Not the 1-D MTF, but the 2-D MTF does.
   //}

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
   }

   fftw_mpi_cleanup();
   MPI_Finalize();

   return EXIT_SUCCESS;
}


