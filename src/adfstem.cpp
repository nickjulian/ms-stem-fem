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
// File: adfstem.cpp
// Purpose:
//
// Simulate the conventional TEM of a thick sample as viewed perpendicular
//  to the third axis (q[2]).
// This file implements the multislice approximation, phase grating 
//  approximation, and maybe even the weak phase approximation. 

#ifndef ADFSTEM_CPP
#define ADFSTEM_CPP

#include "adfstem.hpp"



int TEM_NS::adfstem(
      const unsigned int& input_flag_fem,
      const unsigned int& input_flag_aberration_correction,
      const unsigned int& input_flag_complex_realspace_sum,
      const unsigned int& input_flag_image_output,
      const unsigned int& input_flag_netcdf_images,
      const unsigned int& input_flag_netcdf_variance,
      const unsigned int& input_flag_debug,
      // parameters taken from simulation of system evolution:
      //const std::list< slice* >& slicesOfScatterers,
      std::list< slice* >& slicesOfScatterers,
      //const size_t& total_population,
      const double& bwcutoff_pap,
      const double& bwcutoff_t,
      const double& azimuthal_binning_size_factor,
      //const double& xperiod, const double& yperiod, const double& zperiod,
      const double& xmin,      // lowest value of x domain
      const double& ymin,      // lowest value of y domain
      //const double& zmin,      // lowest value of z domain
      //const double& raster_spacing,// used for both x & y raster steps
      const std::list<double>& x_p,
      const std::list<double>& y_p,
      const double& xperiod,
      const double& yperiod,
      const double& xperiod_duped,
      const double& yperiod_duped,
      const double* const xx_local,
      const double* const xx_joined,
      const double* const yy,
      const double* const kx_local,
      const double* const kx_joined,
      const ptrdiff_t& Nx,  // samples for each dimensional discretization
      const ptrdiff_t& Nx_local, // samples for each dimensional discretization
      const double* const ky,
      const ptrdiff_t& Ny, 
      // parameters of the tem microscope:
      //const double& VV,       // accelerating voltage
      const double& Cs3,       // third order spherical aberration
      const double& Cs5,       // fifth order spherical aberration
      const double& defocus,  // focal point distance from sample bottom
      const double& defocus_spread,
      const double& alpha_max, // objective aperture limiting angle
      const double& detector_inner_angle, // detector dimension [\AA^{-1}]
      const double& detector_outer_angle, // detector dimension [\AA^{-1}]
      const double& lambda,
      const string& outFileName_prefix,
      const ptrdiff_t& local_alloc_size_fftw,
      const ptrdiff_t& local_0_start_fftw,
      //const ptrdiff_t& probe_idx_local_start_x,
      //const int& threads_ok,
      //const int& nthreads,
      const int& mynode,      // MPI rank of local node
      const int& rootnode,    // MPI rank of root node
      const int& totalnodes,
      MPI_Comm comm
      )
{

   const double lambda_sqr = lambda * lambda;
   const double alpha_max_sqr = alpha_max * alpha_max;
   
   size_t resolutionUnit, resolutionUnit_recip;
   double xResolution, yResolution,             // for real space images
          xResolution_stem, yResolution_stem,   // for real space images
          xResolution_recip, yResolution_recip; // for diffraction images

   double diffraction_scale_factor;

   //if ( input_flag_image_output || input_flag_netcdf_output ) 
   //{
   // intialize tif output class
   //writeImage imageWriter;  // not used by this function

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

   //MPI_Bcast(&xResolution, 1, MPI_DOUBLE, rootnode, comm);
   //MPI_Bcast(&yResolution, 1, MPI_DOUBLE, rootnode, comm);
   //MPI_Bcast(&xResolution_stem, 1, MPI_DOUBLE, rootnode, comm);
   //MPI_Bcast(&yResolution_stem, 1, MPI_DOUBLE, rootnode, comm);
   //MPI_Bcast(&xResolution_recip, 1, MPI_DOUBLE, rootnode, comm);
   //MPI_Bcast(&yResolution_recip, 1, MPI_DOUBLE, rootnode, comm);

   MPI_Bcast(bcast_resolutions, 6, MPI_DOUBLE, rootnode, comm);
   //}


   ///////////////////////////////////////////////////////////////////
   // Split the sample along the transmission axis into slices
   ///////////////////////////////////////////////////////////////////
   // Notes: 
   // "If the slice thickness in the multislice calculation matches 
   // the natural periodicity in the specimen (i.e., an integer number
   // of slices in the repeat length of the specimen) then the 
   // multislice simulation will reproduce the higher order Laue 
   // zones. If the slice thickness does not match the specimen 
   // periodicity then beating can occur between the slice thickness 
   // and the specimen periodicity to produce artifacts in the image.
   // The slice thickness can produce an artificial periodicity in the
   // specimen if it does not match the natural periodicity of the 
   // specimen. If the multislice slice thickness is much larger than
   // the natural periodicity of the speciment then false HOLZ lines 
   // can be created at: 
   // k = \sqrt(\frac{2}{\Delta z \lambda}) eqn 6.95, 
   // where \Delta z is the multislice slice thickness, \lambda 
   // is the electron wavelength and \alpha = \lambda k is the 
   // electron scattering angle. The slice thickness effectively takes
   // the place of the normal crystal lattice spacing along z. An 
   // overly large slice thickness can produce serious artifacts in a
   // simulated image and 
   // should be avoided." 
   // ...
   // "The standard multislice method [(6.92) and (6.93)] is only 
   // accurate to order \Delta z." ... "As long as the slice thickness
   // is large 
   // compared to the effective range of the atomic potential this 
   // approximation is valid. However, if the slice thickness \Delta z
   // is made less than the range of the atomic potential (about 1 
   // \AA) then this approximation is no longer valid. Therefore, if 
   // the total projected atomic potential is used (as is typical) the
   // minimum slice thickness is about 1 \AA, which sets a limit on
   // the maximum achievable accuracy of the multislice calculation.
   // There is some question if 1 \AA is thin enough for heavy atoms,
   // such as gold, at 100 keV, however higher voltage or lower atomic
   // number should be all right (Watanabe [372]). The electron
   // wavelength is typically of order 0.03 \AA, which is much smaller
   // than the minimum slice thickness. If the full wave function were
   // used then the slice thickness would have to be much smaller than
   // the electron wavelength and the multislice method would not work
   // at all. However, only the slowly varying part of the wave 
   // function is calculated so a slice thickness of several 
   // Angstroms is acceptable."
   // - Kirkland (2010)


   // NOTE: locations of slicing are passed by parent function

   const double NxNy = Nx * Ny;
   const double NxNy_sqr = NxNy * NxNy;
   double sqrtNxNy; 
   if ( Nx == Ny ) sqrtNxNy =  Nx ;   
   else sqrtNxNy = sqrt( NxNy );


   //size_t sliceNumber = 0; //debug

   // TODO: make these commands specific to STEM imaging conditional 
   //       (not needed for fem)
   ///////////////////////////////////////////////////////////////////
   // Instantiate the STEM image
   ///////////////////////////////////////////////////////////////////
   double* stem_image;
   if ( mynode == rootnode )
      stem_image = new double[ x_p.size() * y_p.size() ];

   ///////////////////////////////////////////////////////////////////
   // Instantiate the electron wave
   ///////////////////////////////////////////////////////////////////
   fftw_complex* psi;
   psi = fftw_alloc_complex( local_alloc_size_fftw );

   ///////////////////////////////////////////////////////////////////
   // Instantiate the large STEM probe
   ///////////////////////////////////////////////////////////////////

   //////////////////////////////////////////////////////////////////
   // Allocate local variables and domains required for the large probe 
   const ptrdiff_t large_probe_factor = 1;

   const ptrdiff_t Nx_large = large_probe_factor * Nx;
   const ptrdiff_t Ny_large = large_probe_factor * Ny;
   ptrdiff_t Nx_large_local;
   ptrdiff_t local_alloc_size_fftw_large;
   ptrdiff_t local_0_start_fftw_large;

   local_alloc_size_fftw_large
      = fftw_mpi_local_size_2d(     // fftw will be 2-D here
            Nx_large, Ny_large,
            MPI_COMM_WORLD,
            &Nx_large_local,
            &local_0_start_fftw_large );

   const double xperiod_duped_large 
                  = large_probe_factor * xperiod_duped;
   const double yperiod_duped_large 
                  = large_probe_factor * yperiod_duped;

   double* kx_large_joined; kx_large_joined = new double[ Nx_large ];
   double* xx_large_joined; xx_large_joined = new double[ Nx_large ];
   double* kx_large_local; kx_large_local = new double[ Nx_large_local ];
   double* xx_large_local; xx_large_local = new double[ Nx_large_local ];
   double* ky_large; ky_large = new double[ Ny_large ];
   double* yy_large; yy_large = new double[ Ny_large ];

   double delta_x, delta_y, delta_kx, delta_ky;
   delta_x = xperiod_duped_large / Nx; 
   delta_y = yperiod_duped_large / Ny; 
   delta_kx = 1.0/xperiod_duped_large;
   delta_ky = 1.0/yperiod_duped_large;

   //const double kxperiod_large = Nx_large / xperiod_duped; 
   //const double kyperiod_large = Ny_large / yperiod_duped;
   const double kxperiod_large = Nx_large / xperiod_duped_large; 
   const double kyperiod_large = Ny_large / yperiod_duped_large;

   // larger probe than used in stem
   fftw_complex* large_probe;
   large_probe = fftw_alloc_complex( local_alloc_size_fftw_large );

   // fftw plan is required to transform the cached probe into realspace
   fftw_plan pb_c2c_large_probe_split;

   // c2c in-place reverse fft
   pb_c2c_large_probe_split = fftw_mpi_plan_dft_2d( 
                           Nx_large, Ny_large,
                           large_probe, large_probe,
                           //comm, FFTW_BACKWARD, FFTW_ESTIMATE );
                           //comm, FFTW_BACKWARD, FFTW_PATIENT );
                           //comm, FFTW_BACKWARD, FFTW_EXHAUSTIVE );
                           comm, FFTW_BACKWARD, FFTW_MEASURE );
   ////////////////////////////////////////////////////////////////////

   // TODO: Compensate for periodic boundary conditions on the probe
   // TODO: Evaluate the probe on a much larger domain to reduce
   //       boundary effects, then transform it into realspace, 
   //       copy a subsection of the probe into a smaller domain and 
   //       transform it back into reciprocal space.

   ////////////////////////////////////////////////////////////////////
   // evaluate and scatter large domains needed to create large probe
   if ( mynode == rootnode )
   {
      // The reciprocal space domain my be restricted to 2-D, since the
      //  Fourier projection theorem is being used.

      domain_2D_periodic_recip( Nx_large, Ny_large, 
            xperiod_duped_large, yperiod_duped_large, 
            kxperiod_large, kyperiod_large, 
            delta_kx, delta_ky,
            kx_large_joined, ky_large
            );

      domain_2D_periodic( Nx_large, Ny_large, 
            xperiod_duped_large, yperiod_duped_large, 
            xmin, ymin, 
            delta_x, delta_y,
            xx_large_joined, yy_large
            );
   }

   MPI_Scatter( xx_large_joined, Nx_large_local, MPI_DOUBLE, // 2-D
               xx_large_local, Nx_large_local, MPI_DOUBLE,
               rootnode, MPI_COMM_WORLD);
   MPI_Scatter( kx_large_joined, Nx_large_local, MPI_DOUBLE, // 2-D
               kx_large_local, Nx_large_local, MPI_DOUBLE,
               rootnode, MPI_COMM_WORLD);
   MPI_Bcast( ky_large, Ny_large, MPI_DOUBLE, rootnode, MPI_COMM_WORLD);
   MPI_Bcast( yy_large, Ny_large, MPI_DOUBLE, rootnode, MPI_COMM_WORLD);
   MPI_Bcast( kx_large_joined, Nx_large, 
         MPI_DOUBLE, rootnode, MPI_COMM_WORLD);
   MPI_Bcast( xx_large_joined, Nx_large, 
         MPI_DOUBLE, rootnode, MPI_COMM_WORLD);
   ////////////////////////////////////////////////////////////////////

   ////////////////////////////////////////////////////////////////////
   // evaluate the probe values in reciprocal space
   if ( input_flag_aberration_correction )
   {
      if ( mynode == rootnode 
            && input_flag_debug )
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
            << xx_large_joined[0] << ", " << yy_large[0] << ")" 
            << endl << "   Nx_large_local, Ny_large : "
            << Nx_large_local << ", " << Ny_large
            << endl;

      probe_wavefunction_correctedtoCs5_unnormalized(
            xx_large_joined[0], yy_large[0],
            kx_large_local, Nx_large_local, ky_large, Ny_large, 
            Cs3, Cs5, defocus, alpha_max_sqr, 
            lambda, lambda_sqr,
            large_probe
            );

      if ( input_flag_debug )
         cout << "Evaluated probe ..." << endl;// debug
   }
   else
   {
      if ( mynode == rootnode 
            && input_flag_debug )
         cout << "calling"
            << " probe_wavefunction_uncorrected_unnormalized"
            << endl << "   Cs3 : " << Cs3
            << endl << "   defocus : " << defocus
            << endl << "   alpha_max_sqr : "
            << alpha_max_sqr 
            << endl << "   lambda, lambda_sqr : "
            << lambda << ", " << lambda_sqr
            << endl << "   position : ("
            << xx_large_joined[0] << ", " << yy_large[0] << ")" 
            << endl << "   Nx_large_local, Ny_large : "
            << Nx_large_local << ", " << Ny_large
            << endl;

      probe_wavefunction_uncorrected_unnormalized(
            xx_large_joined[0], yy_large[0],
            kx_large_local, Nx_large_local, ky_large, Ny_large, 
            Cs3, defocus, alpha_max_sqr, 
            lambda, lambda_sqr,
            large_probe
            );

      if ( input_flag_debug )
         cout << "Evaluated probe ..." << endl;// debug
   }

   delete[] kx_large_joined;
   delete[] xx_large_joined;
   delete[] kx_large_local;
   delete[] xx_large_local;
   delete[] ky_large;
   delete[] yy_large;
   ////////////////////////////////////////////////////////////////////

   ///////////////////////////////////////////////////////////////////
   // do this further below ... - Normalize the probe
   ///////////////////////////////////////////////////////////////////
   // Periodic boundary conditions ensure that the probe norm
   //  is invariant with probe position.
   //if ( mynode == rootnode && input_flag_debug )
   //   cout << "Calculating probe norm in reciprocal space ... " << endl;

   //double Ap; Ap = 0.0; // probe normalization factor
   //probe_wavefunction_norm(
   //      Ap, 
   //      Nx_local, Ny, 
   //      psi,
   //      comm
   //      );

   //if ( mynode == rootnode && input_flag_debug )
   //   cout << "Normalizing probe ..." << endl;

   //for (ptrdiff_t i=0; i < local_alloc_size_fftw; ++i)
   //{  // reassign probe values for MPI_Allgather while removing fftw scale
   //   //  factor
   //   psi[i][0] = Ap * psi[i][0];
   //   psi[i][1] = Ap * psi[i][1];
   //}
   //if (mynode == rootnode ) // debug
   //   cout << "Probe norm in reciprocal space : " << Ap << endl; // debug

   ///////////////////////////////////////////////////////////////////
   // Transform the split probe into realspace
   if ( mynode == rootnode && input_flag_debug )
      cout << "Transforming probe to realspace ..." << endl;

   fftw_execute( pb_c2c_large_probe_split );
   ///////////////////////////////////////////////////////////////////

   ///////////////////////////////////////////////////////////////////
   // Gather the large probe realspace wavefunction onto each node
   ///////////////////////////////////////////////////////////////////

   double* large_probe_split_re;
   large_probe_split_re = new double[ local_alloc_size_fftw_large ];
   double* large_probe_split_im;
   large_probe_split_im = new double[ local_alloc_size_fftw_large ];


   // TODO: limit the probe in realspace, to avoid aliasing 
   //       artifacts in reciprocal space image. Since this is multiplied
   //       against the transmission function in realspace, there is no
   //       need to limit the transmission function in realspace.
   ///////////////////////////////////////////////////////////////////
   // Limit the STEM probe extent in realspace
   ///////////////////////////////////////////////////////////////////
   // TODO: the following is just a test and should be deleted
   //  Also try realspace limiting the propagator and 
   //   transmission functions.
   //if ( mynode == rootnode && input_flag_debug )
   //   cout << "Cutting probe in realspace, centering at" 
   //      << "  (xx_local[0], yy[0]): ("
   //         << xx_local[0] << ", " << yy[0] << ")" << endl;

   //double probecutoff; probecutoff = 0.125;
   //for ( ptrdiff_t i=0; i<Nx_local; ++i)
   //   for ( ptrdiff_t j=0; j<Ny; ++j)
   //   {
   //      if (
   //            (pow(xx_joined[0] - xx_local[i],2) 
   //                + pow(yy[0] - yy[j],2)
   //               > pow(probecutoff * xperiod_duped,2) 
   //                  + pow(probecutoff * yperiod_duped,2)
   //            )
   //            &&
   //            (pow(xx_joined[0] + xperiod_duped - xx_local[i],2) 
   //                + pow(yy[0] - yy[j],2)
   //               > pow(probecutoff * xperiod_duped,2) 
   //                  + pow(probecutoff *yperiod_duped,2)
   //            )
   //            &&
   //            (pow(xx_joined[0] - xx_local[i],2) 
   //                + pow(yy[0] + yperiod_duped - yy[j],2)
   //               > pow(probecutoff * xperiod_duped,2) 
   //                  + pow(probecutoff * yperiod_duped,2)
   //            )
   //            &&
   //            (pow(xx_joined[0] + xperiod_duped - xx_local[i],2) 
   //                + pow(yy[0] + yperiod_duped - yy[j],2)
   //               > pow(probecutoff * xperiod_duped,2) 
   //                  + pow(probecutoff * yperiod_duped,2)
   //            )
   //            &&
   //            (pow(xx_joined[0] - xperiod_duped - xx_local[i],2) 
   //                + pow(yy[0] - yy[j],2)
   //               > pow(probecutoff * xperiod_duped,2) 
   //                  + pow(probecutoff * yperiod_duped,2)
   //            )
   //            &&
   //            (pow(xx_joined[0] - xx_local[i],2) 
   //                + pow(yy[0] - yperiod_duped - yy[j],2)
   //               > pow(probecutoff * xperiod_duped,2) 
   //                  + pow(probecutoff * yperiod_duped,2)
   //            )
   //            &&
   //            (pow(xx_joined[0] - xperiod_duped - xx_local[i],2) 
   //                + pow(yy[0] - yperiod_duped - yy[j],2)
   //               > pow(probecutoff * xperiod_duped,2) 
   //                  + pow(probecutoff * yperiod_duped,2)
   //            )
   //            &&
   //            (pow(xx_joined[0] + xperiod_duped - xx_local[i],2) 
   //                + pow(yy[0] - yperiod_duped - yy[j],2)
   //               > pow(probecutoff * xperiod_duped,2) 
   //                  + pow(probecutoff * yperiod_duped,2)
   //            )
   //            &&
   //            (pow(xx_joined[0] - xperiod_duped - xx_local[i],2) 
   //                + pow(yy[0] + yperiod_duped - yy[j],2)
   //               > pow(probecutoff * xperiod_duped,2) 
   //                  + pow(probecutoff * yperiod_duped,2)
   //            )
   //         )
   //      {
   //         psi[j + i*Ny][0] = 0.0;
   //         psi[j + i*Ny][1] = 0.0;
   //      }
   //   }
   // TODO: the above is just a test and should be deleted


   for (ptrdiff_t i=0; i < local_alloc_size_fftw_large; ++i)
   {  // reassign probe values for MPI_Allgather while removing fftw scale
      //  factor
      large_probe_split_re[i] = large_probe[i][0] / sqrtNxNy;
      large_probe_split_im[i] = large_probe[i][1] / sqrtNxNy;
   }

   if ( mynode == rootnode && input_flag_debug )
      cout << "Allgathering enlarged probe ..." << endl;

   double* large_probe_joined_re;
   large_probe_joined_re = new double[Nx_large * Ny_large];
   double* large_probe_joined_im;
   large_probe_joined_im = new double[Nx_large * Ny_large];

   MPI_Allgather(
         large_probe_split_re,
         local_alloc_size_fftw_large, MPI_DOUBLE,
         large_probe_joined_re,
         local_alloc_size_fftw_large, MPI_DOUBLE,
         comm);

   MPI_Allgather(
         large_probe_split_im,
         local_alloc_size_fftw_large, MPI_DOUBLE,
         large_probe_joined_im,
         local_alloc_size_fftw_large, MPI_DOUBLE,
         comm);

   delete[] large_probe_split_re;
   delete[] large_probe_split_im;

   ///////////////////////////////////////////////////////////////////
   // Copy a reduced portion of large probe onto the wavefunction
   ///////////////////////////////////////////////////////////////////

   // stem_probe_joined_re, stem_probe_joined_im : will hold cached
   //    realspace values to be shifted and assigned to stem_probe_split
   double* stem_probe_joined_re; 
   stem_probe_joined_re = new double[Nx * Ny];
   double* stem_probe_joined_im; 
   stem_probe_joined_im = new double[Nx * Ny];

   for ( ptrdiff_t i=0; i<Nx/2; ++i)
   {
      for ( ptrdiff_t j=0; j < Ny/2; ++j)
      {
         stem_probe_joined_re[j + i*Ny] 
            = large_probe_joined_re[j + i*Ny_large];
         stem_probe_joined_im[j + i*Ny] 
            = large_probe_joined_im[j + i*Ny_large];
      }
      for ( ptrdiff_t j=Ny/2; j < Ny; ++j)
      {
         stem_probe_joined_re[j + i*Ny] 
            = large_probe_joined_re[(Ny_large - Ny + j) + i*Ny_large];
         stem_probe_joined_im[j + i*Ny] 
            = large_probe_joined_im[(Ny_large - Ny + j) + i*Ny_large];
      }
   }
   for ( ptrdiff_t i=Nx/2; i<Nx; ++i)
   {
      for ( ptrdiff_t j=0; j < Ny/2; ++j)
      {
         stem_probe_joined_re[j + i*Ny] 
            = large_probe_joined_re[j + (Nx_large - Nx + i)*Ny_large];
         stem_probe_joined_im[j + i*Ny] 
            = large_probe_joined_im[j + (Nx_large - Nx + i)*Ny_large];
      }
      for ( ptrdiff_t j=Ny/2; j < Ny; ++j)
      {
         stem_probe_joined_re[j + i*Ny] 
            = large_probe_joined_re[
                  (Ny_large - Ny + j) + (Nx_large - Nx + i)*Ny_large
               ];
         stem_probe_joined_im[j + i*Ny] 
            = large_probe_joined_im[
                  (Ny_large - Ny + j) + (Nx_large - Nx + i)*Ny_large
               ];
      }
   }

   delete[] large_probe_joined_re;
   delete[] large_probe_joined_im;

   fftw_destroy_plan( pb_c2c_large_probe_split );
   fftw_free( large_probe );

   ///////////////////////////////////////////////////////////////////
   // Normalize the probe in real space
   ///////////////////////////////////////////////////////////////////
   // Periodic boundary conditions ensure that the probe norm
   //  is invariant with probe position.
   if ( mynode == rootnode && input_flag_debug )
      cout << "Calculating probe norm in real space ... " << endl;

   double Ap; Ap = 0.0; // probe normalization factor
   probe_wavefunction_norm(
         Ap, 
         Nx_local, Ny, 
         //psi,
         stem_probe_joined_re,
         stem_probe_joined_im,
         comm
         );
   if (mynode == rootnode ) // debug
      cout << "Probe norm in real space : " << Ap << endl;//debug
   for ( ptrdiff_t i=0; i< Nx * Ny; ++i)
   {
      stem_probe_joined_re[i] = Ap * stem_probe_joined_re[i];
      stem_probe_joined_im[i] = Ap * stem_probe_joined_im[i];
   }


   ///////////////////////////////////////////////////////////////////
   // Instantiate the first and second moment of diffracted intensity
   //  for FTEM calculation.
   ///////////////////////////////////////////////////////////////////
   double* diffracted_wave_mag_sum;       // Nx*Ny / number of nodes
   double* diffracted_wave_mag_sqr_sum;   // Nx*Ny / number of nodes
   double* diffracted_wave_re_sum;        // Nx*Ny / number of nodes
   double* diffracted_wave_im_sum;

   double* diffracted_wave_radial_intensity_sum_local;
   double* diffracted_wave_radial_intensity_sum;
   double* diffracted_wave_radial_intensity_sqr_sum;
   double* diffracted_wave_radial_intensity_sum_tmp;

   double* ftem_variance;

   unsigned int number_of_bins;  // azimuthal integration subintervals

   int* bin_counts_local;  // not size_t since MPI datatype required
   int* bin_counts_aggregated;

   double delta_k;   // length of aziumthal integration subintervals
   std::vector<double> binning_boundaries;

   if ( input_flag_fem )
   {
      // used for 2-D variance image
      diffracted_wave_mag_sum 
         = new double[ local_alloc_size_fftw ];
      diffracted_wave_mag_sqr_sum 
         = new double[ local_alloc_size_fftw ];

      if (input_flag_complex_realspace_sum) 
      {
         // used 2-D variance image when adding phase in the image plane
         diffracted_wave_re_sum 
            = new double[ local_alloc_size_fftw ];
         diffracted_wave_im_sum 
            = new double[ local_alloc_size_fftw ];
      }
      
      for ( ptrdiff_t ii=0; ii < local_alloc_size_fftw; ++ii )
      {// ensure that they are initialized to 0.0
         // 2-D image
         diffracted_wave_mag_sum[ii] = 0.0;
         diffracted_wave_mag_sqr_sum[ii] = 0.0;

         if (input_flag_complex_realspace_sum) 
         {
            diffracted_wave_re_sum[ii] = 0.0;
            diffracted_wave_im_sum[ii] = 0.0;
         }
      }

      // Initalize binning_boundaries for azimuthal integration.
      //
      // TODO: determine a good separation between binning boundaries,
      //        delta_k. The current value is just a guess.
      //if ( kx_local[1] - kx_local[0] > ky[1] - ky[0] ) 

      delta_k =   
         azimuthal_binning_size_factor * sqrt( 
            (kx_local[1] - kx_local[0]) * (kx_local[1] - kx_local[0])
             + (ky[1] - ky[0]) * (ky[1] - ky[0])
            );

      if ( mynode == rootnode ) number_of_bins = bwcutoff_t / delta_k;
      MPI_Bcast( &number_of_bins, 1, MPI_UNSIGNED, rootnode, comm);

      if ( mynode == rootnode )
         for ( unsigned int i=0; i< number_of_bins + 1; ++i)
            binning_boundaries.push_back( i * delta_k );
      else
         binning_boundaries.resize(number_of_bins + 1);

      MPI_Bcast( &binning_boundaries[0], number_of_bins + 1, 
                  MPI_DOUBLE, rootnode, comm);

      // 1-D variance curve variables allocated after number_of_bins
      diffracted_wave_radial_intensity_sum_local
         = new double[ number_of_bins ];
      bin_counts_local = new int[number_of_bins];

      for ( ptrdiff_t ii=0; ii < number_of_bins; ++ii )
      {
         diffracted_wave_radial_intensity_sum_local[ii] = 0.0;
         bin_counts_local[ii] = 0;
      }

      if ( mynode == rootnode )
      {
         diffracted_wave_radial_intensity_sqr_sum
            = new double[ number_of_bins ];
         diffracted_wave_radial_intensity_sum
            = new double[ number_of_bins ];
         diffracted_wave_radial_intensity_sum_tmp // temporary storage
            = new double[ number_of_bins ];       // to allow squaring
         bin_counts_aggregated = new int[ number_of_bins ];
         for ( ptrdiff_t ii=0; ii < number_of_bins; ++ii )
         {
            diffracted_wave_radial_intensity_sqr_sum[ii] = 0.0;
            diffracted_wave_radial_intensity_sum[ii] = 0.0;
            diffracted_wave_radial_intensity_sum_tmp[ii] = 0.0;
            bin_counts_aggregated[ii] = 0;
         }
      }
 
   }


   ///////////////////////////////////////////////////////////////////
   // Instantiate the fftw plan for psi
   ///////////////////////////////////////////////////////////////////
   //if ( threads_ok ) fftw_plan_with_nthreads( nthreads );
   fftw_plan pf_c2c_psi, pb_c2c_psi, pf_c2c_t, pb_c2c_t;//,pb_c2c_pap;


   // Notes regarding the multislice method: 
   //       eqn 6.92:
   //       \psi_{n+1}(x,y) = FT^{-1}[ 
   //             P_{n}(k_{x}, k_{y}, \Delta z_{n})
   //             FT[ t_{n}(x,y) \psi_{n}(x,y) ]
   //          ]  
   //
   // Each slice requires a transmission function
   //       t_{n}(\vec{x} = \exp[i \sigma v_{zn}(\vec{x})]
   // The propagator function is 
   //    eqn 6.90
   //    p(x,y,\Deltaz) = \frac{1}{i \lambda \Delta z}
   //                     \exp\left[
   //                       \frac{i \pi}{\lambda \Delta z} (x^{2} + y^{2})
   //                      \right]
   //       
   //   and might be identical for each slice having the same thickness
   //
   // Condition to eliminate aliasing:
   //
   // bwcutoff_t + bwcutoff_propagator + bwcutoff_psi < 2*k_max
   //
   // Sequence of operations:
   //
   //  for each slice n {
   //    - evaluate \sigma V_{zn}(k_{x}, k_{y}) 
   //    - transform \sigma V_{zn}(k_{x}, k_{y}) into \sigma v_{zn}(x,y) 
   //    - evaluate eqn 6.71
   //      P_{n}(k_{x}, k_{y}, \Delta z)
   //         = \exp[ -i \pi \lambda (k^{2}_{x} + k^{2}_{y}) \Delta z_{n}]
   //            ( k_{x} is split over nodes, and thus P_{n} is too )
   //      Note: the only benefit of making P_{n} a stored member of each 
   //      slice is to allow its reuse when \Delta z values are 
   //      duplicates.
   //      If \Delta z values are all unique, you might as well just
   //        evaluate P_{n} during wave propagation and save the memory.
   //  }
   //
   //  - initialize \psi_{0}(x,y) 
   //    * size of \psi on local node is the same as pap (kx_local * ky)
   //    * alpha_max functions as a reciprocal space bandwidth limit of 
   //       the stem probe
   //
   //  evaluate in order of slice number n=0,1, ..., {
   //    - evaluate t_{n}(x,y) \psi_{n}(x,y) 
   //                  = \exp[i \sigma v_{zn}(x,y)] \psi_{n}(x,y)
   //          ( v_{zn}(x,y) is split over nodes, and thus \psi_{n}(x,y) 
   //           needs to be split over nodes too, but for n=0 it's just 1,
   //           so no x_local domain values are necessary )
   //    - transform into FT[ t_{n}(x,y) \psi_{n}(x,y) ]
   //            ( FT[ ... ] is split over nodes since \psi is )
   //    - evaluate P_{n}(k_{x}, k_{y}, \Delta z) 
   //               * FT[ t_{n}(x,y) \psi_{n}(x,y) ]
   //    - transform into \psi_{n+1}(x,y) = 
   //         FT^{-1}[ 
   //            P(k_{x}, k_{y}, \Delta z) 
   //            FT[ 
   //               t_{n}(x,y) \psi_{n}(x,y)
   //            ]
   //         ]
   //  }
   //
   //  - transform \psi(x,y) into \Psi(k_{x},k_{y}) (for BFCTEM)
   //  - evaluate H_0(k_{x},k_{y}) \Psi(k_{x},k_{y}) (for BFCTEM)
   //  - transform into the image wave function (for BFCTEM)
   //       \psi_{i}(x,y) =  FT^{-1}[ H_0(k_{x},k_{y}) \Psi(k_{x},k_{y}) ]
   
   //=================================================================
   // For each point p of the sample surface, translate the 
   //  probe wavefunction and propagate it through every slice
   //=================================================================
   size_t probe_idx_shift_x; probe_idx_shift_x = 0;
   size_t probe_idx_shift_y; probe_idx_shift_y = 0;
   size_t probe_idx_x; probe_idx_x = 0;
   size_t probe_idx_y; probe_idx_y = 0;

   bool first_probe_flag; first_probe_flag  = true;
   double detected_intensity, detected_intensity_reduced;
   size_t pixel_number_x = 0;
   size_t pixel_number_y = 0;
   size_t number_of_raster_points //= N_pixels_x * N_pixels_y;
                            = x_p.size() * y_p.size();

   double half_delta_x = 0.5 * (xx_joined[1] - xx_joined[0]);
   double half_delta_y = 0.5 * (yy[1] - yy[0]);


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
      {
         // TODO: Consider using the stem raster spacing to increment 
         //       beam position instead of using x_p and y_p. 

         // debug
         //if ( mynode == rootnode )
         //{
         ////cout << "node : " << mynode << ", ";
         //   cout << "probing position : " 
         //      //<< pixel_location_x << ", " 
         //      << *x_itr << ", " 
         //      //<< pixel_location_y
         //      << *y_itr
         //      //<< y_p[j] 
         //      << endl;
         //}
         // end debug
         
         /////////////////////////////////////////////////////////////
         // Create forward and reverse fftw plans for psi
         /////////////////////////////////////////////////////////////
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

         //////////////////////////////////////////////////////////////
         // Translate the STEM probe in realspace 
         //////////////////////////////////////////////////////////////
         // Determine indices of xx_joined[] and yy[] which represent
         //  the point closest to (*x_itr, *y_itr).
         // The following assumes that xx_joined[0] == min(xx_joined[]),
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

         // copy stem_probe_joined_re,im to psi[][0,1] with 
         //  shifted x and y indices 

         // Those indices should be the same as the number of positions
         //  which the cached probe array should be shifted.

         // Iterate over the domains, assign probe locally to
         //  the variable which will be used in the fft.
         for ( size_t i=0; i < Nx_local; ++i)
         {
            // enforce periodic boundary condition, x direction
            //if ( i + idx_local_start_x < probe_idx_shift_x )
            if ( i + local_0_start_fftw < probe_idx_shift_x )
               probe_idx_x = i + Nx + local_0_start_fftw 
                              - probe_idx_shift_x;
               //idx_x = i + Nx + idx_local_start_x - probe_idx_shift_x;
            else
               probe_idx_x = i + local_0_start_fftw - probe_idx_shift_x;
               //idx_x = i + idx_local_start_x - probe_idx_shift_x;

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
                     << ", " << probe_idx_y + probe_idx_x * Ny << ", " 
                     << local_alloc_size_fftw
                     << ") exceeds bounds (Nx,Ny)"
                     << Nx << ", " << Ny << ")" << endl;
               // end debug

               psi[j + i * Ny][0] = 
                  stem_probe_joined_re[ probe_idx_y + probe_idx_x * Ny ];

               psi[j + i * Ny][1] = 
                   stem_probe_joined_im[ probe_idx_y + probe_idx_x * Ny ];
            }
         }

         //TODO: uncomment this to see if evaluating the probe at each
         //       STEM image point gives different results than copying it
         //if ( input_flag_aberration_correction )
         //{
         //   if ( mynode == rootnode && first_probe_flag 
         //         && input_flag_debug )
         //      cout 
         //         << "calling"
         //         << " probe_wavefunction_correctedtoCs5_unnormalized"
         //         << endl << "   Cs3 : " << Cs3
         //         << endl << "   Cs5 : " << Cs5
         //         << endl << "   defocus : " << defocus
         //         << endl << "   alpha_max_sqr : "
         //         << alpha_max_sqr 
         //         << endl << "   *x_itr: " << *x_itr
         //         << endl << "   *y_itr: " << *y_itr
         //         << endl;

         //   probe_wavefunction_correctedtoCs5_unnormalized(
         //         *x_itr, *y_itr,
         //         kx_local, Nx_local, ky, Ny, 
         //         Cs3, Cs5, defocus, alpha_max_sqr, 
         //         lambda, lambda_sqr,
         //         psi
         //         );
         //}
         //else
         //{
         //   if ( mynode == rootnode && first_probe_flag 
         //         && input_flag_debug )
         //      cout << "calling"
         //         << " probe_wavefunction_uncorrected_unnormalized"
         //         << endl << "   Cs3 : " << Cs3
         //         << endl << "   defocus : " << defocus
         //         << endl << "   alpha_max_sqr : "
         //         << alpha_max_sqr 
         //         << endl << "   *x_itr: " << *x_itr
         //         << endl << "   *y_itr: " << *y_itr
         //         << endl;

         //   probe_wavefunction_uncorrected_unnormalized(
         //         *x_itr, *y_itr,
         //         kx_local, Nx_local, ky, Ny, 
         //         Cs3, defocus, alpha_max_sqr, 
         //         lambda, lambda_sqr,
         //         psi
         //         );
         //}

         //if ( first_probe_flag ) 
         //{  // Only evaluate the probe norm once. 
         //   // Periodic boundary conditions ensure that the probe norm
         //   //  is invariant with probe position.
         //   //cout << "evaluating probe norm Ap" << endl; // debug
         //   probe_wavefunction_norm(
         //         Ap, 
         //         Nx_local, Ny, 
         //         psi, 
         //         comm
         //         );
         //}

         //cout << "Ap : " << Ap << endl; // debug

         // Multiply psi by Ap ( == (\int |psi|^2)^{-1} ) and fftw scale
         //cout << "Probe norm Ap: " << Ap << endl; // debug
         //for ( ptrdiff_t i=0; i<local_alloc_size_fftw; ++i)
         //{
         //   psi[i][0] = Ap * psi[i][0];// / sqrtNxNy;
         //   psi[i][1] = Ap * psi[i][1];// / sqrtNxNy;
         //}
   
         // Note: there is no need to bandwidth limit the stem probe in
         //    reciprocal space, alpha_max (probe forming aperture) 
         //    already limits the probe in reciprocal space.

         // TODO : consider elliminating the following block
         if( //first_probe_flag 
             //  && 
               input_flag_debug 
               &&
               (input_flag_image_output || input_flag_netcdf_images))
         {
            //if ( mynode == rootnode && input_flag_debug 
            //      && input_flag_image_output )
            //   cout << "saving initial probe in reciprocal space" << endl;

            //if ( input_flag_image_output )
            //{
            //   diffraction_scale_factor = 1.0e30;
            //   output_diffraction_with_renormalization(
            //         psi,
            //         diffraction_scale_factor,
            //         local_alloc_size_fftw,
            //         Nx_local, Nx, Ny,
            //         resolutionUnit_recip,
            //         xResolution_recip, yResolution_recip,
            //         outFileName_prefix + "_initial_probe",
            //         mynode, rootnode, comm
            //         );
            //}

            // TODO: uncomment if evaluating the probe rather than copying
            //        it
            //fftw_execute( pb_c2c_psi ); // debug

            //for ( ptrdiff_t i=0; i<local_alloc_size_fftw; ++i)
            //{
            //   psi[i][0] = psi[i][0] / sqrtNxNy;
            //   psi[i][1] = psi[i][1] / sqrtNxNy;
            //}

            if ( input_flag_image_output )
            {
               if ( mynode == rootnode && input_flag_debug )  // debug
                  cout << "saving initial probe to tiff" << endl; // debug

               debug_output_complex_fftw_operand(//TODO: rename this function
                  psi,
                  0, 
                  local_alloc_size_fftw,
                  Nx_local, Nx, Ny,
                  resolutionUnit,
                  xResolution, yResolution,
                  outFileName_prefix 
                  + "_initial_probe_realspace",
                  mynode, rootnode, comm);
            }

            if ( input_flag_netcdf_images )
            {
               if ( mynode == rootnode && input_flag_debug)
                  cout << "saving initial probe to netCDF" << endl;

               if( output_psi_realspace_to_netcdf(
                     psi,
                     local_alloc_size_fftw,
                     Nx_local, kx_joined, Nx, ky, Ny,
                     outFileName_prefix 
                        + "_initial_probe_realspace",
                     mynode, rootnode, comm
                     ) == EXIT_FAILURE)
                  cout << "output_psi_realspace_to_netcdf() failed" 
                     << endl;
            }
            // TODO: uncomment if evaluating the probe rather than copying
            //        it
            //fftw_execute( pf_c2c_psi ); // debug
         }

         // transform translated beam into reciprocal space for slice loop
         fftw_execute( pf_c2c_psi ); 

         for ( ptrdiff_t i=0; i<local_alloc_size_fftw; ++i)
         {
            psi[i][0] = psi[i][0] / sqrtNxNy;
            psi[i][1] = psi[i][1] / sqrtNxNy;
         }
         // Bandwidth limit the probe again to elliminate any extra
         //  frequencies that might have been introduced by the FFT
         //bw_limit(
         //      psi,
         //      Nx_local, kx_local, Ny, ky,
         //      bwcutoff_t   
         //      );
         if ( input_flag_image_output && input_flag_debug )
         {
            diffraction_scale_factor = 1.0e+0;
            output_diffraction(
                  psi,
                  diffraction_scale_factor,
                  local_alloc_size_fftw,
                  Nx_local, Nx, Ny,
                  resolutionUnit_recip,
                  xResolution_recip, yResolution_recip,
                  outFileName_prefix + "_initial_probe_recip_",
                  mynode, rootnode, comm
                  );
         }

         /////////////////////////////////////////////////////////////
         // Propagate the probe through the slices
         /////////////////////////////////////////////////////////////
         size_t sliceNumber = 0; //debug
         for(
             std::list< slice* >::iterator sliceList_itr
               = slicesOfScatterers.begin();
             sliceList_itr != slicesOfScatterers.end(); 
             ++sliceList_itr )
         {
            ++sliceNumber;//debug

            //cout << " slice " << sliceNumber << endl;//debug

            //////////////////////////////////////////////////////////
            // Multiply psi by the transmission function
            //    t_{n}(x,y) \psi_{n}(x,y) = 
            //       \exp[i \sigma v_{zn}(x,y)] \psi_{n}(x,y) 
            //////////////////////////////////////////////////////////

            // bw_limit psi in reciprocal space before multiplying 
            //  against it in realspace
            //bw_limit(
            //      psi,
            //      Nx_local, kx_local, Ny, ky,
            //      bwcutoff_t   // TODO: this is a problem
            //      );

            // bring psi back into realspace
            fftw_execute( pb_c2c_psi );

            //TODO: make a flag for this output since it's sometimes nice,
            //       but usually a waste of disk space
            //if ( first_probe_flag && input_flag_image_output )  // debug
            //   debug_output_complex_fftw_operand(
            //         psi,
            //         2, 
            //         local_alloc_size_fftw,
            //         Nx_local, Nx, Ny,
            //         resolutionUnit,
            //         xResolution, yResolution,
            //         outFileName_prefix 
            //           + "_probe_" + TEM_NS::to_string(sliceNumber),
            //         mynode, rootnode, comm
            //         );

            if ( mynode == rootnode && input_flag_debug )
            {
               cout << "multiplying psi by "
                  << "transmission function, slice : "
                  << sliceNumber
                  << ", node : " << mynode << endl;
            }

            // limit the probe in realspace, to prevent aliasing artifacts 
            //  from appearing in reciprocal space
            //double probecutoff; probecutoff = 0.25;
            //double probecutoff; probecutoff = 0.125;
            for ( ptrdiff_t i=0; i < Nx_local; i++)
            {
               for ( ptrdiff_t j=0; j < Ny; j++)
               {
                  //if (
                  //      (pow(xx_joined[0] - xx_local[i],2) 
                  //          + pow(yy[0] - yy[j],2)
                  //         > pow(probecutoff * xperiod_duped,2) 
                  //            + pow(probecutoff * yperiod_duped,2)
                  //      )
                  //      &&
                  //      (pow(xx_joined[0] + xperiod_duped - xx_local[i],2) 
                  //          + pow(yy[0] - yy[j],2)
                  //         > pow(probecutoff * xperiod_duped,2) 
                  //            + pow(probecutoff *yperiod_duped,2)
                  //      )
                  //      &&
                  //      (pow(xx_joined[0] - xx_local[i],2) 
                  //          + pow(yy[0] + yperiod_duped - yy[j],2)
                  //         > pow(probecutoff * xperiod_duped,2) 
                  //            + pow(probecutoff * yperiod_duped,2)
                  //      )
                  //      &&
                  //      (pow(xx_joined[0] + xperiod_duped - xx_local[i],2) 
                  //          + pow(yy[0] + yperiod_duped - yy[j],2)
                  //         > pow(probecutoff * xperiod_duped,2) 
                  //            + pow(probecutoff * yperiod_duped,2)
                  //      )
                  //      &&
                  //      (pow(xx_joined[0] - xperiod_duped - xx_local[i],2) 
                  //          + pow(yy[0] - yy[j],2)
                  //         > pow(probecutoff * xperiod_duped,2) 
                  //            + pow(probecutoff * yperiod_duped,2)
                  //      )
                  //      &&
                  //      (pow(xx_joined[0] - xx_local[i],2) 
                  //          + pow(yy[0] - yperiod_duped - yy[j],2)
                  //         > pow(probecutoff * xperiod_duped,2) 
                  //            + pow(probecutoff * yperiod_duped,2)
                  //      )
                  //      &&
                  //      (pow(xx_joined[0] - xperiod_duped - xx_local[i],2) 
                  //          + pow(yy[0] - yperiod_duped - yy[j],2)
                  //         > pow(probecutoff * xperiod_duped,2) 
                  //            + pow(probecutoff * yperiod_duped,2)
                  //      )
                  //      &&
                  //      (pow(xx_joined[0] + xperiod_duped - xx_local[i],2) 
                  //          + pow(yy[0] - yperiod_duped - yy[j],2)
                  //         > pow(probecutoff * xperiod_duped,2) 
                  //            + pow(probecutoff * yperiod_duped,2)
                  //      )
                  //      &&
                  //      (pow(xx_joined[0] - xperiod_duped - xx_local[i],2) 
                  //          + pow(yy[0] + yperiod_duped - yy[j],2)
                  //         > pow(probecutoff * xperiod_duped,2) 
                  //            + pow(probecutoff * yperiod_duped,2)
                  //      )
                  //   )
                  //{
                  //   psi[j + i*Ny][0] = 0.0;
                  //   psi[j + i*Ny][1] = 0.0;
                  //}
                  //else
                  //{

                     // complex multiplication of the transmission 
                     //  function and psi
                     // (a + ib)(c + id) = ac - bd + i(ad + bc)
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
                  //}
               }
            }

            //if ( input_flag_debug )
            //{
            //   cout << "first_probe_flag " << first_probe_flag << endl
            //      << "sliceNumber " << sliceNumber << endl;
            //}

            // scale factor for logarithmic scaling of diffraction images
            //double diffraction_scale_factor;
            // TODO: comment this out or give it its own flag
            //if ( input_flag_image_output && input_flag_debug
            //      && first_probe_flag )  // debug
            //{
            //   debug_output_complex_fftw_operand(
            //         psi,
            //         3, 
            //         local_alloc_size_fftw,
            //         Nx_local, Nx, Ny,
            //         resolutionUnit,
            //         xResolution, yResolution,
            //         outFileName_prefix 
            //            + "_slice_" + TEM_NS::to_string(sliceNumber),
            //         mynode, rootnode, comm
            //         );
            //}

            //////////////////////////////////////////////////////////
            // Transform psi into reciprocal space: 
            //    FT[ t_{n}(x,y) \psi_{n}(x,y) ]
            //////////////////////////////////////////////////////////

            fftw_execute( pf_c2c_psi );

            for ( ptrdiff_t i=0; i < local_alloc_size_fftw; i++)
            {
               psi[i][0] = psi[i][0] /sqrtNxNy;
               psi[i][1] = psi[i][1] /sqrtNxNy;
            }

            //////////////////////////////////////////////////////////
            // (skip) Evaluate P_{n}(k_{x},k_{y},\Delta z)
            //       = \exp[ 
            //         -i \pi \lambda (k^{2}_{x} + k^{2}_{y}) \Delta z_{n}
            //         ]
            //////////////////////////////////////////////////////////
            // NOTE: saving propagators as slice members might only be 
            //       useful to conserve cpu at the expense of memory when 
            //       slice thicknesses and thus propagators are duplicated
            //       or re-used, as in STEM.
            //
            //       Saving the propagator as a slice member would also be
            //       useful if the tem calculation were performed many 
            //       times without deleting the slices.
            //
            // TODO: Considering that reverse monte carlo only alters the 
            //        position of one atom at a time, is there any way to 
            //        preserve all slice data in memory between tem() 
            //        calls and only update those which change upon atom 
            //        movement?
            //       Is there a way to compare slices to see if they 
            //        differ?
            //
            //evaluate_propagator(
            //      lambda,
            //      kx_local, Nx_local, ky, Ny,
            //      (*sliceList_itr)->thickness,
            //      (*sliceList_itr)->propagator
            //      );
  
            //////////////////////////////////////////////////////////
            // Evaluate P_{n}(k_{x},k_{y},\Delta z) 
            //             * FT[t_{n}(x,y) \psi_{n}(x,y)]
            //////////////////////////////////////////////////////////

            //if ( mynode == rootnode && input_flag_debug )//debug
            //{
            //   cout << "propagating, slice : " // debug
            //      << sliceNumber 
            //      //<< ", node : " << mynode 
            //      << endl;// debug
            //}

            bw_limit(
                  psi,
                  Nx_local, kx_local, Ny, ky,
                  bwcutoff_t   
                  );

            // TODO: double check that propagate() bandwidth limits psi
            // // .... it appears to.

            (*sliceList_itr)->propagate(
                  lambda,
                  kx_local, Nx_local, ky, Ny, 
                  bwcutoff_t,
                  psi
                  );

            // debug
//            fftw_execute( pb_c2c_psi );
//            for ( ptrdiff_t i=0; i < local_alloc_size_fftw; i++)
//            {
//               psi[i][0] = psi[i][0] /sqrtNxNy;
//               psi[i][1] = psi[i][1] /sqrtNxNy;
//            }
//
//            //if ( first_probe_flag && input_flag_debug)
//            //{
//               //if( output_psi_realspace_to_netcdf(
//               //         psi,
//               //         local_alloc_size_fftw,
//               //         Nx_local, kx_joined, Nx, ky, Ny,
//               //         outFileName_prefix 
//               //            + "_psi_slice_" + to_string(sliceNumber),
//               //         mynode, rootnode, comm
//               //      ) == EXIT_FAILURE )
//               //  cout << "output_psi_realspace_to_netcdf() failed" << endl;
//            //}
//
//            fftw_execute( pf_c2c_psi );
//
//            for ( ptrdiff_t i=0; i < local_alloc_size_fftw; i++)
//            {
//               psi[i][0] = psi[i][0] /sqrtNxNy;
//               psi[i][1] = psi[i][1] /sqrtNxNy;
//            }

            // end debug

            //debug
            //fftw_execute( pb_c2c_psi );
            //debug_output_complex_fftw_operand(
            //      psi,
            //      4, 
            //      local_alloc_size_fftw,
            //      Nx_local, Nx, Ny,
            //      outFileName_prefix 
            //         + "_slice_" + TEM_NS::to_string(sliceNumber),
            //      mynode, rootnode, comm
            //      );
            //fftw_execute( pf_c2c_psi );
            //for ( ptrdiff_t i=0; i < Nx_local; i++)
            //   for ( ptrdiff_t j=0; j < Ny; j++)
            //   {
            //      //if (  
            //      // ( pow( kx_local[j + i*Ny], 2) + pow(ky[j + i*Ny], 2) )
            //      //      < bwcutoff_t_sqr 
            //      //)// TODO: should the output wave function be bw 
            //         //          limited 
            //          //       as it is here???
            //          //       No. It should be bw limited, but not here. Do
            //          //       it at the top of the loop so that the output 
            //          //       has the higher frequencies for the last 
            //          //       slice.
            //          //
            //          //       If you take the approach in Kirkland (2010),
            //          //       then you'd need 
            //          //         bwcutoff_t + bwcutoff_propagator = k_max
            //          //       and no cutoff for psi, limiting frequency 
            //          //       contributions to 1/2 k_max.
            //          //       Yet here we can have 2/3 k_max for all 
            //          //       transmission function, propagator, and psi.
            //          //
            //          //       The condition to remove aliasing is:
            //          //     bwcutoff_psi + bwcutoff_t + bwcutoff_propagator
            //          //        < 2 k_max
            //      //{
            //         psi[i*Ny + j][0] = psi[i*Ny + j][0] / NxNy; 
            //         psi[i*Ny + j][1] = psi[i*Ny + j][1] / NxNy;
            //      //}
            //      //else
            //      //{
            //      //   psi[i*Ny + j][0] = 0.0; 
            //      //   psi[i*Ny + j][1] = 0.0; 
            //      //}
            //   }
            //end debug remove the above normalization when removing the debug tif

         } // end of iteration over slices

         if ( mynode == rootnode && input_flag_debug )
            cout << "finished propagating through slices,"
              << " probe centered at : (x,y) = ( " 
              << *x_itr << ", " << *y_itr << ")" << endl;

         if ( first_probe_flag 
               && (input_flag_image_output || input_flag_netcdf_images) )
         {
            if ( mynode == rootnode && input_flag_debug 
                  && input_flag_netcdf_images )
            {
               cout << "writing first diffraction to netcdf" << endl;
               if( output_psi_reciprocalspace_to_netcdf(
                     psi,
                     local_alloc_size_fftw,
                     Nx_local, kx_joined, Nx, ky, Ny,
                     outFileName_prefix 
                        + "_first_diffraction_" 
                        + TEM_NS::to_string(*x_itr) 
                        + "_" + TEM_NS::to_string(*y_itr),
                     mynode, rootnode, comm
                     ) == EXIT_FAILURE )
               {
                  cout << "output_psi_realspace_to_netcdf() failed" 
                     << endl;
               }
            }

            // scale factor for logarithmic scaling of diffraction images
            
            diffraction_scale_factor = 1.0e-5;
            output_diffraction(
                  psi,
                  diffraction_scale_factor,
                  local_alloc_size_fftw,
                  Nx_local, Nx, Ny,
                  resolutionUnit_recip,
                  xResolution_recip, yResolution_recip,
                  outFileName_prefix + "_first_diffraction_",
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
                  outFileName_prefix + "_first_diffraction_",
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
                  outFileName_prefix + "_first_diffraction_",
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
                  outFileName_prefix + "_first_diffraction_",
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
                  outFileName_prefix + "_first_diffraction_",
                  mynode, rootnode, comm
                  );
               }

         first_probe_flag = false;

         //-----------------------------------------------------------
         // End of beam transmission through individual slices
         //-----------------------------------------------------------

         /////////////////////////////////////////////////////////////
         // FTEM: accumulate the first and second moments of 
         //       diffracted intensity over real-space
         //       \langle I(|\vec{k}|, K_{aperture}) \rangle
         //       \langle I^{2}(|\vec{k}|, K_{aperture}) \rangle
         /////////////////////////////////////////////////////////////
         //  For 1-D variance:
         //   I^{2}(|\vec{k}|) is held on rootnode, while 
         //   I(|\vec{k}|) is still split among the nodes and needs to
         //   be reduced during final variance calculation.

         if ( input_flag_fem )
         {
            // Accumulate values for 2-D variance image.
            for ( ptrdiff_t ii=0; ii < local_alloc_size_fftw; ++ii )
            { 
               diffracted_wave_mag_sum[ii] 
                  += 
                  sqrt(psi[ii][0] * psi[ii][0] 
                        + psi[ii][1] * psi[ii][1]);

               diffracted_wave_mag_sqr_sum[ii] 
                  += psi[ii][0] * psi[ii][0] 
                     + psi[ii][1] * psi[ii][1];
                  //(psi[ii][0] * psi[ii][0] 
                  //   + psi[ii][1] * psi[ii][1])
                  //   * (psi[ii][0] * psi[ii][0] 
                  //   + psi[ii][1] * psi[ii][1]);
               
               if (input_flag_complex_realspace_sum) 
               {
                  diffracted_wave_re_sum[ii] += psi[ii][0] ;
                  diffracted_wave_im_sum[ii] += psi[ii][1] ;
               }
            }

            // Azimuthally integrate diffraction for 1-D variance.
            integrate_out_theta_fftw(
                  psi,
                  kx_local, Nx_local,
                  ky, Ny,
                  binning_boundaries, 
                  bin_counts_local,
                  diffracted_wave_radial_intensity_sum_local
                  );

            // TODO: is the following initialization necessary?
            // No.
            //if ( mynode == rootnode )
            //   for (unsigned int ii=0; ii < number_of_bins; ++ii)
            //      diffracted_wave_radial_intensity_sum_tmp[ii] =0.0;
            if ( number_of_bins != binning_boundaries.size() - 1)//debug
            {//debug
               cerr << "Error - number_of_bins != binning_boundaries.size() - 1" << endl;//debug
               // TODO: make a function to release memory for a clean exit
            }//debug

            MPI_Reduce(
                  diffracted_wave_radial_intensity_sum_local, // local
                  diffracted_wave_radial_intensity_sum_tmp, //rootnode
                  number_of_bins,
                  MPI_DOUBLE, MPI_SUM,
                  rootnode, comm
                  );

            MPI_Reduce(
                  bin_counts_local, // local
                  bin_counts_aggregated, //rootnode
                  number_of_bins,
                  MPI_INT, MPI_SUM,
                  rootnode, comm
                  );

            if ( mynode == rootnode )
            {
               for ( unsigned int i=0; i<number_of_bins; ++i)
               {
                  // division by bin_counts_aggregated implements 
                  //  azimuthal averaging 
                  diffracted_wave_radial_intensity_sum[i]
                     += diffracted_wave_radial_intensity_sum_tmp[i] 
                        / bin_counts_aggregated[i];

                  diffracted_wave_radial_intensity_sqr_sum[i]
                     += (diffracted_wave_radial_intensity_sum_tmp[i]
                        * diffracted_wave_radial_intensity_sum_tmp[i])
                           / (bin_counts_aggregated[i] 
                                 * bin_counts_aggregated[i]);
               }
            }
         } // end of input_flag_fem dependent block

         /////////////////////////////////////////////////////////////
         // STEM: integrate the diffracted intensity impingent upon 
         //       the annular detector
         /////////////////////////////////////////////////////////////

         if ( input_flag_image_output )
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
               stem_image[ pixel_number_y + pixel_number_x * y_p.size() ]
                  = detected_intensity_reduced;
            }
         }

         fftw_destroy_plan( pb_c2c_psi );
         fftw_destroy_plan( pf_c2c_psi );

         ++pixel_number_y; 
      }
      ++pixel_number_x; 
   }

   //=================================================================
   // End of STEM beam position rastering
   //=================================================================
   if ( input_flag_image_output )
   {
      //output_diffraction_with_renormalization(
      //diffraction_scale_factor = 1.0e13;//1e-20;
      diffraction_scale_factor = 1.0e+40;//1e-20;
      output_diffraction(
            psi,
            diffraction_scale_factor,
            local_alloc_size_fftw,
            Nx_local, Nx, Ny,
            resolutionUnit_recip,
            xResolution_recip, yResolution_recip,
            outFileName_prefix + "_final_0_",
            mynode, rootnode, comm
            );
   }



   ////////////////////////////////////////////////////////////////
   // Normalize psi by a factor of ??? remaining from fftw
   ////////////////////////////////////////////////////////////////
   //for ( ptrdiff_t i=0; i < local_alloc_size_fftw; ++i )
   //{// TODO: is this the correct factor?
   //   psi[i][0] = psi[i][0] / NxNy;
   //   psi[i][1] = psi[i][1] / NxNy;
   //}

   double* diffracted_wave_mag_variance; 
   fftw_complex* diffracted_wave_fftw_complex_sum;// only used for 1 image
   double* diffracted_wave_complex_sum_mag;
   // TODO: modify supporting functions so that the above variable can be 
   //       removed.
   double diffracted_wave_mag_sqr_sum_max = 0.0;
   double diffracted_wave_mag_sum_max = 0.0;
   if ( input_flag_fem )
   {// fem specific commands
      //double diffracted_wave_mag_sqr_sum_max = 0.0;
      //double diffracted_wave_mag_sum_max = 0.0;
      // TODO: put this division elsewhere or eliminate when 
      //        unnecessary

      ///////////////////////////////////////////////////////////////
      // Divide the sum by the number of points to obtain the average
      // 2-D diffracted intensity.
      ///////////////////////////////////////////////////////////////

      for ( ptrdiff_t i=0; i < local_alloc_size_fftw; ++i )
      {
         diffracted_wave_mag_sqr_sum[i] 
            = diffracted_wave_mag_sqr_sum[i] / number_of_raster_points;

         diffracted_wave_mag_sum[i] 
            = diffracted_wave_mag_sum[i] / number_of_raster_points; 

         if (input_flag_complex_realspace_sum)
         {
            diffracted_wave_re_sum[i] 
               = diffracted_wave_re_sum[i] / number_of_raster_points;

            diffracted_wave_im_sum[i] 
               = diffracted_wave_im_sum[i] / number_of_raster_points;
         }

         if (diffracted_wave_mag_sum_max < diffracted_wave_mag_sum[i])
            diffracted_wave_mag_sum_max =  diffracted_wave_mag_sum[i];

         if (
               diffracted_wave_mag_sqr_sum_max 
               < diffracted_wave_mag_sqr_sum[i] 
            )
            diffracted_wave_mag_sqr_sum_max 
               =  diffracted_wave_mag_sqr_sum[i];
         // end debug 
      }
      if( input_flag_debug )
      {
         cout << "node " << mynode 
            << ", E[X^2] diffracted_wave_mag_sqr_sum_max : "
            << diffracted_wave_mag_sqr_sum_max << endl; // debug
         cout << "node " << mynode 
            << ", (E[X])^2 diffracted_wave_mag_sum_max : "
            << (diffracted_wave_mag_sum_max 
                  * diffracted_wave_mag_sum_max )
                  << endl; // debug
         cout << "node " << mynode << ",  E[X^2] - (E[X])^2 : "
            << diffracted_wave_mag_sqr_sum_max -
            (diffracted_wave_mag_sum_max * diffracted_wave_mag_sum_max ) 
            << endl; // debug
      }

      ////////////////////////////////////////////////////////////////
      // Combine re_sum and im_sum to obtain the magnitudes
      //  for 2-D variance images
      ////////////////////////////////////////////////////////////////

      // TODO: modify to obtain 1-D variance curve using the complex sum
      if (input_flag_complex_realspace_sum) 
      {
         diffracted_wave_complex_sum_mag 
            = new double[ local_alloc_size_fftw ];

         diffracted_wave_fftw_complex_sum 
            = fftw_alloc_complex( local_alloc_size_fftw );

         for ( size_t i=0; i < local_alloc_size_fftw; ++i )
         {
            diffracted_wave_complex_sum_mag[i] 
               = sqrt(
                     diffracted_wave_re_sum[i] * diffracted_wave_re_sum[i]
                   + diffracted_wave_im_sum[i] * diffracted_wave_im_sum[i]
                 );

            //cout << "diffracted_wave_im_sum[" << i << "] : " // debug
            //   << diffracted_wave_im_sum[i] << endl; // debug
            //cout << "diffracted_wave_complex_sum[" << i << "] : " // debug
            //   << diffracted_wave_complex_sum[i] << endl; // debug

            // Store the sums as fftw_complex variables for 
            //    input into output_psi_reciprocalspace_to_netcdf
            // TODO : check to see if this is necessary, or even useful
            // TODO : do these need to be divided by number_of_raster_points?
            diffracted_wave_fftw_complex_sum[i][0] 
               = diffracted_wave_re_sum[i];
            diffracted_wave_fftw_complex_sum[i][1] 
               = diffracted_wave_im_sum[i];

            // TODO: remove this averaging when performed by 
            //        variance_2D_STEM ... WHY?
            diffracted_wave_complex_sum_mag[i]
               = diffracted_wave_complex_sum_mag[i] 
                  / number_of_raster_points;
         }
      }

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
      
      // TODO: covariance 


      ////////////////////////////////////////////////////////////////
      // Save average of 2-D diffraction magnitude taken from all 
      //  raster points
      ////////////////////////////////////////////////////////////////
      if ( input_flag_image_output )
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
               outFileName_prefix + "_diffracted_wave_avg",
               mynode, rootnode, comm
               );

         diffraction_scale_factor = 1.0e13 ;//1e-53;
         output_diffraction(
               diffracted_wave_mag_sqr_sum,
               diffraction_scale_factor,
               local_alloc_size_fftw,
               Nx_local, Nx, Ny,
               resolutionUnit_recip,
               xResolution_recip, yResolution_recip,
               outFileName_prefix + "_diffracted_wave_sqr_avg",
               mynode, rootnode, comm
               );
      }

      if (input_flag_image_output && input_flag_complex_realspace_sum) 
      {
         diffraction_scale_factor = 1.0e13 ;//1e-53;
         output_diffraction(
               diffracted_wave_re_sum,
               diffraction_scale_factor,
               local_alloc_size_fftw,
               Nx_local, Nx, Ny,
               resolutionUnit_recip,
               xResolution_recip, yResolution_recip,
               outFileName_prefix + "_diffracted_wave_re_sum",
               mynode, rootnode, comm
               );


         output_diffraction(
               diffracted_wave_im_sum,
               diffraction_scale_factor,
               local_alloc_size_fftw,
               Nx_local, Nx, Ny,
               resolutionUnit_recip,
               xResolution_recip, yResolution_recip,
               outFileName_prefix + "_diffracted_wave_im_sum",
               mynode, rootnode, comm
               );

         output_diffraction(
               diffracted_wave_complex_sum_mag,
               diffraction_scale_factor,
               local_alloc_size_fftw,
               Nx_local, Nx, Ny,
               resolutionUnit_recip,
               xResolution_recip, yResolution_recip,
               outFileName_prefix + "_diffracted_wave_complex_sum_mag_avg",
               mynode, rootnode, comm
               );

      }
      if ( input_flag_complex_realspace_sum && input_flag_netcdf_images)
        if(
           output_psi_reciprocalspace_to_netcdf(
                 diffracted_wave_fftw_complex_sum,
                 local_alloc_size_fftw,
                 Nx_local,
                 kx_joined, Nx,
                 ky, Ny,
                 outFileName_prefix 
                    + "_diffracted_wave_fftw_complex_sum",
                 mynode, rootnode, comm
                 ) == EXIT_FAILURE
        )
        {
           cout << "output_psi_reciprocalspace_to_netcdf() failed" 
              << endl;
            fftw_free( diffracted_wave_fftw_complex_sum ); // debug

         }

      ////////////////////////////////////////////////////////////////
      // create an image of the variance in 2D reciprocoal space 
      ////////////////////////////////////////////////////////////////
      //diffracted_wave_mag_variance = new double[ local_alloc_size_fftw ];
      //for ( ptrdiff_t i=0; i<local_alloc_size_fftw; ++i )
      //{
      //   if (
      //         ( diffracted_wave_mag_sum[i] != 0.0 )
      //         //&& ( diffracted_wave_mag_sum != inf )
      //      )
      //   {
      //      diffracted_wave_mag_variance[i]
      //       = (
      //           diffracted_wave_mag_sqr_sum[i] 
      //           -
      //          (diffracted_wave_mag_sum[i] * diffracted_wave_mag_sum[i])
      //         )
      //         /
      //         (diffracted_wave_mag_sum[i] * diffracted_wave_mag_sum[i]);
      //      //diffracted_wave_mag_variance[i]
      //      //   = (
      //      //      diffracted_wave_mag_sqr_sum[i] 
      //      //      /
      //      //     ( diffracted_wave_mag_sum[i] 
      //      //        * diffracted_wave_mag_sum[i] )
      //      //   )
      //      //- 1.0;
      //   }
      //   else
      //      diffracted_wave_mag_variance[i] = 0.0;
      //}
      //cout << "writing image of variance" << endl; // debug
      //output_diffraction(
      //      diffracted_wave_mag_variance,
      //      1.0e0,
      //      //diffraction_scale_factor,
      //      local_alloc_size_fftw,
      //      Nx_local, Nx, Ny,
      //      outFileName_prefix + "_diffracted_wave_mag_variance",
      //      mynode, rootnode, comm
      //      );

      //if(
      //   output_psi_mag_reciprocalspace_to_netcdf(
      //      diffracted_wave_mag_variance,
      //      local_alloc_size_fftw,
      //      Nx_local,
      //      kx_joined, Nx,
      //      ky, Ny,
      //      outFileName_prefix + "_diffracted_wave_mag_variance",
      //      mynode, rootnode, comm
      //      ) == EXIT_FAILURE
      //)
      //   cout << "output_psi_reciprocalspace_to_netcdf() failed" << endl;
      //delete[] diffracted_wave_mag_variance ;

      ////////////////////////////////////////////////////////////////
      // FTEM : Calculate the variance of the diffracted intensities 
      //          and save it to a netCDF format file
      ////////////////////////////////////////////////////////////////
      //cout << "variance calculation, node " << mynode // debug
      //   << ", kx_local[0], kx_local[1], ky[0], ky[1], delta_k : " // debug
      //   << kx_local[0] << ", " << kx_local[1] << ", " // debug
      //   << ky[0] << ", " << ky[1] << ", " // debug
      //   << delta_k << endl; // debug

      if ( variance_2D_STEM(
            diffracted_wave_mag_sum,         // denominator
            diffracted_wave_mag_sqr_sum,     // numerator
            number_of_raster_points,
            //binning_boundaries,
            kx_local, Nx_local,
            kx_joined, 
            Nx,
            ky, Ny,
            resolutionUnit_recip,
            xResolution_recip, yResolution_recip,
            delta_k, // binning separation in number of pixels
            bwcutoff_t,// data was set to 0 at magnitudes > bwcutoff_t
            outFileName_prefix + "_mag_sqr_sum_over_mag_sum",
            mynode, rootnode, comm
            ) == EXIT_FAILURE )
         cout << "variance_2D_STEM() failed" << endl;

      if (input_flag_complex_realspace_sum) 
      {
         if( variance_2D_STEM(
               diffracted_wave_complex_sum_mag,    // denominator
               diffracted_wave_mag_sqr_sum,        // numerator
               number_of_raster_points,
               //binning_boundaries,
               kx_local, Nx_local,
               kx_joined, 
               Nx,
               ky, Ny,
               resolutionUnit_recip,
               xResolution_recip, yResolution_recip,
               delta_k, // binning separation in number of pixels
               bwcutoff_t, //data was set to 0 at magnitudes > bwcutoff_t
               outFileName_prefix + "_variance_using_complex_sum",
               mynode, rootnode, comm
               ) == EXIT_FAILURE)
            cout << "variance_2D_STEM() failed" << endl;
      }

      if ( input_flag_debug && mynode == rootnode )
         cout << "finished variance_2D_STEM() calls" << endl; // debug 

      ////////////////////////////////////////////////////////////////
      // Calculate the FTEM 1-D variance V(|\vec{k}|,K_{aperture})
      ////////////////////////////////////////////////////////////////
      // V(|\vec{k}|,K_{aperture})
      //  = \frac{ \langle I^{2}(|\vec{k}|,K_{aperture}) \rangle }
      //          { \langle I(|\vec{k}|,K_{aperture}) \rangle^{2} } -1
      // where I() and I^{2}() are 
      // diffracted_wave_radial_intensity_sum and 
      // diffracted_wave_radial_intensity_sqr_sum, respectively.

      if ( mynode == rootnode )
      {
         ftem_variance = new double[ number_of_bins ];
         if ( 
               variance_1D_STEM(
                  diffracted_wave_radial_intensity_sum,
                  diffracted_wave_radial_intensity_sqr_sum,
                  binning_boundaries,
                  number_of_raster_points,
                  input_flag_netcdf_variance,
                  outFileName_prefix,
                  ftem_variance
                  )
               == EXIT_FAILURE )
            cout << "variance_1D_STEM() failed" << endl;
      }

   }// end fem specific commands

   ////////////////////////////////////////////////////////////////
   // Save stem_image as a tif
   ////////////////////////////////////////////////////////////////
   if ( mynode == rootnode && input_flag_image_output )
   {
      output_stem_image(
            stem_image,
            x_p.size(), y_p.size(),
            resolutionUnit,
            xResolution_stem, yResolution_stem,
            outFileName_prefix
            //mynode, rootnode, comm
            );
   }

   ////////////////////////////////////////////////////////////////
   // Save diffraction as tif
   ////////////////////////////////////////////////////////////////
   
   if ( input_flag_image_output )
   {
      //output_diffraction_with_renormalization(
      //diffraction_scale_factor = 1.0e13;//1e-20;
      diffraction_scale_factor = 1.0e+40;//1e-20;
      output_diffraction(
            psi,
            diffraction_scale_factor,
            local_alloc_size_fftw,
            Nx_local, Nx, Ny,
            resolutionUnit_recip,
            xResolution_recip, yResolution_recip,
            outFileName_prefix + "_final_1",
            mynode, rootnode, comm
            );

      //debug_output_complex_fftw_operand(
      //      psi,
      //      5, 
      //      local_alloc_size_fftw,
      //      Nx_local, Nx, Ny,
      //      resolutionUnit,
      //      xResolution, yResolution,
      //      outFileName_prefix + "_diffraction_" 
      //         + TEM_NS::to_string(sliceNumber),
      //      mynode, rootnode, comm
      //      );
   }


   /*
   ///////////////////////////////////////////////////////////////////
   // CTEM: apply the objective lens transfer function
   ///////////////////////////////////////////////////////////////////
   
   double H_0_re;
   double H_0_im;

   for ( ptrdiff_t i=0; i < Nx_local; i++)
      for ( ptrdiff_t j=0; j < Ny; j++)
      {
         //objective_lens_transfer_function_coher(
         objective_lens_transfer_function_wp(
               pow(kx_local[i], 2) + pow(ky[j], 2),
               Cs3, defocus, alpha_max_sqr,
               lambda, lambda_sqr,
               H_0_re, H_0_im
               );
         //H_0_mag = sqrt( pow( H_0_re, 2) + pow( H_0_im, 2) );
         //multiply the working image by H_0 (complex number multiplication)
         psi[i*Ny + j][0]
            = psi[i*Ny + j][0] * H_0_re - psi[i*Ny + j][1] * H_0_im;

         psi[i*Ny + j][1]
            = psi[i*Ny + j][0] * H_0_im + psi[i*Ny + j][1] * H_0_re;
      }


   // Execute the reverse fftw plan
   fftw_execute( pb_c2c_psi );

   // psi needs normalization by NxNy at this point
   
   debug_output_complex_fftw_operand(
         psi,
         5, 
         local_alloc_size_fftw,
         Nx_local, Nx, Ny,
         resolutionUnit,
         xResolution, yResolution,
         outFileName_prefix + "_slice_" + TEM_NS::to_string(sliceNumber),
         mynode, rootnode, comm
         );
 
   // using Kirkland (2016) eq 11 for image intensity 
   for ( ptrdiff_t i=0; i < local_alloc_size_fftw; i++)
   {// Normalize the output of the forward then reverse Fourier transforms
      psi[i][0]
         = 1 + 2 * psi[i][0] / NxNy;
      psi[i][1]
         = 1 + 2 * psi[i][1] / NxNy; 
   }

   debug_output_complex_fftw_operand(
         psi,
         6, 
         local_alloc_size_fftw,
         Nx_local, Nx, Ny,
         resolutionUnit,
         xResolution, yResolution,
         outFileName_prefix + "_slice_" + TEM_NS::to_string(sliceNumber),
         mynode, rootnode, comm
         );

   */
   
   // debug output
   //if ( mynode == rootnode )
   //{
   //   cout << "VV : " << VV << " [V]" << endl; // debug
   //   cout << "m_e0 : " << m_e0 << " [kg]" << endl; // debug
   //   cout << "cc : " << cc << " [m/s]" << endl; // debug
   //   cout << "ee : " << ee << " [C]" << endl; // debug
   //   cout << "hh : " << hh << " [m^{2} kg/s]" << endl; // debug
   //   cout << "lambda : " << lambda << " [AA]" << endl; // debug
   //   cout << "lambda_inv : " << lambda_inv << " [AA^{-1}]" << endl;//debug
   //   cout << "sigma : " << sigma << " [1/(V AA)]" << endl; // debug
   //   cout << "gamma : " << gamma << endl; // debug
   //   cout << "Cs3 : " << Cs3 << endl; // debug

   //   cout << endl;
   //   cout << "Optimal parameters without aberration compensation "
   //      << " with fixed Cs3 (Kirkland 2016):" << endl;
   //   cout << "defocus : " << scherzer_defocus_adfstem_uncorrected(
   //   //cout << "defocus : " << scherzer_defocus_bfctem_uncorrected(
   //         lambda, Cs3) << endl;
   //   cout << "alpha_max : " << scherzer_alphamax_adfstem_uncorrected(
   //   //cout << "alpha_max : " << scherzer_alphamax_bfctem_uncorrected(
   //         lambda, Cs3) << endl;
   //   cout << "approx. min. resolution : " 
   //      << scherzer_approxresolution_adfstem_uncorrected(
   //      //<< scherzer_approxresolution_bfctem_uncorrected(
   //         lambda, Cs3) << endl;
   //}
   // end debug output
   
   ///////////////////////////////////////////////////////////////////
   // clean up
   ///////////////////////////////////////////////////////////////////

   //if ( mynode == rootnode )
   //   if ( ! flag_wisdomFile )
   //   {
   //      cout << "Exporting wisdom to file: " << wisdom_filename << endl;
   //      if ( ! fftw_export_wisdom_to_filename( wisdom_filename.c_str()))
   //      {
   //         cerr << "Error - could not export wisdom to file: "
   //           << wisdom_filename << endl;
   //      }
   //   }

   if ( mynode == rootnode )
   {
      //delete[] kx_joined;
      delete[] stem_image;
      //delete[] x_p;
      //delete[] y_p;
   }

   delete stem_probe_joined_re, stem_probe_joined_im;

   //fftw_destroy_plan( pb_c2c_psi );
   //fftw_destroy_plan( pf_c2c_psi );

   
   fftw_free( psi );

   //if ( input_flag_fem )
   //{
   //   fftw_free( diffracted_wave_fftw_complex_sum ); // debug
   //}
   
   //fftw_mpi_cleanup();

   //delete[] ky;
   //delete[] kx_local;

   if ( input_flag_fem )
   {
      // 2-D variance
      delete[] diffracted_wave_mag_sum;
      delete[] diffracted_wave_mag_sqr_sum;
      
      if (input_flag_complex_realspace_sum) 
      {
         delete[] diffracted_wave_complex_sum_mag; // debug

         delete[] diffracted_wave_re_sum;
         delete[] diffracted_wave_im_sum;
      }

      // 1-D variance
      delete[] bin_counts_local;
      if ( mynode == rootnode )
      {
         delete[] bin_counts_aggregated;
         delete[] diffracted_wave_radial_intensity_sqr_sum;
         delete[] diffracted_wave_radial_intensity_sum;
         delete[] diffracted_wave_radial_intensity_sum_tmp;

         delete[] ftem_variance;
      }

      delete[] diffracted_wave_radial_intensity_sum_local;
   }

   //delete_slices_from_list_stem( slicesOfScatterers );
     // ^^ calls fftw_free() on paps, exp_i_sigma_v, and deletes 
     //  propagator_x_re, propagator_x_im, propagator_y_re,
     //  and propagator_y_im
   
   return EXIT_SUCCESS;
}

#endif
