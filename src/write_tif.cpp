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
// File: write_tif.cpp
// Purpose: to define tools for writing data to 8, 16, & 32 bit tif images 
//    using libtiff

#ifndef WRITE_TIF_CPP
#define WRITE_TIF_CPP

#include "write_tif.hpp"

////////////////////////////////////////////////////////////////////////////
// public member functions /////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////

//-----------------------------------------------------------------------
int TEM_NS::writeImage::write_tif_grayscale8( 
      const double* const outDataRaw, 
      const double& min, const double& max,// bounds of scale for input data
      const size_t& height, const size_t& width,// image dimensions [pixels]
      const size_t& resolutionUnit,
      const double& xResolution, const double& yResolution,
      const string& outTifPrefix )
{
   /////////////////////////////////////////////////////////////////////////
   // write the image
   /////////////////////////////////////////////////////////////////////////
   if ( write_tif_grayscale<unsigned char>(
            outDataRaw,
            min, max,
            height, width, 
            resolutionUnit, 
            xResolution, yResolution,
            outTifPrefix,
            false )  == EXIT_FAILURE)
   {
      cerr << "Error : write_tif_grayscale8() failed" << endl;
      return EXIT_FAILURE;
   }
   return EXIT_SUCCESS;
} // write_tif_grayscale8

//-----------------------------------------------------------------------
int TEM_NS::writeImage::write_tif_grayscale16(
      const double* const outDataRaw, 
      const double& min, const double& max,// bounds of scale for input data
      const size_t& height, const size_t& width,// image dimensions [pixels]
      const size_t& resolutionUnit,
      const double& xResolution, const double& yResolution,
      const string& outTifPrefix )
{
   /////////////////////////////////////////////////////////////////////////
   // write the image
   /////////////////////////////////////////////////////////////////////////
   if ( write_tif_grayscale<uint16_t>(
            outDataRaw,
            min, max,
            height, width, 
            resolutionUnit,
            xResolution, yResolution,
            outTifPrefix,
            false )  == EXIT_FAILURE)
   {
      cerr << "Error : write_tif_grayscale16() failed" << endl;
      return EXIT_FAILURE;
   }
   return EXIT_SUCCESS;
} // write_tif_grayscale16

//-----------------------------------------------------------------------
int TEM_NS::writeImage::write_tif_grayscale32( 
      const double* const outDataRaw, 
      const double& min, const double& max,// bounds of scale for input data
      const size_t& height, const size_t& width,// image dimensions [pixels]
      const size_t& resolutionUnit,
      const double& xResolution, const double& yResolution,
      const string& outTifPrefix )
{
   /////////////////////////////////////////////////////////////////////////
   // write the image
   /////////////////////////////////////////////////////////////////////////
   if ( write_tif_grayscale<uint32_t>( 
            outDataRaw,
            min, max,
            height, width, 
            resolutionUnit,
            xResolution, yResolution,
            outTifPrefix, 
            false )  == EXIT_FAILURE)
   {
      cerr << "Error : write_tif_grayscale32() failed" << endl;
      return EXIT_FAILURE;
   }
   return EXIT_SUCCESS;
} // write_tif_grayscale32


//-----------------------------------------------------------------------
int TEM_NS::writeImage::write_tif_grayscale8_logscale( 
      const double* const outDataRaw, 
      const double& min, const double& max,// bounds of scale for input data
      const size_t& height, const size_t& width,// image dimensions [pixels]
      const size_t& resolutionUnit,
      const double& xResolution, const double& yResolution,
      const string& outTifPrefix )
{
   /////////////////////////////////////////////////////////////////////////
   // write the image
   /////////////////////////////////////////////////////////////////////////
   if ( write_tif_grayscale<unsigned char>(
            outDataRaw,
            min, max,
            height, width, 
            resolutionUnit,
            xResolution, yResolution,
            outTifPrefix,
            true )  == EXIT_FAILURE)
   {
      cerr << "Error : write_tif_grayscale8_logscale() failed" << endl;
      return EXIT_FAILURE;
   }
   return EXIT_SUCCESS;
} // write_tif_grayscale8

//-----------------------------------------------------------------------
int TEM_NS::writeImage::write_tif_grayscale16_logscale(
      const double* const outDataRaw, 
      const double& min, const double& max,// bounds of scale for input data
      const size_t& height, const size_t& width,// image dimensions [pixels]
      const size_t& resolutionUnit,
      const double& xResolution, const double& yResolution,
      const string& outTifPrefix )
{
   /////////////////////////////////////////////////////////////////////////
   // write the image
   /////////////////////////////////////////////////////////////////////////
   if ( write_tif_grayscale<uint16_t>(
            outDataRaw,
            min, max,
            height, width, 
            resolutionUnit,
            xResolution, yResolution,
            outTifPrefix,
            true )  == EXIT_FAILURE)
   {
      cerr << "Error : write_tif_grayscale16_logscale() failed" << endl;
      return EXIT_FAILURE;
   }
   return EXIT_SUCCESS;
} // write_tif_grayscale16

//-----------------------------------------------------------------------
int TEM_NS::writeImage::write_tif_grayscale32_logscale( 
      const double* const outDataRaw, 
      const double& min, const double& max,// bounds of scale for input data
      const size_t& height, const size_t& width,// image dimensions [pixels]
      const size_t& resolutionUnit,
      const double& xResolution, const double& yResolution,
      const string& outTifPrefix)
{
   /////////////////////////////////////////////////////////////////////////
   // write the image
   /////////////////////////////////////////////////////////////////////////
   if ( write_tif_grayscale<uint32_t>( 
            outDataRaw,
            min, max,
            height, width, 
            resolutionUnit,
            xResolution, yResolution,
            outTifPrefix, 
            true )  == EXIT_FAILURE)
   {
      cerr << "Error : write_tif_grayscale32_logscale() failed" << endl;
      return EXIT_FAILURE;
   }
   return EXIT_SUCCESS;
} // write_tif_grayscale32

//-----------------------------------------------------------------------
int TEM_NS::output_diffraction_with_renormalization(
      const fftw_complex* const psi,
      const double& scale_factor,
      const ptrdiff_t& local_alloc_size_fftw,
      const ptrdiff_t& Nx_local,
      const ptrdiff_t& Nx,
      const ptrdiff_t& Ny,
      const size_t& resolutionUnit,
      const double& xResolution, const double& yResolution,
      const string& outFileName_prefix,
      const int& mynode,
      const int& rootnode,
      MPI_Comm comm
      )
{
   // Precondition:
   // - psi represents the data in reciprocal space, requiring rearrangement
   //    and also requiring division by sqrt(Nx * Ny) due to fftw.

   writeImage imageWriter;

   double* psi_mag;
   psi_mag = new double[ local_alloc_size_fftw ];

   const double NxNy = Nx * Ny ;
   const double NxNyNxNy = NxNy * NxNy ;

   double max_psi_mag;
   double min_psi_mag;
   double max_psi_mag_joined;
   double min_psi_mag_joined;

   // Calculate intensities, scaling them onto a log scale, and determine 
   //  the max and min of the resulting log scale intensities.
   // Kirkland (2010) eq 4.28:
   // D(k_{x}, k_{y}) = log( 1 + c | F(k_{x}, k_{y}) |^{2} ) 

   //cout << "scale_factor : " << setprecision(20) << scale_factor << endl;//debug
   // set the center of the direct beam to zero
   //if ( local_alloc_size_fftw > 0 )
   //{
   //   if ( mynode == rootnode ) psi_mag[0] = 0.0;
   //   else
   //      psi_mag[0] //  = 0.0;
   //         = log(
   //               1.0 + (
   //                     scale_factor 
   //                     * ( psi[0][0] * psi[0][0] + psi[0][1] * psi[0][1]  )
   //                     / NxNyNxNy
   //                   )
   //            );
   //}

   if ( local_alloc_size_fftw > 0 )
   {
      max_psi_mag = psi_mag[0];
      min_psi_mag = max_psi_mag;
   }
   for ( ptrdiff_t i=0; i < local_alloc_size_fftw; ++i)
   //for ( ptrdiff_t i=1; i < local_alloc_size_fftw; ++i)
   {
      psi_mag[i]
         = log(
               1.0 + (
                     scale_factor 
                     * ( psi[i][0] * psi[i][0] 
                         +  psi[i][1] * psi[i][1] 
                       ) / NxNyNxNy
                   )
            );
         if ( psi_mag[i] > max_psi_mag) max_psi_mag = psi_mag[i];
         if ( psi_mag[i] < min_psi_mag) min_psi_mag = psi_mag[i];
   }

   // debug printing real and imaginary parts separately
   //double* psi_re;
   //double* psi_im;
   //double max_psi_re;
   //double min_psi_re;
   //double max_psi_re_joined;
   //double min_psi_re_joined;
   //double max_psi_im;
   //double min_psi_im;
   //double max_psi_im_joined;
   //double min_psi_im_joined;
   //psi_re = new double[ local_alloc_size_fftw ];
   //psi_im = new double[ local_alloc_size_fftw ];
   //if ( local_alloc_size_fftw > 0 )
   //{
   //   //if ( mynode == rootnode ) psi_re[0] = 0.0;
   //   //else
   //   psi_re[0] //  = 0.0;
   //      = log(
   //            1.0 + ( scale_factor * ( psi[0][0] ) / NxNy )
   //         );
   //   max_psi_re = psi_re[0];
   //   min_psi_re = max_psi_re;
   //   psi_im[0] //  = 0.0;
   //      = log(
   //            1.0 + ( scale_factor * ( psi[0][1] ) / NxNy )
   //         );
   //   max_psi_im = psi_im[0];
   //   min_psi_im = max_psi_im;
   //}

   //for ( ptrdiff_t i=1; i < local_alloc_size_fftw; ++i)
   //{
   //   psi_re[i]
   //      = log(
   //            1.0 + ( scale_factor * ( psi[i][0] ) / NxNy )
   //         );
   //      if ( psi_re[i] > max_psi_re) max_psi_re = psi_re[i];
   //      if ( psi_re[i] < min_psi_re) min_psi_re = psi_re[i];
   //   psi_im[i]
   //      = log(
   //            1.0 + ( scale_factor * ( psi[i][1] ) / NxNy )
   //         );
   //      if ( psi_im[i] > max_psi_im) max_psi_im = psi_im[i];
   //      if ( psi_im[i] < min_psi_im) min_psi_im = psi_im[i];
   //}
   // end debug

   double* psi_mag_joined;
   //double* psi_re_joined; // debug
   //double* psi_im_joined; // debug
   if( mynode == rootnode )
   {
      psi_mag_joined = new double[ Nx * Ny ];
      //psi_re_joined = new double[ Nx * Ny ]; // debug
      //psi_im_joined = new double[ Nx * Ny ]; // debug
   }

   MPI_Gather( psi_mag, local_alloc_size_fftw, MPI_DOUBLE, 
         psi_mag_joined, local_alloc_size_fftw, MPI_DOUBLE,
         rootnode, comm);

   MPI_Reduce( &max_psi_mag, &max_psi_mag_joined, 
         1, MPI_DOUBLE, 
         MPI_MAX,
         rootnode, comm);

   MPI_Reduce( &min_psi_mag, &min_psi_mag_joined, 
         1, MPI_DOUBLE, 
         MPI_MIN,
         rootnode, comm);

   //MPI_Gather( psi_re, local_alloc_size_fftw, MPI_DOUBLE, 
   //      psi_re_joined, local_alloc_size_fftw, MPI_DOUBLE,
   //      rootnode, comm);

   //MPI_Reduce( &max_psi_re, &max_psi_re_joined, 
   //      1, MPI_DOUBLE, 
   //      MPI_MAX,
   //      rootnode, comm);

   //MPI_Reduce( &min_psi_re, &min_psi_re_joined, 
   //      1, MPI_DOUBLE, 
   //      MPI_MIN,
   //      rootnode, comm);

   //MPI_Gather( psi_im, local_alloc_size_fftw, MPI_DOUBLE, 
   //      psi_im_joined, local_alloc_size_fftw, MPI_DOUBLE,
   //      rootnode, comm);

   //MPI_Reduce( &max_psi_im, &max_psi_im_joined, 
   //      1, MPI_DOUBLE, 
   //      MPI_MAX,
   //      rootnode, comm);

   //MPI_Reduce( &min_psi_im, &min_psi_im_joined, 
   //      1, MPI_DOUBLE, 
   //      MPI_MIN,
   //      rootnode, comm);

   if ( mynode == rootnode )
   {
      //psi_mag_joined[0] = 0.0;
      //cout << "Writing diffraction to tiff" << endl // debug
      //   << " max_diffration_psi_mag_joined : " // debug
      //   << setprecision(20) // debug
      //   << max_psi_mag_joined << endl; // debug
      //cout << " min_diffraction_psi_mag_joined : " // debug
      //   << setprecision(20) // debug
      //   << min_psi_mag_joined << endl; // debug

      //cout << " max_diffration_psi_re_joined : " // debug
      //   << setprecision(20) // debug
      //   << max_psi_re_joined << endl; // debug
      //cout << " min_diffration_psi_re_joined : " // debug
      //   << setprecision(20) // debug
      //   << min_psi_re_joined << endl; // debug
      //cout << " max_diffraction_psi_im_joined : " // debug
      //   << setprecision(20) // debug
      //   << max_psi_im_joined << endl; // debug
      //cout << " min_diffraction_psi_im_joined : " // debug
      //   << setprecision(20) // debug
      //   << min_psi_im_joined << endl; // debug
      
      // Rearrange so that the origin is in the center
      double tmp_mag;
      //double tmp_re; // debug
      //double tmp_im; // debug

      for ( ptrdiff_t i=0; i < Nx/2; ++i)
      {
         for ( ptrdiff_t j=0; j < Ny/2; ++j)
         {
            tmp_mag = psi_mag_joined[j + i * Ny];
            psi_mag_joined[j + i * Ny] 
               = psi_mag_joined[(j + Ny/2) + (i + Nx/2) * Ny];
            psi_mag_joined[(j + Ny/2) + (i + Nx/2) * Ny] = tmp_mag;

            //tmp_re = psi_re_joined[j + i * Ny]; // debug
            //psi_re_joined[j + i * Ny]  // debug
            //   = psi_re_joined[(j + Ny/2) + (i + Nx/2) * Ny]; // debug
            //psi_re_joined[(j + Ny/2) + (i + Nx/2) * Ny] = tmp_re; // debug

            //tmp_im = psi_im_joined[j + i * Ny]; // debug
            //psi_im_joined[j + i * Ny]  // debug
            //   = psi_im_joined[(j + Ny/2) + (i + Nx/2) * Ny]; // debug
            //psi_im_joined[(j + Ny/2) + (i + Nx/2) * Ny] = tmp_im; // debug
         }

         for ( ptrdiff_t j=Ny/2; j < Ny; ++j)
         {
            tmp_mag = psi_mag_joined[j + i * Ny];
            psi_mag_joined[j + i * Ny] 
               = psi_mag_joined[(j - Ny/2) + (i + Nx/2) * Ny];
            psi_mag_joined[(j - Ny/2) + (i + Nx/2) * Ny] = tmp_mag;

            //tmp_re = psi_re_joined[j + i * Ny]; // debug
            //psi_re_joined[j + i * Ny]  // debug
            //   = psi_re_joined[(j - Ny/2) + (i + Nx/2) * Ny]; // debug
            //psi_re_joined[(j - Ny/2) + (i + Nx/2) * Ny] = tmp_re; // debug

            //tmp_im = psi_im_joined[j + i * Ny]; // debug
            //psi_im_joined[j + i * Ny]  // debug
            //   = psi_im_joined[(j - Ny/2) + (i + Nx/2) * Ny]; // debug
            //psi_im_joined[(j - Ny/2) + (i + Nx/2) * Ny] = tmp_im; // debug
         }
      }

      imageWriter.write_tif_grayscale16(
            psi_mag_joined,
            min_psi_mag_joined,
            max_psi_mag_joined,
            Nx, Ny,
            resolutionUnit,
            xResolution, yResolution,
            outFileName_prefix + "_diffraction_scalefactor_" 
               + to_string(scale_factor)
            );

      //imageWriter.write_tif_grayscale16( // debug
      //      psi_re_joined, // debug
      //      min_psi_re_joined, // debug
      //      max_psi_re_joined, // debug
      //      Nx, Ny, // debug
      //      resolutionUnit,
      //      xResolution, yResolution,
      //      outFileName_prefix + "_diffraction_re_scalefactor_"  // debug
      //         + to_string(scale_factor) // debug
      //      ); // debug

      //imageWriter.write_tif_grayscale16( // debug
      //      psi_im_joined, // debug
      //      min_psi_im_joined, // debug
      //      max_psi_im_joined, // debug
      //      Nx, Ny, // debug
      //      resolutionUnit,
      //      xResolution, yResolution,
      //      outFileName_prefix + "_diffraction_im_scalefactor_"  // debug
      //         + to_string(scale_factor) // debug
      //      ); // debug

      delete[] psi_mag_joined;
      //delete[] psi_re_joined; // debug
      //delete[] psi_im_joined; // debug
   }

   delete[] psi_mag;
   //delete[] psi_re; // debug
   //delete[] psi_im; // debug
   return EXIT_SUCCESS;
}

//-----------------------------------------------------------------------
int TEM_NS::output_diffraction(
      const fftw_complex* const psi,
      const double& scale_factor,
      const ptrdiff_t& local_alloc_size_fftw,
      const ptrdiff_t& Nx_local,
      const ptrdiff_t& Nx,
      const ptrdiff_t& Ny,
      const size_t& resolutionUnit,
      const double& xResolution, const double& yResolution,
      const string& outFileName_prefix,
      const int& mynode,
      const int& rootnode,
      MPI_Comm comm
      )
{
   // Precondition:
   // - psi represents the data in reciprocal space, requiring rearrangement

   writeImage imageWriter;

   double* psi_mag;
   psi_mag = new double[ local_alloc_size_fftw ];

   double max_psi_mag;
   double min_psi_mag;
   double max_psi_mag_joined;
   double min_psi_mag_joined;

   // Calculate intensities, scaling them onto a log scale, and determine 
   //  the max and min of the resulting log scale intensities.
   // Kirkland (2010) eq 4.28:
   // D(k_{x}, k_{y}) = log( 1 + c | F(k_{x}, k_{y}) |^{2} ) 

   //cout << "scale_factor : " << setprecision(20) << scale_factor << endl;//debug
   if ( local_alloc_size_fftw > 0 )
   {
      //if ( mynode == rootnode ) psi_mag[0] = 0.0;   // TODO: this removes
      //else                                          //  the direct beam
         psi_mag[0] //  = 0.0;
            = log(
                  1.0 + (
                        scale_factor 
                        * ( psi[0][0] * psi[0][0] + psi[0][1] * psi[0][1]  )
                      )
               );
      max_psi_mag = psi_mag[0];
      min_psi_mag = max_psi_mag;
   }

   for ( ptrdiff_t i=1; i < local_alloc_size_fftw; ++i)
   {
      psi_mag[i]
         = log(
               1.0 + (
                     scale_factor 
                     * ( psi[i][0] * psi[i][0] +  psi[i][1] * psi[i][1] )
                   )
            );
         if ( psi_mag[i] > max_psi_mag) max_psi_mag = psi_mag[i];
         if ( psi_mag[i] < min_psi_mag) min_psi_mag = psi_mag[i];
   }

   // debug printing real and imaginary parts separately
   double* psi_re;
   double* psi_im;
   double max_psi_re;
   double min_psi_re;
   double max_psi_re_joined;
   double min_psi_re_joined;
   double max_psi_im;
   double min_psi_im;
   double max_psi_im_joined;
   double min_psi_im_joined;
   psi_re = new double[ local_alloc_size_fftw ];
   psi_im = new double[ local_alloc_size_fftw ];
   if ( local_alloc_size_fftw > 0 )
   {
      //if ( mynode == rootnode ) psi_re[0] = 0.0;
      //else
      psi_re[0] //  = 0.0;
         = log(
               1.0 + ( scale_factor * ( psi[0][0] ) )  
            );
      max_psi_re = psi_re[0];
      min_psi_re = max_psi_re;
      psi_im[0] //  = 0.0;
         = log(
               1.0 + ( scale_factor * ( psi[0][1] ) )
            );
      max_psi_im = psi_im[0];
      min_psi_im = max_psi_im;
   }

   for ( ptrdiff_t i=1; i < local_alloc_size_fftw; ++i)
   {
      psi_re[i]
         = log(
               1.0 + ( scale_factor * ( psi[i][0] ) )
            );
         if ( psi_re[i] > max_psi_re) max_psi_re = psi_re[i];
         if ( psi_re[i] < min_psi_re) min_psi_re = psi_re[i];
      psi_im[i]
         = log(
               1.0 + ( scale_factor * ( psi[i][1] ) )
            );
         if ( psi_im[i] > max_psi_im) max_psi_im = psi_im[i];
         if ( psi_im[i] < min_psi_im) min_psi_im = psi_im[i];
   }
   // end debug

   double* psi_mag_joined;
   //double* psi_re_joined; // debug
   //double* psi_im_joined; // debug
   if( mynode == rootnode )
   {
      psi_mag_joined = new double[ Nx * Ny ];
      //psi_re_joined = new double[ Nx * Ny ]; // debug
      //psi_im_joined = new double[ Nx * Ny ]; // debug
   }

   MPI_Gather( psi_mag, local_alloc_size_fftw, MPI_DOUBLE, 
         psi_mag_joined, local_alloc_size_fftw, MPI_DOUBLE,
         rootnode, comm);

   MPI_Reduce( &max_psi_mag, &max_psi_mag_joined, 
         1, MPI_DOUBLE, 
         MPI_MAX,
         rootnode, comm);

   MPI_Reduce( &min_psi_mag, &min_psi_mag_joined, 
         1, MPI_DOUBLE, 
         MPI_MIN,
         rootnode, comm);

   //MPI_Gather( psi_re, local_alloc_size_fftw, MPI_DOUBLE, 
   //      psi_re_joined, local_alloc_size_fftw, MPI_DOUBLE,
   //      rootnode, comm);

   //MPI_Reduce( &max_psi_re, &max_psi_re_joined, 
   //      1, MPI_DOUBLE, 
   //      MPI_MAX,
   //      rootnode, comm);

   //MPI_Reduce( &min_psi_re, &min_psi_re_joined, 
   //      1, MPI_DOUBLE, 
   //      MPI_MIN,
   //      rootnode, comm);

   //MPI_Gather( psi_im, local_alloc_size_fftw, MPI_DOUBLE, 
   //      psi_im_joined, local_alloc_size_fftw, MPI_DOUBLE,
   //      rootnode, comm);

   //MPI_Reduce( &max_psi_im, &max_psi_im_joined, 
   //      1, MPI_DOUBLE, 
   //      MPI_MAX,
   //      rootnode, comm);

   //MPI_Reduce( &min_psi_im, &min_psi_im_joined, 
   //      1, MPI_DOUBLE, 
   //      MPI_MIN,
   //      rootnode, comm);

   if ( mynode == rootnode )
   {
      ////psi_mag_joined[0] = 0.0;
      //cout << "Writing diffraction to tiff" << endl // debug
      //   << " max_diffration_psi_mag_joined : " // debug
      //   << setprecision(20) // debug
      //   << max_psi_mag_joined << endl; // debug
      //cout << " min_diffraction_psi_mag_joined : " // debug
      //   << setprecision(20) // debug
      //   << min_psi_mag_joined << endl; // debug

      //cout << " max_diffration_psi_re_joined : " // debug
      //   << setprecision(20) // debug
      //   << max_psi_re_joined << endl; // debug
      //cout << " min_diffration_psi_re_joined : " // debug
      //   << setprecision(20) // debug
      //   << min_psi_re_joined << endl; // debug
      //cout << " max_diffraction_psi_im_joined : " // debug
      //   << setprecision(20) // debug
      //   << max_psi_im_joined << endl; // debug
      //cout << " min_diffraction_psi_im_joined : " // debug
      //   << setprecision(20) // debug
      //   << min_psi_im_joined << endl; // debug
      
      // Rearrange so that the origin is in the center
      double tmp_mag;
      //double tmp_re; // debug
      //double tmp_im; // debug

      for ( ptrdiff_t i=0; i < Nx/2; ++i)
      {
         for ( ptrdiff_t j=0; j < Ny/2; ++j)
         {
            tmp_mag = psi_mag_joined[j + i * Ny];
            psi_mag_joined[j + i * Ny] 
               = psi_mag_joined[(j + Ny/2) + (i + Nx/2) * Ny];
            psi_mag_joined[(j + Ny/2) + (i + Nx/2) * Ny] = tmp_mag;

            //tmp_re = psi_re_joined[j + i * Ny]; // debug
            //psi_re_joined[j + i * Ny]  // debug
            //   = psi_re_joined[(j + Ny/2) + (i + Nx/2) * Ny]; // debug
            //psi_re_joined[(j + Ny/2) + (i + Nx/2) * Ny] = tmp_re; // debug

            //tmp_im = psi_im_joined[j + i * Ny]; // debug
            //psi_im_joined[j + i * Ny]  // debug
            //   = psi_im_joined[(j + Ny/2) + (i + Nx/2) * Ny]; // debug
            //psi_im_joined[(j + Ny/2) + (i + Nx/2) * Ny] = tmp_im; // debug
         }

         for ( ptrdiff_t j=Ny/2; j < Ny; ++j)
         {
            tmp_mag = psi_mag_joined[j + i * Ny];
            psi_mag_joined[j + i * Ny] 
               = psi_mag_joined[(j - Ny/2) + (i + Nx/2) * Ny];
            psi_mag_joined[(j - Ny/2) + (i + Nx/2) * Ny] = tmp_mag;

            //tmp_re = psi_re_joined[j + i * Ny]; // debug
            //psi_re_joined[j + i * Ny]  // debug
            //   = psi_re_joined[(j - Ny/2) + (i + Nx/2) * Ny]; // debug
            //psi_re_joined[(j - Ny/2) + (i + Nx/2) * Ny] = tmp_re; // debug

            //tmp_im = psi_im_joined[j + i * Ny]; // debug
            //psi_im_joined[j + i * Ny]  // debug
            //   = psi_im_joined[(j - Ny/2) + (i + Nx/2) * Ny]; // debug
            //psi_im_joined[(j - Ny/2) + (i + Nx/2) * Ny] = tmp_im; // debug
         }
      }

      imageWriter.write_tif_grayscale16(
            psi_mag_joined,
            min_psi_mag_joined,
            max_psi_mag_joined,
            Nx, Ny,
            resolutionUnit,
            xResolution, yResolution,
            outFileName_prefix + "_diffraction_scalefactor_" 
               + to_string(scale_factor)
            );

      //imageWriter.write_tif_grayscale16( // debug
      //      psi_re_joined, // debug
      //      min_psi_re_joined, // debug
      //      max_psi_re_joined, // debug
      //      Nx, Ny, // debug
      //      resolutionUnit,
      //      xResolution, yResolution,
      //      outFileName_prefix + "_diffraction_re_scalefactor_"  // debug
      //         + to_string(scale_factor) // debug
      //      ); // debug

      //imageWriter.write_tif_grayscale16( // debug
      //      psi_im_joined, // debug
      //      min_psi_im_joined, // debug
      //      max_psi_im_joined, // debug
      //      Nx, Ny, // debug
      //      resolutionUnit,
      //      xResolution, yResolution,
      //      outFileName_prefix + "_diffraction_im_scalefactor_"  // debug
      //         + to_string(scale_factor) // debug
      //      ); // debug

      delete[] psi_mag_joined;
      //delete[] psi_re_joined; // debug
      //delete[] psi_im_joined; // debug
   }

   delete[] psi_mag;
   //delete[] psi_re; // debug
   //delete[] psi_im; // debug
   return EXIT_SUCCESS;
}

//-----------------------------------------------------------------------
int TEM_NS::output_diffraction(
      const double* const psi,
      const double& scale_factor,
      const ptrdiff_t& local_alloc_size_fftw,
      const ptrdiff_t& Nx_local,
      const ptrdiff_t& Nx,
      const ptrdiff_t& Ny,
      const size_t& resolutionUnit,
      const double& xResolution, const double& yResolution,
      const string& outFileName_prefix,
      const int& mynode,
      const int& rootnode,
      MPI_Comm comm
      )
{
   // Precondition:
   // - psi represents the data in reciprocal space, requiring rearrangement
   //    and also requiring division by sqrt(Nx * Ny) due to fftw.

   writeImage imageWriter;

   double* psi_mag;
   psi_mag = new double[ local_alloc_size_fftw ];

   const double NxNyNxNy = Nx * Ny * Nx * Ny ;

   double max_psi_mag;
   double min_psi_mag;
   double max_psi_mag_joined;
   double min_psi_mag_joined;

   // Calculate intensities, scaling them onto a log scale, and determine 
   //  the max and min of the resulting log scale intensities.
   // Kirkland (2010) eq 4.28:
   // D(k_{x}, k_{y}) = log( 1 + c | F(k_{x}, k_{y}) |^{2} ) 

   //cout << "scale_factor : " << setprecision(20) << scale_factor << endl;//debug
   if ( local_alloc_size_fftw > 0 )
   {
      //if ( mynode == rootnode ) psi_mag[0] = 0.0;
      //else
         psi_mag[0] //  = 0.0;
            = log(
                  1.0 + ( scale_factor * ( psi[0] * psi[0]  ) )
               );
      max_psi_mag = psi_mag[0];
      min_psi_mag = max_psi_mag;
   }

   for ( ptrdiff_t i=1; i < local_alloc_size_fftw; ++i)
   {
      psi_mag[i]
         = log(
               1.0 + ( scale_factor * ( psi[i] * psi[i] ) )
            );
         if ( psi_mag[i] > max_psi_mag) max_psi_mag = psi_mag[i];
         if ( psi_mag[i] < min_psi_mag) min_psi_mag = psi_mag[i];
   }


   double* psi_mag_joined;
   if( mynode == rootnode )
   {
      psi_mag_joined = new double[ Nx * Ny ];
   }

   MPI_Gather( psi_mag, local_alloc_size_fftw, MPI_DOUBLE, 
         psi_mag_joined, local_alloc_size_fftw, MPI_DOUBLE,
         rootnode, comm);

   MPI_Reduce( &max_psi_mag, &max_psi_mag_joined, 
         1, MPI_DOUBLE, 
         MPI_MAX,
         rootnode, comm);

   MPI_Reduce( &min_psi_mag, &min_psi_mag_joined, 
         1, MPI_DOUBLE, 
         MPI_MIN,
         rootnode, comm);

   if ( mynode == rootnode )
   {
      //psi_mag_joined[0] = 0.0;
      //cout << "Writing diffraction to tiff" << endl // debug
      //   << " max_diffraction_psi_mag_joined : " // debug
      //   << setprecision(20) // debug
      //   << max_psi_mag_joined << endl; // debug
      //cout << " min_diffraction_psi_mag_joined : " // debug
      //   << setprecision(20) // debug
      //   << min_psi_mag_joined << endl; // debug

      // Rearrange so that the origin is in the center
      double tmp_mag;

      for ( ptrdiff_t i=0; i < Nx/2; ++i)
      {
         for ( ptrdiff_t j=0; j < Ny/2; ++j)
         {
            tmp_mag = psi_mag_joined[j + i * Ny];
            psi_mag_joined[j + i * Ny] 
               = psi_mag_joined[(j + Ny/2) + (i + Nx/2) * Ny];
            psi_mag_joined[(j + Ny/2) + (i + Nx/2) * Ny] = tmp_mag;
         }

         for ( ptrdiff_t j=Ny/2; j < Ny; ++j)
         {
            tmp_mag = psi_mag_joined[j + i * Ny];
            psi_mag_joined[j + i * Ny] 
               = psi_mag_joined[(j - Ny/2) + (i + Nx/2) * Ny];
            psi_mag_joined[(j - Ny/2) + (i + Nx/2) * Ny] = tmp_mag;
         }
      }

      imageWriter.write_tif_grayscale16(
            psi_mag_joined,
            min_psi_mag_joined,
            max_psi_mag_joined,
            Nx, Ny,
            resolutionUnit,
            xResolution, yResolution,
            outFileName_prefix + "_diffraction_scalefactor_" 
               + to_string(scale_factor)
            );

      delete[] psi_mag_joined;
   }

   delete[] psi_mag;
   return EXIT_SUCCESS;
}

void TEM_NS::output_stem_image(
      const double* const stem_image,
      const ptrdiff_t& Nx,
      const ptrdiff_t& Ny,
      const size_t& resolutionUnit,
      const double& xResolution, 
      const double& yResolution, 
      const string& outFileName_prefix
      )
{

   writeImage imageWriter;

   double min_stem_image;
   double max_stem_image;
   
   // Determine the max and min intensities
   min_stem_image = stem_image[0];
   max_stem_image = min_stem_image;
   for ( ptrdiff_t i = 0; i < Nx * Ny; i++)
   {  
      if( stem_image[i] > max_stem_image ) max_stem_image = stem_image[i]; 
      if( stem_image[i] < min_stem_image ) min_stem_image = stem_image[i];
   }
 
   //cout << " max_stem_image : "
   //   << setprecision(20) << max_stem_image << endl;
   //cout <<" min_stem_image : "
   //   << setprecision(20) << min_stem_image << endl;

   imageWriter.write_tif_grayscale16(
         stem_image,
         min_stem_image,
         max_stem_image,
         Nx, Ny,
         resolutionUnit,
         xResolution, yResolution,
         outFileName_prefix  
            +"_stem_image_output");
 
   imageWriter.write_tif_grayscale16_logscale(
         stem_image,
         min_stem_image,
         max_stem_image,
         Nx, Ny,
         resolutionUnit,
         xResolution, yResolution,
         outFileName_prefix 
            +  "_stem_image_output"
            + "_logscale" );
 

   return;
} // output_stem_image

void TEM_NS::output_correlograph_image(
      const double* const correlograph,
      const ptrdiff_t& Nx,
      const ptrdiff_t& Ny,
      const size_t& resolutionUnit,
      const double& xResolution, 
      const double& yResolution, 
      const string& outFileName_prefix,
      const unsigned int& logflag
      )
{

   writeImage imageWriter;

   double min_image;
   double max_image;
   
   // Determine the max and min intensities
   min_image = correlograph[0];
   max_image = min_image;
   for ( ptrdiff_t i = 0; i < Nx * Ny; i++)
   {  
      if( correlograph[i] > max_image ) max_image = correlograph[i]; 
      if( correlograph[i] < min_image ) min_image = correlograph[i];
   }

   if ( ! logflag )
      imageWriter.write_tif_grayscale16(
            correlograph,
            min_image,
            max_image,
            Nx, Ny,
            resolutionUnit,
            xResolution, yResolution,
            outFileName_prefix);
 
   if ( logflag )
      imageWriter.write_tif_grayscale16_logscale(
            correlograph,
            min_image,
            max_image,
            Nx, Ny,
            resolutionUnit,
            xResolution, yResolution,
            outFileName_prefix 
               + "_logscale" );
   return;
} // output_correlograph_image
#endif
