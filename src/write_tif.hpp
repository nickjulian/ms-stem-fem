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
// File: write_tif.hpp
// Purpose: to define tools for writing data to 8, 16, & 32 bit tif images 
//    using libtiff

#ifndef WRITE_TIF_HPP
#define WRITE_TIF_HPP

#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>   // EXIT_SUCCESS & EXIT_FAILURE 
#include <cstdio>   
#include <limits>    // std::numeric_limits::digits to find bits per sample
#include <string>
#include <string.h>  // required for memcpy
#include <iomanip>   // setprecision()

#include <stdint.h>

#include <mpi.h>

#include <tiffio.h>

#include "to_string.hpp"

#include <fftw3-mpi.h>
//#include "../include/fftw3-mpi.h"

using std::cout;
using std::endl;
using std::string;
using std::cerr;
using std::setprecision;

namespace TEM_NS
{
class writeImage
{
   //TODO: use c++ std::vector<> instead of arrays 
   public:
      int write_tif_grayscale8( const double* const outDataRaw,
            const double& min, const double& max,
            //const size_t& width, const size_t& height, 
            const size_t& height, const size_t& width,
            const size_t& resolutionUnit,
            const double& xResolution, const double& yResolution,
            const string& outTifPrefix);
      // Precondition:
      //    - outDataRaw[] contains width*height double values in 
      //       row-major format (contiguous rows in memory)
      //       format, with maximum and minimum values max and min 
      //       respectively.
      // Postcondition:
      //    - an 8-bit tif image file has been written to outTifPrefix.tif
      
      int write_tif_grayscale16( const double* const outDataRaw,
            const double& min, const double& max,
            //const size_t& width, const size_t& height, 
            const size_t& height, const size_t& width,
            const size_t& resolutionUnit,
            const double& xResolution, const double& yResolution,
            const string& outTifPrefix);
      
      int write_tif_grayscale16_append( const double* const outDataRaw,
            const double& min, const double& max,
            //const size_t& width, const size_t& height, 
            const size_t& height, const size_t& width,
            const size_t& resolutionUnit,
            const double& xResolution, const double& yResolution,
            const string& outTifPrefix);
      
      int write_tif_grayscale32( const double* const outDataRaw,
            const double& min, const double& max,
            //const size_t& width, const size_t& height, 
            const size_t& height, const size_t& width,
            const size_t& resolutionUnit,
            const double& xResolution, const double& yResolution,
            const string& outTifPrefix);

      int write_tif_grayscale8_logscale( const double* const outDataRaw,
            const double& min, const double& max,
            //const size_t& width, const size_t& height, 
            const size_t& height, const size_t& width,
            const size_t& resolutionUnit,
            const double& xResolution, const double& yResolution,
            const string& outTifPrefix);

      int write_tif_grayscale16_logscale( const double* const outDataRaw,
            const double& min, const double& max,
            //const size_t& width, const size_t& height, 
            const size_t& height, const size_t& width,
            const size_t& resolutionUnit,
            const double& xResolution, const double& yResolution,
            const string& outTifPrefix);

      int write_tif_grayscale32_logscale( const double* const outDataRaw,
            const double& min, const double& max,
            const size_t& height, const size_t& width,
            const size_t& resolutionUnit,
            const double& xResolution, const double& yResolution,
            const string& outTifPrefix);
   private:
      template<typename SampleType> 
         int write_tif_grayscale( const double* const outDataRaw,
            const double& min, const double& max,
            const size_t& height, const size_t& width,
            const size_t& resolutionUnit,
            const double& xResolution, const double& yResolution,
            const string& outTifPrefix,
            const bool& logarithmTag );

      template<typename SampleType> 
         int write_tif_grayscale_append( const double* const outDataRaw,
            const double& min, const double& max,
            const size_t& height, const size_t& width,
            const size_t& resolutionUnit,
            const double& xResolution, const double& yResolution,
            const string& outTifPrefix,
            const bool& logarithmTag );
}; // writeTif class

template<typename SampleType> 
   int writeImage::write_tif_grayscale( const double* const outDataRaw,
      const double& min, const double& max,
      const size_t& height, const size_t& width,// height == Nx, width == Ny
      const size_t& resolutionUnit,
      const double& xResolution, const double& yResolution,
      const string& outTifPrefix, 
      const bool& logarithmTag )
{
   /////////////////////////////////////////////////////////////////////////
   // debug, input review
   /////////////////////////////////////////////////////////////////////////
   //cout << "sizeof(outImage[0]) : " << sizeof(outImage[0]) << endl;
   //cout << "width : " << width << endl;
   //cout << "height : " << height << endl;
   //cout << "outTifPrefix : " << outTifPrefix << endl;
   
   /////////////////////////////////////////////////////////////////////////
   // debug, system characterization
   /////////////////////////////////////////////////////////////////////////
   //cout << "numeric_limits< char >::digits : " // debug 
   //   << std::numeric_limits< char >::digits //debug            7
   //   << endl;  // debug
   // cout << "numeric_limits< char >::radix : " // debug         
   //   << std::numeric_limits< char >::radix  //debug            2 
   //   << endl;  // debug
   //cout << "numeric_limits< unsigned char >::digits : " // debug 
   //   << std::numeric_limits< unsigned char >::digits //debug   8
   //   << endl;  // debug
   //cout << "numeric_limits< unsigned char >::radix : " // debug         
   //   << std::numeric_limits< unsigned char >::radix  //debug   2
   //   << endl;  // debug
   //cout << "numeric_limits< uint16_t >::digits : " // debug
   //   << std::numeric_limits< uint16_t >::digits //debug        16
   //   << endl;  // debug
   // cout << "numeric_limits< uint16_t >::radix : " // debug
   //   << std::numeric_limits< uint16_t >::radix  //debug        2
   //   << endl;  // debug
   //cout << "numeric_limits< uint32_t >::digits : " // debug
   //   << std::numeric_limits< uint32_t >::digits //debug        32
   //   << endl;  // debug
   //cout << "numeric_limits< uint32_t >::radix : " // debug
   //   << std::numeric_limits< uint32_t >::radix  //debug        2
   //   << endl;  // debug
   //cout << "numeric_limits< float >::digits : " // debug
   //   << std::numeric_limits< float >::digits //debug           24
   //   << endl;  // debug
   //cout << "numeric_limits< float >::radix : " // debug
   //   << std::numeric_limits< float >::radix  //debug           2
   //   << endl;  // debug
   //cout << "numeric_limits< double >::digits : " // debug
   //   << std::numeric_limits< double >::digits //debug          53
   //   << endl;  // debug
   //cout << "numeric_limits< double >::radix : " // debug
   //   << std::numeric_limits< double >::radix  //debug          2
   //   << endl;  // debug

   /////////////////////////////////////////////////////////////////////////
   // input validation 
   /////////////////////////////////////////////////////////////////////////
   int samplesPerPixel = 1; // grayscale
   // int samplesPerPixel = 4; // RGBA
   int bitsPerSample = std::numeric_limits<SampleType>::digits ;
   //cout << "bitsPerSample : " << bitsPerSample << endl; //debug
   if (  bitsPerSample != 8 // && bitsPerSample !=4 //4 isn't divisible by 8
         && bitsPerSample != 16 && bitsPerSample != 32 )
   {
      cerr << "Error - TIFFTAG_BITSPERSAMPLE currently restricted to "
         << "4, 8, 16, or 32. Current value of bitsPerSample : " 
         << bitsPerSample << endl;
      return EXIT_FAILURE;
   }
   
   /////////////////////////////////////////////////////////////////////////
   // scale the input data to fit the output image intensity range
   /////////////////////////////////////////////////////////////////////////

   // Note:
   // Nx == number of rows == height
   // Ny == number of columns == width
   // row-major <-> rows contiguous in memory
   // column-major <-> columns contiguous in memory
   
   // outDataRaw is a row-major array (contiguous rows)
   // outImage is also a row-major array (contiguous rows)
   double dataIntervalLength = max - min;
   SampleType* outImage; outImage = new SampleType[width * height];
   double SampleTypeDigitsFactor = 
      (pow(2, std::numeric_limits< SampleType >::digits )-1);
   
   if ( logarithmTag == true )
   {
      // TODO: assert for edge case of max == min
      double data_tmp;
      double dataIntervalLength_log = log( max - min + 1);
      for ( size_t i = 0; i < height; i++) // rows
         for ( size_t j = 0; j < width; j++) // columns
         {
            data_tmp = ( log( outDataRaw[ j + width* i ] - min + 1) ) 
               / dataIntervalLength_log; 
            if ( data_tmp > 0.0 )
            {
               outImage[ j + width* i ] = // outImage is row-major array
                  static_cast< SampleType >(  
                        data_tmp * SampleTypeDigitsFactor + 0.5 );
            } else outImage[ j + width* i ] = static_cast<SampleType>(0.0);
         }
   } else
   {
      double dataIntervalLength = max - min;
      double commonScalingFactor =
         SampleTypeDigitsFactor / dataIntervalLength;
      for ( size_t i = 0; i < height; i++)
         for ( size_t j = 0; j < width; j++)
         {
            // (curr_value - min)*255/(dataIntervalLength)
            outImage[ j + width* i ] = static_cast< SampleType >( 
                  (outDataRaw[ j + width* i ] - min) 
                  * commonScalingFactor
                  + 0.5 
                  ); 
         }
   }

   /////////////////////////////////////////////////////////////////////////
   // open tiff file
   /////////////////////////////////////////////////////////////////////////
   const string outTifName = outTifPrefix + ".tif";
   TIFF* outTif = TIFFOpen( outTifName.c_str(), "w" );
   if ( ! outTif )
   {
      cerr << "Output file could no be opened: " << outTifName << endl;
      return EXIT_FAILURE;
   }

   /////////////////////////////////////////////////////////////////////////
   // assign TIFF tags
   /////////////////////////////////////////////////////////////////////////

   TIFFSetField( outTif, TIFFTAG_IMAGEWIDTH, width );
   TIFFSetField( outTif, TIFFTAG_IMAGELENGTH, height );
   TIFFSetField( outTif, TIFFTAG_SAMPLESPERPIXEL, samplesPerPixel );
   TIFFSetField( outTif, TIFFTAG_BITSPERSAMPLE, bitsPerSample );

   // RESOLUTIONUNIT = 1 (no absolute unit), 2 (inch), 3 (cm)
   TIFFSetField( outTif, TIFFTAG_RESOLUTIONUNIT, resolutionUnit);
   // XRESOLUTION = The number of pixels per RESOLUTIONUNIT along IMAGEWIDTH
   TIFFSetField( outTif, TIFFTAG_XRESOLUTION, xResolution);
   TIFFSetField( outTif, TIFFTAG_YRESOLUTION, yResolution);

   // set the origin of the image to the top left
   TIFFSetField( outTif, TIFFTAG_ORIENTATION, ORIENTATION_TOPLEFT );
   TIFFSetField( outTif, TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG);
   //TIFFSetField( outTif, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_RGB ); // color
   TIFFSetField( outTif, TIFFTAG_PHOTOMETRIC, 1); // black is zero
   // TIFFSetField( outTif, TIFFTAG_SAMPLEFORMAT, 3);//floating point; fails
   TIFFSetField( outTif, TIFFTAG_ROWSPERSTRIP, 
         TIFFDefaultStripSize( outTif, width * samplesPerPixel) );

   /////////////////////////////////////////////////////////////////////////
   // write the image data 
   /////////////////////////////////////////////////////////////////////////
   tsize_t linebytes = width * samplesPerPixel * bitsPerSample/8; // TODO
   // bitsPerSample is taken from std::numeric_limits<SampleType>::digits 
   //  which appears to have possible values of: 7,8,16,24,32,53
   //  and so is restricted in the input validation section to 8,16, or 32

   if (TIFFScanlineSize( outTif ) != linebytes) 
   {
      cerr << "TIFFScanlineSize ("<< TIFFScanlineSize( outTif ) 
         << ") != linebytes (" << linebytes << ")" << endl;
      if ( outTif ) TIFFClose( outTif );
      return EXIT_FAILURE;
   }
   
   SampleType* buf; buf = ( SampleType* ) _TIFFmalloc( linebytes );

   for ( uint32 row = 0; row < height; row++ )
   {
      memcpy(buf, &outImage[(height - row - 1)*( width )], linebytes);
      if ( TIFFWriteScanline( outTif, buf, row, 0) < 0) break;
   }

   /////////////////////////////////////////////////////////////////////////
   // clean up
   /////////////////////////////////////////////////////////////////////////
   if ( buf ) _TIFFfree( buf );
   if ( outTif ) TIFFClose( outTif );

   delete[] outImage;

   return EXIT_SUCCESS;
} // write_tif_grayscale

template<typename SampleType> 
   int writeImage::write_tif_grayscale_append( 
      const double* const outDataRaw,
      const double& min, const double& max,
      const size_t& height, const size_t& width,// height == Nx, width == Ny
      const size_t& resolutionUnit,
      const double& xResolution, const double& yResolution,
      const string& outTifPrefix, 
      const bool& logarithmTag )
{
   /////////////////////////////////////////////////////////////////////////
   // input validation 
   /////////////////////////////////////////////////////////////////////////
   int samplesPerPixel = 1; // grayscale
   // int samplesPerPixel = 4; // RGBA
   int bitsPerSample = std::numeric_limits<SampleType>::digits ;
   //cout << "bitsPerSample : " << bitsPerSample << endl; //debug
   if (  bitsPerSample != 8 // && bitsPerSample !=4 //4 isn't divisible by 8
         && bitsPerSample != 16 && bitsPerSample != 32 )
   {
      cerr << "Error - TIFFTAG_BITSPERSAMPLE currently restricted to "
         << "4, 8, 16, or 32. Current value of bitsPerSample : " 
         << bitsPerSample << endl;
      return EXIT_FAILURE;
   }
   
   /////////////////////////////////////////////////////////////////////////
   // scale the input data to fit the output image intensity range
   /////////////////////////////////////////////////////////////////////////

   // Note:
   // Nx == number of rows == height
   // Ny == number of columns == width
   // row-major <-> rows contiguous in memory
   // column-major <-> columns contiguous in memory
   
   // outDataRaw is a row-major array (contiguous rows)
   // outImage is also a row-major array (contiguous rows)
   double dataIntervalLength = max - min;
   SampleType* outImage; outImage = new SampleType[width * height];
   double SampleTypeDigitsFactor = 
      (pow(2, std::numeric_limits< SampleType >::digits )-1);
   
   if ( logarithmTag == true )
   {
      // TODO: assert for edge case of max == min
      double data_tmp;
      double dataIntervalLength_log = log( max - min + 1);
      for ( size_t i = 0; i < height; i++) // rows
         for ( size_t j = 0; j < width; j++) // columns
         {
            data_tmp = ( log( outDataRaw[ j + width* i ] - min + 1) ) 
               / dataIntervalLength_log; 
            if ( data_tmp > 0.0 )
            {
               outImage[ j + width* i ] = // outImage is row-major array
                  static_cast< SampleType >(  
                        data_tmp * SampleTypeDigitsFactor + 0.5 );
            } else outImage[ j + width* i ] = static_cast<SampleType>(0.0);
         }
   } else
   {
      double dataIntervalLength = max - min;
      double commonScalingFactor =
         SampleTypeDigitsFactor / dataIntervalLength;
      for ( size_t i = 0; i < height; i++)
         for ( size_t j = 0; j < width; j++)
         {
            // (curr_value - min)*255/(dataIntervalLength)
            outImage[ j + width* i ] = static_cast< SampleType >( 
                  (outDataRaw[ j + width* i ] - min) 
                  * commonScalingFactor
                  + 0.5 
                  ); 
         }
   }

   /////////////////////////////////////////////////////////////////////////
   // open tiff file
   /////////////////////////////////////////////////////////////////////////
   const string outTifName = outTifPrefix + ".tif";
   TIFF* outTif = TIFFOpen( outTifName.c_str(), "a" );
   if ( ! outTif )
   {
      cerr << "Output file could no be opened: " << outTifName << endl;
      return EXIT_FAILURE;
   }

   /////////////////////////////////////////////////////////////////////////
   // assign TIFF tags
   /////////////////////////////////////////////////////////////////////////

   TIFFSetField( outTif, TIFFTAG_IMAGEWIDTH, width );
   TIFFSetField( outTif, TIFFTAG_IMAGELENGTH, height );
   TIFFSetField( outTif, TIFFTAG_SAMPLESPERPIXEL, samplesPerPixel );
   TIFFSetField( outTif, TIFFTAG_BITSPERSAMPLE, bitsPerSample );

   // RESOLUTIONUNIT = 1 (no absolute unit), 2 (inch), 3 (cm)
   TIFFSetField( outTif, TIFFTAG_RESOLUTIONUNIT, resolutionUnit);
   // XRESOLUTION = The number of pixels per RESOLUTIONUNIT along IMAGEWIDTH
   TIFFSetField( outTif, TIFFTAG_XRESOLUTION, xResolution);
   TIFFSetField( outTif, TIFFTAG_YRESOLUTION, yResolution);

   // set the origin of the image to the top left
   TIFFSetField( outTif, TIFFTAG_ORIENTATION, ORIENTATION_TOPLEFT );
   TIFFSetField( outTif, TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG);
   //TIFFSetField( outTif, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_RGB ); // color
   TIFFSetField( outTif, TIFFTAG_PHOTOMETRIC, 1); // black is zero
   // TIFFSetField( outTif, TIFFTAG_SAMPLEFORMAT, 3);//floating point; fails
   TIFFSetField( outTif, TIFFTAG_ROWSPERSTRIP, 
         TIFFDefaultStripSize( outTif, width * samplesPerPixel) );

   /////////////////////////////////////////////////////////////////////////
   // write the image data 
   /////////////////////////////////////////////////////////////////////////
   tsize_t linebytes = width * samplesPerPixel * bitsPerSample/8; // TODO
   // bitsPerSample is taken from std::numeric_limits<SampleType>::digits 
   //  which appears to have possible values of: 7,8,16,24,32,53
   //  and so is restricted in the input validation section to 8,16, or 32

   if (TIFFScanlineSize( outTif ) != linebytes) 
   {
      cerr << "TIFFScanlineSize ("<< TIFFScanlineSize( outTif ) 
         << ") != linebytes (" << linebytes << ")" << endl;
      if ( outTif ) TIFFClose( outTif );
      return EXIT_FAILURE;
   }
   
   SampleType* buf; buf = ( SampleType* ) _TIFFmalloc( linebytes );

   for ( uint32 row = 0; row < height; row++ )
   {
      memcpy(buf, &outImage[(height - row - 1)*( width )], linebytes);
      if ( TIFFWriteScanline( outTif, buf, row, 0) < 0) break;
   }

   /////////////////////////////////////////////////////////////////////////
   // clean up
   /////////////////////////////////////////////////////////////////////////
   if ( buf ) _TIFFfree( buf );
   if ( outTif ) TIFFClose( outTif );

   delete[] outImage;

   return EXIT_SUCCESS;
} // write_tif_grayscale


///////////////////////////////////////////////////////////////////////////
// public member functions ////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////
//
//
//int writeImage::write_tif_grayscale8( const double* const outDataRaw, 
//      const double& max, const double& min,// bounds of scale for input data
//      const size_t& width, const size_t& height,// image dimensions [pixels]
//      const string& outTifPrefix )
//{
//   /////////////////////////////////////////////////////////////////////////
//   // write the image
//   /////////////////////////////////////////////////////////////////////////
//   if ( write_tif_grayscale<unsigned char>(
//            outDataRaw,
//            min, max,
//            width, height, 
//            outTifPrefix )  == EXIT_FAILURE)
//   {
//      cerr << "Error : write_tif_grayscale8() failed" << endl;
//      return EXIT_FAILURE;
//   }
//   return EXIT_SUCCESS;
//} // write_tif_grayscale8
//
//int writeImage::write_tif_grayscale16( const double* const outDataRaw, 
//      const double& min, const double& max,// bounds of scale for input data
//      const size_t& width, const size_t& height,// image dimensions [pixels]
//      const string& outTifPrefix )
//{
//   /////////////////////////////////////////////////////////////////////////
//   // write the image
//   /////////////////////////////////////////////////////////////////////////
//   if ( write_tif_grayscale<uint16_t>(
//            outDataRaw,
//            min, max,
//            width, height, 
//            outTifPrefix )  == EXIT_FAILURE)
//   {
//      cerr << "Error : write_tif_grayscale16() failed" << endl;
//      return EXIT_FAILURE;
//   }
//   return EXIT_SUCCESS;
//} // write_tif_grayscale16
//
//int writeImage::write_tif_grayscale32( const double* const outDataRaw, 
//      const double& min, const double& max,// bounds of scale for input data
//      const size_t& width, const size_t& height,// image dimensions [pixels]
//      const string& outTifPrefix )
//{
//   /////////////////////////////////////////////////////////////////////////
//   // write the image
//   /////////////////////////////////////////////////////////////////////////
//   if ( write_tif_grayscale<uint32_t>(
//            outDataRaw,
//            min, max,
//            width, height, 
//            outTifPrefix )  == EXIT_FAILURE)
//   {
//      cerr << "Error : write_tif_grayscale32() failed" << endl;
//      return EXIT_FAILURE;
//   }
//   return EXIT_SUCCESS;
//} // write_tif_grayscale16

int output_diffraction_with_renormalization(
      const fftw_complex* const psi,
      const double& scale_factor,
      const ptrdiff_t& local_alloc_size_fftw,
      const ptrdiff_t& Nx_local,
      const ptrdiff_t& Nx,
      const ptrdiff_t& Ny,
      const size_t& resolutionUnit,
      const double& xResolution, const double& yResolution,
      const string& outFileName_prefix,
      const int* const psi_mag_strides,
      const int* const psi_mag_displacements,
      const int& mynode,
      const int& rootnode,
      MPI_Comm comm 
      );

int output_diffraction_with_renormalization_append(
      const fftw_complex* const psi,
      const double& scale_factor,
      const ptrdiff_t& local_alloc_size_fftw,
      const ptrdiff_t& Nx_local,
      const ptrdiff_t& Nx,
      const ptrdiff_t& Ny,
      const size_t& resolutionUnit,
      const double& xResolution, const double& yResolution,
      const string& outFileName_prefix,
      const int* const psi_mag_strides,
      const int* const psi_mag_displacements,
      const int& mynode,
      const int& rootnode,
      MPI_Comm comm 
      );

int output_diffraction(
      const fftw_complex* const psi,
      const double& scale_factor,
      const ptrdiff_t& local_alloc_size_fftw,
      const ptrdiff_t& Nx_local,
      const ptrdiff_t& Nx,
      const ptrdiff_t& Ny,
      const size_t& resolutionUnit,
      const double& xResolution, const double& yResolution,
      const string& outFileName_prefix,
      const int* const psi_mag_strides,
      const int* const psi_mag_displacements,
      const int& mynode,
      const int& rootnode,
      MPI_Comm comm 
      );

int output_diffraction(
      const double* const psi,
      const double& scale_factor,
      const ptrdiff_t& local_alloc_size_fftw,
      const ptrdiff_t& Nx_local,
      const ptrdiff_t& Nx,
      const ptrdiff_t& Ny,
      const size_t& resolutionUnit,
      const double& xResolution, const double& yResolution,
      const string& outFileName_prefix,
      const int* const psi_mag_strides,
      const int* const psi_mag_displacements,
      const int& mynode,
      const int& rootnode,
      MPI_Comm comm 
      );

int output_diffraction_append(
      const fftw_complex* const psi,
      const double& scale_factor,
      const double& lambda_sqr,
      const double& alpha_max_sqr,
      const ptrdiff_t& local_alloc_size_fftw,
      const ptrdiff_t& Nx_local,
      const double* const kx_local,   // kx domain elements
      const ptrdiff_t& Nx,
      //const double* const kx_joined,   // kx domain elements
      const ptrdiff_t& Ny,
      const double* const ky, // ky domain elements
      const size_t& resolutionUnit,
      const double& xResolution, const double& yResolution,
      const string& outFileName_prefix,
      const int* const psi_mag_strides,
      const int* const psi_mag_displacements,
      const int& mynode,
      const int& rootnode,
      MPI_Comm comm 
      );

void output_stem_image(
      const double* const stem_image,
      const ptrdiff_t& Nx,
      const ptrdiff_t& Ny,
      const size_t& resolutionUnit,
      const double& xResolution, 
      const double& yResolution, 
      const string& outFileName_prefix
      );

void output_correlograph_image(
      const double* const correlograph,
      const ptrdiff_t& Nx,
      const ptrdiff_t& Ny,
      const size_t& resolutionUnit,
      const double& xResolution, 
      const double& yResolution, 
      const string& outFileName_prefix,
      const unsigned int& logflag
      );
} // TEM_NS
#endif
