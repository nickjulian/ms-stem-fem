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
// File: variance.cpp
// Purpose:

#ifndef VARIANCE_CPP
#define VARIANCE_CPP


#include <cstdlib>
#include <iostream>

#include "variance.hpp"

//using namespace std;
using std::sqrt;
using std::cerr;
using std::endl;
using std::setw;
using std::setprecision;

using std::vector;
using std::list;

// need to send a size_t via MPI, but it's a derived type and we need to 
//  choose the proper MPI_ datatype using a macro
//  from stackoverflow user 5239503/Gilles
#include <stdint.h>
#include <limits.h>
//#if SIZE_MAX == UCHAR_MAX
//   #define my_MPI_SIZE_T MPI_UNSIGNED_CHAR
//#elif SIZE_MAX == USHRT_MAX
//   #define my_MPI_SIZE_T MPI_UNSIGNED_SHORT
//#elif SIZE_MAX == UINT_MAX
//   #define my_MPI_SIZE_T MPI_UNSIGNED
//#elif SIZE_MAX == ULONG_MAX
//   #define my_MPI_SIZE_T MPI_UNSIGNED_LONG
//#elif SIZE_MAX == ULLONG_MAX
//   #define my_MPI_SIZE_T MPI_UNSIGNED_LONG_LONG
//#else
//   #error "Unable to define custom MPI data type corresponding to size_t"
//#endif


int TEM_NS::integrate_out_theta_fftw( 
      // integrate to make f( \vec{k} ) into f( | \vec{k} | )
      //const double** const psi,// data to be binned,[][0] real,[][1] im
      const fftw_complex* const psi,//data to be integrated azimuthally
      const double* const kx_local, const size_t& Nx_local,
      const double* const ky, const size_t& Ny,
      const std::vector<double>& binning_boundaries,// bin boundaries
      int* bin_counts,
      double* data1D_local // output, 
             // having binning_boundaries.size() - 1 allocated elements
      )
{
   // Precondition : 
   // // - data1D_local has been allocated to have 
   // //    number_of_bins == bwcutoff_t / delta_k  
   // //       == binning_boundaries.size() - 1
   // //       == number of elements of bin_counts[]
   // //    elements
   // - binning boundaries do not exceed the maximum or minimum of the
   //    domain kx, ky of the data to be binned, but the non-zero data may 
   //    lie entirely within the bounds of the binning
   
   // Plan:
   // - create a vector of indexed_vector_magnitude elements, each of
   //    which contains array indices and a magnitude equal to the sum of 
   //    squares of those indices
   // - sort them by their magnitude member (v_mag_sqr)
   // - starting with the lowest bin (upper/lower boundary pair), iterate 
   //    through the sorted indexed_vector_magnitude elements, if 
   //    their v_mag_sqr exceeds the bin's upper boundary, move to the next
   //    bin.
   //    If v_mag_sqr is between the upper and lower boundaries 
   //    (squared) then increment data1D_local[i] by psi[j+i*Ny] 

   vector<indexed_vector_magnitude_sqr> indexed_magnitudes;

   for ( ptrdiff_t i=0; i<Nx_local; ++i)
      for ( ptrdiff_t j=0; j<Ny; ++j)
      { //if ( kx[i] * kx[i] + ky[j] * ky[j] == 16.0 ) // debug
         //{ // debug
         //   cout << "kx[" << i << "]^2 + ky[" << j << "]^2 == 16 == " 
         //      << kx[i] << "^2 + " << ky[j] << "^2 " << endl;
         //} // debug
         indexed_magnitudes.push_back(
               indexed_vector_magnitude_sqr( 
                     i, j,
                     kx_local[i] * kx_local[i] + ky[j] * ky[j] // |k|^2
                     // |k|^{2}, not |k|
                  )
               );
      }  // size of indexed_vector_magnitude : Nx_local * Ny

   // sort the indexed_magnitudes by the magnitude of their |k|^2 values
   std::sort(
         indexed_magnitudes.begin(),
         indexed_magnitudes.end(),
         indexed_vector_magnitude_sqr_lt   // pointer to "<" function 
         );
   //cout << "min |k|^2 of indexed_magnitudes : "  // debug
   //   << indexed_magnitudes.front().v_mag_sqr << endl; // debug
   //cout << "max |k|^2 of indexed_magnitudes : "  // debug
   //   << indexed_magnitudes.back().v_mag_sqr << endl; // debug

   //cout << "max of binning_boundaries : " 
   //   << binning_boundaries.back() << endl; // debug

   // Iterate over psi using indices sorted by corresponding domain 
   //  magnitudes,
   // sum the values of psi into bins having domain spacings delta_k
   //

   // - for all indexed_vector_magnitude elements,
   //    binning lower boundaries, 
   //    - while mag_itr->mag_sqr remains less than 
   //       (boundary + delta_k)^{2}
   //       - add psi[ mag_itr->j + mag_itr->i * Ny ] 
   //       - increment to the next indexed_vector_magnitude


   // zero the data1D_local[]
   for ( size_t i=0; i < binning_boundaries.size() - 1; ++i) 
   {
      data1D_local[i] = 0.0;
      bin_counts[i] = 0;
   }

   // - starting with the lowest bin (upper/lower boundary pair), iterate 
   //    through the sorted indexed_vector_magnitude elements, if 
   //    their v_mag_sqr exceeds the bin's upper boundary, move to the next
   //    bin.
   //    If v_mag_sqr is between the upper and lower boundaries 
   //    (squared) then increment data1D_local[i] by psi[j+i*Ny] 

   std::vector<indexed_vector_magnitude_sqr>::iterator 
         mag_itr = indexed_magnitudes.begin(); // the values to be binned

   size_t number_of_points_binned = 0;
   //double bin_element_count; // number of elements encountered in one bin

   //double lower_bound = binning_boundaries.front();
   double upper_bound;
   size_t ii=0; // index of bin ( data1D_local[ii] )

   double tmp_re, tmp_im;

   for ( std::vector<double>::const_iterator
        binning_boundary_itr = ++(binning_boundaries.begin());
        binning_boundary_itr != binning_boundaries.end();
        ++binning_boundary_itr )
   {
      upper_bound = *binning_boundary_itr;

      //cout << "binning between : "  // debug
      //   << lower_bound << ", " << upper_bound << endl; // debug

      if ( ii >= binning_boundaries.size() )
      {
         cerr << "Error: integrate_out_theta() oob for data1D_local[]"
            << endl;
         return EXIT_FAILURE;
      }

      if ( mag_itr == indexed_magnitudes.end() )
      {
         // This situation might occur if psi is split over nodes by MPI,
         //  in which case kx_local will be smaller than Nx.
         cerr << "Warning : integrate_out_theta() hit end of "  // debug
            << "indexed_magnitudes before hitting the end of " // debug
            << "binning_boundaries" // debug
            << endl; // debug
         cerr << "indexed_magnitudes.size(), binning_boundaries.size(), "
           << "number_of_points_binned : " 
            << indexed_magnitudes.size() << ", "
            << binning_boundaries.size() << ", "
            << number_of_points_binned
            << endl;
         cerr << "kx_local[Nx_local/2 -1], ky[Ny/2 -1] : "
            << kx_local[Nx_local/2 -1] << ", " <<  ky[Ny/2 -1] << endl;
         cerr << "current upper_bound : " << upper_bound << endl;
         cerr << "final upper_bound : " 
            << binning_boundaries.back() << endl;
         break;
      }

      //bin_element_count = 0;

      for ( ; 
            mag_itr->v_mag_sqr <= (upper_bound  * upper_bound)
            && 
            mag_itr != indexed_magnitudes.end();
            ++mag_itr )
      {
         //cout << "ii, mag_itr->v_mag_sqr, upper_bound^2, " // debug
         //   << "lower_bound^2 : " // debug
         //   << ii << ", " << mag_itr->v_mag_sqr << ", " 
         //   << upper_bound * upper_bound // debug
         //   << ", " << lower_bound * lower_bound << endl;  // debug


         tmp_re = psi[ (mag_itr->j) + Ny * (mag_itr->i) ][0];
         tmp_im = psi[ (mag_itr->j) + Ny * (mag_itr->i) ][1];
         data1D_local[ii] 
            += sqrt( (tmp_re * tmp_re) + (tmp_im * tmp_im) );
         ++number_of_points_binned; // debug
         bin_counts[ii] += 1;
         //++bin_element_count; 

      }

      ++ii;
      //lower_bound = upper_bound;// debug lower_bound isn't useful otherwise
   }

//   // Calculate the mean and standard deviation of the incremental 
//   //  difference of data1D_local[]
//   double tmp_prev=data1D_local[0];
//   double avg= 0.0;//tmp_prev;
//   double avg_sqr = 0.0;// = avg * avg;
//   for ( size_t i=1; i < binning_boundaries.size() - 1; ++i)
//   {
//     avg += (data1D_local[i] - tmp_prev); 
//     // cout << "data1D_local[" << i << "] : " 
//     //    << data1D_local[i]  ;//debug
//     // cout << " avg, tmp_prev : " << avg << ", " << tmp_prev // debug
//     //    << endl; // debug
//     if ( tmp_prev > data1D_local[i] )
//     {
//        cout << "decrease in number of points binned" << endl; // debug
//     }
//     avg_sqr += (data1D_local[i] - tmp_prev) 
//                  * (data1D_local[i] - tmp_prev);
//     tmp_prev = data1D_local[i]; 
//   }
//
//   cout << "sum : " << avg << endl;//debug
//   avg = avg / (binning_boundaries.size() - 1);
//   cout << "avg : " << avg << endl;//debug
//   avg_sqr = avg_sqr / (binning_boundaries.size() - 1);
//   cout << "avg_sqr : " << avg_sqr << endl;//debug
//
//
//   double deviation = sqrt( avg_sqr - avg * avg );
//
//   cout << "Average change : " 
//      << avg // /(binning_boundaries.size() - 1) // debug
//      << " +/- " << deviation // /(binning_boundaries.size() - 1) // debug
//      << endl;// debug
//
//   cout << "Average change relative to number of points binned : " 
//      << avg / number_of_points_binned 
//      // /(binning_boundaries.size() - 1) // debug
//      << " +/- " << deviation / number_of_points_binned // debug
//      << endl;// debug

   return EXIT_SUCCESS;
}


// TODO: change the name to portray that it averages instead of integrates
int TEM_NS::integrate_out_theta_double( 
      // integrate to make f( \vec{k} ) into f( | \vec{k} | )
      //const double** const psi,// data to be binned,[][0] real,[][1] im
      const double* const psi,//data to be integrated azimuthally
      const double* const kx_local, const size_t& Nx_local,
      const double* const ky, const size_t& Ny,
      const std::vector<double>& binning_boundaries,// bin boundaries
      double* data1D_local // output, 
             // having binning_boundaries.size() - 1 allocated elements
      )
{
   // Precondition : 
   // // - data1D_local has been allocated to have 
   // //    number_of_bins == bwcutoff_t / delta_k  
   // //       == binning_boundaries.size() - 1
   // //    elements
   // - binning boundaries do not exceed the maximum or minimum of the
   //    domain kx, ky of the data to be binned, but the non-zero data may 
   //    lie entirely within the bounds of the binning
   
   // Plan:
   // - create a vector of indexed_vector_magnitude elements, each of
   //    which contains array indices and a magnitude equal to the sum of 
   //    squares of those indices
   // - sort them by their magnitude member (v_mag_sqr)
   // - starting with the lowest bin (upper/lower boundary pair), iterate 
   //    through the sorted indexed_vector_magnitude elements, if 
   //    their v_mag_sqr exceeds the bin's upper boundary, move to the next
   //    bin.
   //    If v_mag_sqr is between the upper and lower boundaries 
   //    (squared) then increment data1D_local[i] by psi[j+i*Ny] 

   // TODO: Revision: established FTEM methods average azimuthally instead 
   //       of integrate azimuthally, so this function must be changed.

   vector<indexed_vector_magnitude_sqr> indexed_magnitudes;
   for ( ptrdiff_t i=0; i<Nx_local; ++i)
      for ( ptrdiff_t j=0; j<Ny; ++j)
      {
         //if ( kx[i] * kx[i] + ky[j] * ky[j] == 16.0 ) // debug
         //{ // debug
         //   cout << "kx[" << i << "]^2 + ky[" << j << "]^2 == 16 == " 
         //      << kx[i] << "^2 + " << ky[j] << "^2 " << endl;
         //} // debug
         indexed_magnitudes.push_back(
               indexed_vector_magnitude_sqr( 
                     i, j,
                     kx_local[i] * kx_local[i] + ky[j] * ky[j] // |k|^2
                     // |k|^{2}, not |k|
                  )
               );
      }  // size of indexed_vector_magnitude : Nx_local * Ny

   // sort the indexed_magnitudes by the magnitude of their |k|^2 values
   std::sort(
         indexed_magnitudes.begin(),
         indexed_magnitudes.end(),
         indexed_vector_magnitude_sqr_lt   // pointer to "<" function 
         );
   //cout << "min |k|^2 of indexed_magnitudes : "  // debug
   //   << indexed_magnitudes.front().v_mag_sqr << endl; // debug
   //cout << "max |k|^2 of indexed_magnitudes : "  // debug
   //   << indexed_magnitudes.back().v_mag_sqr << endl; // debug

   //cout << "max of binning_boundaries : " 
   //   << binning_boundaries.back() << endl; // debug

   // Iterate over psi using indices sorted by corresponding domain 
   //  magnitudes,
   // sum the values of psi into bins having domain spacings delta_k
   //

   // - for all indexed_vector_magnitude elements,
   //    binning lower boundaries, 
   //    - while mag_itr->mag_sqr remains less than 
   //       (boundary + delta_k)^{2}
   //       - add psi[ mag_itr->j + mag_itr->i * Ny ] 
   //       - increment to the next indexed_vector_magnitude


   // zero the data1D_local[]
   for ( size_t i=0; i < binning_boundaries.size() - 1; ++i) 
      data1D_local[i] = 0.0;

   // - starting with the lowest bin (upper/lower boundary pair), iterate 
   //    through the sorted indexed_vector_magnitude elements, if 
   //    their v_mag_sqr exceeds the bin's upper boundary, move to the next
   //    bin.
   //    If v_mag_sqr is between the upper and lower boundaries 
   //    (squared) then increment data1D_local[i] by psi[j+i*Ny] 

   std::vector<indexed_vector_magnitude_sqr>::iterator 
         mag_itr = indexed_magnitudes.begin(); // the values to be binned

   size_t number_of_points_binned = 0;
   double bin_element_count;

   //double lower_bound = binning_boundaries.front();
   double upper_bound;
   size_t ii=0; // index of bin ( data1D_local[ii] )

   double tmp_re, tmp_im;


   for ( std::vector<double>::const_iterator
        binning_boundary_itr = ++(binning_boundaries.begin());
        binning_boundary_itr != binning_boundaries.end();
        ++binning_boundary_itr )
   {
      upper_bound = *binning_boundary_itr;

      //cout << "binning between : "  // debug
      //   << lower_bound << ", " << upper_bound << endl; // debug

      if ( ii >= binning_boundaries.size() )
      {
         cerr << "Error: integrate_out_theta() oob for data1D_local[]"
            << endl;
         return EXIT_FAILURE;
      }

      if ( mag_itr == indexed_magnitudes.end() )
      {
         // This situation might occur if psi is split over nodes by MPI,
         //  in which case kx_local will be smaller than Nx.
         cerr << "Warning : integrate_out_theta() hit end of "  // debug
            << "indexed_magnitudes before hitting the end of " // debug
            << "binning_boundaries" // debug
            << endl; // debug
         cerr << "indexed_magnitudes.size(), binning_boundaries.size(), "
           << "number_of_points_binned : " 
            << indexed_magnitudes.size() << ", "
            << binning_boundaries.size() << ", "
            << number_of_points_binned
            << endl;
         cerr << "kx_local[Nx_local/2 -1], ky[Ny/2 -1] : "
            << kx_local[Nx_local/2 -1] << ", " <<  ky[Ny/2 -1] << endl;
         cerr << "current upper_bound : " << upper_bound << endl;
         cerr << "final upper_bound : " 
            << binning_boundaries.back() << endl;
         break;
      }

      bin_element_count = 0.0; 

      for ( ; 
            mag_itr->v_mag_sqr <= (upper_bound  * upper_bound)
            && 
            mag_itr != indexed_magnitudes.end();
            ++mag_itr )
      {
         //cout << "ii, mag_itr->v_mag_sqr, upper_bound^2, " // debug
         //   << "lower_bound^2 : " // debug
         //   << ii << ", " << mag_itr->v_mag_sqr << ", " 
         //   << upper_bound * upper_bound // debug
         //   << ", " << lower_bound * lower_bound << endl;  // debug

         ++number_of_points_binned; // debug
         ++bin_element_count;

         data1D_local[ii] += psi[ (mag_itr->j) + Ny * (mag_itr->i) ];
      }

      // implement averaging within the current bin
      if ( bin_element_count != 0)
         data1D_local[ii] = data1D_local[ii] / bin_element_count;

      ++ii;
      //lower_bound = upper_bound;// debug lower_bound isn't useful otherwise
   }

//   // Calculate the mean and standard deviation of the incremental 
//   //  difference of data1D_local[]
//   double tmp_prev=data1D_local[0];
//   double avg= 0.0;//tmp_prev;
//   double avg_sqr = 0.0;// = avg * avg;
//   for ( size_t i=1; i < binning_boundaries.size() - 1; ++i)
//   {
//     avg += (data1D_local[i] - tmp_prev); 
//     // cout << "data1D_local[" << i << "] : " 
//     //    << data1D_local[i]  ;//debug
//     // cout << " avg, tmp_prev : " << avg << ", " << tmp_prev // debug
//     //    << endl; // debug
//     if ( tmp_prev > data1D_local[i] )
//     {
//        cout << "decrease in number of points binned" << endl; // debug
//     }
//     avg_sqr += (data1D_local[i] - tmp_prev) 
//                  * (data1D_local[i] - tmp_prev);
//     tmp_prev = data1D_local[i]; 
//   }
//
//   cout << "sum : " << avg << endl;//debug
//   avg = avg / (binning_boundaries.size() - 1);
//   cout << "avg : " << avg << endl;//debug
//   avg_sqr = avg_sqr / (binning_boundaries.size() - 1);
//   cout << "avg_sqr : " << avg_sqr << endl;//debug
//
//
//   double deviation = sqrt( avg_sqr - avg * avg );
//
//   cout << "Average change : " 
//      << avg // /(binning_boundaries.size() - 1) // debug
//      << " +/- " << deviation // /(binning_boundaries.size() - 1) // debug
//      << endl;// debug
//
//   cout << "Average change relative to number of points binned : " 
//      << avg / number_of_points_binned 
//      // /(binning_boundaries.size() - 1) // debug
//      << " +/- " << deviation / number_of_points_binned // debug
//      << endl;// debug

   return EXIT_SUCCESS;
}




int TEM_NS::variance_1D_STEM( 
      const double* const diffracted_wave_radial_intensity_sum,
      const double* const diffracted_wave_radial_intensity_sqr_sum,
      //const size_t& number_of_bins,
      const std::vector<double>& binning_boundaries, // for output
      const size_t& number_of_raster_points,
      const unsigned int& input_flag_netcdf_variance,
      const string& outFilePrefix,
      double* variance
      )
{
   size_t number_of_bins; number_of_bins = binning_boundaries.size() - 1;
   // debug
   //cout << "variance_1D_STEM() number_of_raster_points : " 
   //   << number_of_raster_points
   //   << ", number_of_bins : " << number_of_bins
   //   << endl;
   // end debug

   for ( size_t idx=0; idx < number_of_bins; ++idx)
   {
      if (
            diffracted_wave_radial_intensity_sum[idx] == 0.0
            //||
            //diffracted_wave_radial_intensity_sqr_sum[idx] == 0.0
         )
      {
         variance[idx] = 0.0;
      }
      else
      {
         variance[idx] 
            = ( 
                  number_of_raster_points 
                  *  
                  diffracted_wave_radial_intensity_sqr_sum[idx]
                     / ( 
                        diffracted_wave_radial_intensity_sum[idx]
                        * diffracted_wave_radial_intensity_sum[idx]
                       )
               ) - 1.0;
      }

      if ( variance[idx] < 0.0 )
      {
         cout << " WARNING, variance_1D_STEM, variance[" 
            << idx << "] < 0.0 : "
            << variance[idx] << endl;
      }
   }

   // debug
   if( input_flag_netcdf_variance )
   {
      output_variance_to_netcdf(
         diffracted_wave_radial_intensity_sqr_sum,
         binning_boundaries,
         outFilePrefix + "_radial_intensity_sqr_sum"
         );
      output_variance_to_netcdf(
         diffracted_wave_radial_intensity_sum,
         binning_boundaries,
         outFilePrefix + "_radial_intensity_sum"
         );
   } else {
      output_variance_to_txt(
         diffracted_wave_radial_intensity_sqr_sum,
         binning_boundaries,
         outFilePrefix + "_radial_intensity_sqr_sum"
         );
      output_variance_to_txt(
         diffracted_wave_radial_intensity_sum,
         binning_boundaries,
         outFilePrefix + "_radial_intensity_sum"
         );
   }
   //end debug

   // write the variance to a file
   if( input_flag_netcdf_variance )
   {
      if ( 
            output_variance_to_netcdf(
               variance,
               binning_boundaries,
               outFilePrefix + "_radial_intensity_variance"
            ) != EXIT_SUCCESS)
      {
         cout << " failed to write variance data to netCDF file : " 
            << outFilePrefix << endl;
      }
   } else {
      if ( 
         output_variance_to_txt(
               variance,
               binning_boundaries,
               outFilePrefix + "_radial_intensity_variance"
            ) != EXIT_SUCCESS)
      {
         cout << " failed to write variance data to a txt file : " 
            << outFilePrefix << endl;
      }
   }
}

// TODO: rewrite varianc_2D_STEM() to work within revised adfstem()
int TEM_NS::variance_2D_STEM( 
      const double* const data2D_avgs_local,       // denominator
      const double* const data2D_sqr_avgs_local,   // numerator
      const unsigned int number_of_raster_points,
      //const std::vector<double>& binning_boundaries,// bin boundaries
      const double* const kx_local, const size_t& Nx_local, 
      const double* const kx_joined, // non-null address only on rootnode
      const size_t& Nx_joined, 
      const double* const ky, const size_t& Ny, 
      const size_t& resolutionUnit_recip,
      const double& xResolution_recip, const double& yResolution_recip,
      const double& delta_k,
      const double& bwcutoff_t, // data known to be 0 beyond this |k|
      const string& outFilePrefix,
      const int& mynode,
      const int& rootnode,
      MPI_Comm comm
      )
{
   //  data2D_avgs_local : average over raster position of intensity
   //  data2D_sqr_avgs_local : average over raster position of intensity^{2}


   // NOTE: The following is not appropriate for V(|k|), only for 
   //        V( \vec{k} )
   //   - integrate_out_theta() for both data2D_avgs_local and 
   //      data2D_sqr_avgs_local
   //   - assign the variance of each pixel 
   //      variance[i] 
   //         = data2D_avgs_local[i]^2 / data2D_sqr_avgs_local[i] - 1

   // NOTE: appropriate plan
   //    - determine binning boundaries on the root node and broadcast them
   //    - integrate w.r.t. \theta in polar coordinates to yield 1-D 
   //       data1D_sqr_avgs_local and data1D_avgs_local
   //    - MPI_Reduce() to sum each 
   //    - calculate variance from data1D_avgs_reduced and 
   //       data1D_sqr_avgs_reduced

   ///////////////////////////////////////////////////////////////////
   // TODO : pull this block out and put it into tem()
   // Assign binning boundary values for the integration
   ///////////////////////////////////////////////////////////////////
   //size_t number_of_bins = binning_boundaries.size();
   //double* variance;
   //if( mynode == rootnode ) variance 
   //   = new double[ number_of_bins ];
   ///////////////////////////////////////////////////////////////////
   // end of block to pull out and into tem()
   ///////////////////////////////////////////////////////////////////
   
   ///////////////////////////////////////////////////////////////////
   // This function is only used for creating a 2-D image of the
   //  variance of diffraction.
   ///////////////////////////////////////////////////////////////////
   double* diffracted_wave_mag_variance_2D;
   diffracted_wave_mag_variance_2D = new double[ Nx_local * Ny ];

   //debug  check the input to see if the numerator < denominator
   for ( ptrdiff_t i=0; i< Nx_local * Ny; ++i)
   {
      //diffracted_wave_mag_variance_2D[i]
      //   =  data2D_sqr_avgs_local[i] 
      //         - ( data2D_avgs_local[i] * data2D_avgs_local[i]);
      //if ( diffracted_wave_mag_variance_2D[i] < 0.0 )
      //{
      //   cout << " node " << mynode << ", variance_2D() difference : " 
      //       << "diffracted_wave_mag_variance_2D[" << i << "] : "
      //      << diffracted_wave_mag_variance_2D[i]
      //       << " < 0.0 " 
      //       << endl;// debug
      //   diffracted_wave_mag_variance_2D[i] = 0.0;
      //}
      ////if ( diffracted_wave_mag_variance_2D[i] == 0.0 )
      ////{
      ////   diffracted_wave_mag_variance_2D[i] = 1.0;
      ////}
      ////else 
      //if ( diffracted_wave_mag_variance_2D[i] != 0.0 )
      //   if ( data2D_avgs_local[i] != 0.0 )
      //   {
      //      diffracted_wave_mag_variance_2D[i]
      //         = diffracted_wave_mag_variance_2D[i]
      //               / ( data2D_avgs_local[i] * data2D_avgs_local[i]);
      //   }
      //   else
      //   {
      //      diffracted_wave_mag_variance_2D[i] = 0.0;
      //   }
      // rewriting the above
      if ( data2D_avgs_local[i] == 0.0 ) 
         diffracted_wave_mag_variance_2D[i] = 0.0;
      else
         diffracted_wave_mag_variance_2D[i] 
            = (
                  number_of_raster_points * data2D_sqr_avgs_local[i] / 
                  (data2D_avgs_local[i] * data2D_avgs_local[i])
              ) - 1.0;
      if ( diffracted_wave_mag_variance_2D[i] < 0.0 ) 
         cout << " node " << mynode << ", variance_2D_STEM() variance["
            << i << "] < 0.0 : " << diffracted_wave_mag_variance_2D[i] 
            << endl; // debug

   } // end debug
   output_diffraction(
      diffracted_wave_mag_variance_2D,
      1.0e-40,
      Nx_local * Ny,
      Nx_local, Nx_joined, Ny,
      resolutionUnit_recip,
      xResolution_recip, yResolution_recip,
      outFilePrefix + "_variance2D",
      mynode, rootnode, comm
      );

   output_diffraction(
      diffracted_wave_mag_variance_2D,
      1.0e-30,
      Nx_local * Ny,
      Nx_local, Nx_joined, Ny,
      resolutionUnit_recip,
      xResolution_recip, yResolution_recip,
      outFilePrefix + "_variance2D",
      mynode, rootnode, comm
      );

   output_diffraction(
      diffracted_wave_mag_variance_2D,
      1.0e-20,
      Nx_local * Ny,
      Nx_local, Nx_joined, Ny,
      resolutionUnit_recip,
      xResolution_recip, yResolution_recip,
      outFilePrefix + "_variance2D",
      mynode, rootnode, comm
      );

   output_diffraction(
      diffracted_wave_mag_variance_2D,
      1.0e-15,
      Nx_local * Ny,
      Nx_local, Nx_joined, Ny,
      resolutionUnit_recip,
      xResolution_recip, yResolution_recip,
      outFilePrefix + "_variance2D",
      mynode, rootnode, comm
      );

   output_diffraction(
      diffracted_wave_mag_variance_2D,
      1.0e-10,
      Nx_local * Ny,
      Nx_local, Nx_joined, Ny,
      resolutionUnit_recip,
      xResolution_recip, yResolution_recip,
      outFilePrefix + "_variance2D",
      mynode, rootnode, comm
      );

   output_diffraction(
      diffracted_wave_mag_variance_2D,
      1.0e-5,
      Nx_local * Ny,
      Nx_local, Nx_joined, Ny,
      resolutionUnit_recip,
      xResolution_recip, yResolution_recip,
      outFilePrefix + "_variance2D",
      mynode, rootnode, comm
      );

   output_diffraction(
      diffracted_wave_mag_variance_2D,
      1.0e-2,
      Nx_local * Ny,
      Nx_local, Nx_joined, Ny,
      resolutionUnit_recip,
      xResolution_recip, yResolution_recip,
      outFilePrefix + "_variance2D",
      mynode, rootnode, comm
      );

   output_diffraction(
      diffracted_wave_mag_variance_2D,
      1.0e-1,
      Nx_local * Ny,
      Nx_local, Nx_joined, Ny,
      resolutionUnit_recip,
      xResolution_recip, yResolution_recip,
      outFilePrefix + "_variance2D",
      mynode, rootnode, comm
      );

   output_diffraction(
      diffracted_wave_mag_variance_2D,
      1.0,
      Nx_local * Ny,
      Nx_local, Nx_joined, Ny,
      resolutionUnit_recip,
      xResolution_recip, yResolution_recip,
      outFilePrefix + "_variance2D",
      mynode, rootnode, comm
      );

   output_diffraction(
      diffracted_wave_mag_variance_2D,
      10.0,
      Nx_local * Ny,
      Nx_local, Nx_joined, Ny,
      resolutionUnit_recip,
      xResolution_recip, yResolution_recip,
      outFilePrefix + "_variance2D",
      mynode, rootnode, comm
      );
   //delete[] diffracted_wave_mag_variance_1D_local; 
   delete[] diffracted_wave_mag_variance_2D;
   ///////////////////////////////////////////////////////////////////
   // TODO : end of portion to eliminate when no longer desiring 2-D 
   //       image of variance of the diffraction.
   ///////////////////////////////////////////////////////////////////
   
   /////////////////////////////////////////////////////////////////////
   //// Accumulate avg values into bins according to the magnitude of
   ////  their domain vectors. 
   ////  (integrating out the angular dependence of |<I(\vec{k})>| to
   ////   get a function of |\vec{k}|
   /////////////////////////////////////////////////////////////////////
   ////
   ////
   ////
   //// TODO: eliminate the following variance_1D calculation
   ////
   ////
   ////double* diffracted_wave_mag_variance_1D_local;
   ////diffracted_wave_mag_variance_1D_local = new double[ number_of_bins ];
   //double* data1D_sqr_avgs_local = new double[number_of_bins];
   //double* data1D_avgs_local = new double[number_of_bins];
   //double* data1D_sqr_avgs_reduced;
   //double* data1D_avgs_reduced;
   ////double* variance;
   ////double* diffracted_wave_mag_variance_1D_reduced;
   //if ( mynode == rootnode )
   //{
   //   //diffracted_wave_mag_variance_1D_reduced = new double[number_of_bins];
   //   data1D_avgs_reduced = new double[number_of_bins];
   //   data1D_sqr_avgs_reduced = new double[number_of_bins];
   //   //variance = new double[number_of_bins];
   //}
   //if (
   //      integrate_out_theta_double(
   //            data2D_avgs_local,
   //            kx_local, Nx_local,
   //            ky, Ny,
   //            binning_boundaries,
   //            data1D_avgs_local//,
   //            //number_of_bins
   //            )
   //            != EXIT_SUCCESS
   //   )
   //{
   //   cerr << "Error, node " << mynode << ": integrate_out_theta() failed" 
   //      << endl;
   //   return EXIT_FAILURE;
   //}
   //if (
   //      integrate_out_theta_double(
   //            data2D_sqr_avgs_local,
   //            kx_local, Nx_local,
   //            ky, Ny,
   //            binning_boundaries,
   //            data1D_sqr_avgs_local//,
   //            //number_of_bins
   //            )
   //            != EXIT_SUCCESS
   //   )
   //{
   //   cerr << "Error, node " << mynode << ": integrate_out_theta() failed" 
   //      << endl;
   //   return EXIT_FAILURE;
   //}

   //MPI_Reduce( 
   //            data1D_avgs_local,
   //            data1D_avgs_reduced,
   //            //diffracted_wave_mag_variance_1D_local, 
   //            //diffracted_wave_mag_variance_1D_reduced,
   //            number_of_bins,
   //            MPI_DOUBLE, MPI_SUM,
   //            rootnode, comm );

   //MPI_Reduce( 
   //            data1D_sqr_avgs_local,
   //            data1D_sqr_avgs_reduced,
   //            //diffracted_wave_mag_variance_1D_local, 
   //            //diffracted_wave_mag_variance_1D_reduced,
   //            number_of_bins,
   //            MPI_DOUBLE, MPI_SUM,
   //            rootnode, comm );

   //if ( mynode == rootnode )
   //{
   //   // calculate variance
   //   for ( unsigned int idx=0; idx < number_of_bins; ++idx)
   //   {
   //      //if (
   //      //      data1D_avgs_reduced[idx]   // denominator
   //      //      >
   //      //      data1D_sqr_avgs_reduced[idx]  // numerator
   //      //   )
   //      //{
   //      //   cout << " WARNING, variance_2D_STEM: negative variance , "
   //      //      << "data1D_avgs_reduced[" << idx 
   //      //      << "] > data1D_sqr_avgs_reduced[" << idx << "]" 
   //      //      << endl
   //      //      << "  setting variance[" << idx << "] = 0.0" << endl;

   //      //   variance[idx] = 0.0;
   //      //}
   //      //else
   //      if ( 
   //            data1D_sqr_avgs_reduced[idx] == 0.0
   //            ||
   //            data1D_avgs_reduced[idx] == 0.0
   //         )
   //      {
   //         variance[idx] = 0.0;
   //      }
   //      else
   //      {
   //         variance[idx] 
   //            = (
   //               data1D_sqr_avgs_reduced[idx]
   //               / ( data1D_avgs_reduced[idx] * data1D_avgs_reduced[idx] )
   //              )
   //               - 1.0;
   //      }
   //      //if ( data1D_avgs_reduced[idx] > 0.0 )
   //      //{
   //      //   variance[idx] 
   //      //      = variance[idx] 
   //      //         / (data1D_avgs_reduced[idx] * data1D_avgs_reduced[idx]);
   //      //}
   //      if ( variance[idx] < 0.0 )
   //      {
   //         cout << " WARNING, variance_2D_STEM, variance[" 
   //            << idx << "] < 0.0 : "
   //            << variance[idx] << endl;
   //      }
   //   }
   //   //debug
   //   output_variance_to_netcdf(
   //      data1D_avgs_reduced,
   //      //diffracted_wave_mag_variance_1D_reduced,
   //      binning_boundaries, // domain binning boundaries 
   //      outFilePrefix + "_data1D_avgs_reduced"
   //      );
   //   output_variance_to_netcdf(
   //      data1D_sqr_avgs_reduced,
   //      //diffracted_wave_mag_variance_1D_reduced,
   //      binning_boundaries, // domain binning boundaries 
   //      outFilePrefix + "_data1D_sqr_avgs_reduced"
   //      );
   //   //end debug
   //   // write the variance to a netcdf file
   //   if ( 
   //         output_variance_to_netcdf(
   //            variance,
   //            //diffracted_wave_mag_variance_1D_reduced,
   //            binning_boundaries, // domain binning boundaries 
   //            outFilePrefix + "_variance1D"
   //         ) != EXIT_SUCCESS
   //      )
   //   {
   //      cout << " failed to write variance data to netCDF file : " 
   //         << outFilePrefix << endl;
   //   }
   //   //delete[] diffracted_wave_mag_variance_1D_reduced;
   //}

   //// end revised approach
   

   //if (mynode == rootnode)
   //{
   //   delete[] data1D_avgs_reduced;
   //   delete[] data1D_sqr_avgs_reduced;

   //   ///////////////////////////////////////////////////////////////////
   //   // TODO: pull out and into tem()
   //   delete[] variance;
   //   ///////////////////////////////////////////////////////////////////
   //}

   //delete[] data1D_avgs_local;
   //delete[] data1D_sqr_avgs_local;

   return EXIT_SUCCESS;
}



int TEM_NS::variance_2D_BFCTEM( 
      double* data2D_avgs_local,
      double* data2D_sqr_avgs_local,
      const radial_discretization& sample_rotations,
      const string& outFilePrefix,
      const int& mynode,
      const int& rootnode,
      MPI_Comm comm
      )
{
   double* variance;
   double* data1D_avgs_reduced;
   double* data1D_sqr_avgs_reduced ;
   size_t number_of_q = sample_rotations.size();

   if ( mynode == rootnode )
   {
      data1D_avgs_reduced = new double[ number_of_q ];
      data1D_sqr_avgs_reduced = new double[ number_of_q ];
      variance = new double[ number_of_q ];
   }
   
   MPI_Reduce( data2D_avgs_local, data1D_avgs_reduced, 
               number_of_q, MPI_DOUBLE, MPI_SUM, 
               rootnode, comm);
   MPI_Reduce( data2D_sqr_avgs_local, data1D_sqr_avgs_reduced, 
               number_of_q, MPI_DOUBLE, MPI_SUM, 
               rootnode, comm);


   if ( mynode == rootnode )
   {

      for( size_t i=0; i<number_of_q; ++i)
      {
         if ( data1D_avgs_reduced[i] != 0.0 )
         {
            variance[i] 
               = 
               (
                data1D_sqr_avgs_reduced[i] 
                / ( data1D_avgs_reduced[i] * data1D_avgs_reduced[i] )
               ) - 1.0;
         }
         else
         {
            variance[i] = 0.0;
         }
      }

      // output_variance_to_netcdf() requires a list<double> of binning
      //  boundaries
      //std::list<double> binning_boundaries;
      //for( radial_discretization::const_iterator 
      //      rad_itr = sample_rotations.begin();
      //      rad_itr != sample_rotations.end();
      //      ++rad_itr)
      //{
      //   binning_boundaries.push_back( rad_itr->radius );
      //}

      if ( output_variance_to_netcdf(
               variance,
               sample_rotations,
               outFilePrefix + "_avgs"
               ) != EXIT_SUCCESS)
      {
         cout << "failed to write variance data to netCDF file : "
            << outFilePrefix << endl;
      }
   }

   if ( mynode == rootnode )
   {
      delete[] data1D_avgs_reduced;
      delete[] data1D_sqr_avgs_reduced;
      //delete[] variance;
   }

   return EXIT_SUCCESS;
}
#endif
