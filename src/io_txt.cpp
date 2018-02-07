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
// File: io_txt.cpp
// Purpose: 

#ifndef IO_TXT_CPP
#define IO_TXT_CPP

#include "io_txt.hpp"

int TEM_NS::output_variance_to_txt(
      const double* const outDataRaw,  // array of range elements
      const std::vector<double>& binning_boundaries, // domain bdy elements
      const string& outFilePrefix
      )
{
   // outDataRaw : f_of_k_magnitudes_reduced[]
   // Nx : number of elements in f_of_k_magnitudes_reduced[]


   //double myDomain[ Nx ];
   //std::vector<double>::const_iterator bdy_itr = binning_boundaries.begin();
   //for ( size_t i=0; i < Nx ; ++i)
   //{
   //   myDomain[i] = *bdy_itr;
   //   ++bdy_itr;
   //}

   // Create the txt file, if the file exists it will be overwritten
   string outFileName = outFilePrefix + ".txt";
   std::vector<double>::const_iterator 
      bdy_itr = binning_boundaries.begin();

   // open the output file
   std::ofstream outputFile;
   outputFile.open( outFileName.c_str(), std::ios::out);
   if ( outputFile.is_open() && outputFile.good()) {

      //outputFile.setf( std::scientific, std::ios::floatfield);
      // alternatively use 
      outputFile << std::scientific;
      //outputFile.precision(10);
      // alternatively use the following with #include <iomanip> :
      outputFile << std::setprecision(10); 

      size_t Nx = binning_boundaries.size() - 1;
      // TODO: consider writing a header line
      for ( size_t i=0; i < Nx ; ++i)
      {
         if( outputFile.good() )
            outputFile << *bdy_itr << " " << outDataRaw[i] << "\n";
         else
         {
            cerr << "Failure while writing variance to file: " 
               << outFileName 
               << ", outDataRaw[" << i << "] : " << outDataRaw[i] 
               << ", *bdy_itr : " << *bdy_itr << endl;
            return EXIT_FAILURE;
         }
         ++bdy_itr;
      }

      outputFile.close();
   }
   else
   {
      cerr << "Unable to write variance to file: " << outFileName << endl;
      return EXIT_FAILURE;
   }

   //cout << "Variance data written to : " << outFilePrefix + ".txt" << endl;

   return EXIT_SUCCESS;
}


//int TEM_NS::output_variance_to_txt(
//      const double* const outDataRaw,  // array of range elements
//      const radial_discretization& sample_rotations, // domain bdy elements
//      const string& outFilePrefix
//      )
//{
//   // outDataRaw : f_of_k_magnitudes_reduced[]
//   // Nx : number of elements in f_of_k_magnitudes_reduced[]
//   ///////////////////////////////////////////////////////////////////
//   // validate input
//   ///////////////////////////////////////////////////////////////////
//
//   size_t Nr = sample_rotations.size() ;
//
//   ///////////////////////////////////////////////////////////////////
//   // The netCDF libraries require the domain to be an array,
//   //  so copy the domain radius elements to an array
//   ///////////////////////////////////////////////////////////////////
//
//   double myDomain[ Nr ];
//   radial_discretization::const_iterator rad_itr = sample_rotations.begin();
//   for ( size_t i=0; i < Nr ; ++i)
//   {
//      myDomain[i] = rad_itr->radius;
//      ++rad_itr;
//   }
//
//   // Create the netCDF file, if the file exists it will be overwritten
//   //  due to the use of the Replace parameter
//   string outFileName = outFilePrefix + ".nc";
//   NcFile dataFile( outFileName.c_str(), NcFile::Replace );
//   if ( ! dataFile.is_valid() ) 
//   {
//      //cerr << "failed to open file : " << outFilePrefix + ".nc" << endl;
//      cerr << "failed to open file : " << outFileName << endl;
//      return EXIT_FAILURE;
//   }
//
//   // domain dimensions each require an NcDim and NcVar, in order to 
//   //  describe the domain values, units, and perhaps other attributes
//   NcDim* kxDimension = dataFile.add_dim("k", Nr);
//   NcVar* kxVariable = dataFile.add_var("k", ncDouble, kxDimension);
//   kxVariable->add_att("units", "\305^-1");
//
//   // Values which are functions of a domain each require an NcVar 
//   //  specifying the NcDim elements of their domain
//   NcVar* varianceVariable = dataFile.add_var("V", ncDouble, kxDimension );
//   varianceVariable->add_att("units", "(arb.)");
//
//   // write the variables to file
//   kxVariable->put( &myDomain[0], Nr );
//   varianceVariable->put( &outDataRaw[0], Nr );
//   //TODO: note that this direct output won't work for row/column -major
//   //       arrays
//
//   cout << "Variance data written to : " << outFilePrefix + ".nc" << endl;
//   cout << "Nr : " << Nr << endl; // debug
//
//   return EXIT_SUCCESS;
//}



//int TEM_NS::output_psi_realspace_to_txt(
//      const fftw_complex* const psi,   // local array of range elements
//      const ptrdiff_t& local_alloc_size_fftw, // psi size after MPI_Gather
//      const ptrdiff_t Nx_local,  // local size of psi in the x-dimension
//      const double* const kx_joined,   // kx domain elements
//      const ptrdiff_t Nx,     // number of elements in kx_joined
//      const double* const ky, // ky domain elements
//      const ptrdiff_t Ny,     // number of elements in ky
//      const string& outFilePrefix,
//      const int& mynode,
//      const int& rootnode,
//      MPI_Comm comm
//      )
//{
//
//   //double* psi_mag_joined; double* psi_re_joined; double* psi_im_joined;
//   //double* psi_mag_local; double* psi_re_local; double* psi_im_local;
//   //psi_mag_local = new double[ Nx_local ][ Ny ];
//   //psi_re_local = new double[ Nx_local ][ Ny ];
//   //psi_im_local = new double[ Nx_local ][ Ny ];
//
//   //double psi_mag_local[ Nx_local ][ Ny ]; 
//   //double psi_re_local[ Nx_local ][ Ny ]; 
//   //double psi_im_local[ Nx_local ][ Ny ];
//
//   double* psi_mag_local; 
//   psi_mag_local = new double[ local_alloc_size_fftw];
//   double* psi_re_local; 
//   psi_re_local = new double[ local_alloc_size_fftw];
//   double* psi_im_local; 
//   psi_im_local = new double[ local_alloc_size_fftw];
//
//
//   //cout << "output_psi_realspace_to_netcdf() line 159" << endl; // debug
//   //for (ptrdiff_t i=0; i < Nx_local; ++i)
//   //   for (ptrdiff_t j=0; j < Ny; ++j)
//   //   {
//   //      psi_mag_local[i][j] 
//   //         = sqrt(
//   //               psi[j + i * Ny][0] * psi[j + i * Ny][0]
//   //               + psi[j + i * Ny][1] * psi[j + i * Ny][1]
//   //            );
//   //      psi_re_local[i][j] = psi[j + i * Ny][0];
//   //      psi_im_local[i][j] = psi[j + i * Ny][1];
//   //   }
//   for (ptrdiff_t i=0; i < local_alloc_size_fftw; ++i)
//   {
//      psi_mag_local[i]
//         = sqrt(
//               psi[i][0] * psi[i][0] + psi[i][1] * psi[i][1]
//            );
//      psi_re_local[i] = psi[i][0];
//      psi_im_local[i] = psi[i][1];
//   }
//
//
//   //cout << "output_psi_realspace_to_netcdf() line 182" << endl; // debug
//   //double psi_mag_joined[ Nx ][ Ny ]; 
//   //double psi_re_joined[ Nx ][ Ny ]; 
//   //double psi_im_joined[ Nx ][ Ny ];
//   double* psi_mag_joined;
//   double* psi_re_joined; 
//   double* psi_im_joined;
//
//   if ( mynode == rootnode )
//   {
//      psi_mag_joined = new double[ Nx * Ny ];
//      psi_re_joined = new double[ Nx * Ny ];
//      psi_im_joined = new double[ Nx * Ny ];
//      //for ( ptrdiff_t i=0; i<Nx; ++i)
//      //{
//      //   psi_mag_joined[i] = new double[ Ny ];
//      //   psi_re_joined[i] = new double[ Ny ];
//      //   psi_im_joined[i] = new double[ Ny ];
//      //}
//   }
//
//   //cout << "output_psi_realspace_to_netcdf() line 186" << endl; // debug
//   MPI_Gather( psi_mag_local, local_alloc_size_fftw, MPI_DOUBLE,
//         psi_mag_joined, local_alloc_size_fftw, MPI_DOUBLE,
//         rootnode, comm);
//   MPI_Gather( psi_re_local, local_alloc_size_fftw, MPI_DOUBLE,
//         psi_re_joined, local_alloc_size_fftw, MPI_DOUBLE,
//         rootnode, comm);
//   MPI_Gather( psi_im_local, local_alloc_size_fftw, MPI_DOUBLE,
//         psi_im_joined, local_alloc_size_fftw, MPI_DOUBLE,
//         rootnode, comm);
//   //cout << "output_psi_realspace_to_netcdf() line 196" << endl; // debug
//
//
//   if ( mynode == rootnode )
//   {
//      string outFileName = outFilePrefix + ".nc";
//      NcFile dataFile( outFileName.c_str(), NcFile::Replace );
//      if ( ! dataFile.is_valid() ) 
//      {
//         cerr << "failed to open file : " << outFileName << endl;
//         return EXIT_FAILURE;
//      }
//
//      // Domain dimensions each require an NcDim and NcVar, in order to
//      //  describe the domain values, units, and perhaps other attributes.
//      NcDim* kxDimension = dataFile.add_dim("k_x", Nx);
//      NcVar* kxVariable = dataFile.add_var("k_x", ncDouble, kxDimension);
//      kxVariable->add_att("units", "Angstrom");
//
//      NcDim* kyDimension = dataFile.add_dim("k_y", Ny);
//      NcVar* kyVariable = dataFile.add_var("k_y", ncDouble, kyDimension);
//      kyVariable->add_att("units", "Angstrom");
//
//      // Values which are functions of a domain each require an NcVar
//      //  specifying the NcDim elements of their domain
//      NcVar* psi_magVariable = dataFile.add_var("psi_mag", ncDouble, 
//            kxDimension, kyDimension);
//      psi_magVariable->add_att("units", "(arb.)");
//
//      NcVar* psi_reVariable = dataFile.add_var("psi_re", ncDouble, 
//            kxDimension, kyDimension);
//      psi_reVariable->add_att("units", "(arb.)");
//
//      NcVar* psi_imVariable = dataFile.add_var("psi_im", ncDouble, 
//            kxDimension, kyDimension);
//      psi_imVariable->add_att("units", "(arb.)");
//
//      // write the variables to file
//      kxVariable->put( &kx_joined[0], Nx);
//      kyVariable->put( &ky[0], Nx);
//
//      double** psi_mag_joined_2Darray;
//      double** psi_re_joined_2Darray;
//      double** psi_im_joined_2Darray;
//      psi_mag_joined_2Darray = new double*[Nx];
//      psi_re_joined_2Darray = new double*[Nx];
//      psi_im_joined_2Darray = new double*[Nx];
//      for ( ptrdiff_t i=0; i<Nx; ++i )
//      {
//         psi_mag_joined_2Darray[i] = new double[Ny];
//         psi_re_joined_2Darray[i] = new double[Ny];
//         psi_im_joined_2Darray[i] = new double[Ny];
//
//         for ( ptrdiff_t j=0; j<Ny; ++j)
//         {
//            psi_mag_joined_2Darray[i][j] = psi_mag_joined[j + i*Ny];
//            psi_re_joined_2Darray[i][j] = psi_re_joined[j + i*Ny];
//            psi_im_joined_2Darray[i][j] = psi_im_joined[j + i*Ny];
//         }
//      }
//
//      psi_magVariable->put( &psi_mag_joined_2Darray[0][0], Nx, Ny );
//      psi_reVariable->put( &psi_re_joined_2Darray[0][0], Nx, Ny );
//      psi_imVariable->put( &psi_im_joined_2Darray[0][0], Nx, Ny );
//
//      cout << "wave function written to : " << outFileName << endl; // debug
//
//      for ( ptrdiff_t i=0; i<Nx; ++i)
//      {
//         delete[] psi_mag_joined_2Darray[i];
//         delete[] psi_re_joined_2Darray[i];
//         delete[] psi_im_joined_2Darray[i];
//      }
//      delete[] psi_mag_joined_2Darray;
//      delete[] psi_re_joined_2Darray;
//      delete[] psi_im_joined_2Darray;
//
//      delete[] psi_mag_joined;
//      delete[] psi_re_joined;
//      delete[] psi_im_joined;
//   }
//
//   delete[] psi_mag_local;
//   delete[] psi_re_local;
//   delete[] psi_im_local;
//
//   return EXIT_SUCCESS;
//}


//int TEM_NS::output_psi_reciprocalspace_to_txt(
//      const fftw_complex* const psi,   // local array of range elements
//      const ptrdiff_t& local_alloc_size_fftw, // psi size after MPI_Gather
//      const ptrdiff_t Nx_local,  // local size of psi in the x-dimension
//      const double* const kx_joined,   // kx domain elements
//      const ptrdiff_t Nx,     // number of elements in kx_joined
//      const double* const ky, // ky domain elements
//      const ptrdiff_t Ny,     // number of elements in ky
//      const string& outFilePrefix,
//      const int& mynode,
//      const int& rootnode,
//      MPI_Comm comm
//      )
//{
//   // This function differs from its realspace counterpart only in that it
//   //  rearranges psi to put the corners in the center.
//
//   //double* psi_mag_joined; double* psi_re_joined; double* psi_im_joined;
//   //double* psi_mag_local; double* psi_re_local; double* psi_im_local;
//   //psi_mag_local = new double[ Nx_local ][ Ny ];
//   //psi_re_local = new double[ Nx_local ][ Ny ];
//   //psi_im_local = new double[ Nx_local ][ Ny ];
//
//   double* psi_mag_local; 
//   psi_mag_local = new double[ local_alloc_size_fftw];
//   double* psi_re_local; 
//   psi_re_local = new double[ local_alloc_size_fftw];
//   double* psi_im_local; 
//   psi_im_local = new double[ local_alloc_size_fftw];
//
//   //double psi_mag_local[ Nx_local ][ Ny ]; 
//   //double psi_re_local[ Nx_local ][ Ny ]; 
//   //double psi_im_local[ Nx_local ][ Ny ];
//
//   //for (ptrdiff_t i=0; i < Nx_local; ++i)
//   //   for (ptrdiff_t j=0; j < Ny; ++j)
//   //   {
//   //      psi_mag_local[i][j] 
//   //         = sqrt(
//   //               psi[j + i * Ny][0] * psi[j + i * Ny][0]
//   //               + psi[j + i * Ny][1] * psi[j + i * Ny][1]
//   //            );
//   //      psi_re_local[i][j] = psi[j + i * Ny][0];
//   //      psi_im_local[i][j] = psi[j + i * Ny][1];
//   //   }
//
//   for ( ptrdiff_t i=0; i<local_alloc_size_fftw; ++i )
//   {
//      psi_mag_local[i]  
//            = sqrt(
//                  psi[i][0] * psi[i][0]
//                  + psi[i][1] * psi[i][1]
//               );
//      psi_re_local[i] = psi[i][0];
//      psi_im_local[i] = psi[i][1];
//   }
//
//   //if ( mynode == rootnode )
//   //{
//   // TODO: this allocation of arrays is wastefull; these arrays are only
//   //       being allocated on all nodes because I'm not allowed to 
//   //       use 'new double[ Nx ][ Ny ]', and I can't limit their scope
//   //       because they're used by MPI_Gather.
//   //
//   //   double psi_mag_joined[ Nx ][ Ny ]; 
//   //   double psi_re_joined[ Nx ][ Ny ]; 
//   //   double psi_im_joined[ Nx ][ Ny ];
//   double* psi_mag_joined;
//   double* psi_re_joined; 
//   double* psi_im_joined;
//
//   if ( mynode == rootnode )
//   {
//      psi_mag_joined = new double[ Nx * Ny ];
//      psi_re_joined = new double[ Nx * Ny ];
//      psi_im_joined = new double[ Nx * Ny ];
//   }
//   //psi_mag_joined = new double[ Nx][ Ny ];
//   //psi_re_joined = new double[ Nx ][ Ny ];
//   //psi_im_joined = new double[ Nx ][ Ny ];
//
//   MPI_Gather( psi_mag_local, local_alloc_size_fftw, MPI_DOUBLE,
//         psi_mag_joined, local_alloc_size_fftw, MPI_DOUBLE,
//         rootnode, comm);
//   MPI_Gather( psi_re_local, local_alloc_size_fftw, MPI_DOUBLE,
//         psi_re_joined, local_alloc_size_fftw, MPI_DOUBLE,
//         rootnode, comm);
//   MPI_Gather( psi_im_local, local_alloc_size_fftw, MPI_DOUBLE,
//         psi_im_joined, local_alloc_size_fftw, MPI_DOUBLE,
//         rootnode, comm);
//
//   if ( mynode == rootnode )
//   {
//      // Rearrange so that the origin is in the center
//      double tmp_mag, tmp_re, tmp_im;
//      for ( ptrdiff_t i=0; i < Nx/2; ++i)
//      {
//         for ( ptrdiff_t j=0; j < Ny/2; ++j)
//         {
//            tmp_mag = psi_mag_joined[ j + i*Ny ];
//            psi_mag_joined[ j + i*Ny ] 
//               = psi_mag_joined[ (j + Ny/2) + (i + Nx/2) *Ny ];
//            psi_mag_joined[ (j + Ny/2) + (i + Nx/2) *Ny ] = tmp_mag;
//
//            tmp_re = psi_re_joined[ j + i*Ny ];
//            psi_re_joined[ j + i*Ny ] 
//               = psi_re_joined[ (j + Ny/2) + (i + Nx/2) *Ny ];
//            psi_re_joined[ (j + Ny/2) + (i + Nx/2) *Ny ] = tmp_re;
//
//            tmp_im = psi_im_joined[ j + i*Ny ];
//            psi_im_joined[ j + i*Ny ] 
//               = psi_im_joined[ (j + Ny/2) + (i + Nx/2) *Ny ];
//            psi_im_joined[ (j + Ny/2) + (i + Nx/2) *Ny ] = tmp_im;
//
//         }
//
//         for ( ptrdiff_t j=Ny/2; j < Ny; ++j)
//         {
//            //tmp_mag = psi_mag_joined[ i ][ j ];
//            //psi_mag_joined[ i ][ j ] = psi_mag_joined[ i + Nx/2 ][ j - Ny/2 ];
//            //psi_mag_joined[ i + Nx/2 ][ j - Ny/2 ] = tmp_mag;
//
//            //tmp_re = psi_re_joined[ i ][ j ];
//            //psi_re_joined[ i ][ j ] = psi_re_joined[ i + Nx/2 ][ j - Ny/2 ];
//            //psi_re_joined[ i + Nx/2 ][ j - Ny/2 ] = tmp_re;
//
//            //tmp_im = psi_im_joined[ i ][ j ];
//            //psi_im_joined[ i ][ j ] = psi_im_joined[ i + Nx/2 ][ j - Ny/2 ];
//            //psi_im_joined[ i + Nx/2 ][ j - Ny/2 ] = tmp_im;
//
//
//            tmp_mag = psi_mag_joined[ j + i*Ny ];
//            psi_mag_joined[ j + i*Ny ] 
//               = psi_mag_joined[ (j - Ny/2) + (i + Nx/2) *Ny ];
//            psi_mag_joined[ (j - Ny/2) + (i + Nx/2) *Ny ] = tmp_mag;
//
//            tmp_re = psi_re_joined[ j + i*Ny ];
//            psi_re_joined[ j + i*Ny ] 
//               = psi_re_joined[ (j - Ny/2) + (i + Nx/2) *Ny ];
//            psi_re_joined[ (j - Ny/2) + (i + Nx/2) *Ny ] = tmp_re;
//
//            tmp_im = psi_im_joined[ j + i*Ny ];
//            psi_im_joined[ j + i*Ny ] 
//               = psi_im_joined[ (j - Ny/2) + (i + Nx/2) *Ny ];
//            psi_im_joined[ (j - Ny/2) + (i + Nx/2) *Ny ] = tmp_im;
//
//         }
//      }
//
//      double** psi_mag_joined_2Darray;
//      double** psi_re_joined_2Darray;
//      double** psi_im_joined_2Darray;
//      psi_mag_joined_2Darray = new double*[Nx];
//      psi_re_joined_2Darray = new double*[Nx];
//      psi_im_joined_2Darray = new double*[Nx];
//      for ( ptrdiff_t i=0; i<Nx; ++i )
//      {
//         psi_mag_joined_2Darray[i] = new double[Ny];
//         psi_re_joined_2Darray[i] = new double[Ny];
//         psi_im_joined_2Darray[i] = new double[Ny];
//
//         for ( ptrdiff_t j=0; j<Ny; ++j)
//         {
//            psi_mag_joined_2Darray[i][j] = psi_mag_joined[j + i*Ny];
//            psi_re_joined_2Darray[i][j] = psi_re_joined[j + i*Ny];
//            psi_im_joined_2Darray[i][j] = psi_im_joined[j + i*Ny];
//         }
//      }
//      
//      // Create and write to the netCDF file
//      string outFileName = outFilePrefix + ".nc";
//      NcFile dataFile( outFileName.c_str(), NcFile::Replace );
//      if ( ! dataFile.is_valid() ) 
//      {
//         cerr << "failed to open file : " << outFileName << endl;
//         return EXIT_FAILURE;
//      }
//
//      // Domain dimensions each require an NcDim and NcVar, in order to
//      //  describe the domain values, units, and perhaps other attributes.
//      NcDim* kxDimension = dataFile.add_dim("k_x", Nx);
//      NcVar* kxVariable = dataFile.add_var("k_x", ncDouble, kxDimension);
//      kxVariable->add_att("units", "Angstrom");
//
//      NcDim* kyDimension = dataFile.add_dim("k_y", Ny);
//      NcVar* kyVariable = dataFile.add_var("k_y", ncDouble, kyDimension);
//      kyVariable->add_att("units", "Angstrom");
//
//      // Values which are functions of a domain each require an NcVar
//      //  specifying the NcDim elements of their domain
//      NcVar* psi_magVariable = dataFile.add_var("psi_mag", ncDouble, 
//            kxDimension, kyDimension);
//      psi_magVariable->add_att("units", "(arb.)");
//
//      NcVar* psi_reVariable = dataFile.add_var("psi_re", ncDouble, 
//            kxDimension, kyDimension);
//      psi_reVariable->add_att("units", "(arb.)");
//
//      NcVar* psi_imVariable = dataFile.add_var("psi_im", ncDouble, 
//            kxDimension, kyDimension);
//      psi_imVariable->add_att("units", "(arb.)");
//
//      ptrdiff_t counts[2]; counts[0] = Nx; counts[1] = Ny;
//      // write the variables to file
//      kxVariable->put( &kx_joined[0], Nx);
//      kyVariable->put( &ky[0], Nx);
//      psi_magVariable->put( &psi_mag_joined_2Darray[0][0], Nx, Ny );
//      psi_reVariable->put( &psi_re_joined_2Darray[0][0], Nx, Ny );
//      psi_imVariable->put( &psi_im_joined_2Darray[0][0], Nx, Ny );
//
//      cout << "wave function written to : " << outFileName << endl; // debug
//
//      for ( ptrdiff_t i=0; i<Nx; ++i)
//      {
//         delete[] psi_mag_joined_2Darray[i];
//         delete[] psi_re_joined_2Darray[i];
//         delete[] psi_im_joined_2Darray[i];
//      }
//      delete[] psi_mag_joined_2Darray;
//      delete[] psi_re_joined_2Darray;
//      delete[] psi_im_joined_2Darray;
//
//      delete[] psi_mag_joined;
//      delete[] psi_re_joined;
//      delete[] psi_im_joined;
//   }
//
//   delete[] psi_mag_local;
//   delete[] psi_re_local;
//   delete[] psi_im_local;
//
//   return EXIT_SUCCESS;
//}

//int TEM_NS::output_psi_mag_reciprocalspace_to_txt(
//      const double* const psi_mag_local,   // local array of range elements
//      const ptrdiff_t& local_alloc_size_fftw, // psi size after MPI_Gather
//      const ptrdiff_t Nx_local,  // local size of psi in the x-dimension
//      const double* const kx_joined,   // kx domain elements
//      const ptrdiff_t Nx,     // number of elements in kx_joined
//      const double* const ky, // ky domain elements
//      const ptrdiff_t Ny,     // number of elements in ky
//      const string& outFilePrefix,
//      const int& mynode,
//      const int& rootnode,
//      MPI_Comm comm
//      )
//{
//   // This function differs from its realspace counterpart only in that it
//   //  rearranges psi to put the corners in the center.
//
//   //double* psi_mag_joined; double* psi_re_joined; double* psi_im_joined;
//   //double* psi_mag_local; double* psi_re_local; double* psi_im_local;
//   //psi_mag_local = new double[ Nx_local ][ Ny ];
//   //psi_re_local = new double[ Nx_local ][ Ny ];
//   //psi_im_local = new double[ Nx_local ][ Ny ];
//
//   //double** psi_mag_local;
//   //psi_mag_local = new double*[ Nx_local ];
//
//   double* psi_mag_local_nonconst;
//   psi_mag_local_nonconst = new double[ local_alloc_size_fftw ];
//   for (size_t i=0; i<local_alloc_size_fftw; ++i)
//   {
//      psi_mag_local_nonconst[i] = psi_mag_local[i];
//   }
//
//   double* psi_mag_joined;
//   if ( mynode == rootnode )
//   {
//      psi_mag_joined = new double[ Nx * Ny ];
//   }
// 
//   MPI_Gather( psi_mag_local_nonconst, local_alloc_size_fftw, MPI_DOUBLE,
//         psi_mag_joined, local_alloc_size_fftw, MPI_DOUBLE,
//         rootnode, comm);
//
//   delete[] psi_mag_local_nonconst;
//
//   double** psi_mag_joined_2Darray;
//   if ( mynode == rootnode )
//   {
//      psi_mag_joined_2Darray = new double*[Nx];
//      for ( size_t i=0; i<Nx; ++i)
//      {
//         psi_mag_joined_2Darray[i] = new double[Ny];
//         for ( size_t j=0; j<Ny; ++j)
//            psi_mag_joined_2Darray[i][j] = psi_mag_joined[j + i*Ny];
//      }
//
//      delete[] psi_mag_joined;
//      
//      // Rearrange so that the origin is in the center
//      double tmp_mag, tmp_re, tmp_im;
//      for ( ptrdiff_t i=0; i < Nx/2; ++i)
//      {
//         for ( ptrdiff_t j=0; j < Ny/2; ++j)
//         {
//            tmp_mag = psi_mag_joined_2Darray[ i ][ j ];
//            psi_mag_joined_2Darray[ i ][ j ]
//               = psi_mag_joined_2Darray[ i + Nx/2 ][ j + Ny/2 ];
//            psi_mag_joined_2Darray[ i + Nx/2 ][ j + Ny/2 ] = tmp_mag;
//         }
//
//         for ( ptrdiff_t j=Ny/2; j < Ny; ++j)
//         {
//            tmp_mag = psi_mag_joined_2Darray[ i ][ j ];
//            psi_mag_joined_2Darray[ i ][ j ]
//               = psi_mag_joined_2Darray[ i + Nx/2 ][ j - Ny/2 ];
//            psi_mag_joined_2Darray[ i + Nx/2 ][ j - Ny/2 ] = tmp_mag;
//         }
//      }
//      
//      // Create and write to the netCDF file
//      string outFileName = outFilePrefix + ".nc";
//      NcFile dataFile( outFileName.c_str(), NcFile::Replace );
//      if ( ! dataFile.is_valid() ) 
//      {
//         cerr << "failed to open file : " << outFileName << endl;
//         return EXIT_FAILURE;
//      }
//
//      // Domain dimensions each require an NcDim and NcVar, in order to
//      //  describe the domain values, units, and perhaps other attributes.
//      NcDim* kxDimension = dataFile.add_dim("k_x", Nx);
//      NcVar* kxVariable = dataFile.add_var("k_x", ncDouble, kxDimension);
//      kxVariable->add_att("units", "Angstrom");
//
//      NcDim* kyDimension = dataFile.add_dim("k_y", Ny);
//      NcVar* kyVariable = dataFile.add_var("k_y", ncDouble, kyDimension);
//      kyVariable->add_att("units", "Angstrom");
//
//      // Values which are functions of a domain each require an NcVar
//      //  specifying the NcDim elements of their domain
//      NcVar* psi_magVariable = dataFile.add_var("psi_mag", ncDouble, 
//            kxDimension, kyDimension);
//      psi_magVariable->add_att("units", "(arb.)");
//
//      //NcVar* psi_reVariable = dataFile.add_var("psi_re", ncDouble, 
//      //      kxDimension, kyDimension);
//      //psi_reVariable->add_att("units", "(arb.)");
//
//      //NcVar* psi_imVariable = dataFile.add_var("psi_im", ncDouble, 
//      //      kxDimension, kyDimension);
//      //psi_imVariable->add_att("units", "(arb.)");
//
//      // write the variables to file
//      kxVariable->put( &kx_joined[0], Nx);
//      kyVariable->put( &ky[0], Nx);
//      psi_magVariable->put( &psi_mag_joined_2Darray[0][0], Nx, Ny );
//      //psi_reVariable->put( &psi_re_joined[0][0], Nx, Ny );
//      //psi_imVariable->put( &psi_im_joined[0][0], Nx, Ny );
//
//      cout << "wave function written to : " << outFileName << endl; // debug
//   }
//
//   if ( mynode == rootnode )
//   {
//
//      for ( size_t i=0; i<Nx; ++i)
//      {
//         delete[] psi_mag_joined_2Darray[i];
//      }
//
//      delete[] psi_mag_joined_2Darray;
//   }
//
//
//   return EXIT_SUCCESS;
//}

#endif
