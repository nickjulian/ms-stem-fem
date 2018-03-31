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
// File: io_netcdf.cpp
// Purpose: 

#ifndef IO_NETCDF_CPP
#define IO_NETCDF_CPP

#include "io_netcdf.hpp"

int TEM_NS::output_variance_to_netcdf(
      const double* const outDataRaw,  // array of range elements
      const std::vector<double>& binning_boundaries, // domain bdy elements
      const string& outFilePrefix
      )
{
   // outDataRaw : f_of_k_magnitudes_reduced[]
   // Nx : number of elements in f_of_k_magnitudes_reduced[]
   ///////////////////////////////////////////////////////////////////
   // validate input
   ///////////////////////////////////////////////////////////////////

   size_t Nx = binning_boundaries.size() - 1;

   ///////////////////////////////////////////////////////////////////
   // The netCDF libraries require the domain to be an array,
   //  so copy all but the final element of binning_boundaries to an array
   ///////////////////////////////////////////////////////////////////

   double myDomain[ Nx ];
   std::vector<double>::const_iterator bdy_itr = binning_boundaries.begin();
   for ( size_t i=0; i < Nx ; ++i)
   {
      myDomain[i] = *bdy_itr;
      ++bdy_itr;
   }

   // Create the netCDF file, if the file exists it will be overwritten
   //  due to the use of the Replace parameter
   string outFileName = outFilePrefix + ".nc";
   NcFile dataFile( outFileName.c_str(), NcFile::Replace );
   if ( ! dataFile.is_valid() ) 
   {
      //cerr << "failed to open file : " << outFilePrefix + ".nc" << endl;
      cerr << "failed to open file : " << outFileName << endl;
      return EXIT_FAILURE;
   }

   // domain dimensions each require an NcDim and NcVar, in order to 
   //  describe the domain values, units, and perhaps other attributes
   NcDim* kxDimension = dataFile.add_dim("k", Nx);
   NcVar* kxVariable = dataFile.add_var("k", ncDouble, kxDimension);
   kxVariable->add_att("units", "\305^-1");

   // Values which are functions of a domain each require an NcVar 
   //  specifying the NcDim elements of their domain
   NcVar* varianceVariable = dataFile.add_var("V", ncDouble, kxDimension );
   varianceVariable->add_att("units", "arb.");

   // write the variables to file
   kxVariable->put( &myDomain[0], Nx );
   varianceVariable->put( &outDataRaw[0], Nx );
   //TODO: note that this direct output won't work for row/column -major
   //       arrays

   //cout << "Variance data written to : " << outFilePrefix + ".nc" << endl;
   //cout << "Nx : " << Nx << endl; // debug

   return EXIT_SUCCESS;
}


int TEM_NS::output_variance_to_netcdf(
      const double* const outDataRaw,  // array of range elements
      const radial_discretization& sample_rotations, // domain bdy elements
      const string& outFilePrefix
      )
{
   // outDataRaw : f_of_k_magnitudes_reduced[]
   // Nx : number of elements in f_of_k_magnitudes_reduced[]
   ///////////////////////////////////////////////////////////////////
   // validate input
   ///////////////////////////////////////////////////////////////////

   size_t Nr = sample_rotations.size() ;

   ///////////////////////////////////////////////////////////////////
   // The netCDF libraries require the domain to be an array,
   //  so copy the domain radius elements to an array
   ///////////////////////////////////////////////////////////////////

   double myDomain[ Nr ];
   radial_discretization::const_iterator rad_itr = sample_rotations.begin();
   for ( size_t i=0; i < Nr ; ++i)
   {
      myDomain[i] = rad_itr->radius;
      ++rad_itr;
   }

   // Create the netCDF file, if the file exists it will be overwritten
   //  due to the use of the Replace parameter
   string outFileName = outFilePrefix + ".nc";
   NcFile dataFile( outFileName.c_str(), NcFile::Replace );
   if ( ! dataFile.is_valid() ) 
   {
      //cerr << "failed to open file : " << outFilePrefix + ".nc" << endl;
      cerr << "failed to open file : " << outFileName << endl;
      return EXIT_FAILURE;
   }

   // domain dimensions each require an NcDim and NcVar, in order to 
   //  describe the domain values, units, and perhaps other attributes
   NcDim* kxDimension = dataFile.add_dim("k", Nr);
   NcVar* kxVariable = dataFile.add_var("k", ncDouble, kxDimension);
   kxVariable->add_att("units", "\305^-1");

   // Values which are functions of a domain each require an NcVar 
   //  specifying the NcDim elements of their domain
   NcVar* varianceVariable = dataFile.add_var("V", ncDouble, kxDimension );
   varianceVariable->add_att("units", "(arb.)");

   // write the variables to file
   kxVariable->put( &myDomain[0], Nr );
   varianceVariable->put( &outDataRaw[0], Nr );
   //TODO: note that this direct output won't work for row/column -major
   //       arrays

   //cout << "Variance data written to : " << outFilePrefix + ".nc" << endl;
   //cout << "Nr : " << Nr << endl; // debug

   return EXIT_SUCCESS;
}



int TEM_NS::output_psi_realspace_to_netcdf(
      const fftw_complex* const psi,   // local array of range elements
      const ptrdiff_t& local_alloc_size_fftw, // psi size after MPI_Gather
      const ptrdiff_t Nx_local,  // local size of psi in the x-dimension
      const double* const kx_joined,   // kx domain elements
      const ptrdiff_t Nx,     // number of elements in kx_joined
      const double* const ky, // ky domain elements
      const ptrdiff_t Ny,     // number of elements in ky
      const string& outFilePrefix,
      const int* const psi_mag_strides,
      const int* const psi_mag_displacements,
      const int& mynode,
      const int& rootnode,
      MPI_Comm comm
      )
{

   //double* psi_mag_joined; double* psi_re_joined; double* psi_im_joined;
   //double* psi_mag_local; double* psi_re_local; double* psi_im_local;
   //psi_mag_local = new double[ Nx_local ][ Ny ];
   //psi_re_local = new double[ Nx_local ][ Ny ];
   //psi_im_local = new double[ Nx_local ][ Ny ];

   //double psi_mag_local[ Nx_local ][ Ny ]; 
   //double psi_re_local[ Nx_local ][ Ny ]; 
   //double psi_im_local[ Nx_local ][ Ny ];

   double* psi_mag_local; 
   psi_mag_local = new double[ local_alloc_size_fftw];
   double* psi_re_local; 
   psi_re_local = new double[ local_alloc_size_fftw];
   double* psi_im_local; 
   psi_im_local = new double[ local_alloc_size_fftw];


   //cout << "output_psi_realspace_to_netcdf() line 159" << endl; // debug
   //for (ptrdiff_t i=0; i < Nx_local; ++i)
   //   for (ptrdiff_t j=0; j < Ny; ++j)
   //   {
   //      psi_mag_local[i][j] 
   //         = sqrt(
   //               psi[j + i * Ny][0] * psi[j + i * Ny][0]
   //               + psi[j + i * Ny][1] * psi[j + i * Ny][1]
   //            );
   //      psi_re_local[i][j] = psi[j + i * Ny][0];
   //      psi_im_local[i][j] = psi[j + i * Ny][1];
   //   }
   for (ptrdiff_t i=0; i < local_alloc_size_fftw; ++i)
   {
      psi_mag_local[i]
         = sqrt(
               psi[i][0] * psi[i][0] + psi[i][1] * psi[i][1]
            );
      psi_re_local[i] = psi[i][0];
      psi_im_local[i] = psi[i][1];
   }


   //cout << "output_psi_realspace_to_netcdf() line 182" << endl; // debug
   //double psi_mag_joined[ Nx ][ Ny ]; 
   //double psi_re_joined[ Nx ][ Ny ]; 
   //double psi_im_joined[ Nx ][ Ny ];
   double* psi_mag_joined;
   double* psi_re_joined; 
   double* psi_im_joined;

   if ( mynode == rootnode )
   {
      psi_mag_joined = new double[ Nx * Ny ];
      psi_re_joined = new double[ Nx * Ny ];
      psi_im_joined = new double[ Nx * Ny ];
      //for ( ptrdiff_t i=0; i<Nx; ++i)
      //{
      //   psi_mag_joined[i] = new double[ Ny ];
      //   psi_re_joined[i] = new double[ Ny ];
      //   psi_im_joined[i] = new double[ Ny ];
      //}
   }

   //cout << "output_psi_realspace_to_netcdf() line 186" << endl; // debug
   MPI_Gatherv( 
         psi_mag_local, 
         Nx_local* Ny, //local_alloc_size_fftw, 
         MPI_DOUBLE,
         psi_mag_joined, 
         psi_mag_strides,       // recvcount[]
         psi_mag_displacements, // displs[]
         MPI_DOUBLE,
         rootnode, comm);
   MPI_Gatherv( 
         psi_re_local, 
         Nx_local* Ny, //local_alloc_size_fftw, 
         MPI_DOUBLE,
         psi_re_joined, 
         psi_mag_strides,       // recvcount[]
         psi_mag_displacements, // displs[]
         MPI_DOUBLE,
         rootnode, comm);
   MPI_Gatherv( 
         psi_im_local, 
         Nx_local* Ny, //local_alloc_size_fftw, 
         MPI_DOUBLE,
         psi_im_joined, 
         psi_mag_strides,       // recvcount[]
         psi_mag_displacements, // displs[]
         MPI_DOUBLE,
         rootnode, comm);
   //cout << "output_psi_realspace_to_netcdf() line 196" << endl; // debug


   if ( mynode == rootnode )
   {
      string outFileName = outFilePrefix + ".nc";
      NcFile dataFile( outFileName.c_str(), NcFile::Replace );
      if ( ! dataFile.is_valid() ) 
      {
         cerr << "failed to open file : " << outFileName << endl;
         return EXIT_FAILURE;
      }

      // Domain dimensions each require an NcDim and NcVar, in order to
      //  describe the domain values, units, and perhaps other attributes.
      NcDim* kxDimension = dataFile.add_dim("k_x", Nx);
      NcVar* kxVariable = dataFile.add_var("k_x", ncDouble, kxDimension);
      kxVariable->add_att("units", "Angstrom");

      NcDim* kyDimension = dataFile.add_dim("k_y", Ny);
      NcVar* kyVariable = dataFile.add_var("k_y", ncDouble, kyDimension);
      kyVariable->add_att("units", "Angstrom");

      // Values which are functions of a domain each require an NcVar
      //  specifying the NcDim elements of their domain
      NcVar* psi_magVariable = dataFile.add_var("psi_mag", ncDouble, 
            kxDimension, kyDimension);
      psi_magVariable->add_att("units", "(arb.)");

      NcVar* psi_reVariable = dataFile.add_var("psi_re", ncDouble, 
            kxDimension, kyDimension);
      psi_reVariable->add_att("units", "(arb.)");

      NcVar* psi_imVariable = dataFile.add_var("psi_im", ncDouble, 
            kxDimension, kyDimension);
      psi_imVariable->add_att("units", "(arb.)");

      // write the variables to file
      kxVariable->put( &kx_joined[0], Nx);
      kyVariable->put( &ky[0], Nx);

      double** psi_mag_joined_2Darray;
      double** psi_re_joined_2Darray;
      double** psi_im_joined_2Darray;
      psi_mag_joined_2Darray = new double*[Nx];
      psi_re_joined_2Darray = new double*[Nx];
      psi_im_joined_2Darray = new double*[Nx];
      for ( ptrdiff_t i=0; i<Nx; ++i )
      {
         psi_mag_joined_2Darray[i] = new double[Ny];
         psi_re_joined_2Darray[i] = new double[Ny];
         psi_im_joined_2Darray[i] = new double[Ny];

         for ( ptrdiff_t j=0; j<Ny; ++j)
         {
            psi_mag_joined_2Darray[i][j] = psi_mag_joined[j + i*Ny];
            psi_re_joined_2Darray[i][j] = psi_re_joined[j + i*Ny];
            psi_im_joined_2Darray[i][j] = psi_im_joined[j + i*Ny];
         }
      }

      psi_magVariable->put( &psi_mag_joined_2Darray[0][0], Nx, Ny );
      psi_reVariable->put( &psi_re_joined_2Darray[0][0], Nx, Ny );
      psi_imVariable->put( &psi_im_joined_2Darray[0][0], Nx, Ny );

      //cout << "wave function written to : " << outFileName << endl; // debug

      for ( ptrdiff_t i=0; i<Nx; ++i)
      {
         delete[] psi_mag_joined_2Darray[i];
         delete[] psi_re_joined_2Darray[i];
         delete[] psi_im_joined_2Darray[i];
      }
      delete[] psi_mag_joined_2Darray;
      delete[] psi_re_joined_2Darray;
      delete[] psi_im_joined_2Darray;

      delete[] psi_mag_joined;
      delete[] psi_re_joined;
      delete[] psi_im_joined;
   }

   delete[] psi_mag_local;
   delete[] psi_re_local;
   delete[] psi_im_local;

   return EXIT_SUCCESS;
}


int TEM_NS::output_psi_reciprocalspace_to_netcdf(
      const fftw_complex* const psi,   // local array of range elements
      const ptrdiff_t& local_alloc_size_fftw, // psi size after MPI_Gather
      const ptrdiff_t Nx_local,  // local size of psi in the x-dimension
      const double* const kx_joined,   // kx domain elements
      const ptrdiff_t Nx,     // number of elements in kx_joined
      const double* const ky, // ky domain elements
      const ptrdiff_t Ny,     // number of elements in ky
      const string& outFilePrefix,
      const int* const psi_mag_strides,
      const int* const psi_mag_displacements,
      const int& mynode,
      const int& rootnode,
      MPI_Comm comm
      )
{
   // This function differs from its realspace counterpart only in that it
   //  rearranges psi to put the corners in the center.

   //double* psi_mag_joined; double* psi_re_joined; double* psi_im_joined;
   //double* psi_mag_local; double* psi_re_local; double* psi_im_local;
   //psi_mag_local = new double[ Nx_local ][ Ny ];
   //psi_re_local = new double[ Nx_local ][ Ny ];
   //psi_im_local = new double[ Nx_local ][ Ny ];

   double* psi_mag_local; 
   psi_mag_local = new double[ local_alloc_size_fftw];
   double* psi_re_local; 
   psi_re_local = new double[ local_alloc_size_fftw];
   double* psi_im_local; 
   psi_im_local = new double[ local_alloc_size_fftw];

   //double psi_mag_local[ Nx_local ][ Ny ]; 
   //double psi_re_local[ Nx_local ][ Ny ]; 
   //double psi_im_local[ Nx_local ][ Ny ];

   //for (ptrdiff_t i=0; i < Nx_local; ++i)
   //   for (ptrdiff_t j=0; j < Ny; ++j)
   //   {
   //      psi_mag_local[i][j] 
   //         = sqrt(
   //               psi[j + i * Ny][0] * psi[j + i * Ny][0]
   //               + psi[j + i * Ny][1] * psi[j + i * Ny][1]
   //            );
   //      psi_re_local[i][j] = psi[j + i * Ny][0];
   //      psi_im_local[i][j] = psi[j + i * Ny][1];
   //   }

   for ( ptrdiff_t i=0; i<local_alloc_size_fftw; ++i )
   {
      psi_mag_local[i]  
            = sqrt(
                  psi[i][0] * psi[i][0]
                  + psi[i][1] * psi[i][1]
               );
      psi_re_local[i] = psi[i][0];
      psi_im_local[i] = psi[i][1];
   }

   //if ( mynode == rootnode )
   //{
   // TODO: this allocation of arrays is wastefull; these arrays are only
   //       being allocated on all nodes because I'm not allowed to 
   //       use 'new double[ Nx ][ Ny ]', and I can't limit their scope
   //       because they're used by MPI_Gather.
   //
   //   double psi_mag_joined[ Nx ][ Ny ]; 
   //   double psi_re_joined[ Nx ][ Ny ]; 
   //   double psi_im_joined[ Nx ][ Ny ];
   double* psi_mag_joined;
   double* psi_re_joined; 
   double* psi_im_joined;

   if ( mynode == rootnode )
   {
      psi_mag_joined = new double[ Nx * Ny ];
      psi_re_joined = new double[ Nx * Ny ];
      psi_im_joined = new double[ Nx * Ny ];
   }
   //psi_mag_joined = new double[ Nx][ Ny ];
   //psi_re_joined = new double[ Nx ][ Ny ];
   //psi_im_joined = new double[ Nx ][ Ny ];

   MPI_Gatherv( 
         psi_mag_local, 
         Nx_local* Ny, //local_alloc_size_fftw, 
         MPI_DOUBLE,
         psi_mag_joined, 
         psi_mag_strides,       // recvcount[]
         psi_mag_displacements, // displs[]
         MPI_DOUBLE,
         rootnode, comm);
   MPI_Gatherv( 
         psi_re_local, 
         Nx_local* Ny, //local_alloc_size_fftw, 
         MPI_DOUBLE,
         psi_re_joined, 
         psi_mag_strides,       // recvcount[]
         psi_mag_displacements, // displs[]
         MPI_DOUBLE,
         rootnode, comm);
   MPI_Gatherv( 
         psi_im_local, 
         Nx_local* Ny, //local_alloc_size_fftw, 
         MPI_DOUBLE,
         psi_im_joined, 
         psi_mag_strides,       // recvcount[]
         psi_mag_displacements, // displs[]
         MPI_DOUBLE,
         rootnode, comm);

   if ( mynode == rootnode )
   {
      // Rearrange so that the origin is in the center
      double tmp_mag, tmp_re, tmp_im;
      for ( ptrdiff_t i=0; i < Nx/2; ++i)
      {
         for ( ptrdiff_t j=0; j < Ny/2; ++j)
         {
            tmp_mag = psi_mag_joined[ j + i*Ny ];
            psi_mag_joined[ j + i*Ny ] 
               = psi_mag_joined[ (j + Ny/2) + (i + Nx/2) *Ny ];
            psi_mag_joined[ (j + Ny/2) + (i + Nx/2) *Ny ] = tmp_mag;

            tmp_re = psi_re_joined[ j + i*Ny ];
            psi_re_joined[ j + i*Ny ] 
               = psi_re_joined[ (j + Ny/2) + (i + Nx/2) *Ny ];
            psi_re_joined[ (j + Ny/2) + (i + Nx/2) *Ny ] = tmp_re;

            tmp_im = psi_im_joined[ j + i*Ny ];
            psi_im_joined[ j + i*Ny ] 
               = psi_im_joined[ (j + Ny/2) + (i + Nx/2) *Ny ];
            psi_im_joined[ (j + Ny/2) + (i + Nx/2) *Ny ] = tmp_im;

         }

         for ( ptrdiff_t j=Ny/2; j < Ny; ++j)
         {
            //tmp_mag = psi_mag_joined[ i ][ j ];
            //psi_mag_joined[ i ][ j ] = psi_mag_joined[ i + Nx/2 ][ j - Ny/2 ];
            //psi_mag_joined[ i + Nx/2 ][ j - Ny/2 ] = tmp_mag;

            //tmp_re = psi_re_joined[ i ][ j ];
            //psi_re_joined[ i ][ j ] = psi_re_joined[ i + Nx/2 ][ j - Ny/2 ];
            //psi_re_joined[ i + Nx/2 ][ j - Ny/2 ] = tmp_re;

            //tmp_im = psi_im_joined[ i ][ j ];
            //psi_im_joined[ i ][ j ] = psi_im_joined[ i + Nx/2 ][ j - Ny/2 ];
            //psi_im_joined[ i + Nx/2 ][ j - Ny/2 ] = tmp_im;


            tmp_mag = psi_mag_joined[ j + i*Ny ];
            psi_mag_joined[ j + i*Ny ] 
               = psi_mag_joined[ (j - Ny/2) + (i + Nx/2) *Ny ];
            psi_mag_joined[ (j - Ny/2) + (i + Nx/2) *Ny ] = tmp_mag;

            tmp_re = psi_re_joined[ j + i*Ny ];
            psi_re_joined[ j + i*Ny ] 
               = psi_re_joined[ (j - Ny/2) + (i + Nx/2) *Ny ];
            psi_re_joined[ (j - Ny/2) + (i + Nx/2) *Ny ] = tmp_re;

            tmp_im = psi_im_joined[ j + i*Ny ];
            psi_im_joined[ j + i*Ny ] 
               = psi_im_joined[ (j - Ny/2) + (i + Nx/2) *Ny ];
            psi_im_joined[ (j - Ny/2) + (i + Nx/2) *Ny ] = tmp_im;

         }
      }

      double** psi_mag_joined_2Darray;
      double** psi_re_joined_2Darray;
      double** psi_im_joined_2Darray;
      psi_mag_joined_2Darray = new double*[Nx];
      psi_re_joined_2Darray = new double*[Nx];
      psi_im_joined_2Darray = new double*[Nx];
      for ( ptrdiff_t i=0; i<Nx; ++i )
      {
         psi_mag_joined_2Darray[i] = new double[Ny];
         psi_re_joined_2Darray[i] = new double[Ny];
         psi_im_joined_2Darray[i] = new double[Ny];

         for ( ptrdiff_t j=0; j<Ny; ++j)
         {
            psi_mag_joined_2Darray[i][j] = psi_mag_joined[j + i*Ny];
            psi_re_joined_2Darray[i][j] = psi_re_joined[j + i*Ny];
            psi_im_joined_2Darray[i][j] = psi_im_joined[j + i*Ny];
         }
      }
      
      // Create and write to the netCDF file
      string outFileName = outFilePrefix + ".nc";
      NcFile dataFile( outFileName.c_str(), NcFile::Replace );
      if ( ! dataFile.is_valid() ) 
      {
         cerr << "failed to open file : " << outFileName << endl;
         return EXIT_FAILURE;
      }

      // Domain dimensions each require an NcDim and NcVar, in order to
      //  describe the domain values, units, and perhaps other attributes.
      NcDim* kxDimension = dataFile.add_dim("k_x", Nx);
      NcVar* kxVariable = dataFile.add_var("k_x", ncDouble, kxDimension);
      kxVariable->add_att("units", "Angstrom");

      NcDim* kyDimension = dataFile.add_dim("k_y", Ny);
      NcVar* kyVariable = dataFile.add_var("k_y", ncDouble, kyDimension);
      kyVariable->add_att("units", "Angstrom");

      // Values which are functions of a domain each require an NcVar
      //  specifying the NcDim elements of their domain
      NcVar* psi_magVariable = dataFile.add_var("psi_mag", ncDouble, 
            kxDimension, kyDimension);
      psi_magVariable->add_att("units", "(arb.)");

      NcVar* psi_reVariable = dataFile.add_var("psi_re", ncDouble, 
            kxDimension, kyDimension);
      psi_reVariable->add_att("units", "(arb.)");

      NcVar* psi_imVariable = dataFile.add_var("psi_im", ncDouble, 
            kxDimension, kyDimension);
      psi_imVariable->add_att("units", "(arb.)");

      ptrdiff_t counts[2]; counts[0] = Nx; counts[1] = Ny;
      // write the variables to file
      kxVariable->put( &kx_joined[0], Nx);
      kyVariable->put( &ky[0], Nx);
      psi_magVariable->put( &psi_mag_joined_2Darray[0][0], Nx, Ny );
      psi_reVariable->put( &psi_re_joined_2Darray[0][0], Nx, Ny );
      psi_imVariable->put( &psi_im_joined_2Darray[0][0], Nx, Ny );

      //cout << "wave function written to : " << outFileName << endl; // debug

      for ( ptrdiff_t i=0; i<Nx; ++i)
      {
         delete[] psi_mag_joined_2Darray[i];
         delete[] psi_re_joined_2Darray[i];
         delete[] psi_im_joined_2Darray[i];
      }
      delete[] psi_mag_joined_2Darray;
      delete[] psi_re_joined_2Darray;
      delete[] psi_im_joined_2Darray;

      delete[] psi_mag_joined;
      delete[] psi_re_joined;
      delete[] psi_im_joined;
   }

   delete[] psi_mag_local;
   delete[] psi_re_local;
   delete[] psi_im_local;

   return EXIT_SUCCESS;
}

int TEM_NS::output_psi_mag_reciprocalspace_to_netcdf(
      const double* const psi_mag_local,   // local array of range elements
      const ptrdiff_t& local_alloc_size_fftw, // psi size after MPI_Gather
      const ptrdiff_t Nx_local,  // local size of psi in the x-dimension
      const double* const kx_joined,   // kx domain elements
      const ptrdiff_t Nx,     // number of elements in kx_joined
      const double* const ky, // ky domain elements
      const ptrdiff_t Ny,     // number of elements in ky
      const string& outFilePrefix,
      const int* const psi_mag_strides,
      const int* const psi_mag_displacements,
      const int& mynode,
      const int& rootnode,
      MPI_Comm comm
      )
{
   // This function differs from its realspace counterpart only in that it
   //  rearranges psi to put the corners in the center.

   //double* psi_mag_joined; double* psi_re_joined; double* psi_im_joined;
   //double* psi_mag_local; double* psi_re_local; double* psi_im_local;
   //psi_mag_local = new double[ Nx_local ][ Ny ];
   //psi_re_local = new double[ Nx_local ][ Ny ];
   //psi_im_local = new double[ Nx_local ][ Ny ];

   //double** psi_mag_local;
   //psi_mag_local = new double*[ Nx_local ];

   double* psi_mag_local_nonconst;
   psi_mag_local_nonconst = new double[ local_alloc_size_fftw ];
   for (size_t i=0; i<local_alloc_size_fftw; ++i)
   {
      psi_mag_local_nonconst[i] = psi_mag_local[i];
   }

   double* psi_mag_joined;
   if ( mynode == rootnode )
   {
      psi_mag_joined = new double[ Nx * Ny ];
   }
 
   MPI_Gatherv( 
         psi_mag_local_nonconst, 
         Nx_local* Ny, //local_alloc_size_fftw, 
         MPI_DOUBLE,
         psi_mag_joined, 
         psi_mag_strides,
         psi_mag_displacements,
         MPI_DOUBLE,
         rootnode, comm);

   delete[] psi_mag_local_nonconst;

   double** psi_mag_joined_2Darray;
   if ( mynode == rootnode )
   {
      psi_mag_joined_2Darray = new double*[Nx];
      for ( size_t i=0; i<Nx; ++i)
      {
         psi_mag_joined_2Darray[i] = new double[Ny];
         for ( size_t j=0; j<Ny; ++j)
            psi_mag_joined_2Darray[i][j] = psi_mag_joined[j + i*Ny];
      }

      delete[] psi_mag_joined;
      
      // Rearrange so that the origin is in the center
      double tmp_mag, tmp_re, tmp_im;
      for ( ptrdiff_t i=0; i < Nx/2; ++i)
      {
         for ( ptrdiff_t j=0; j < Ny/2; ++j)
         {
            tmp_mag = psi_mag_joined_2Darray[ i ][ j ];
            psi_mag_joined_2Darray[ i ][ j ]
               = psi_mag_joined_2Darray[ i + Nx/2 ][ j + Ny/2 ];
            psi_mag_joined_2Darray[ i + Nx/2 ][ j + Ny/2 ] = tmp_mag;
         }

         for ( ptrdiff_t j=Ny/2; j < Ny; ++j)
         {
            tmp_mag = psi_mag_joined_2Darray[ i ][ j ];
            psi_mag_joined_2Darray[ i ][ j ]
               = psi_mag_joined_2Darray[ i + Nx/2 ][ j - Ny/2 ];
            psi_mag_joined_2Darray[ i + Nx/2 ][ j - Ny/2 ] = tmp_mag;
         }
      }
      
      // Create and write to the netCDF file
      string outFileName = outFilePrefix + ".nc";
      NcFile dataFile( outFileName.c_str(), NcFile::Replace );
      if ( ! dataFile.is_valid() ) 
      {
         cerr << "failed to open file : " << outFileName << endl;
         return EXIT_FAILURE;
      }

      // Domain dimensions each require an NcDim and NcVar, in order to
      //  describe the domain values, units, and perhaps other attributes.
      NcDim* kxDimension = dataFile.add_dim("k_x", Nx);
      NcVar* kxVariable = dataFile.add_var("k_x", ncDouble, kxDimension);
      kxVariable->add_att("units", "Angstrom");

      NcDim* kyDimension = dataFile.add_dim("k_y", Ny);
      NcVar* kyVariable = dataFile.add_var("k_y", ncDouble, kyDimension);
      kyVariable->add_att("units", "Angstrom");

      // Values which are functions of a domain each require an NcVar
      //  specifying the NcDim elements of their domain
      NcVar* psi_magVariable = dataFile.add_var("psi_mag", ncDouble, 
            kxDimension, kyDimension);
      psi_magVariable->add_att("units", "(arb.)");

      //NcVar* psi_reVariable = dataFile.add_var("psi_re", ncDouble, 
      //      kxDimension, kyDimension);
      //psi_reVariable->add_att("units", "(arb.)");

      //NcVar* psi_imVariable = dataFile.add_var("psi_im", ncDouble, 
      //      kxDimension, kyDimension);
      //psi_imVariable->add_att("units", "(arb.)");

      // write the variables to file
      kxVariable->put( &kx_joined[0], Nx);
      kyVariable->put( &ky[0], Nx);
      psi_magVariable->put( &psi_mag_joined_2Darray[0][0], Nx, Ny );
      //psi_reVariable->put( &psi_re_joined[0][0], Nx, Ny );
      //psi_imVariable->put( &psi_im_joined[0][0], Nx, Ny );

      //cout << "wave function written to : " << outFileName << endl; // debug
   }

   if ( mynode == rootnode )
   {

      for ( size_t i=0; i<Nx; ++i)
      {
         delete[] psi_mag_joined_2Darray[i];
      }

      delete[] psi_mag_joined_2Darray;
   }


   return EXIT_SUCCESS;
}

// TODO: fix the following function
//int TEM_NS::append_correlograph_to_netcdf(
//   const double* const outDataRaw,
//   const std::vector<double>& k_binning_boundaries,
//   const std::vector<double>& phi_binning_boundaries,
//   const string& outFilePrefix
//   )
//{
//   size_t Nx = k_binning_boundaries.size() - 1;
//   size_t Ny = phi_binning_boundaries.size() - 1;
//   double myXDomain[ Nx ];
//   double myYDomain[ Ny ];
//   std::vector<double>::const_iterator 
//      k_bdy_itr = k_binning_boundaries.begin();
//   std::vector<double>::const_iterator 
//      phi_bdy_itr = phi_binning_boundaries.begin();
//   for ( size_t i=0; i < Nx ; ++i)
//   {
//      myXDomain[i] = *k_bdy_itr;
//      ++k_bdy_itr;
//   }
//   for ( size_t i=0; i < Ny ; ++i)
//   {
//      myYDomain[i] = *phi_bdy_itr;
//      ++phi_bdy_itr;
//   }
//
//   // Since row/column major data doesn't work with netcdf put methods,
//   //  copy the output to a 2-D array
//   double twoDimData[ Nx ][ Ny ];
//   for (size_t i=0; i<Nx; ++i)
//      for (size_t j=0; j<Ny; ++j)
//         twoDimeData[i][j] = outDataRaw[j + i * Ny];
//
//   string outFileName = outFilePrefix + ".nc";
//   NcFile dataFile( outFileName.c_str(), NcFile::Write );
//   if ( ! dataFile.is_valid() ) 
//   {
//      cerr << "failed to open file : " << outFileName << endl;
//      return EXIT_FAILURE;
//   }
//
//   NcDim* kxDimension = dataFile.add_dim("k", Nx);
//   NcVar* kxVariable = dataFile.add_var("k", ncDouble, kxDimension);
//   kxVariable->add_att("units", "\305^-1");
//
//   NcDim* phiDimension = dataFile.add_dim("phi", Ny);
//   NcVar* phiVariable = dataFile.add_var("phi", ncDouble, phiDimension);
//   phiVariable->add_att("units", "radian");
//
//
//   kxVariable->put( &myXDomain[0], Nx );
//   phiVariable->put( &myYDomain[0], Ny );
//
//   Variable->put( &outDataRaw[0], Nx );
//
//   return EXIT_SUCCESS;
//}

#endif
