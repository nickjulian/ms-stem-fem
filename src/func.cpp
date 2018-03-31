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
// File: func.cpp
// Purpose:
//    - implement functions used to generate data for the test.

#ifndef FUNC_CPP
#define FUNC_CPP


#include <cstdlib>
#include <iostream>

#include "func.hpp"

using namespace std;
using std::sqrt;
using std::cerr;
using std::endl;
using std::setw;
using std::setprecision;

using std::vector;
using std::list;

using TEM_NS::to_string;

void TEM_NS::domain_1D(
         const int N,
         const double P,
         double *x
         )
{
   //const double pi = 3.14159265359;
   const double NOverP = ((double) N) / P; 
   for(int i=0; i<N; i++)
      x[i] = ( 1 / NOverP ) * (2 * PI) * i;

   return ;
}

void TEM_NS::domain_2D(
         const int& Nx, const int& Ny,
         const double& qmin_x, const double& qmax_x, 
         const double& qmin_y, const double& qmax_y,
         double *qx, double *qy
         )
{
   // periodicity: qmax_x + dx == qmin_x
   const double dx = (qmax_x - qmin_x)/(Nx-1);
   const double dy = (qmax_y - qmin_y)/(Ny-1);
   //cout << "dx, dy : " << dx << ", " << dy << endl; // debug

   for(int i=0; i<Nx; i++) qx[i] = qmin_x + ( i * dx) ;
   for(int i=0; i<Ny; i++) qy[i] = qmin_y + ( i * dy) ;

   return ;
}

void TEM_NS::domain_2D_periodic(
         const int& Nx, const int& Ny,
         const double& qperiod_x, const double& qperiod_y, 
         const double& qmin_x, const double& qmin_y,
         const double& dx, const double& dy,
         double *qx, double *qy
         )
{
   // periodicity: qperiod_x == (dx + qmax_x - qmin_x);
   //dx = qperiod_x/Nx;
   //dy = qperiod_y/Ny;
   //cout << "dx, dy : " << dx << ", " << dy << endl; // debug

   for(int i=0; i<Nx; i++) qx[i] = qmin_x + ( i * dx) ;
   for(int i=0; i<Ny; i++) qy[i] = qmin_y + ( i * dy) ;

   //cout << "xdomain: [" << qx[0] << ": " << dx << ": "<< qx[Nx-1] << "]" << endl;//debug
   //cout << "ydomain: [" << qy[0] << ": " << dy << ": "<< qy[Ny-1] << "]" << endl;//debug

   return ;
}

void TEM_NS::assign_stem_raster_points(
         const int& Nx, const int& Ny,
         const double& qperiod_x, const double& qmin_x, 
         const double& qperiod_y, const double& qmin_y,
         //double* qx, double* qy
         list<double>& qx, list<double>& qy
         )
{
   // periodicity: qmax_x + dx == qmin_x
   const double dx = qperiod_x/Nx;
   const double dy = qperiod_y/Ny;
   //cout << "dx, dy : " << dx << ", " << dy << endl; // debug

   for(int i=0; i<Nx; i++) qx.push_back( qmin_x + ( i * dx) );
   for(int j=0; j<Ny; j++) qy.push_back( qmin_y + ( j * dy) );
   //for(int i=0; i<Nx; i++) qx[i] = qmin_x + ( i * dx) ;
   //for(int j=0; j<Ny; j++) qy[j] = qmin_y + ( j * dy) ;

   return ;
}


void TEM_NS::domain_2D_periodic_recip(
         const int& Nx, const int& Ny,
         const double& xperiod, const double& yperiod,
         const double& kxperiod, const double& kyperiod,
         const double& dkx, const double& dky,
         double *kx, double *ky
         )
{
   // periodicity: Nx/xperiod == 0, Nx/(2 xperiod) = -Nx/(2 xperiod)
   //const double kxperiod = Nx / xperiod; 
   //const double kyperiod = Ny / yperiod;
   //dkx = 1/xperiod; // == kxperiod / Nx == (Nx/xperiod)/Nx
   //dky = 1/yperiod; 
   //cout << "dx, dy : " << dx << ", " << dy << endl; // debug

   for(int i=0; i<Nx/2; i++) kx[i] = i * dkx ;
   for(int i=Nx/2; i<Nx; i++) kx[i] = i * dkx - kxperiod ;
   for(int i=0; i<Ny/2; i++) ky[i] = i * dky ;
   for(int i=Ny/2; i<Ny; i++) ky[i] = i * dky - kyperiod ;

   //int precision = 7;//debug
   //int width = precision + 2;//debug
   //cout << "(kxperiod, xperiod): (" 
   //   << setwidth(width) << setprecision(precision)//debug
   //   << kxperiod << ", " //debug
   //   << setwidth(width) << setprecision(precision)//debug
   //   << xperiod << ")" << endl;//debug
   //cout << "kxdomain: [" //debug
   //   << setwidth(width) << setprecision(precision)//debug
   //   << kx[0] << ": " //debug
   //   << setwidth(width) << setprecision(precision)//debug
   //   << dkx << ": "//debug
   //   << setwidth(width) << setprecision(precision)//debug
   //   << kx[Nx-1] << "]" << endl;//debug
   //cout << "kydomain: [" //debug
   //   << setwidth(width) << setprecision(precision)//debug
   //   << ky[0] << ": " //debug
   //   << setwidth(width) << setprecision(precision)//debug
   //   << dky << ": "//debug
   //   << setwidth(width) << setprecision(precision)//debug
   //   << ky[Ny-1] << "]" << endl;//debug

   return ;
}

void TEM_NS::domain_2D_periodic_recip_pos(
         const int& Nx, const int& Ny,
         //const double& kxperiod,  // ? should it be 2 \pi Nx/xperiod
         //const double& kyperiod,  // or just Nx / xperiod ?
         const double& xperiod,  // x period in real space
         const double& yperiod,  // y period in real space
         double *kx, double *ky
         )
{
   // periodicity: Nx/xperiod == 0, Nx/(2 xperiod) = -Nx/(2 xperiod)
   //const double kxperiod = Nx / xperiod; 
   //const double kyperiod = Ny / yperiod;
   //const double kxperiod = 2 * PI * Nx / xperiod; 
   //const double kyperiod = 2 * PI * Ny / yperiod;
   const double dkx = 1 /xperiod; // == kxperiod / Nx == (Nx/xperiod)/Nx
   const double dky = 1 /yperiod; 
   //const double dkx = 2 * PI /xperiod; // == kxperiod / Nx == (Nx/xperiod)/Nx
   //const double dky = 2 * PI /yperiod; 
   //cout << "dx, dy : " << dx << ", " << dy << endl; // debug

   //for(int i=0; i<Nx; i++) kx[i] = -( PI/xperiod) + i * dkx ;
   //for(int i=0; i<Ny; i++) ky[i] = -( PI/yperiod) + i * dky ;
   for(int i=0; i<Nx; i++) kx[i] = i * dkx ;
   for(int i=0; i<Ny; i++) ky[i] = i * dky ;

   return ;
}

void TEM_NS::domain_3D_periodic(
         const int& Nx, const int& Ny, const int& Nz,// pixels per dimension
         const double& qperiod_x, const double& qmin_x, 
         const double& qperiod_y, const double& qmin_y,
         const double& qperiod_z, const double& qmin_z, 
         double* qx, double* qy, double* qz
         )
{
   // periodicity: qmax_x + dx == qmin_x
   const double dx = qperiod_x/Nx;
   const double dy = qperiod_y/Ny;
   const double dz = qperiod_z/Nz;
   //cout << "dx, dy : " << dx << ", " << dy << endl; // debug

   for(int i=0; i<Nx; i++) qx[i] = qmin_x + ( i * dx);
   for(int i=0; i<Ny; i++) qy[i] = qmin_y + ( i * dy);
   for(int i=0; i<Nz; i++) qz[i] = qmin_z + ( i * dz);

   return ;
}

void TEM_NS::domain_3D_periodic_recip(
         const int& Nx, const int& Ny, const int& Nz,
         const double& xperiod,  // x period in real space
         const double& yperiod,  // y period in real space
         const double& zperiod,  // z period in real space
         double* kx, double* ky, double* kz
         )
{
   // periodicity: Nx/xperiod == 0, Nx/(2 xperiod) = -Nx/(2 xperiod)
   const double kxperiod = Nx / xperiod; 
   const double kyperiod = Ny / yperiod;
   const double kzperiod = Nz / zperiod;
   const double dkx = 1/xperiod; // == kxperiod / Nx == (Nx/xperiod)/Nx
   const double dky = 1/yperiod; 
   const double dkz = 1/zperiod; 

   for(int i=0; i<Nx/2; i++) kx[i] = i * dkx ;
   for(int i=Nx/2; i<Nx; i++) kx[i] = i * dkx - kxperiod ;
   for(int i=0; i<Ny/2; i++) ky[i] = i * dky ;
   for(int i=Ny/2; i<Ny; i++) ky[i] = i * dky - kyperiod ;
   for(int i=0; i<Nz/2; i++) kz[i] = i * dkz ;
   for(int i=Nz/2; i<Nz; i++) kz[i] = i * dkz - kzperiod ;

   return ;
}


void TEM_NS::cos_1D(
      const int N,
      const double A, 
      const double W,
      const double q0,
      double *q,
      double *f
      )
{
   for (int i=0; i< N; i++)
      f[i] += A* cos( W * (q[i] - q0) );
   return ; 
}

void TEM_NS::gaussian_1D(
         const int N,
         const double A,
         const double mu,
         const double sigma,
         double *q,
         double *f
         )
{

   for (int i=0; i< N/2 ; i++)
      f[i] += A * exp( -1* pow( (q[i] - mu), 2) / (pow(2 * sigma, 2)) );

   for (int i=N/2; i<N; i++)
   //TODO:this shouldn't be implemented as a reflection,  but as an addition
      f[i] += A * exp( -1* pow( (q[N - 1 - i] - mu), 2) 
               / (pow(2 * sigma, 2)) );

   return ;
}

void TEM_NS::step_2D(
      const int& Nx, const int& Ny,
      const double& q0x, const double& q0y,
      const double* const qx, const double* const qy,
      double *f
      )
{

   if ( 
         ( q0x < qx[0] ) || ( q0x >  qx[size_t(Nx/2 - 1)] ) 
         || 
         ( q0y < qy[0] ) || ( q0y >  qy[size_t(Ny/2 - 1)] )
      )
   {
      cerr << "q0x : " << q0x << ", q0y : " << q0y 
         << " out of range" << endl;
      cerr << "qx[0] : " << qx[0] << ", qx[ Nx/2 -1 ] : " 
         << qx[size_t(Nx/2 -1)] 
         << "qy[0] : " << qy[0] << ", qy[ Ny/2 -1 ] : " 
         << qy[size_t(Ny/2 -1)] << endl;
      exit( 1 );
      //return EXIT_FAILURE;
   }

   // 2-D row major array access: f(i,j) == f[j+Ny*i]

   size_t idx = 0; 
   double qmax_x = qx[Nx-1];
   double qmax_y = qy[Ny-1];

   for (int i=0; i<Nx; i++)
      for(int j=0; j<Ny; j++)
      {
         idx = j + Ny * i;
         if( 
            ( (qx[i] < q0x) || (qx[i] > (qmax_x - q0x) ) )
            &&
            ( (qy[j] < q0y) || (qy[j] > (qmax_y - q0y) ) )
           )
            f[idx] = 1; 
         else 
            f[idx] = 0; 
      }
   return ;
}

void TEM_NS::gaussian_2D(
      // Precondition: 
      //    - q[] & f[] are row-major arrays of 2-D data
      //    - the domain q[Nx*Ny] has been initialized with values.
      //    - the range f[Nx*Ny] has been initialized with values.
      // Postcondition: 
      //    - f[Nx*Ny] values have been incremented by a gaussian curve
      //       centered at (mux,muy), with stddev sigma and magnitude A.
      const int& Nx,
      const int& Ny,
      const double& A,
      const double& mux, const double& muy,
      const double& sigma,
      const double* const qx,
      const double* const qy,
      double *f
      )
{
   int idx = 0;
   for (int i=0; i< Nx ; i++)
      for (int j=0; j< Ny ; j++)
      {
         idx = j + Ny * i; // (i,j)th element of a 2-D row-major array
         f[idx] 
            += A * exp( -1* 
                  (pow( (qx[i] - mux), 2) + pow( (qy[j] - muy), 2) )
                  / (pow(2 * sigma, 2)) );
      }

// //TODO:this shouldn't be implemented as a reflection,  but as an addition
//   idx = (Ny/2) + (Ny/2)*(Nx/2);
//   for (int i=Nx/2; i<Nx; i++)
//      for (int j=Ny/2; j<Ny; j++)
//      {
//         idx = (Ny - 1 -j) + Ny * (Nx - 1 - i); 
//         f[idx] 
//            += A * exp( -1* 
//                  ( pow( (qx[Nx - 1 - i] - mux), 2) 
//                     + pow( (qy[Ny - 1 - j] - muy), 2) )
//                  / (pow(2 * sigma, 2)) );
//      }
//
   return;
}

void TEM_NS::gaussian_2D_lattice_square(
      // Precondition: 
      //    - q[] & f[] are row-major arrays of 2-D data
      //    - the domain q[Nx*Ny] has been initialized with values.
      //    - the range f[Nx*Ny] has been initialized with values.
      // Postcondition: 
      //    - f[Nx*Ny] values have been incremented by a gaussian curve
      //       centered at (mux,muy), with stddev sigma and magnitude A.
      const int& Nx,
      const int& Ny,
      const double& A,
      const int& Nkx, const int& Nky,
      const double& sigma,
      const double* const qx,
      const double* const qy,
      double *f
      )
{

   for ( int i = 0; i < Nkx; i++ )
      for ( int j = 0; j < Nky; j++ )
      {

         gaussian_2D( Nx, Ny, A,
               qx[ size_t( i * Nx / Nkx ) ], qy[ size_t( j * Ny / Nky ) ], 
               sigma, qx, qy, f);
      }
   for ( int i = 0; i < Nkx; i++ )
      gaussian_2D( Nx, Ny, A,
            qx[ size_t(  i * Nx / Nkx ) ], qy[ size_t( Ny - 1 ) ], 
            sigma, qx, qy, f);

   for ( int j = 0; j < Nky; j++ )
      gaussian_2D( Nx, Ny, A,
            qx[ size_t(  Nx - 1  ) ], qy[ size_t( j * Ny / Nky ) ], 
            sigma, qx, qy, f);

   gaussian_2D( Nx, Ny, A,
         qx[ size_t(  Nx - 1 ) ], qy[ size_t( Ny - 1 ) ], 
         sigma, qx, qy, f);

   return; 
}

void TEM_NS::positions_2D_lattice_square(
   const size_t& Nx,    // number of atomic lattice positions, x direction
   const size_t& Ny,    // number of atomic lattice positions, y direction
   const size_t& Nz,    // number of atomic lattice positions, z direction
   const double* const x,  // x domain
   const double& x_size,   // number of points in the discretized x domain
   const double* const y,  // y domain
   const double& y_size,   // number of points in the discretized y domain
   //const double* const z,  // z domain
   //const double& z_size,   // number of points in the discretized z domain
   double* const q         // lattice point positions (x,y,z), output
   )
{
   for( size_t i = 0; i < Nx; i++)
      for( size_t j = 0; j < Ny; j++)
      //   for( size_t k = 0; k < Nz; k++) // 3-D
      {
         q[ 2*(i * Ny + j) + 0 ] = x[ size_t( i * x_size / Nx )];// 2-D
         q[ 2*(i * Ny + j) + 1 ] = y[ size_t( j * y_size / Ny )];// 2-D
         //q[2*((i * Ny + j)*Nz + k) + 0 ] = x[ size_t( i * x_size / Nx )];
         //q[2*((i * Ny + j)*Nz + k) + 1 ] = y[ size_t( j * y_size / Ny )];
         //q[2*((i * Ny + j)*Nz + k) + 2 ] = z[ size_t( k * z_size / Nz )];
         //population++;
      }

   return;
}

void TEM_NS::positions_3D_lattice_diamond_001(
   const size_t& x_cells,  // number of unit cells, x direction
   const size_t& y_cells,  // number of unit cells, x direction
   const size_t& z_cells,  // number of unit cells, x direction
   const double* const x,  // x domain
   const double& x_pixels, //number of points in discretized x domain 
   const double* const y,  // y domain
   const double& y_pixels, //number of points in discretized y domain 
   const double* const z,  // z domain
   const double& z_pixels, //number of points in discretized z domain 
   double* const q         // lattice point positions (x,y,z), output
   //std::list<double[3]> q       // lattice point positions (x,y,z), output
   )
{
   for( size_t i = 0; i < x_cells; i++)
      for( size_t j = 0; j < y_cells; j++)
         for( size_t k = 0; k < z_cells; k++) // 3-D
      {
         // (i,j,k) unit cell index
         // 8 atoms * 3 coordinates = 24 doubles per unit cell
         
         // atom at unit cell origin
         q[ 24*((i * y_cells + j) * z_cells + k) + 0 ] 
            = x[ size_t( i * x_pixels / x_cells )];
         q[ 24*((i * y_cells + j) * z_cells + k) + 1 ] 
            = y[ size_t( j * y_pixels / y_cells )];
         q[ 24*((i * y_cells + j) * z_cells + k) + 2 ] 
            = z[ size_t( k * z_pixels / z_cells )];

         // atom at (1/2,0,1/2)
         q[ 24*((i * y_cells + j) * z_cells + k) + 3 ] 
            = x[ size_t( (i + 0.5) * (x_pixels / x_cells) )];
         q[ 24*((i * y_cells + j) * z_cells + k) + 4 ] 
            = y[ size_t( j * (y_pixels / y_cells) )];
         q[ 24*((i * y_cells + j) * z_cells + k) + 5 ] 
            = z[ size_t( (k + 0.5) * (z_pixels / z_cells) )];

         // atom at (1/4,1/4,1/4)
         q[ 24*((i * y_cells + j) * z_cells + k) + 6 ] 
            = x[ size_t( (i + 0.25) * (x_pixels / x_cells) )];
         q[ 24*((i * y_cells + j) * z_cells + k) + 7 ] 
            = y[ size_t( (j + 0.25) * (y_pixels / y_cells) )];
         q[ 24*((i * y_cells + j) * z_cells + k) + 8 ] 
            = z[ size_t( (k + 0.25) * (z_pixels / z_cells) )];

         // atom at (0,1/2,1/2)
         q[ 24*((i * y_cells + j) * z_cells + k) + 9 ] 
            = x[ size_t( i * (x_pixels / x_cells) )];
         q[ 24*((i * y_cells + j) * z_cells + k) + 10 ] 
            = y[ size_t( (j + 0.5) * (y_pixels / y_cells) )];
         q[ 24*((i * y_cells + j) * z_cells + k) + 11 ] 
            = z[ size_t( (k + 0.5) * (z_pixels / z_cells) )];

         // atom at (1/2,1/2,0)
         q[ 24*((i * y_cells + j) * z_cells + k) + 12 ] 
            = x[ size_t( (i + 0.5) * (x_pixels / x_cells) )];
         q[ 24*((i * y_cells + j) * z_cells + k) + 13 ] 
            = y[ size_t( (j + 0.5) * (y_pixels / y_cells) )];
         q[ 24*((i * y_cells + j) * z_cells + k) + 14 ] 
            = z[ size_t( k * (z_pixels / z_cells) )];
         
         // atom at (3/4,3/4,1/4)
         q[ 24*((i * y_cells + j) * z_cells + k) + 15 ] 
            = x[ size_t( (i + 0.75) * (x_pixels / x_cells) )];
         q[ 24*((i * y_cells + j) * z_cells + k) + 16 ] 
            = y[ size_t( (j + 0.75) * (y_pixels / y_cells) )];
         q[ 24*((i * y_cells + j) * z_cells + k) + 17 ] 
            = z[ size_t( (k + 0.25) * (z_pixels / z_cells) )];

         // atom at (3/4,1/4,3/4)
         q[ 24*((i * y_cells + j) * z_cells + k) + 18 ] 
            = x[ size_t( (i + 0.75) * (x_pixels / x_cells) )];
         q[ 24*((i * y_cells + j) * z_cells + k) + 19 ] 
            = y[ size_t( (j + 0.25) * (y_pixels / y_cells) )];
         q[ 24*((i * y_cells + j) * z_cells + k) + 20 ] 
            = z[ size_t( (k + 0.75) * (z_pixels / z_cells) )];

         // atom at (1/4,3/4,3/4)
         q[ 24*((i * y_cells + j) * z_cells + k) + 21 ] 
            = x[ size_t( (i + 0.25) * (x_pixels / x_cells) )];
         q[ 24*((i * y_cells + j) * z_cells + k) + 22 ] 
            = y[ size_t( (j + 0.75) * (y_pixels / y_cells) )];
         q[ 24*((i * y_cells + j) * z_cells + k) + 23 ] 
            = z[ size_t( (k + 0.75) * (z_pixels / z_cells) )];
      }

   return;
}

void TEM_NS::positions_3D_lattice_diamond_001(
   const size_t& x_cells,  // number of unit cells, x direction
   const size_t& y_cells,  // number of unit cells, x direction
   const size_t& z_cells,  // number of unit cells, x direction
   const double* const x,  // x domain
   const double& x_pixels, //number of points in discretized x domain 
   const double* const y,  // y domain
   const double& y_pixels, //number of points in discretized y domain 
   const double* const z,  // z domain
   const double& z_pixels, //number of points in discretized z domain 
   //double* const q         // lattice point positions (x,y,z), output
   std::list<double> q       // lattice point positions (x,y,z), output
   )
{
   double q_tmp[3];
   for( size_t i = 0; i < x_cells; i++)
      for( size_t j = 0; j < y_cells; j++)
         for( size_t k = 0; k < z_cells; k++) // 3-D
      {
         // (i,j,k) unit cell index
         // 8 atoms * 3 coordinates = 24 doubles per unit cell
         
         // atom at unit cell origin
         q.push_back( x[ size_t( i * x_pixels / x_cells )] );
         q.push_back( y[ size_t( j * y_pixels / y_cells )] );
         q.push_back( z[ size_t( k * z_pixels / z_cells )] );
         //q[ 24*((i * y_cells + j) * z_cells + k) + 0 ] 
         //   = 
         //   x[ size_t( i * x_pixels / x_cells )];
         //q[ 24*((i * y_cells + j) * z_cells + k) + 1 ] 
         //   = y[ size_t( j * y_pixels / y_cells )];
         //q[ 24*((i * y_cells + j) * z_cells + k) + 2 ] 
         //   = z[ size_t( k * z_pixels / z_cells )];

         // atom at (1/2,0,1/2)
         q.push_back( x[ size_t( (i+0.5) * x_pixels / x_cells )] );
         q.push_back( y[ size_t( (j) * y_pixels / y_cells )] );
         q.push_back( z[ size_t( (k+0.5) * z_pixels / z_cells )] );
         //q[ 24*((i * y_cells + j) * z_cells + k) + 3 ] 
         //   = x[ size_t( (i + 0.5) * (x_pixels / x_cells) )];
         //q[ 24*((i * y_cells + j) * z_cells + k) + 4 ] 
         //   = y[ size_t( j * (y_pixels / y_cells) )];
         //q[ 24*((i * y_cells + j) * z_cells + k) + 5 ] 
         //   = z[ size_t( (k + 0.5) * (z_pixels / z_cells) )];

         // atom at (1/4,1/4,1/4)
         q.push_back( x[ size_t( (i+0.25) * x_pixels / x_cells )] );
         q.push_back( y[ size_t( (j+0.25) * y_pixels / y_cells )] );
         q.push_back( z[ size_t( (k+0.25) * z_pixels / z_cells )] );
         //q[ 24*((i * y_cells + j) * z_cells + k) + 6 ] 
         //   = x[ size_t( (i + 0.25) * (x_pixels / x_cells) )];
         //q[ 24*((i * y_cells + j) * z_cells + k) + 7 ] 
         //   = y[ size_t( (j + 0.25) * (y_pixels / y_cells) )];
         //q[ 24*((i * y_cells + j) * z_cells + k) + 8 ] 
         //   = z[ size_t( (k + 0.25) * (z_pixels / z_cells) )];

         // atom at (0,1/2,1/2)
         q.push_back( x[ size_t( (i+0.0) * x_pixels / x_cells )] );
         q.push_back( y[ size_t( (j+0.5) * y_pixels / y_cells )] );
         q.push_back( z[ size_t( (k+0.5) * z_pixels / z_cells )] );
         //q[ 24*((i * y_cells + j) * z_cells + k) + 9 ] 
         //   = x[ size_t( i * (x_pixels / x_cells) )];
         //q[ 24*((i * y_cells + j) * z_cells + k) + 10 ] 
         //   = y[ size_t( (j + 0.5) * (y_pixels / y_cells) )];
         //q[ 24*((i * y_cells + j) * z_cells + k) + 11 ] 
         //   = z[ size_t( (k + 0.5) * (z_pixels / z_cells) )];

         // atom at (1/2,1/2,0)
         q.push_back( x[ size_t( (i+0.5) * x_pixels / x_cells )] );
         q.push_back( y[ size_t( (j+0.5) * y_pixels / y_cells )] );
         q.push_back( z[ size_t( (k+0.0) * z_pixels / z_cells )] );
         //q[ 24*((i * y_cells + j) * z_cells + k) + 12 ] 
         //   = x[ size_t( (i + 0.5) * (x_pixels / x_cells) )];
         //q[ 24*((i * y_cells + j) * z_cells + k) + 13 ] 
         //   = y[ size_t( (j + 0.5) * (y_pixels / y_cells) )];
         //q[ 24*((i * y_cells + j) * z_cells + k) + 14 ] 
         //   = z[ size_t( k * (z_pixels / z_cells) )];
         
         // atom at (3/4,3/4,1/4)
         q.push_back( x[ size_t( (i+0.75) * x_pixels / x_cells )] );
         q.push_back( y[ size_t( (j+0.75) * y_pixels / y_cells )] );
         q.push_back( z[ size_t( (k+0.25) * z_pixels / z_cells )] );
         //q[ 24*((i * y_cells + j) * z_cells + k) + 15 ] 
         //   = x[ size_t( (i + 0.75) * (x_pixels / x_cells) )];
         //q[ 24*((i * y_cells + j) * z_cells + k) + 16 ] 
         //   = y[ size_t( (j + 0.75) * (y_pixels / y_cells) )];
         //q[ 24*((i * y_cells + j) * z_cells + k) + 17 ] 
         //   = z[ size_t( (k + 0.25) * (z_pixels / z_cells) )];

         // atom at (3/4,1/4,3/4)
         q.push_back( x[ size_t( (i+0.75) * x_pixels / x_cells )] );
         q.push_back( y[ size_t( (j+0.25) * y_pixels / y_cells )] );
         q.push_back( z[ size_t( (k+0.75) * z_pixels / z_cells )] );
         //q[ 24*((i * y_cells + j) * z_cells + k) + 18 ] 
         //   = x[ size_t( (i + 0.75) * (x_pixels / x_cells) )];
         //q[ 24*((i * y_cells + j) * z_cells + k) + 19 ] 
         //   = y[ size_t( (j + 0.25) * (y_pixels / y_cells) )];
         //q[ 24*((i * y_cells + j) * z_cells + k) + 20 ] 
         //   = z[ size_t( (k + 0.75) * (z_pixels / z_cells) )];

         // atom at (1/4,3/4,3/4)
         q.push_back( x[ size_t( (i+0.25) * x_pixels / x_cells )] );
         q.push_back( y[ size_t( (j+0.75) * y_pixels / y_cells )] );
         q.push_back( z[ size_t( (k+0.75) * z_pixels / z_cells )] );
         //q[ 24*((i * y_cells + j) * z_cells + k) + 21 ] 
         //   = x[ size_t( (i + 0.25) * (x_pixels / x_cells) )];
         //q[ 24*((i * y_cells + j) * z_cells + k) + 22 ] 
         //   = y[ size_t( (j + 0.75) * (y_pixels / y_cells) )];
         //q[ 24*((i * y_cells + j) * z_cells + k) + 23 ] 
         //   = z[ size_t( (k + 0.75) * (z_pixels / z_cells) )];
      }

   return;
}


void TEM_NS::positions_3D_lattice_diamond_011(
      // 001: x=0 y=0 z=1; 101: x=1 y=0 z=1
   const size_t& x_cells,  // number of unit cells, x direction
   const size_t& y_cells,  // number of unit cells, x direction
   const size_t& z_cells,  // number of unit cells, x direction
   const double* const x,  // x domain
   const double& x_pixels, //number of points in discretized x domain 
   const double* const y,  // y domain
   const double& y_pixels, //number of points in discretized y domain 
   const double* const z,  // z domain
   const double& z_pixels, //number of points in discretized z domain 
   double* const q         // lattice point positions (x,y,z), output
   )
{
   for( size_t i = 0; i < x_cells; i++)
      for( size_t j = 0; j < y_cells; j++)
         for( size_t k = 0; k < z_cells; k++) // 3-D
      {
         // (i,j,k) unit cell index
         // 8 atoms * 3 coordinates = 24 doubles per unit cell
         
         // atom at unit cell origin
         q[ 12*((i * y_cells + j) * z_cells + k) + 0 ] 
            = x[ size_t( (i) * x_pixels / x_cells )];
         q[ 12*((i * y_cells + j) * z_cells + k) + 1 ] 
            = y[ size_t( (j) * y_pixels / y_cells )];
         q[ 12*((i * y_cells + j) * z_cells + k) + 2 ] 
            = z[ size_t( (k) * z_pixels / z_cells )];
         // atom at (1/2,1/2,1/2)
         q[ 12*((i * y_cells + j) * z_cells + k) + 3 ] 
            = x[ size_t( (i + 0.5) * (x_pixels / x_cells) )];
         q[ 12*((i * y_cells + j) * z_cells + k) + 4 ] 
            = y[ size_t( (j + 0.5) * (y_pixels / y_cells) )];
         q[ 12*((i * y_cells + j) * z_cells + k) + 5 ] 
            = z[ size_t( (k + 0.5) * (z_pixels / z_cells) )];
         // atom at (1/4,1/2,0)
         q[ 12*((i * y_cells + j) * z_cells + k) + 6 ] 
            = x[ size_t( (i + 0.25) * (x_pixels / x_cells) )];
         q[ 12*((i * y_cells + j) * z_cells + k) + 7 ] 
            = y[ size_t( (j + 0.5) * (y_pixels / y_cells) )];
         q[ 12*((i * y_cells + j) * z_cells + k) + 8 ] 
            = z[ size_t( (k) * (z_pixels / z_cells) )];
         // atom at (3/4,0,1/2)
         q[ 12*((i * y_cells + j) * z_cells + k) + 9 ] 
            = x[ size_t( (i + 0.75) * (x_pixels / x_cells) )];
         q[ 12*((i * y_cells + j) * z_cells + k) + 10 ] 
            = y[ size_t( (j ) * (y_pixels / y_cells) )];
         q[ 12*((i * y_cells + j) * z_cells + k) + 11 ] 
            = z[ size_t( (k + 0.5) * (z_pixels / z_cells) )];
      }

   return;
}

void TEM_NS::positions_3D_lattice_diamond_111_bad(
      // 001: x=0 y=0 z=1; 101: x=1 y=0 z=1
   const size_t& x_cells,  // number of unit cells, x direction
   const size_t& y_cells,  // number of unit cells, x direction
   const size_t& z_cells,  // number of unit cells, x direction
   const double* const x,  // x domain
   const double& x_pixels, //number of points in discretized x domain 
   const double* const y,  // y domain
   const double& y_pixels, //number of points in discretized y domain 
   const double* const z,  // z domain
   const double& z_pixels, //number of points in discretized z domain 
   double* const q         // lattice point positions (x,y,z), output
   )
{
   // First instantiate a 001 system, then rotate it into the 111 
   //  orientation while translating atoms across periodic boundaries.
   positions_3D_lattice_diamond_001(
            x_cells,  // number of unit cells, x direction
            y_cells,  // number of unit cells, x direction
            z_cells,  // number of unit cells, x direction
            x,  // x domain
            x_pixels, //number of points in discretized x domain 
            y,  // y domain
            y_pixels, //number of points in discretized y domain 
            z,  // z domain
            z_pixels, //number of points in discretized z domain 
            q         // lattice point positions (x,y,z), output
         );
   double x_boundary_pos = x[ size_t(x_pixels) - 1];
   double x_boundary_neg = x[ 0 ];
   double x_period = x_boundary_pos - x_boundary_neg;
   cout << "x_period : " << x_period << endl; // debug
   double y_boundary_pos = y[ size_t(y_pixels) - 1];
   double y_boundary_neg = y[ 0 ];
   double y_period = y_boundary_pos - y_boundary_neg;
   cout << "y_period : " << y_period << endl; // debug
   double z_boundary_pos = z[ size_t(z_pixels) - 1];
   double z_boundary_neg = z[ 0 ];
   double z_period = z_boundary_pos - z_boundary_neg;
   cout << "z_period : " << z_period << endl; // debug

   double rotationMatrix[3][3];
   const double sqrt2 = sqrt(2);
   const double atansqrt2 = atan(sqrt2);
   const double sinatansqrt2 = sin(atansqrt2);
   const double cosatansqrt2 = cos(atansqrt2);

   rotationMatrix[0][0] = sinatansqrt2;
   rotationMatrix[0][1] = -cosatansqrt2/sqrt2;
   rotationMatrix[0][2] = -cosatansqrt2/sqrt2;
   rotationMatrix[1][0] = 0;
   rotationMatrix[1][1] = 1/sqrt2;
   rotationMatrix[1][2] = -1/sqrt2;
   rotationMatrix[2][0] = cosatansqrt2;
   rotationMatrix[2][1] = sinatansqrt2/sqrt2;
   rotationMatrix[2][2] = sinatansqrt2/sqrt2;

   // For each atom, multiply its position vector by the rotation matrix
   //  and translate it if the result exceeds the periodic boundaries.
   const size_t total_population = 8 * x_cells * y_cells * z_cells;
   double x_tmp, y_tmp, z_tmp;
   for (size_t i=0; i<total_population; ++i)
   {

      x_tmp = q[3*i + 0];
      y_tmp = q[3*i + 1];
      z_tmp = q[3*i + 2];

      q[3*i + 0] = 
         rotationMatrix[0][0] * x_tmp
         + rotationMatrix[0][1] * y_tmp
         + rotationMatrix[0][2] * z_tmp;

      q[3*i + 1] = 
         rotationMatrix[1][0] * x_tmp
         + rotationMatrix[1][1] * y_tmp
         + rotationMatrix[1][2] * z_tmp;

      q[3*i + 2] = 
         rotationMatrix[2][0] * x_tmp
         + rotationMatrix[2][1] * y_tmp
         + rotationMatrix[2][2] * z_tmp;

      while ( q[3*i + 0] >= x_boundary_pos ) 
         q[3*i + 0] = q[3*i + 0] - x_period;

      while( q[3*i + 0] < x_boundary_neg ) 
         q[3*i + 0] = q[3*i + 0] + x_period;

      while( q[3*i + 1] >= y_boundary_pos ) 
         q[3*i + 1] = q[3*i + 1] - y_period;

      while( q[3*i + 1] < y_boundary_neg ) 
         q[3*i + 1] = q[3*i + 1] + y_period;

      while( q[3*i + 2] >= z_boundary_pos ) 
         q[3*i + 2] = q[3*i + 2] - z_period;

      while( q[3*i + 2] < z_boundary_neg ) 
         q[3*i + 2] = q[3*i + 2] + z_period;
   }

   return;
}

void TEM_NS::positions_3D_lattice_diamond_111(
      // 001: x=0 y=0 z=1; 101: x=1 y=0 z=1
   const size_t& x_cells,  // number of unit cells, x direction
   const size_t& y_cells,  // number of unit cells, x direction
   const size_t& z_cells,  // number of unit cells, x direction
   const double& lattice_constant,
   //const double* const x,  // x domain
   //const double& x_pixels, //number of points in discretized x domain 
   //const double* const y,  // y domain
   //const double& y_pixels, //number of points in discretized y domain 
   //const double* const z,  // z domain
   //const double& z_pixels, //number of points in discretized z domain 
   //std::list<size_t> Z,       // atom type
   std::list<double>& q       // lattice point positions (x,y,z), output
   )
{
   if( q.size() != 0 )//|| Z.size() != 0 )
   {
      cerr << "positions_3D_lattice_diamond_111() : "
         << "input list q is not empty but should have been." << endl;
      return ;
   }
   // First instantiate a 001 system, then rotate it into the 111 
   //  orientation while translating atoms across periodic boundaries.
   //positions_3D_lattice_diamond_001(
   //         x_cells,  // number of unit cells, x direction
   //         y_cells,  // number of unit cells, x direction
   //         z_cells,  // number of unit cells, x direction
   //         x,  // x domain
   //         x_pixels, //number of points in discretized x domain 
   //         y,  // y domain
   //         y_pixels, //number of points in discretized y domain 
   //         z,  // z domain
   //         z_pixels, //number of points in discretized z domain 
   //         q         // lattice point positions (x,y,z), output
   //      );

   //std::list<int> Z_basis;
   std::list<double> q_basis;
   std::list<double>::iterator q_itr;

   //double x_boundary_pos = x[ size_t(x_pixels) - 1];
   //double x_boundary_neg = x[ 0 ];
   //double x_period = x_boundary_pos - x_boundary_neg;
   //cout << "x_period : " << x_period << endl; // debug
   //double y_boundary_pos = y[ size_t(y_pixels) - 1];
   //double y_boundary_neg = y[ 0 ];
   //double y_period = y_boundary_pos - y_boundary_neg;
   //cout << "y_period : " << y_period << endl; // debug
   //double z_boundary_pos = z[ size_t(z_pixels) - 1];
   //double z_boundary_neg = z[ 0 ];
   //double z_period = z_boundary_pos - z_boundary_neg;
   //cout << "z_period : " << z_period << endl; // debug

   // Rotation Matrix:
   // [ 
   //  sin(atan(sqrt(2))) -cos(atan(sqrt(2)))/sqrt(2) -cos(atan(sqrt(2)))/sqrt(2)
   //          0               1/sqrt(2)          -1/sqrt(2)
   //  cos(atan(sqrt(2))) sin(atan(sqrt(2)))/sqrt(2) sin(atan(sqrt(2)))/sqrt(2)
   // ]
   
   double rotationMatrix[3][3];
   const double sqrt2 = sqrt(2);
   const double atansqrt2 = atan(sqrt2);
   const double sinatansqrt2 = sin(atansqrt2);
   const double cosatansqrt2 = cos(atansqrt2);

   rotationMatrix[0][0] = sinatansqrt2;
   rotationMatrix[0][1] = -cosatansqrt2/sqrt2;
   rotationMatrix[0][2] = -cosatansqrt2/sqrt2;
   rotationMatrix[1][0] = 0;
   rotationMatrix[1][1] = 1/sqrt2;
   rotationMatrix[1][2] = -1/sqrt2;
   rotationMatrix[2][0] = cosatansqrt2;
   rotationMatrix[2][1] = sinatansqrt2/sqrt2;
   rotationMatrix[2][2] = sinatansqrt2/sqrt2;

   // 001 basis vectors : 
   //  [0.0 0.0 0.0], [0.0 0.5 0.5], [0.5, 0.0, 0.5], [0.5, 0.5, 0.0]
   //   [0.25 0.25 0.25], [0.25 0.75 0.75], [0.75, 0.25, 0.75], 
   //   [0.75, 0.75, 0.25]
   //
   // 001 coordinates which will rotate into a 111 basis :
   // Ga:
   // [0.0 0.0 0.0], [1.0 0.0 0.0], [0.5, 0.5, 0.0], [1.0 1.0 0.0],
   // [1.5 1.0 -0.5]
   // [0.5 0.0 -0.5], [0.5 1.0 -0.5], [1.0 0.5 -0.5], [0.0 0.5 -0.5]
   // [1.5 0.5 0.0]
   // [0.5 0.5 -1.0],
   // [1.0 0.5 0.5]
   // As:
   // [0.25 0.25 0.25], [1.25 0.25 0.25], [0.75, 0.75, 0.25], [1.25 1.25 0.25],
   // [1.75 0.75 0.25]
   // [0.75 0.25 -0.25], [0.75 1.25 -0.25], [1.25 0.75 -0.25], [0.25 0.75 -0.25]
   // [1.75 1.25 -0.25]
   // [0.75 0.75 -0.75],
   // [1.25 0.75 0.75]
   // Instantiate the atoms listed directly above^
   // type 1
   q_basis.push_back(0.0); q_basis.push_back(0.0); q_basis.push_back(0.0);
   //Z_basis.push_back(1);
   q_basis.push_back(1.0); q_basis.push_back(0.0); q_basis.push_back(0.0);
   //Z_basis.push_back(1);
   q_basis.push_back(0.5); q_basis.push_back(0.5); q_basis.push_back(0.0);
   //Z_basis.push_back(1);
   q_basis.push_back(1.0); q_basis.push_back(1.0); q_basis.push_back(0.0);
   //Z_basis.push_back(1);
   q_basis.push_back(1.5); q_basis.push_back(1.0); q_basis.push_back(-0.5);
   //Z_basis.push_back(1);
   q_basis.push_back(0.5); q_basis.push_back(0.0); q_basis.push_back(-0.5);
   //Z_basis.push_back(1);
   q_basis.push_back(0.5); q_basis.push_back(1.0); q_basis.push_back(-0.5);
   //Z_basis.push_back(1);
   q_basis.push_back(1.0); q_basis.push_back(0.5); q_basis.push_back(-0.5);
   //Z_basis.push_back(1);
   q_basis.push_back(0.0); q_basis.push_back(0.5); q_basis.push_back(-0.5);
   //Z_basis.push_back(1);
   q_basis.push_back(1.5); q_basis.push_back(0.5); q_basis.push_back(0.0);
   //Z_basis.push_back(1);
   q_basis.push_back(0.5); q_basis.push_back(0.5); q_basis.push_back(-1.0);
   //Z_basis.push_back(1);
   q_basis.push_back(1.0); q_basis.push_back(0.5); q_basis.push_back(0.5);
   //Z_basis.push_back(1);
   // type 2
   q_basis.push_back(0.25); q_basis.push_back(0.25); 
      q_basis.push_back(0.25);
   //Z_basis.push_back(2);
   q_basis.push_back(1.25); q_basis.push_back(0.25); 
      q_basis.push_back(0.25);
   //Z_basis.push_back(2);
   q_basis.push_back(0.75); q_basis.push_back(0.75); 
      q_basis.push_back(0.25);
   //Z_basis.push_back(2);
   q_basis.push_back(1.25); q_basis.push_back(1.25);  
      q_basis.push_back(0.25);
   //Z_basis.push_back(2);
   q_basis.push_back(1.75); q_basis.push_back(0.75); 
      q_basis.push_back(0.25);
   //Z_basis.push_back(2);
   q_basis.push_back(0.75); q_basis.push_back(0.25); 
      q_basis.push_back(-0.25);
   //Z_basis.push_back(2);
   q_basis.push_back(0.75); q_basis.push_back(1.25); 
      q_basis.push_back(-0.25);
   //Z_basis.push_back(2);
   q_basis.push_back(1.25); q_basis.push_back(0.75); 
      q_basis.push_back(-0.25);
   //Z_basis.push_back(2);
   q_basis.push_back(0.25); q_basis.push_back(0.75); 
      q_basis.push_back(-0.25);
   //Z_basis.push_back(2);
   
   q_basis.push_back(1.75); q_basis.push_back(1.25); 
      q_basis.push_back(-0.25);
   
   //Z_basis.push_back(2);
   q_basis.push_back(0.75); q_basis.push_back(0.75); 
      q_basis.push_back(-0.75);
   //Z_basis.push_back(2);
   q_basis.push_back(1.25); q_basis.push_back(0.75); 
      q_basis.push_back(0.75);
   //Z_basis.push_back(2);

   // For each basis atom, multiply its position by the 001 oriented 
   //  lattice constant.
   for (q_itr = q_basis.begin();
         q_itr != q_basis.end(); ++q_itr)
   {
      (*q_itr) = (*q_itr) * lattice_constant;
   }
   
   // For each basis atom, multiply its position vector by the rotation 
   //  matrix

   double x_tmp, y_tmp, z_tmp;
   for (q_itr = q_basis.begin();
         q_itr != q_basis.end(); ++q_itr)
   {

      //x_tmp = q[3*i + 0];
      //y_tmp = q[3*i + 1];
      //z_tmp = q[3*i + 2];
      x_tmp = *q_itr; 
      ++q_itr;
      y_tmp = *q_itr;
      ++q_itr;
      z_tmp = *q_itr;
      --q_itr;
      --q_itr;

      *q_itr = 
         rotationMatrix[0][0] * x_tmp
         + rotationMatrix[0][1] * y_tmp
         + rotationMatrix[0][2] * z_tmp;
      ++q_itr;

      //q[3*i + 1] = 
      *q_itr = 
         rotationMatrix[1][0] * x_tmp
         + rotationMatrix[1][1] * y_tmp
         + rotationMatrix[1][2] * z_tmp;
      ++q_itr;

      //q[3*i + 2] = 
      *q_itr = 
         rotationMatrix[2][0] * x_tmp
         + rotationMatrix[2][1] * y_tmp
         + rotationMatrix[2][2] * z_tmp;
   }
   // Find the new rectangular lattice constants
   // x lattice constant is 6* the third atom's x-coordinate
   // y lattice constant is 4* the third atom's x-coordinate
   // z lattice constant is 3* the third atom's x-coordinate
   // TODO
   q_itr = q_basis.begin();
   ++q_itr; ++q_itr; ++q_itr; 
   ++q_itr; ++q_itr; ++q_itr; 
   // x component of the 3rd atom
   double lattice_constant_x;
   double lattice_constant_y ;
   double lattice_constant_z; 
   lattice_constant_x = 6.0 * (*q_itr);
   ++q_itr; // y component of the 3rd atom
   lattice_constant_y = 4.0 * (*q_itr);
   ++q_itr; // z component of the 3rd atom
   lattice_constant_z = 3.0 * (*q_itr);

   // begin debug
   cout << "lattice constants : " 
      << lattice_constant_x << ", "
      << lattice_constant_y << ", "
      << lattice_constant_z << endl;
   //cout << "basis atoms :" << endl;
   //for ( std::list<double>::iterator q_itr = q_basis.begin();
   //     q_itr != q_basis.end();  ++q_itr)
   //{
   //   cout << *q_itr << ", ";
   //   ++q_itr;
   //   cout << *q_itr << ", ";
   //   ++q_itr;
   //   cout << *q_itr << ", " << endl;
   //}
   //cout << endl;//debug
   // end debug

   
   //std::list<int> Z_itr;
   // For each cell requested, copy translated basis atoms to q
   for ( size_t i=0; i < x_cells; ++i)
      for ( size_t j=0; j < y_cells; ++j)
         for ( size_t k=0; k < z_cells; ++k)
         {
            //Z_itr = Z_basis.begin();
            for(  q_itr = q_basis.begin();
                  q_itr != q_basis.end();
                  ++q_itr)
            {
               q.push_back( (*q_itr) + i * lattice_constant_x );
               //cout << *q_itr << ", ";//debug
               ++q_itr;

               q.push_back( (*q_itr) + j * lattice_constant_y );
               //cout << *q_itr << ", ";//debug
               ++q_itr;

               q.push_back( (*q_itr) + k * lattice_constant_z );
               //cout << *q_itr << endl;//debug
               //++q_itr;
               //Z.push_back( *Z_itr );
            }
         }

   return;
}


void TEM_NS::positions_3D_lattice_zincblende_111(
      // 001: x=0 y=0 z=1; 101: x=1 y=0 z=1
   const size_t& x_cells,  // number of unit cells, x direction
   const size_t& y_cells,  // number of unit cells, x direction
   const size_t& z_cells,  // number of unit cells, x direction
   const double& lattice_constant,
   //const double* const x,  // x domain
   //const double& x_pixels, //number of points in discretized x domain 
   //const double* const y,  // y domain
   //const double& y_pixels, //number of points in discretized y domain 
   //const double* const z,  // z domain
   //const double& z_pixels, //number of points in discretized z domain 
   const size_t& Z1, const size_t& Z2, // atom type input
   std::list<size_t>& Z,       // atom types ordered the same as atoms in q
   std::list<double>& q       // lattice point positions (x,y,z), output
   )
{
   if( q.size() != 0 )//|| Z.size() != 0 )
   {
      cerr << "positions_3D_lattice_diamond_111() : "
         << "input list q is not empty but should have been." << endl;
      return ;
   }
   // First instantiate a 001 system, then rotate it into the 111 
   //  orientation while translating atoms across periodic boundaries.
   //positions_3D_lattice_diamond_001(
   //         x_cells,  // number of unit cells, x direction
   //         y_cells,  // number of unit cells, x direction
   //         z_cells,  // number of unit cells, x direction
   //         x,  // x domain
   //         x_pixels, //number of points in discretized x domain 
   //         y,  // y domain
   //         y_pixels, //number of points in discretized y domain 
   //         z,  // z domain
   //         z_pixels, //number of points in discretized z domain 
   //         q         // lattice point positions (x,y,z), output
   //      );

   std::list<size_t> Z_basis;
   std::list<double> q_basis;
   std::list<double>::iterator q_itr;

   //double x_boundary_pos = x[ size_t(x_pixels) - 1];
   //double x_boundary_neg = x[ 0 ];
   //double x_period = x_boundary_pos - x_boundary_neg;
   //cout << "x_period : " << x_period << endl; // debug
   //double y_boundary_pos = y[ size_t(y_pixels) - 1];
   //double y_boundary_neg = y[ 0 ];
   //double y_period = y_boundary_pos - y_boundary_neg;
   //cout << "y_period : " << y_period << endl; // debug
   //double z_boundary_pos = z[ size_t(z_pixels) - 1];
   //double z_boundary_neg = z[ 0 ];
   //double z_period = z_boundary_pos - z_boundary_neg;
   //cout << "z_period : " << z_period << endl; // debug

   // Rotation Matrix:
   // [ 
   //  sin(atan(sqrt(2))) -cos(atan(sqrt(2)))/sqrt(2) -cos(atan(sqrt(2)))/sqrt(2)
   //          0               1/sqrt(2)          -1/sqrt(2)
   //  cos(atan(sqrt(2))) sin(atan(sqrt(2)))/sqrt(2) sin(atan(sqrt(2)))/sqrt(2)
   // ]
   
   double rotationMatrix[3][3];
   const double sqrt2 = sqrt(2);
   const double atansqrt2 = atan(sqrt2);
   const double sinatansqrt2 = sin(atansqrt2);
   const double cosatansqrt2 = cos(atansqrt2);

   rotationMatrix[0][0] = sinatansqrt2;
   rotationMatrix[0][1] = -cosatansqrt2/sqrt2;
   rotationMatrix[0][2] = -cosatansqrt2/sqrt2;
   rotationMatrix[1][0] = 0;
   rotationMatrix[1][1] = 1/sqrt2;
   rotationMatrix[1][2] = -1/sqrt2;
   rotationMatrix[2][0] = cosatansqrt2;
   rotationMatrix[2][1] = sinatansqrt2/sqrt2;
   rotationMatrix[2][2] = sinatansqrt2/sqrt2;

   // 001 basis vectors : 
   //  [0.0 0.0 0.0], [0.0 0.5 0.5], [0.5, 0.0, 0.5], [0.5, 0.5, 0.0]
   //   [0.25 0.25 0.25], [0.25 0.75 0.75], [0.75, 0.25, 0.75], 
   //   [0.75, 0.75, 0.25]
   //
   // 001 coordinates which will rotate into a 111 basis :
   // Ga:
   // [0.0 0.0 0.0], [1.0 0.0 0.0], [0.5, 0.5, 0.0], [1.0 1.0 0.0],
   // [1.5 1.0 -0.5]
   // [0.5 0.0 -0.5], [0.5 1.0 -0.5], [1.0 0.5 -0.5], [0.0 0.5 -0.5]
   // [1.5 0.5 0.0]
   // [0.5 0.5 -1.0],
   // [1.0 0.5 0.5]
   // As:
   // [0.25 0.25 0.25], [1.25 0.25 0.25], [0.75, 0.75, 0.25], [1.25 1.25 0.25],
   // [1.75 0.75 0.25]
   // [0.75 0.25 -0.25], [0.75 1.25 -0.25], [1.25 0.75 -0.25], [0.25 0.75 -0.25]
   // [1.75 1.25 -0.25]
   // [0.75 0.75 -0.75],
   // [1.25 0.75 0.75]
   // Instantiate the atoms listed directly above^
   // type 1
   q_basis.push_back(0.0); q_basis.push_back(0.0); q_basis.push_back(0.0);
   Z_basis.push_back(1);
   q_basis.push_back(1.0); q_basis.push_back(0.0); q_basis.push_back(0.0);
   Z_basis.push_back(1);
   q_basis.push_back(0.5); q_basis.push_back(0.5); q_basis.push_back(0.0);
   Z_basis.push_back(1);
   q_basis.push_back(1.0); q_basis.push_back(1.0); q_basis.push_back(0.0);
   Z_basis.push_back(1);
   q_basis.push_back(1.5); q_basis.push_back(1.0); q_basis.push_back(-0.5);
   Z_basis.push_back(1);
   q_basis.push_back(0.5); q_basis.push_back(0.0); q_basis.push_back(-0.5);
   Z_basis.push_back(1);
   q_basis.push_back(0.5); q_basis.push_back(1.0); q_basis.push_back(-0.5);
   Z_basis.push_back(1);
   q_basis.push_back(1.0); q_basis.push_back(0.5); q_basis.push_back(-0.5);
   Z_basis.push_back(1);
   q_basis.push_back(0.0); q_basis.push_back(0.5); q_basis.push_back(-0.5);
   Z_basis.push_back(1);
   q_basis.push_back(1.5); q_basis.push_back(0.5); q_basis.push_back(0.0);
   Z_basis.push_back(1);
   q_basis.push_back(0.5); q_basis.push_back(0.5); q_basis.push_back(-1.0);
   Z_basis.push_back(1);
   q_basis.push_back(1.0); q_basis.push_back(0.5); q_basis.push_back(0.5);
   Z_basis.push_back(1);
   // type 2
   q_basis.push_back(0.25); q_basis.push_back(0.25); 
      q_basis.push_back(0.25);
   Z_basis.push_back(2);
   q_basis.push_back(1.25); q_basis.push_back(0.25); 
      q_basis.push_back(0.25);
   Z_basis.push_back(2);
   q_basis.push_back(0.75); q_basis.push_back(0.75); 
      q_basis.push_back(0.25);
   Z_basis.push_back(2);
   q_basis.push_back(1.25); q_basis.push_back(1.25);  
      q_basis.push_back(0.25);
   Z_basis.push_back(2);
   q_basis.push_back(1.75); q_basis.push_back(0.75); 
      q_basis.push_back(0.25);
   Z_basis.push_back(2);
   q_basis.push_back(0.75); q_basis.push_back(0.25); 
      q_basis.push_back(-0.25);
   Z_basis.push_back(2);
   q_basis.push_back(0.75); q_basis.push_back(1.25); 
      q_basis.push_back(-0.25);
   Z_basis.push_back(2);
   q_basis.push_back(1.25); q_basis.push_back(0.75); 
      q_basis.push_back(-0.25);
   Z_basis.push_back(2);
   q_basis.push_back(0.25); q_basis.push_back(0.75); 
      q_basis.push_back(-0.25);
   Z_basis.push_back(2);
   q_basis.push_back(1.75); q_basis.push_back(1.25); 
      q_basis.push_back(-0.25);
   Z_basis.push_back(2);
   q_basis.push_back(0.75); q_basis.push_back(0.75); 
      q_basis.push_back(-0.75);
   Z_basis.push_back(2);
   q_basis.push_back(1.25); q_basis.push_back(0.75); 
      q_basis.push_back(0.75);
   Z_basis.push_back(2);

   // For each basis atom, multiply its position by the 001 oriented 
   //  lattice constant.
   for (q_itr = q_basis.begin();
         q_itr != q_basis.end(); ++q_itr)
   {
      (*q_itr) = (*q_itr) * lattice_constant;
   }
   
   // For each basis atom, multiply its position vector by the rotation 
   //  matrix

   double x_tmp, y_tmp, z_tmp;
   for (q_itr = q_basis.begin();
         q_itr != q_basis.end(); ++q_itr)
   {

      //x_tmp = q[3*i + 0];
      //y_tmp = q[3*i + 1];
      //z_tmp = q[3*i + 2];
      x_tmp = *q_itr; 
      ++q_itr;
      y_tmp = *q_itr;
      ++q_itr;
      z_tmp = *q_itr;
      --q_itr;
      --q_itr;

      *q_itr = 
         rotationMatrix[0][0] * x_tmp
         + rotationMatrix[0][1] * y_tmp
         + rotationMatrix[0][2] * z_tmp;
      ++q_itr;

      //q[3*i + 1] = 
      *q_itr = 
         rotationMatrix[1][0] * x_tmp
         + rotationMatrix[1][1] * y_tmp
         + rotationMatrix[1][2] * z_tmp;
      ++q_itr;

      //q[3*i + 2] = 
      *q_itr = 
         rotationMatrix[2][0] * x_tmp
         + rotationMatrix[2][1] * y_tmp
         + rotationMatrix[2][2] * z_tmp;
   }
   // Find the new rectangular lattice constants
   // x lattice constant is 6* the third atom's x-coordinate
   // y lattice constant is 4* the third atom's x-coordinate
   // z lattice constant is 3* the third atom's x-coordinate
   // TODO
   q_itr = q_basis.begin();
   ++q_itr; ++q_itr; ++q_itr; 
   ++q_itr; ++q_itr; ++q_itr; 
   // x component of the 3rd atom
   double lattice_constant_x;
   double lattice_constant_y ;
   double lattice_constant_z; 
   lattice_constant_x = 6.0 * (*q_itr);
   ++q_itr; // y component of the 3rd atom
   lattice_constant_y = 4.0 * (*q_itr);
   ++q_itr; // z component of the 3rd atom
   lattice_constant_z = 3.0 * (*q_itr);

   // begin debug
   //cout << "lattice constants : " 
   //   << lattice_constant_x << ", "
   //   << lattice_constant_y << ", "
   //   << lattice_constant_z << endl;
   //cout << "basis atoms :" << endl;
   //for ( std::list<double>::iterator q_itr = q_basis.begin();
   //     q_itr != q_basis.end();  ++q_itr)
   //{
   //   cout << *q_itr << ", ";
   //   ++q_itr;
   //   cout << *q_itr << ", ";
   //   ++q_itr;
   //   cout << *q_itr << ", " << endl;
   //}
   //cout << endl;//debug
   // end debug

   
   std::list<size_t>::iterator Z_itr;
   // For each cell requested, copy translated basis atoms to q
   for ( size_t i=0; i < x_cells; ++i)
      for ( size_t j=0; j < y_cells; ++j)
         for ( size_t k=0; k < z_cells; ++k)
         {
            Z_itr = Z_basis.begin();
            for(  q_itr = q_basis.begin();
                  q_itr != q_basis.end();
                  ++q_itr)
            {
               q.push_back( (*q_itr) + i * lattice_constant_x );
               ++q_itr;

               q.push_back( (*q_itr) + j * lattice_constant_y );
               ++q_itr;

               q.push_back( (*q_itr) + k * lattice_constant_z );

               Z.push_back( *Z_itr );
               ++Z_itr;
            }
         }

   return;
}


void TEM_NS::positions_3D_lattice_zincblende_011(
      // 001: x=0 y=0 z=1; 101: x=1 y=0 z=1
   const size_t& x_cells,  // number of unit cells, x direction
   const size_t& y_cells,  // number of unit cells, x direction
   const size_t& z_cells,  // number of unit cells, x direction
   const double& lattice_constant,
   //const double* const x,  // x domain
   //const double& x_pixels, //number of points in discretized x domain 
   //const double* const y,  // y domain
   //const double& y_pixels, //number of points in discretized y domain 
   //const double* const z,  // z domain
   //const double& z_pixels, //number of points in discretized z domain 
   const size_t& Z1, const size_t& Z2, // atom type input
   std::list<size_t>& Z,       // atom types ordered the same as atoms in q
   std::list<double>& q       // lattice point positions (x,y,z), output
   )
{
   if( q.size() != 0 )//|| Z.size() != 0 )
   {
      cerr << "positions_3D_lattice_diamond_111() : "
         << "input list q is not empty but should have been." << endl;
      return ;
   }
   // First instantiate a 001 system, then rotate it into the 111 
   //  orientation while translating atoms across periodic boundaries.
   //positions_3D_lattice_diamond_001(
   //         x_cells,  // number of unit cells, x direction
   //         y_cells,  // number of unit cells, x direction
   //         z_cells,  // number of unit cells, x direction
   //         x,  // x domain
   //         x_pixels, //number of points in discretized x domain 
   //         y,  // y domain
   //         y_pixels, //number of points in discretized y domain 
   //         z,  // z domain
   //         z_pixels, //number of points in discretized z domain 
   //         q         // lattice point positions (x,y,z), output
   //      );

   std::list<size_t> Z_basis;
   std::list<double> q_basis;
   std::list<double>::iterator q_itr;

   //double x_boundary_pos = x[ size_t(x_pixels) - 1];
   //double x_boundary_neg = x[ 0 ];
   //double x_period = x_boundary_pos - x_boundary_neg;
   //cout << "x_period : " << x_period << endl; // debug
   //double y_boundary_pos = y[ size_t(y_pixels) - 1];
   //double y_boundary_neg = y[ 0 ];
   //double y_period = y_boundary_pos - y_boundary_neg;
   //cout << "y_period : " << y_period << endl; // debug
   //double z_boundary_pos = z[ size_t(z_pixels) - 1];
   //double z_boundary_neg = z[ 0 ];
   //double z_period = z_boundary_pos - z_boundary_neg;
   //cout << "z_period : " << z_period << endl; // debug

   // Rotation Matrix:
   // [ 
   //    1        0           0
   //    0        1/sqrt(2)   -1/sqrt(2)
   //    0        1/sqrt(2)   1/sqrt(2)
   // ]
   
   double rotationMatrix[3][3];
   const double sqrt2_inv = 1.0/sqrt(2);


   rotationMatrix[0][0] = 1.0;
   rotationMatrix[0][1] = 0.0;
   rotationMatrix[0][2] = 0.0;
   rotationMatrix[1][0] = 0.0;
   rotationMatrix[1][1] = sqrt2_inv;
   rotationMatrix[1][2] = -sqrt2_inv;
   rotationMatrix[2][0] = 0.0;
   rotationMatrix[2][1] = sqrt2_inv;
   rotationMatrix[2][2] = sqrt2_inv;

   // 001 basis vectors : 
   //  [0.0 0.0 0.0], [0.0 0.5 0.5], [0.5, 0.0, 0.5], [0.5, 0.5, 0.0]
   //   [0.25 0.25 0.25], [0.25 0.75 0.75], [0.75, 0.25, 0.75], 
   //   [0.75, 0.75, 0.25]
   //
   // 001 coordinates which will rotate into a 011 basis :
   // Ga:
   // [0.0 0.0 0.0], [0.5 0.5 0.0], [0.0, 0.5, 0.5], [0.5 1.0 0.5],
   // As:
   // [0.25 0.25 0.25], [0.75 0.25 -0.25], [0.25, 0.75, 0.75], [0.75 0.75 0.25],
   // Instantiate the atoms listed directly above^
   // type 1
   q_basis.push_back(0.0); q_basis.push_back(0.0); q_basis.push_back(0.0);
   Z_basis.push_back(1);
   q_basis.push_back(0.5); q_basis.push_back(0.5); q_basis.push_back(0.0);
   Z_basis.push_back(1);
   q_basis.push_back(0.0); q_basis.push_back(0.5); q_basis.push_back(0.5);
   Z_basis.push_back(1);
   q_basis.push_back(0.5); q_basis.push_back(1.0); q_basis.push_back(0.5);
   Z_basis.push_back(1);
   // type 2
   q_basis.push_back(0.25); q_basis.push_back(0.25); 
      q_basis.push_back(0.25);
   Z_basis.push_back(2);
   q_basis.push_back(0.75); q_basis.push_back(0.25); 
      q_basis.push_back(-0.25);
   Z_basis.push_back(2);
   q_basis.push_back(0.25); q_basis.push_back(0.75); 
      q_basis.push_back(0.75);
   Z_basis.push_back(2);
   q_basis.push_back(0.75); q_basis.push_back(0.75);  
      q_basis.push_back(0.25);
   Z_basis.push_back(2);

   // For each basis atom, multiply its position by the 001 oriented 
   //  lattice constant.
   for (q_itr = q_basis.begin();
         q_itr != q_basis.end(); ++q_itr)
   {
      (*q_itr) = (*q_itr) * lattice_constant;
   }
   
   // For each basis atom, multiply its position vector by the rotation 
   //  matrix

   double x_tmp, y_tmp, z_tmp;
   for (q_itr = q_basis.begin();
         q_itr != q_basis.end(); ++q_itr)
   {

      //x_tmp = q[3*i + 0];
      //y_tmp = q[3*i + 1];
      //z_tmp = q[3*i + 2];
      x_tmp = *q_itr; 
      ++q_itr;
      y_tmp = *q_itr;
      ++q_itr;
      z_tmp = *q_itr;
      --q_itr;
      --q_itr;

      *q_itr = 
         rotationMatrix[0][0] * x_tmp
         + rotationMatrix[0][1] * y_tmp
         + rotationMatrix[0][2] * z_tmp;
      ++q_itr;

      //q[3*i + 1] = 
      *q_itr = 
         rotationMatrix[1][0] * x_tmp
         + rotationMatrix[1][1] * y_tmp
         + rotationMatrix[1][2] * z_tmp;
      ++q_itr;

      //q[3*i + 2] = 
      *q_itr = 
         rotationMatrix[2][0] * x_tmp
         + rotationMatrix[2][1] * y_tmp
         + rotationMatrix[2][2] * z_tmp;
   }
   // Find the new rectangular lattice constants
   // x lattice constant is 2* the 2nd (0.5 0.5 0.0) atom's x-coordinate
   // y lattice constant is 2* the 2nd (0.5 0.5 0.0) atom's y-coordinate
   // z lattice constant is 2* the 3rd (0.0 0.5 0.5) atom's z-coordinate

   double lattice_constant_x;
   double lattice_constant_y ;
   double lattice_constant_z; 

   q_itr = q_basis.begin();
   ++q_itr; ++q_itr; ++q_itr; 
   // x component of the 2nd atom
   lattice_constant_x = 2.0 * (*q_itr);
   ++q_itr; // y component of the 2nd atom
   lattice_constant_y = 2.0 * (*q_itr);
   ++q_itr; 
   ++q_itr; ++q_itr; ++q_itr;
   // z component of the 3rd atom
   lattice_constant_z = 2.0 * (*q_itr);

   // begin debug
   //cout << "lattice constants : " 
   //   << lattice_constant_x << ", "
   //   << lattice_constant_y << ", "
   //   << lattice_constant_z << endl;
   //cout << "basis atoms :" << endl;
   //for ( std::list<double>::iterator q_itr = q_basis.begin();
   //     q_itr != q_basis.end();  ++q_itr)
   //{
   //   cout << *q_itr << ", ";
   //   ++q_itr;
   //   cout << *q_itr << ", ";
   //   ++q_itr;
   //   cout << *q_itr << ", " << endl;
   //}
   //cout << endl;//debug
   // end debug

   
   std::list<size_t>::iterator Z_itr;
   // For each cell requested, copy translated basis atoms to q
   for ( size_t i=0; i < x_cells; ++i)
      for ( size_t j=0; j < y_cells; ++j)
         for ( size_t k=0; k < z_cells; ++k)
         {
            Z_itr = Z_basis.begin();
            for(  q_itr = q_basis.begin();
                  q_itr != q_basis.end();
                  ++q_itr)
            {
               q.push_back( (*q_itr) + i * lattice_constant_x );
               ++q_itr;

               q.push_back( (*q_itr) + j * lattice_constant_y );
               ++q_itr;

               q.push_back( (*q_itr) + k * lattice_constant_z );

               Z.push_back( *Z_itr );
               ++Z_itr;
            }
         }

   return;
}

void TEM_NS::positions_3D_lattice_zincblende_001(
      // 001: x=0 y=0 z=1; 101: x=1 y=0 z=1
   const size_t& x_cells,  // number of unit cells, x direction
   const size_t& y_cells,  // number of unit cells, x direction
   const size_t& z_cells,  // number of unit cells, x direction
   const double& lattice_constant,
   //const double* const x,  // x domain
   //const double& x_pixels, //number of points in discretized x domain 
   //const double* const y,  // y domain
   //const double& y_pixels, //number of points in discretized y domain 
   //const double* const z,  // z domain
   //const double& z_pixels, //number of points in discretized z domain 
   const size_t& Z1, const size_t& Z2, // atom type input
   std::list<size_t>& Z,       // atom types ordered the same as atoms in q
   std::list<double>& q       // lattice point positions (x,y,z), output
   )
{
   if( q.size() != 0 )//|| Z.size() != 0 )
   {
      cerr << "positions_3D_lattice_diamond_111() : "
         << "input list q is not empty but should have been." << endl;
      return ;
   }
   // First instantiate a 001 system, then rotate it into the 111 
   //  orientation while translating atoms across periodic boundaries.
   //positions_3D_lattice_diamond_001(
   //         x_cells,  // number of unit cells, x direction
   //         y_cells,  // number of unit cells, x direction
   //         z_cells,  // number of unit cells, x direction
   //         x,  // x domain
   //         x_pixels, //number of points in discretized x domain 
   //         y,  // y domain
   //         y_pixels, //number of points in discretized y domain 
   //         z,  // z domain
   //         z_pixels, //number of points in discretized z domain 
   //         q         // lattice point positions (x,y,z), output
   //      );

   std::list<size_t> Z_basis;
   std::list<double> q_basis;
   std::list<double>::iterator q_itr;

   //double x_boundary_pos = x[ size_t(x_pixels) - 1];
   //double x_boundary_neg = x[ 0 ];
   //double x_period = x_boundary_pos - x_boundary_neg;
   //cout << "x_period : " << x_period << endl; // debug
   //double y_boundary_pos = y[ size_t(y_pixels) - 1];
   //double y_boundary_neg = y[ 0 ];
   //double y_period = y_boundary_pos - y_boundary_neg;
   //cout << "y_period : " << y_period << endl; // debug
   //double z_boundary_pos = z[ size_t(z_pixels) - 1];
   //double z_boundary_neg = z[ 0 ];
   //double z_period = z_boundary_pos - z_boundary_neg;
   //cout << "z_period : " << z_period << endl; // debug

   // Rotation Matrix:
   // [ 
   //  sin(atan(sqrt(2))) -cos(atan(sqrt(2)))/sqrt(2) -cos(atan(sqrt(2)))/sqrt(2)
   //          0               1/sqrt(2)          -1/sqrt(2)
   //  cos(atan(sqrt(2))) sin(atan(sqrt(2)))/sqrt(2) sin(atan(sqrt(2)))/sqrt(2)
   // ]
   
   //double rotationMatrix[3][3];
   //const double sqrt2 = sqrt(2);
   //const double atansqrt2 = atan(sqrt2);
   //const double sinatansqrt2 = sin(atansqrt2);
   //const double cosatansqrt2 = cos(atansqrt2);

   //rotationMatrix[0][0] = sinatansqrt2;
   //rotationMatrix[0][1] = -cosatansqrt2/sqrt2;
   //rotationMatrix[0][2] = -cosatansqrt2/sqrt2;
   //rotationMatrix[1][0] = 0;
   //rotationMatrix[1][1] = 1/sqrt2;
   //rotationMatrix[1][2] = -1/sqrt2;
   //rotationMatrix[2][0] = cosatansqrt2;
   //rotationMatrix[2][1] = sinatansqrt2/sqrt2;
   //rotationMatrix[2][2] = sinatansqrt2/sqrt2;

   // 001 basis vectors : 
   // Ga:  [0.0 0.0 0.0], [0.0 0.5 0.5], [0.5, 0.0, 0.5], [0.5, 0.5, 0.0]
   // As:  [0.25 0.25 0.25], [0.25 0.75 0.75], [0.75, 0.25, 0.75], 
   //   [0.75, 0.75, 0.25]
   // Instantiate the atoms
   // type 1
   q_basis.push_back(0.0); q_basis.push_back(0.0); q_basis.push_back(0.0);
   Z_basis.push_back(1);
   q_basis.push_back(0.0); q_basis.push_back(0.5); q_basis.push_back(0.5);
   Z_basis.push_back(1);
   q_basis.push_back(0.5); q_basis.push_back(0.0); q_basis.push_back(0.5);
   Z_basis.push_back(1);
   q_basis.push_back(0.5); q_basis.push_back(0.5); q_basis.push_back(0.0);
   Z_basis.push_back(1);
   // type 2
   q_basis.push_back(0.25); q_basis.push_back(0.25); 
      q_basis.push_back(0.25);
   Z_basis.push_back(2);
   q_basis.push_back(0.25); q_basis.push_back(0.75); 
      q_basis.push_back(0.75);
   Z_basis.push_back(2);
   q_basis.push_back(0.75); q_basis.push_back(0.25);  
      q_basis.push_back(0.75);
   Z_basis.push_back(2);
   q_basis.push_back(0.75); q_basis.push_back(0.75); 
      q_basis.push_back(0.25);
   Z_basis.push_back(2);

   // For each basis atom, multiply its position by the 001 oriented 
   //  lattice constant.
   for (q_itr = q_basis.begin();
         q_itr != q_basis.end(); ++q_itr)
   {
      (*q_itr) = (*q_itr) * lattice_constant;
   }
   
   // For each basis atom, multiply its position vector by the rotation 
   //  matrix

   //double x_tmp, y_tmp, z_tmp;
   //for (q_itr = q_basis.begin();
   //      q_itr != q_basis.end(); ++q_itr)
   //{

   //   //x_tmp = q[3*i + 0];
   //   //y_tmp = q[3*i + 1];
   //   //z_tmp = q[3*i + 2];
   //   x_tmp = *q_itr; 
   //   ++q_itr;
   //   y_tmp = *q_itr;
   //   ++q_itr;
   //   z_tmp = *q_itr;
   //   --q_itr;
   //   --q_itr;

   //   *q_itr = 
   //      rotationMatrix[0][0] * x_tmp
   //      + rotationMatrix[0][1] * y_tmp
   //      + rotationMatrix[0][2] * z_tmp;
   //   ++q_itr;

   //   //q[3*i + 1] = 
   //   *q_itr = 
   //      rotationMatrix[1][0] * x_tmp
   //      + rotationMatrix[1][1] * y_tmp
   //      + rotationMatrix[1][2] * z_tmp;
   //   ++q_itr;

   //   //q[3*i + 2] = 
   //   *q_itr = 
   //      rotationMatrix[2][0] * x_tmp
   //      + rotationMatrix[2][1] * y_tmp
   //      + rotationMatrix[2][2] * z_tmp;
   //}
   // Find the new rectangular lattice constants
   // x lattice constant is 6* the third atom's x-coordinate
   // y lattice constant is 4* the third atom's x-coordinate
   // z lattice constant is 3* the third atom's x-coordinate
   q_itr = q_basis.begin();
   ++q_itr; ++q_itr; ++q_itr; 
   ++q_itr; ++q_itr; ++q_itr; 
   // x component of the 3rd atom
   //double lattice_constant_x;
   //double lattice_constant_y ;
   //double lattice_constant_z; 
   //lattice_constant_x = 6.0 * (*q_itr);
   //++q_itr; // y component of the 3rd atom
   //lattice_constant_y = 4.0 * (*q_itr);
   //++q_itr; // z component of the 3rd atom
   //lattice_constant_z = 3.0 * (*q_itr);

   // begin debug
   //cout << "lattice constants : " 
   //   << lattice_constant_x << ", "
   //   << lattice_constant_y << ", "
   //   << lattice_constant_z << endl;
   //cout << "basis atoms :" << endl;
   //for ( std::list<double>::iterator q_itr = q_basis.begin();
   //     q_itr != q_basis.end();  ++q_itr)
   //{
   //   cout << *q_itr << ", ";
   //   ++q_itr;
   //   cout << *q_itr << ", ";
   //   ++q_itr;
   //   cout << *q_itr << ", " << endl;
   //}
   //cout << endl;//debug
   // end debug

   
   std::list<size_t>::iterator Z_itr;
   // For each cell requested, copy translated basis atoms to q
   for ( size_t i=0; i < x_cells; ++i)
      for ( size_t j=0; j < y_cells; ++j)
         for ( size_t k=0; k < z_cells; ++k)
         {
            Z_itr = Z_basis.begin();
            for(  q_itr = q_basis.begin();
                  q_itr != q_basis.end();
                  ++q_itr)
            {
               q.push_back( (*q_itr) + i * lattice_constant );
               ++q_itr;

               q.push_back( (*q_itr) + j * lattice_constant );
               ++q_itr;

               q.push_back( (*q_itr) + k * lattice_constant );

               Z.push_back( *Z_itr );
               ++Z_itr;
            }
         }

   return;
}


/*
void TEM_NS::positions_3D_lattice_zincblende_111(
      // 001: x=0 y=0 z=1; 101: x=1 y=0 z=1
   const size_t& x_cells,  // number of unit cells, x direction
   const size_t& y_cells,  // number of unit cells, x direction
   const size_t& z_cells,  // number of unit cells, x direction
   const double& lattice_constant,
   //const double* const x,  // x domain
   //const double& x_pixels, //number of points in discretized x domain 
   //const double* const y,  // y domain
   //const double& y_pixels, //number of points in discretized y domain 
   //const double* const z,  // z domain
   //const double& z_pixels, //number of points in discretized z domain 
   std::list<size_t>& Z,       // atom types
   std::list<double>& q       // lattice point positions (x,y,z), output
   )
{
   if( q.size() != 0 || Z.size() != 0 )
   {
      cerr << "positions_3D_lattice_diamond_111() : "
         << "input list q is not empty but should have been." << endl;
      return;
   }
   // First instantiate a 001 system, then rotate it into the 111 
   //  orientation while translating atoms across periodic boundaries.
   //positions_3D_lattice_diamond_001(
   //         x_cells,  // number of unit cells, x direction
   //         y_cells,  // number of unit cells, x direction
   //         z_cells,  // number of unit cells, x direction
   //         x,  // x domain
   //         x_pixels, //number of points in discretized x domain 
   //         y,  // y domain
   //         y_pixels, //number of points in discretized y domain 
   //         z,  // z domain
   //         z_pixels, //number of points in discretized z domain 
   //         q         // lattice point positions (x,y,z), output
   //      );

   std::list<int> Z_basis;
   std::list<double> q_basis;

   //double x_boundary_pos = x[ size_t(x_pixels) - 1];
   //double x_boundary_neg = x[ 0 ];
   //double x_period = x_boundary_pos - x_boundary_neg;
   //cout << "x_period : " << x_period << endl; // debug
   //double y_boundary_pos = y[ size_t(y_pixels) - 1];
   //double y_boundary_neg = y[ 0 ];
   //double y_period = y_boundary_pos - y_boundary_neg;
   //cout << "y_period : " << y_period << endl; // debug
   //double z_boundary_pos = z[ size_t(z_pixels) - 1];
   //double z_boundary_neg = z[ 0 ];
   //double z_period = z_boundary_pos - z_boundary_neg;
   //cout << "z_period : " << z_period << endl; // debug

   // Rotation Matrix:
   // [ 
   //  sin(atan(sqrt(2))) -cos(atan(sqrt(2)))/sqrt(2) -cos(atan(sqrt(2)))/sqrt(2)
   //          0               1/sqrt(2)          -1/sqrt(2)
   //  cos(atan(sqrt(2))) sin(atan(sqrt(2)))/sqrt(2) sin(atan(sqrt(2)))/sqrt(2)
   // ]
   
   double rotationMatrix[3][3];
   const double sqrt2 = sqrt(2);
   const double atansqrt2 = atan(sqrt2);
   const double sinatansqrt2 = sin(atansqrt2);
   const double cosatansqrt2 = cos(atansqrt2);

   rotationMatrix[0][0] = sinatansqrt2;
   rotationMatrix[0][1] = -cosatansqrt2/sqrt2;
   rotationMatrix[0][2] = -cosatansqrt2/sqrt2;
   rotationMatrix[1][0] = 0;
   rotationMatrix[1][1] = 1/sqrt2;
   rotationMatrix[1][2] = -1/sqrt2;
   rotationMatrix[2][0] = cosatansqrt2;
   rotationMatrix[2][1] = sinatansqrt2/sqrt2;
   rotationMatrix[2][2] = sinatansqrt2/sqrt2;

   // 001 basis vectors : 
   //  [0.0 0.0 0.0], [0.0 0.5 0.5], [0.5, 0.0, 0.5], [0.5, 0.5, 0.0]
   //   [0.25 0.25 0.25], [0.25 0.75 0.75], [0.75, 0.25, 0.75], 
   //   [0.75, 0.75, 0.25]
   //
   // 001 coordinates which will rotate into a 111 basis :
   // Ga:
   // [0.0 0.0 0.0], [1.0 0.0 0.0], [0.5, 0.5, 0.0], [1.0 1.0 0.0],
   // [1.5 1.0 -0.5]
   // [0.5 0.0 -0.5], [0.5 1.0 -0.5], [1.0 0.5 -0.5], [0.0 0.5 -0.5]
   // [1.5 0.5 0.0]
   // [0.5 0.5 -1.0],
   // [1.0 0.5 0.5]
   // As:
   // [0.25 0.25 0.25], [1.25 0.25 0.25], [0.75, 0.75, 0.25], [1.25 1.25 0.25],
   // [1.75 0.75 0.25]
   // [0.75 0.25 -0.25], [0.75 1.25 -0.25], [1.25 0.75 -0.25], [0.25 0.75 -0.25]
   // [1.75 1.25 -0.25]
   // [0.75 0.75 -0.75],
   // [1.25 0.75 0.75]
   // Instantiate the atoms listed directly above^
   // type 1
   q_basis.push_back(0.0); q_basis.push_back(0.0); q_basis.push_back(0.0);
   Z_basis.push_back(1);
   q_basis.push_back(1.0); q_basis.push_back(0.0); q_basis.push_back(0.0);
   Z_basis.push_back(1);
   q_basis.push_back(0.5); q_basis.push_back(0.5); q_basis.push_back(0.0);
   Z_basis.push_back(1);
   q_basis.push_back(1.0); q_basis.push_back(1.0); q_basis.push_back(0.0);
   Z_basis.push_back(1);
   q_basis.push_back(1.5); q_basis.push_back(1.0); q_basis.push_back(-0.5);
   Z_basis.push_back(1);
   q_basis.push_back(0.5); q_basis.push_back(0.0); q_basis.push_back(-0.5);
   Z_basis.push_back(1);
   q_basis.push_back(0.5); q_basis.push_back(1.0); q_basis.push_back(-0.5);
   Z_basis.push_back(1);
   q_basis.push_back(1.0); q_basis.push_back(0.5); q_basis.push_back(-0.5);
   Z_basis.push_back(1);
   q_basis.push_back(0.0); q_basis.push_back(0.5); q_basis.push_back(-0.5);
   Z_basis.push_back(1);
   q_basis.push_back(1.5); q_basis.push_back(0.5); q_basis.push_back(0.0);
   Z_basis.push_back(1);
   q_basis.push_back(0.5); q_basis.push_back(0.5); q_basis.push_back(-1.0);
   Z_basis.push_back(1);
   q_basis.push_back(1.0); q_basis.push_back(0.5); q_basis.push_back(0.5);
   Z_basis.push_back(1);
   // type 2
   q_basis.push_back(0.25); q_basis.push_back(0.25); 
      q_basis.push_back(0.25);
   Z_basis.push_back(2);
   q_basis.push_back(1.25); q_basis.push_back(0.25); 
      q_basis.push_back(0.25);
   Z_basis.push_back(2);
   q_basis.push_back(0.75); q_basis.push_back(0.75); 
      q_basis.push_back(0.25);
   Z_basis.push_back(2);
   q_basis.push_back(1.25); q_basis.push_back(1.25);  
      q_basis.push_back(0.25);
   Z_basis.push_back(2);
   q_basis.push_back(1.75); q_basis.push_back(0.75); 
      q_basis.push_back(0.25);
   Z_basis.push_back(2);
   q_basis.push_back(0.75); q_basis.push_back(0.25); 
      q_basis.push_back(-0.25);
   Z_basis.push_back(2);
   q_basis.push_back(0.75); q_basis.push_back(1.25); 
      q_basis.push_back(-0.25);
   Z_basis.push_back(2);
   q_basis.push_back(1.25); q_basis.push_back(0.75); 
      q_basis.push_back(-0.25);
   Z_basis.push_back(2);
   q_basis.push_back(0.25); q_basis.push_back(0.75); 
      q_basis.push_back(-0.25);
   Z_basis.push_back(2);
   q_basis.push_back(1.75); q_basis.push_back(1.25); 
      q_basis.push_back(-0.25);
   Z_basis.push_back(2);
   q_basis.push_back(0.75); q_basis.push_back(0.75); 
      q_basis.push_back(-0.75);
   Z_basis.push_back(2);
   q_basis.push_back(1.25); q_basis.push_back(0.75); 
      q_basis.push_back(0.75);
   Z_basis.push_back(2);
   // For each basis atom, multiply its position by the 001 oriented 
   //  lattice constant.
   for (std::list<double>::iterator q_itr = q_basis.begin();
         q_itr != q_basis.end(); ++q_itr)
   {
      (*q_itr) = (*q_itr) * lattice_constant;
   }
   
   // For each basis atom, multiply its position vector by the rotation 
   //  matrix
   double x_tmp, y_tmp, z_tmp;
   for (std::list<double>::iterator q_itr = q_basis.begin();
         q_itr != q_basis.end(); ++q_itr)
   {

      //x_tmp = q[3*i + 0];
      //y_tmp = q[3*i + 1];
      //z_tmp = q[3*i + 2];
      x_tmp = *q_itr; 
      ++q_itr;
      y_tmp = *q_itr;
      ++q_itr;
      z_tmp = *q_itr;
      --q_itr;
      --q_itr;

      *q_itr = 
         rotationMatrix[0][0] * x_tmp
         + rotationMatrix[0][1] * y_tmp
         + rotationMatrix[0][2] * z_tmp;
      ++q_itr;

      //q[3*i + 1] = 
      *q_itr = 
         rotationMatrix[1][0] * x_tmp
         + rotationMatrix[1][1] * y_tmp
         + rotationMatrix[1][2] * z_tmp;
      ++q_itr;

      //q[3*i + 2] = 
      *q_itr = 
         rotationMatrix[2][0] * x_tmp
         + rotationMatrix[2][1] * y_tmp
         + rotationMatrix[2][2] * z_tmp;
   }
   
   std::list<int>::iterator Z_itr;
   // For each cell requested, copy translated basis atoms to q
   for ( size_t i=0; i < x_cells; ++i)
      for ( size_t j=0; j < y_cells; ++j)
         for ( size_t k=0; k < z_cells; ++k)
         {
            Z_itr = Z_basis.begin();
            for(  std::list<double>::iterator q_itr = q_basis.begin();
                  q_itr != q_basis.end();
                  ++q_itr)
            {
               q.push_back( (*q_itr) + i * x_cells );
               ++q_itr;
               q.push_back( (*q_itr) + j * y_cells );
               ++q_itr;
               q.push_back( (*q_itr) + k * z_cells );
               Z.push_back( *Z_itr );
            }
         }

   return;
}
*/

void TEM_NS::debug_output_complex_fftw_operand(
      const fftw_complex* const complex_fftw_operand,
      const int& step,
      const ptrdiff_t& local_alloc_size_fftw,
      const ptrdiff_t& Nx_local,
      const ptrdiff_t& Nx,
      const ptrdiff_t& Ny,
      const size_t& resolutionUnit,
      const double& xResolution,
      const double& yResolution,
      const string& outFileName_prefix,
      const int* const psi_mag_strides,       // const int sendcounts[]
      const int* const psi_mag_displacements, // const int displacements[]
      const int& mynode,
      const int& rootnode,
      MPI_Comm comm )
{

   writeImage imageWriter;

   double* complex_fftw_operand_mag;
   double* complex_fftw_operand_re;
   double* complex_fftw_operand_im;
   complex_fftw_operand_mag = new double[ local_alloc_size_fftw ];
   complex_fftw_operand_re = new double[ local_alloc_size_fftw ];
   complex_fftw_operand_im = new double[ local_alloc_size_fftw ];
 
 
   double min_complex_fftw_operand_mag, max_complex_fftw_operand_mag;
   min_complex_fftw_operand_mag
      = sqrt( pow( complex_fftw_operand[0][0], 2)
            + pow( complex_fftw_operand[0][1], 2) );
   max_complex_fftw_operand_mag = min_complex_fftw_operand_mag;


   double min_complex_fftw_operand_re, max_complex_fftw_operand_re;
   min_complex_fftw_operand_re = complex_fftw_operand[0][0];
   max_complex_fftw_operand_re = min_complex_fftw_operand_re;
 
   double min_complex_fftw_operand_im, max_complex_fftw_operand_im;
   min_complex_fftw_operand_re = complex_fftw_operand[0][0];
   max_complex_fftw_operand_re = min_complex_fftw_operand_re;
 
 
   // Determine the max and min intensities
   for ( ptrdiff_t i = 0; i < local_alloc_size_fftw ; i++)
   {  
      complex_fftw_operand_re[i] = complex_fftw_operand[i][0];
      complex_fftw_operand_im[i] = complex_fftw_operand[i][1];
      
      complex_fftw_operand_mag[i] 
         = sqrt( pow( complex_fftw_operand[i][0], 2) 
               + pow( complex_fftw_operand[i][1], 2) );
      if( complex_fftw_operand_mag[i] > max_complex_fftw_operand_mag )
         max_complex_fftw_operand_mag = complex_fftw_operand_mag[i]; 
      if( complex_fftw_operand_mag[i] < min_complex_fftw_operand_mag )
         min_complex_fftw_operand_mag = complex_fftw_operand_mag[i];

      if( complex_fftw_operand_re[i] > max_complex_fftw_operand_re )
         max_complex_fftw_operand_re = complex_fftw_operand_re[i]; 
      if( complex_fftw_operand_re[i] < min_complex_fftw_operand_re )
         min_complex_fftw_operand_re = complex_fftw_operand_re[i];

      if( complex_fftw_operand_im[i] > max_complex_fftw_operand_im ) 
         max_complex_fftw_operand_im = complex_fftw_operand_im[i]; 
      if( complex_fftw_operand_im[i] < min_complex_fftw_operand_im ) 
         min_complex_fftw_operand_im = complex_fftw_operand_im[i];
   }
 
   //imageWriter.write_tif_grayscale16(
   //      complex_fftw_operand_mag,
   //      min_complex_fftw_operand_mag, max_complex_fftw_operand_mag,
   //      Nx_local, Ny,
   //      outFileName_prefix + "_" + TEM_NS::to_string(step) 
   //         + "_complex_fftw_operand_mag_output_node"
   //         + TEM_NS::to_string(mynode) );
 
   //imageWriter.write_tif_grayscale16_logscale(
   //      complex_fftw_operand_mag,
   //      min_complex_fftw_operand_mag, max_complex_fftw_operand_mag,
   //      Nx_local, Ny,
   //      outFileName_prefix + "_" + TEM_NS::to_string(step) 
   //         + "_complex_fftw_operand_mag_output_node"
   //         + TEM_NS::to_string(mynode)
   //         + "_logscale" );
 
   double* complex_fftw_operand_mag_joined;
   double* complex_fftw_operand_re_joined;
   double* complex_fftw_operand_im_joined;
   if ( mynode == rootnode )
   {
      complex_fftw_operand_mag_joined = new double[Nx * Ny];
      complex_fftw_operand_re_joined = new double[Nx * Ny];
      complex_fftw_operand_im_joined = new double[Nx * Ny];
   }

 
   MPI_Gatherv( 
         complex_fftw_operand_mag, 
         Nx_local* Ny, //local_alloc_size_fftw, 
         MPI_DOUBLE,
         complex_fftw_operand_mag_joined, 
         psi_mag_strides,       // recvcount[]
         psi_mag_displacements, // displs[]
         MPI_DOUBLE,
         rootnode, comm);
 
   MPI_Gatherv( 
         complex_fftw_operand_re, 
         Nx_local* Ny, //local_alloc_size_fftw, 
         MPI_DOUBLE,
         complex_fftw_operand_re_joined, 
         psi_mag_strides,       // recvcount[]
         psi_mag_displacements, // displs[]
         MPI_DOUBLE,
         rootnode, comm);
 
   MPI_Gatherv( 
         complex_fftw_operand_im, 
         Nx_local* Ny, //local_alloc_size_fftw, 
         MPI_DOUBLE,
         complex_fftw_operand_im_joined, 
         psi_mag_strides,       // recvcount[]
         psi_mag_displacements, // displs[]
         MPI_DOUBLE,
         rootnode, comm);
   //MPI_Gather( complex_fftw_operand_mag, local_alloc_size_fftw, MPI_DOUBLE,
   //      complex_fftw_operand_mag_joined, local_alloc_size_fftw, MPI_DOUBLE,
   //      rootnode, comm);
 
   //MPI_Gather( complex_fftw_operand_re, local_alloc_size_fftw, MPI_DOUBLE,
   //      complex_fftw_operand_re_joined, local_alloc_size_fftw, MPI_DOUBLE,
   //      rootnode, comm);
 
   //MPI_Gather( complex_fftw_operand_im, local_alloc_size_fftw, MPI_DOUBLE,
   //      complex_fftw_operand_im_joined, local_alloc_size_fftw, MPI_DOUBLE,
   //      rootnode, comm);
 
   double max_complex_fftw_operand_mag_joined;
   double min_complex_fftw_operand_mag_joined;
   double max_complex_fftw_operand_re_joined;
   double min_complex_fftw_operand_re_joined;
   double max_complex_fftw_operand_im_joined;
   double min_complex_fftw_operand_im_joined;

   MPI_Reduce( &max_complex_fftw_operand_mag, 
               &max_complex_fftw_operand_mag_joined,
               1, MPI_DOUBLE,
               MPI_MAX,
               rootnode, comm);

   MPI_Reduce( &min_complex_fftw_operand_mag, 
               &min_complex_fftw_operand_mag_joined,
               1, MPI_DOUBLE,
               MPI_MIN,
               rootnode, comm);

   MPI_Reduce( &max_complex_fftw_operand_re, 
               &max_complex_fftw_operand_re_joined,
               1, MPI_DOUBLE,
               MPI_MAX,
               rootnode, comm);

   MPI_Reduce( &min_complex_fftw_operand_re, 
               &min_complex_fftw_operand_re_joined,
               1, MPI_DOUBLE,
               MPI_MIN,
               rootnode, comm);

   MPI_Reduce( &max_complex_fftw_operand_im, 
               &max_complex_fftw_operand_im_joined,
               1, MPI_DOUBLE,
               MPI_MAX,
               rootnode, comm);

   MPI_Reduce( &min_complex_fftw_operand_im, 
               &min_complex_fftw_operand_im_joined,
               1, MPI_DOUBLE,
               MPI_MIN,
               rootnode, comm);

   if ( mynode == rootnode )
   {
//      double min_complex_fftw_operand_mag_joined;
//      double max_complex_fftw_operand_mag_joined;
//      double min_complex_fftw_operand_im_joined;
//      double max_complex_fftw_operand_im_joined;
//      double min_complex_fftw_operand_re_joined;
//      double max_complex_fftw_operand_re_joined;
// 
//      min_complex_fftw_operand_mag_joined
//         =  complex_fftw_operand_mag_joined[0];
//      max_complex_fftw_operand_mag_joined
//         = min_complex_fftw_operand_mag_joined;
// 
//      min_complex_fftw_operand_re_joined
//         =  complex_fftw_operand_re_joined[0];
//      max_complex_fftw_operand_re_joined
//         = min_complex_fftw_operand_re_joined;
// 
//      min_complex_fftw_operand_im_joined
//         =  complex_fftw_operand_im_joined[0];
//      max_complex_fftw_operand_im_joined
//         = min_complex_fftw_operand_im_joined;
// 
//      for ( ptrdiff_t i=0; i < Nx * Ny; i++)
//      {
//         if ( complex_fftw_operand_mag_joined[i]
//               <
//               min_complex_fftw_operand_mag_joined )
//            min_complex_fftw_operand_mag_joined
//               = complex_fftw_operand_mag_joined[i];
// 
//         if ( complex_fftw_operand_mag_joined[i]
//               >
//               max_complex_fftw_operand_mag_joined )
//            max_complex_fftw_operand_mag_joined
//               = complex_fftw_operand_mag_joined[i];
// 
//         if ( complex_fftw_operand_re_joined[i]
//               <
//               min_complex_fftw_operand_re_joined )
//            min_complex_fftw_operand_re_joined
//               = complex_fftw_operand_re_joined[i];
// 
//         if ( complex_fftw_operand_re_joined[i]
//               >
//               max_complex_fftw_operand_re_joined )
//            max_complex_fftw_operand_re_joined
//               = complex_fftw_operand_re_joined[i];
// 
//         if ( complex_fftw_operand_im_joined[i]
//               <
//               min_complex_fftw_operand_im_joined )
//            min_complex_fftw_operand_im_joined
//               = complex_fftw_operand_im_joined[i];
// 
//         if ( complex_fftw_operand_im_joined[i]
//               >
//               max_complex_fftw_operand_im_joined )
//            max_complex_fftw_operand_im_joined
//               = complex_fftw_operand_im_joined[i];
// 
//      }
      cout << "Writing realspace image to tiff" << endl;
      cout << " " << step << ". max_complex_fftw_operand_mag_joined : "
         << setprecision(20)
         << max_complex_fftw_operand_mag_joined << endl;
      cout << " " << step << ". min_complex_fftw_operand_mag_joined : "
         << setprecision(20)
         << min_complex_fftw_operand_mag_joined << endl;
      cout << " " << step << ". max_complex_fftw_operand_re_joined : "
         << setprecision(20)
         << max_complex_fftw_operand_re_joined << endl;
      cout << " " << step << ". min_complex_fftw_operand_re_joined : "
         << setprecision(20)
         << min_complex_fftw_operand_re_joined << endl;
      cout << " " << step << ". max_complex_fftw_operand_im_joined : "
         << setprecision(20)
         << max_complex_fftw_operand_im_joined << endl;
      cout << " " << step << ". min_complex_fftw_operand_im_joined : "
         << setprecision(20)
         << min_complex_fftw_operand_im_joined << endl;
 
      imageWriter.write_tif_grayscale16(
            complex_fftw_operand_mag_joined,
            min_complex_fftw_operand_mag_joined,
            max_complex_fftw_operand_mag_joined,
            Nx, Ny,
            resolutionUnit,
            xResolution, yResolution,
            outFileName_prefix + "_" + to_string(step) 
               +"_complex_fftw_operand_mag_output");
 
      imageWriter.write_tif_grayscale16_logscale(
            complex_fftw_operand_mag_joined,
            min_complex_fftw_operand_mag_joined,
            max_complex_fftw_operand_mag_joined,
            Nx, Ny,
            resolutionUnit,
            xResolution, yResolution,
            outFileName_prefix + "_" + to_string(step) 
               +  "_complex_fftw_operand_mag_output"
               + "_logscale" );
 
      imageWriter.write_tif_grayscale16(
            complex_fftw_operand_re_joined,
            min_complex_fftw_operand_re_joined,
            max_complex_fftw_operand_re_joined,
            Nx, Ny,
            resolutionUnit,
            xResolution, yResolution,
            outFileName_prefix + "_" + to_string(step) 
               + "_complex_fftw_operand_re"
            );
      imageWriter.write_tif_grayscale16_logscale(
            complex_fftw_operand_re_joined,
            min_complex_fftw_operand_re_joined,
            max_complex_fftw_operand_re_joined,
            Nx, Ny,
            resolutionUnit,
            xResolution, yResolution,
            outFileName_prefix + "_" + to_string(step) 
               + "_complex_fftw_operand_re"
               + "_logscale"
            );
 
      imageWriter.write_tif_grayscale16_logscale(
            complex_fftw_operand_im_joined,
            min_complex_fftw_operand_im_joined,
            max_complex_fftw_operand_im_joined,
            Nx, Ny,
            resolutionUnit,
            xResolution, yResolution,
            outFileName_prefix + "_" + to_string(step) 
               + "_complex_fftw_operand_im"
            );
      imageWriter.write_tif_grayscale16_logscale(
            complex_fftw_operand_im_joined,
            min_complex_fftw_operand_im_joined,
            max_complex_fftw_operand_im_joined,
            Nx, Ny,
            resolutionUnit,
            xResolution, yResolution,
            outFileName_prefix + "_" + to_string(step) 
               + "_complex_fftw_operand_im"
            +   "_logscale"
            );
      delete[] complex_fftw_operand_mag_joined;
      delete[] complex_fftw_operand_re_joined;
      delete[] complex_fftw_operand_im_joined;
   }
 
   delete[] complex_fftw_operand_mag;
   delete[] complex_fftw_operand_re;
   delete[] complex_fftw_operand_im;

   return;
} // debug_output_complex_fftw_operand

//void TEM_NS::output_stem_image(
//      const double* const stem_image,
//      const ptrdiff_t& Nx,
//      const ptrdiff_t& Ny,
//      const string& outFileName_prefix
//      )
//{
//
//   writeImage imageWriter;
//
//   double min_stem_image;
//   double max_stem_image;
//   
//   // Determine the max and min intensities
//   min_stem_image = stem_image[0];
//   max_stem_image = min_stem_image;
//   for ( ptrdiff_t i = 0; i < Nx * Ny; i++)
//   {  
//      if( stem_image[i] > max_stem_image ) max_stem_image = stem_image[i]; 
//      if( stem_image[i] < min_stem_image ) min_stem_image = stem_image[i];
//   }
// 
//   cout << " max_stem_image : "
//      << setprecision(20) << max_stem_image << endl;
//   cout <<" min_stem_image : "
//      << setprecision(20) << min_stem_image << endl;
//
//   imageWriter.write_tif_grayscale16(
//         stem_image,
//         min_stem_image,
//         max_stem_image,
//         Nx, Ny,
//         outFileName_prefix  
//            +"_stem_image_output");
// 
//   imageWriter.write_tif_grayscale16_logscale(
//         stem_image,
//         min_stem_image,
//         max_stem_image,
//         Nx, Ny,
//         outFileName_prefix 
//            +  "_stem_image_output"
//            + "_logscale" );
// 
//
//   return;
//} // output_stem_image


//double TEM_NS::scherzer_defocus_bfctem_uncorrected(
//      const double& lambda,
//      const double& Cs3
//      )
//{
//   // Kirkland (2016) Table 1
//   // defocus = (1.5 C_{S3} \lambda) ^{1/2}
//   return sqrt(1.5 * Cs3 * lambda);
//}
//
//double TEM_NS::scherzer_alphamax_bfctem_uncorrected(
//      const double& lambda,
//      const double& Cs3
//      )
//{
//   // Kirkland (2016) Table 1
//   // alphamax = (6 \lambda / C_{S3}) ^{1/4}
//   return pow( 6.0 * lambda / Cs3, 0.25);
//}
//
//double TEM_NS::scherzer_approxresolution_bfctem_uncorrected(
//      const double& lambda,
//      const double& Cs3
//      )
//{
//   // Kirkland (2016) Table 1
//   // dmin = 0.67 (C_{S3} \lambda^{3}) ^{1/4}
//   return 0.67 * pow( Cs3 * pow( lambda, 3.0), 0.25);
//}
//
//double TEM_NS::scherzer_defocus_adfstem_uncorrected(
//      const double& lambda,
//      const double& Cs3
//      )
//{
//   // Kirkland (2016) Table 1
//   // defocus = 0.87 (C_{S3} \lambda)^{1/2}
//   return 0.87 * sqrt( Cs3 * lambda );
//}
//
//double TEM_NS::scherzer_alphamax_adfstem_uncorrected(
//      const double& lambda,
//      const double& Cs3
//      )
//{
//   // Kirkland (2016) Table 1
//   // alphamax = 1.34 ( \lambda / C_{S3}) ^{1/4}
//   return 1.34 * pow( lambda / Cs3, 0.25);
//}
//
//double TEM_NS::scherzer_approxresolution_adfstem_uncorrected(
//      const double& lambda,
//      const double& Cs3
//      )
//{
//   // Kirkland (2016) Table 1
//   // dmin = 0.67 (C_{S3} \lambda^{3}) ^{1/4}
//   return 0.43 * pow( Cs3 * pow( lambda, 3.0), 0.25);
//}
//
//double TEM_NS::scherzer_Cs3_bfctem_correctedtoCs5(
//      const double& lambda,
//      const double& Cs5
//      )
//{
//   // Kirkland (2016) Table 2
//   // C_{S3} = -3.2 ( \lambda C_{S5}^{2} )^{1/3}
//   return -3.2 * pow( lambda * pow(Cs5, 2), 1.0/3.0);
//}
//
//double TEM_NS::scherzer_defocus_bfctem_correctedtoCs5(
//      const double& lambda,
//      const double& Cs5
//      )
//{
//   // Kirkland (2016) Table 2
//   // defocus =  -2 ( \lambda^{2} C_{S5} )^{1/3}
//   return -2.0 * pow( pow( lambda, 2) * Cs5, 1.0/3.0);
//}
//
//double TEM_NS::scherzer_alphamax_bfctem_correctedtoCs5(
//      const double& lambda,
//      const double& Cs5
//      )
//{
//   // Kirkland (2016) Table 2
//   // alphamax = \frac{7}{4} ( \lambda / C_{S5} )^{1/6}
//   return (7.0/4.0) * pow( lambda / Cs5, 1.0/6.0);
//}
//
//double TEM_NS::scherzer_approxresolution_bfctem_correctedtoCs5(
//      const double& lambda,
//      const double& Cs5
//      )
//{
//   // Kirkland (2016) Table 2
//   // dmin = \frac{4}{7} ( C_{S5} \lambda^{5} )^{1/6}
//   return (4.0/7.0) * pow( Cs5 * pow( lambda, 5), 1.0/6.0);
//}
//
//double TEM_NS::scherzer_Cs3_adfstem_correctedtoCs5(
//      const double& lambda,
//      const double& Cs5
//      )
//{
//   // Kirkland (2016) Table 2
//   // C_{S3} = -2.289 ( \lambda C_{S5}^{2} )^{1/3}
//   return -2.289 * pow( lambda * pow(Cs5, 2), 1.0/3.0);
//}
//
//double TEM_NS::scherzer_defocus_adfstem_correctedtoCs5(
//      const double& lambda,
//      const double& Cs5
//      )
//{
//   // defocus =  -0.983 ( \lambda^{2} C_{S5} )^{1/3}
//   return -0.983 * pow( pow( lambda, 2) * Cs5, 1.0/3.0);
//}
//
//double TEM_NS::scherzer_alphamax_adfstem_correctedtoCs5(
//      const double& lambda,
//      const double& Cs5
//      )
//{
//   // Kirkland (2016) Table 2
//   // alphamax = 1.513 ( \lambda / C_{S5} )^{1/6}
//   return 1.513 * pow( lambda / Cs5, 1.0/6.0);
//}
//
//double TEM_NS::scherzer_approxresolution_adfstem_correctedtoCs5(
//      const double& lambda,
//      const double& Cs5
//      )
//{
//   // Kirkland (2016) Table 2
//   // dmin = 0.403 ( C_{S5} \lambda^{5} )^{1/6}
//   return 0.403 * pow( Cs5 * pow( lambda, 5), 1.0/6.0);
//}

int TEM_NS::avg_stddev_of_change_in_array_members(
      const double* const myArray,
      const size_t& number_of_elements,
      double& avg,
      double& stddev
      )
{
   double tmp_prev = myArray[0];
   avg = 0.0;
   double avg_sqr = 0.0;

   for ( size_t i=1; i < number_of_elements; ++i)
   {
      avg += ( myArray[i]  - tmp_prev );
      if ( tmp_prev > myArray[i] ) // debug
      { // debug
         cout << "decrease in number of points binned" << endl; // debug
      } // debug

      avg_sqr += (myArray[i] - tmp_prev) * (myArray[i] - tmp_prev);

      tmp_prev = myArray[i];
   }
   cout << "sum of all array elements : " << avg << endl; // debug
   //cout << "sum of squares of all array elements : "  // debug
   //   << avg_sqr << endl; // debug
   avg = avg / number_of_elements;
   cout << "avg : " << avg << endl; // debug
   avg_sqr = avg_sqr / number_of_elements;
   //cout << "avg_sqr : " << avg_sqr << endl; // debug

   stddev = sqrt( avg_sqr - avg * avg );

   //cout << "Average change : "  // debug
   //   << avg // debug
   //   << " +/- " << deviation // debug
   //   << endl; // debug

   return EXIT_SUCCESS;
}

double TEM_NS::sum_detected_intensity( 
      const fftw_complex* const psi,
      const double* const kx, const double* const ky, 
      const double& k1, const double& k2,
      const size_t& Nx, const size_t& Ny
      )
{
   double detected_intensity = 0.0;
   // TODO: alter algorithm to reduce the number of evaluations of 
   //       kx^2 + ky^2 by assuming that the valid k-points are all 
   //       contiguous in the kx, ky domain
   //       - scan along kx 
   //       begin loop
   //       - when |k| becomes valid continue scanning in 
   //          that direction while accumulating the integrand
   //       - accumulate the integrand while scanning kx in the most
   //          recent direction until |k| is invalid
   //       - increment ky 
   //       - scan kx in two opposite directions, alternating at each 
   //          increment
   //       end loop


   // TODO: modify this algorithm to use a vector of indices sorted by 
   //  the magnitude of their corresponding kx, ky domain vector to 
   //  remove the need to iterate over the entire kx, ky domain Nx * Ny 
   //  times. 
   //  This algorithm has been implemented in the variation 
   //  calculation using a vector<indexed_vector_magnitude_sqr> 
   //  in the integrate_out_phi() function. 
   //  Modify the current function to work like integrate_out_phi() .
   //
   // Until a more efficient algorithm is implemented, just scan the
   //  entire kx, ky domain and accumulate the detected intensity when
   //  k1 < |kx^2 + ky^2| < k2
   double k_sqr;
   for (size_t i=0; i < Nx; ++i)
      for (size_t j=0; j < Ny; ++j)
      {
         k_sqr = kx[i] * kx[i] + ky[j] * ky[j];
         if ( k_sqr > k1 && k_sqr <= k2 )// accepting values in (k1, k2]
            detected_intensity += 
               psi[j + i*Ny][0] * psi[j + i*Ny][0]
               +
               psi[j + i*Ny][1] * psi[j + i*Ny][1];
      }
   return detected_intensity;
}


#endif
