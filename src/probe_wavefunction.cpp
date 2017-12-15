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
// File: probe_wavefunction.cpp
// Purpose:
//
// Simulate the scanning TEM of a thick sample as viewed perpendicular
//  to the third axis (q[2]).
// This file implements the multislice approximation, phase grating 
//  approximation, and the weak phase approximation. 

#ifndef PROBE_WAVEFUNCTION_CPP
#define PROBE_WAVEFUNCTION_CPP

#include "probe_wavefunction.hpp"



int TEM_NS::probe_wavefunction_uncorrected_unnormalized(
      const double& x_p,   // x component of \vec{x}_{p}
      const double& y_p,   // y component of \vec{x}_{p}
      const double* const kx, // x component of the reciprocal space domain
      const ptrdiff_t& Nx,
      const double* const ky, // y component of the reciprocal space domain
      const ptrdiff_t& Ny,
      const double& Cs3,       // spherical aberration
      const double& defocus,  // focal point distance from sample bottom
      const double& alpha_max_sqr, // objective aperture limiting angle
      const double& lambda,
      const double& lambda_sqr,
      fftw_complex* psi
      )
{
   // Postcondition:
   // - the wave function in reciprocal space is stored in psi
   //    which still needs an inverse Fourier Transform to reach realspace
   
   // The wave function in realspace calculated via an inverse Fourier 
   //  transform:
   // Kirkland (2010) eq 5.47
   // \psi_{p}(\vec{x}, \vec{x}_{p}) = A_{p} FT^{-1}[
   //       \exp[ -i \Chi( \vec{k}) + 2 \pi i \vec{k} \cdot \vec{x}_{p}]
   //       A(\vec{k})
   //    ]
   //
   // where A(\vec{k}) is the aperture function 
   //    (eq 5.28, =1 inside, else 0)
   // and A_{p} is a normalization constant chosen to yield
   // \int |\psi_{p}(\vec{x}, \vec{x}_{p}|^{2} d^{2}\vec{x} = 1

   // Kirkland (2016) eq 16 is identical to Kirkland (2010) eq 5.45, from 
   //  which Kirkland (2010) eq 5.47 is derived.
   // Kirkland (2016) eq 16
   // \psi_{p}(\vec{x},\vec{x}_{p}) = A_{p} \int^{k_{max}}_{0}
   //    \exp[ -i \Chi(\vec{k}) - 2 \pi i 
   //          \vec{k} \cdot ( \vec{x} - \vec{x}_{p} ) ] 
   //    d^{2}\vec{k}
   // 

   
   //TODO: Do you really need to calculate the normalization constant A_{p}?
   //      If so, then do you really need to calculate it for every p ?
   // If the x & y boundaries are periodic, then the integral of impinging 
   //  wavefunction amplitued over the entire sample should be invariant of
   //  position p.

   // TODO: evaluate A_p and multiply psi[] by it.
   //       - performed outside of this function

   double ksqr;
   double trig_operand; // value to be inserted into exp[]

   for( ptrdiff_t i=0; i < Nx; ++i)
      for( ptrdiff_t j=0; j < Ny; ++j)
      {
         //////////////////////////////////////////////////////////////////
         // if |\vec{k}|^{2} >= alpha_max_sqr 
         // then \psi = 0,
         // else
         //
         // Evaluate the operand of the exponential
         //
         // \exp[-i \Chi(\vec{k}) + 2 \pi i \vec{k} \cdot \vec{x}_{p}]
         // = \exp[-i (\Chi(ksqr) - 2 \pi (kx*x_p + ky*y_p)) ]
         // = cos( trig_operand ) - i sin( trig_operand )
         // where
         // trig_operand = \Chi(ksqr) - 2 \pi (kx * x_p + ky * y_p)
         //////////////////////////////////////////////////////////////////
         //ksqr = kx[j + i * Ny] * kx[j + i * Ny] 
         //         + ky[j + i * Ny] * ky[j + i * Ny];
         ksqr = kx[i] * kx[i] + ky[j] * ky[j];

         if ( lambda_sqr * ksqr < alpha_max_sqr )
         {
            // \Chi(\vec{k}), uncorrected
            trig_operand 
               = aberration_function_uncorrected(
                  ksqr,
                  Cs3,
                  defocus,
                  lambda,
                  lambda_sqr
                  )
                  + 2 * PI *( ((kx[i]) * x_p) + ((ky[j]) * y_p));
                  //- 2 * PI *( ((kx[i]) * x_p) + ((ky[j]) * y_p));

            psi[j + i * Ny][0] = cos( trig_operand );
            psi[j + i * Ny][1] = -1.0 * sin( trig_operand );
         }
         else
         {
            psi[j + i * Ny][0] = 0.0;
            psi[j + i * Ny][1] = 0.0;
         }
      }
   cout << "Evaluated probe ..." << endl;// debug

   return EXIT_SUCCESS;
}

int TEM_NS::probe_wavefunction_correctedtoCs5_unnormalized(
      const double& x_p,   // x component of \vec{x}_{p}
      const double& y_p,   // y component of \vec{x}_{p}
      const double* const kx, // x component of the reciprocal space domain
      const ptrdiff_t& Nx,
      const double* const ky, // y component of the reciprocal space domain
      const ptrdiff_t& Ny,
      const double& Cs3,       // spherical aberration
      const double& Cs5,
      const double& defocus,  // focal point distance from sample bottom
      //const double& defocus_spread,  // \Delta_{0}
      const double& alpha_max_sqr, // objective aperture limiting angle
      const double& lambda,
      const double& lambda_sqr,
      fftw_complex* psi
      )
{
   // Postcondition:
   // - the wave function in reciprocal space is stored in psi
   //    which still needs an inverse Fourier Transform to reach realspace
   
   // The wave function in realspace calculated via an inverse Fourier 
   //  transform:
   // Kirkland (2010) eq 5.47
   // \psi_{p}(\vec{x}, \vec{x}_{p}) = A_{p} FT^{-1}[
   //       \exp[ -i \Chi( \vec{k}) + 2 \pi i \vec{k} \cdot \vec{x}_{p}]
   //       A(\vec{k})
   //    ]
   //
   // where A(\vec{k}) is the aperture function (eq 5.28, =1 inside, else 0)
   // and A_{p} is a normalization constant chosen to yield
   // \int |\psi_{p}(\vec{x}, \vec{x}_{p}|^{2} d^{2}\vec{x} = 1
   // 
   // For the aberration corrected probe \Chi differs from the 
   //  uncorrected probe, making use of Cs5.
   //
   // Kirkland (2016) eq 16 is identical to Kirkland (2010) eq 5.45, from 
   //  which Kirkland (2010) eq 5.47 is derived.

   
   //TODO: Do you really need to calculate the normalization constant A_{p}?
   //      If so, then do you really need to calculate it for every p ?
   // If the x & y boundaries are periodic, then the integral of impinging 
   //  wavefunction amplitued over the entire sample should be invariant of
   //  position p.

   // TODO: evaluate A_p and multiply psi[] by it.

   double ksqr;
   double trig_operand; // value to be inserted into exp[]

   for( ptrdiff_t i=0; i < Nx; ++i)
      for( ptrdiff_t j=0; j < Ny; ++j)
      {
         //////////////////////////////////////////////////////////////////
         // if |\vec{k}|^{2} >= alpha_max_sqr 
         // then \psi = 0,
         // else
         //
         // Evaluate the operand of the exponential
         //
         // \exp[-i \Chi(\vec{k}) + 2 \pi i \vec{k} \cdot \vec{x}_{p}]
         // = \exp[-i (\Chi(ksqr) - 2 \pi (kx*x_p + ky*y_p)) ]
         // = cos( trig_operand ) - i sin( trig_operand )
         // where
         // trig_operand = \Chi(ksqr) - 2 \pi (kx * x_p + ky * y_p)
         //////////////////////////////////////////////////////////////////
         //ksqr = kx[j + i * Ny] * kx[j + i * Ny] 
         //         + ky[j + i * Ny] * ky[j + i * Ny];
         ksqr = kx[i] * kx[i] + ky[j] * ky[j];

         if ( lambda_sqr * ksqr < alpha_max_sqr )
         {
            // \Chi(\vec{k}), corrected
            trig_operand 
               = aberration_function_correctedtoCs5(
                  ksqr,
                  Cs3,
                  Cs5,
                  defocus,
                  lambda,
                  lambda_sqr
                  )
                  - 2 * PI *( ((kx[i]) * x_p) + ((ky[j]) * y_p));
                  //- 2 * PI * (kx[j + i * Ny] * x_p
                  //            + ky[j + i * Ny] * y_p);

            psi[j + i * Ny][0] = cos( trig_operand );
            psi[j + i * Ny][1] = -1.0 * sin( trig_operand );
            // NOTE: probe norm A_p is evaluated and applied outside
            //       of this function.
         }
         else
         {
            psi[j + i * Ny][0] = 0.0;
            psi[j + i * Ny][1] = 0.0;
         }
      }

   return EXIT_SUCCESS;
}

int TEM_NS::probe_wavefunction_norm(
      double& Ap,
      const ptrdiff_t& Nx,
      const ptrdiff_t& Ny,
      fftw_complex* psi,
      //const int& mynode,
      //const int& rootnode,
      //const int& totalnodes,
      MPI_Comm comm
      )
{
   // A_{p} = 1/sqrt(\int |\psi_{p}(\vec{x}, \vec{x}_{p})|^{2} d^{2}\vec{x})
   //
   // The "Plancherel" theorem states that
   // \int_{-Inf}^{Inf} |f(x)|^{2} dx = \int_{-Inf}^{Inf} |FT[f(x)]|^{2} dx
   //
   // so
   // A_{p} 
   //  = 1/sqrt(
   //       \int | FT[ \psi_{p}(\vec{x}, \vec{x}_{p})] |^{2} d^{2}\vec{x}
   //    )
   // and the norm may thus be calculated without having to transform
   // the probe wavefunction from reciprocal space to real space.
   
   //cout << " probe_wavefunction_norm( " << Ap 
   //   << ", " << Nx << ", " << Ny
   //   << ", " << psi << ", " << comm << ")" << endl;

   const ptrdiff_t NxNy = Nx * Ny;
   Ap = 0.0;
   double Ap_local; Ap_local = 0.0;

   for ( ptrdiff_t i=0; i< NxNy; ++i )
   {
      Ap_local += psi[i][0] * psi[i][0] + psi[i][1] * psi[i][1];
   }

   MPI_Allreduce( &Ap_local, &Ap, 1, MPI_DOUBLE, MPI_SUM, comm);

   //cout << "Ap_local : " << setprecision(5) << Ap_local << endl; // debug

   Ap = 1.0 / sqrt( Ap );

   return EXIT_SUCCESS;
}

int TEM_NS::probe_wavefunction_norm(
      double& Ap,
      const ptrdiff_t& Nx,
      const ptrdiff_t& Ny,
      const double* const psi_re,
      const double* const psi_im,
      //const int& mynode,
      //const int& rootnode,
      //const int& totalnodes,
      MPI_Comm comm
      )
{
   const ptrdiff_t NxNy = Nx * Ny;
   Ap = 0.0;
   double Ap_local; Ap_local = 0.0;

   for ( ptrdiff_t i=0; i< NxNy; ++i )
   {
      Ap_local += psi_re[i] * psi_re[i] + psi_im[i] * psi_im[i];
   }

   MPI_Allreduce( &Ap_local, &Ap, 1, MPI_DOUBLE, MPI_SUM, comm);

   //cout << "Ap_local : " << setprecision(5) << Ap_local << endl; // debug

   Ap = 1.0 / sqrt( Ap );

   return EXIT_SUCCESS;
}


//int TEM_NS::parallel_probe_wavefunction_uncorrected_unnormalized(
//      const double& x_p,   // x component of \vec{x}_{p}
//      const double& y_p,   // y component of \vec{x}_{p}
//      const double* const kx, // x component of the reciprocal space domain
//      const ptrdiff_t& Nx,
//      const double* const ky, // y component of the reciprocal space domain
//      const ptrdiff_t& Ny,
//      const double& alpha_max_sqr, // objective aperture limiting angle
//      const double& lambda_sqr,
//      fftw_complex* psi
//      )
//{
//
//   double ksqr;
//   double trig_operand; // value to be inserted into exp[]
//
//   for( ptrdiff_t i=0; i < Nx; ++i)
//      for( ptrdiff_t j=0; j < Ny; ++j)
//      {
//         ///////////////////////////////////////////////////////////////////
//         // if |\vec{k}|^{2} >= alpha_max_sqr 
//         // then \psi = 0,
//         // else
//         //
//         // Evaluate the operand of the exponential
//         //
//         // \exp[-i \Chi(\vec{k}) + 2 \pi i \vec{k} \cdot \vec{x}_{p}]
//         // = \exp[-i (\Chi(ksqr) - 2 \pi (kx*x_p + ky*y_p)) ]
//         // = cos( trig_operand ) - i sin( trig_operand )
//         // where
//         // trig_operand = \Chi(ksqr) - 2 \pi (kx * x_p + ky * y_p)
//         ///////////////////////////////////////////////////////////////////
//         //ksqr = kx[j + i * Ny] * kx[j + i * Ny] 
//         //         + ky[j + i * Ny] * ky[j + i * Ny];
//         ksqr = kx[i] * kx[i] + ky[j] * ky[j];
//
//         if ( lambda_sqr * ksqr < alpha_max_sqr )
//         {
//            // \Chi(\vec{k}), uncorrected
//            trig_operand 
//               = 
//               // TODO : Is it correct to elliminate the obj lens for 
//               //          nanodiffraction (FTEM) ?
//               //aberration_function_uncorrected(
//               //   ksqr,
//               //   Cs3,
//               //   defocus,
//               //   lambda,
//               //   lambda_sqr
//               //   )
//                  - 2 * PI *( ((kx[i]) * x_p) + ((ky[j]) * y_p));
//                  //- 2 * PI * (kx[j + i * Ny] * x_p
//                  //            + ky[j + i * Ny] * y_p);
//
//            psi[j + i * Ny][0] = cos( trig_operand );
//            psi[j + i * Ny][1] = -1.0 * sin( trig_operand );
//            // TODO: evaluate A_p and multiply psi[] by it.
//         }
//         else
//         {
//            psi[j + i * Ny][0] = 0.0;
//            psi[j + i * Ny][1] = 0.0;
//         }
//      }
//   return EXIT_SUCCESS;
//}

#endif
