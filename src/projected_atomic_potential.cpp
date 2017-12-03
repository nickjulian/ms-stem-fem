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
// File: projected_atomic_potential.cpp
// Purpose:
//

// For the test, use a fake projected_atomic_potential function to retrieve
//  LUT parameters.

#ifndef PROJECTED_ATOMIC_POTENTIAL_CPP
#define PROJECTED_ATOMIC_POTENTIAL_CPP

// TODO: remove these dependencies used for debugging
#include <iomanip>
using std::setw;
using std::setprecision;
//

//#include <cstdlib>

//#include <limits>   // numeric_limits::epsilon()

#include "projected_atomic_potential.hpp"
//#include "bessel_K0.hpp"

// Excerpt from Kirkland (2010) Appendix C
// "The three-dimensional atomic potential $V( x, y, z )$ is (C.19)
// $V( x, y, z) = V( \vec{r} ) = 2 \pi a_{0} e 
//    \int f_{e}(q) \exp(-2  \pi i \vec{q} \cdot \vec{r} ) d^{3}
//    = 2 \pi^{2} a_{0} e 
//       \sum_{i} \frac{a_{i}}{r} \exp(-2 \pi r \sqrt(b_{i})
//       + 2 \pi^{5/2} a_{0} e \sum_{i} c_{i} d_{i}^{-3/2}
//          \exp(-pi^{2} r^{2} /d_{i} )$
//  with $r^{2} = x^{2} + y^{2} + z^{2}$
//  and the projected atomic potential is    (C.20)
//  $V_{z}( x, y ) = \int_{-\infty}^{+\infty} V( x, y, z )dz
//    = 4 \pi^{2} a_{0} e 
//       \sum_{i} a_{i} K_{0}(2 \pi r \sqrt(b_{i}) )
//       + 2 \pi^{2} a_{0} e 
//          \sum_{i} \frac{c_{i}}{d_{i}} \exp( -\pi^{2} r^{2} / d_{i} )$
// where $K_{0}(x)$ is the modified Bessel function. The integral for the
// first summation was helped by expression 3.387.6 of Gradshteyn and
// Ryzhik [ref 128 in Kirkland]. Abramowitz and Stegun [ref 1 in Kirkland]
// (Sect. 9.8) give a convenient numerical expression for evaluating 
// K_{0}(x). The right hand side of the expression for V_{z}( x, y ) (C.20)
// has units of the electron charge $e$. By combining the Rydberg constant
// $Ry = 0.5 e^{2}/a_0}$ ($Ry/e = 13.6$ volts) and the Bohr radius 
// $a_{0} = 0.529 \AA$ the electron charge can be written in an 
// unconventional set of units as $e = 14.4$ Volt-Angstroms which is more 
// convenient for evaluating $V_{z}( x, y )$."
//
// Use of the projected atomic potential is referred to as the weak phase
// object approximation (Cowley and Iijima 1972, electron microscope image 
//  contrast for thin crystals. Z. Naturforsch 27a:445-451).

// Excerpts from Abramowitz and Stegun (1965)
// "
// 9.6.21
// $K_{0}(x) = \int_{0}^{\infty} \cos(x \sinh t)dt 
//    = \int_{0}^{\infty} \frac{\cos( x t )}{\sqrt(t^{2} + 1)} dt (x > 0)$
// "
// "
// 9.8.5 $ 0 < x <= 2 $
// $K_{0}(x) = - \ln (x/2) I_{0}(x) - 0.57721 566
//    + 0.42278 420 (x/2)^{2} + 0.23069 756 (x/2)^{4} 
//    + 0.03488 590 (x/2)^{6} + 0.00262 698 (x/2)^{8}
//    + 0.00010 750 (x/2)^{10} + 0.00000 740 (x/2)^{12} 
//    + \epsilon$
//    $|\epsilon| < 1E-8  $
// "
// "
// 9.8.6 $ 2 <= x < \infty $
// $ x^{1/2} e^{x} K_{0}(x) 
//    = 1.25331 414 - 0.07832 358 (2/x) 
//    + 0.02189 568 (2/x)^{2} + 0.01062 446 (2/x)^{3} 
//    + 0.00587 872 (2/x)^{4} - 0.00251 540 (2/x)^{5}
//    + 0.00053 208 (2/x)^{6} 
//    + \epsilon$
//    $|\epsilon| < 1.9E-7  $
// "
// "
// In equations 9.8.1 to 9.8.4, $t = x/3.75$.
// 9.8.1 $-3.75 <= x <= 3.75$
// $I_{0}(x) = 1 + 3.51562 29 t^{2} + 3.08994 24 t^{4} + 1.20674 92 t^{6}
//    + 0.26597 32 t^{8} + 0.03607 68 t^{10} + 0.00458 13 t^{12} 
//    + \epsilon $
//    $|\epsilon| < 1.6E-7  $
// 9.8.2 $3.75 <= x < \infty$
// $ x^{1/2} e^{-x} I_{0}(x) =
//    0.39894 228 + 0.01328 592 t^{-1}
//    + 0.00225 319 t^{-2} - 0.00157 565 t^{-3}
//    + 0.00916 281 t^{-4} - 0.02057 706 t^{-5}
//    + 0.02635 537 t^{-6} - 0.01647 633 t^{-7}
//    + 0.00392 377 t^{-8} 
//    + \epsilon $
//    $|\epsilon| < 1.9E-7  $
// "

//int TEM_NS::sigmaVz(
//         // projected atomic potential in reciprocal space
//      const std::vector<const scatterer*>& myScatterers,
//      const double& lambda, const double& gamma, const double& ab_inv,
//      const double& cutoff, 
//      const double* const kxdomain, 
//      const ptrdiff_t& Nx,
//      const double* const kydomain, 
//      const ptrdiff_t& Ny,
//      fftw_complex* pap
//      )
//{
//   // TODO: - create pap for each element simulated
//   //       - when sigmaVz is called, translate the existing paps to the 
//   //          appropriate atom positions
//   
//   // Accumulate sigmaVz_single_scatterer() over all of myScatterers over
//   //  all points in the reciprocal space domain.
//   // This version limits bandwidth frequency squared to cutoffsq, all
//   //  points having ||(kx,ky)||^{2} beyond cutoffsq will not be evaluated.
//   
//   //cout << "using sigmaVz from line 883" << endl;//debug
//   const double cutoffsq = pow(cutoff, 2);
//   // TODO: could the following be separated into x & y factors so that
//   //  only Nx + Ny operations are needed for each scatterer, rather than
//   //  Nx * Ny operations?
//   //  Only the calculation of sigmaVz is seperable. The Nx*Ny evaluations
//   //   of sigmaVz_single_scatterer may be reduced to Nx+Ny, but the outer
//   //   loop must still hit all Nx*Ny points of pap.
//   for ( ptrdiff_t i=0; i < Nx; i++) 
//      for ( ptrdiff_t j=0; j < Ny; j++) 
//      {
//         pap[(j+ i *Ny)][0] = 0.0;
//         pap[(j+ i *Ny)][1] = 0.0;
//         if ( pow(kxdomain[i], 2) + pow(kydomain[j],2) < cutoffsq ) 
//            for ( int k=0; k < myScatterers.size(); k++)
//                  sigmaVz_single_scatterer(
//                        lambda, gamma, ab_inv,
//                        kxdomain[i], kydomain[j],
//                        myScatterers[k]->q[0], myScatterers[k]->q[1],
//                        //myScatterers[k]->Z,  // atomic species
//                        myScatterers[k]->a, myScatterers[k]->b,
//                        myScatterers[k]->c, myScatterers[k]->d,
//                        pap[(j+ i *Ny)][0],   // re
//                        pap[(j+ i *Ny)][1]   // im
//                        );
//      }
//
//   return EXIT_SUCCESS;
//}

int TEM_NS::sigmaVz_single_scatterer(
      const double& lambda, const double& gamma, const double& ab_inv,
      const double& kx, const double& ky, // kz =0, Fourier projection thm
      const double& x, const double& y,   // z integrated out
      const double* const a, const double* const b,
      const double* const c, const double* const d,
      double& sigmaVz_re,
      double& sigmaVz_im
      )
{
   // \sigma V_{z} (k_{x}, k_{y}, 0) 
   //    = \lambda \frac{m}{m_{0}} \frac{1}{ab}
   //       \sum_{j} 
   //          f_{ej}( k_{x}, k_{y}, 0) 
   //          \exp{2 \pi i (k_{x} x_{j} + k_{y} y_{j})}
   //
   // \frac{m}{m_{0}} = \gamma 
   // \sigma = \frac{2 pi}{\lambda V} 
   //    \left( \frac{m_{0} c^{2} + e V}{2 m_{0} c^{2} + e V} \right)
   //    = \frac{2 pi m e \lambda}{ h^{2} }
   // \implies
   // m = \frac{h^{2}}{\lambda^{2} e V}
   //    \left( \frac{m_{0} c^{2} + e V}{2 m_{0} c^{2} + e V} \right)
   //    = \frac{ \sigma h^{2} }{ 2 pi e \lambda }
   // \implies 
   // gamma = \frac{m}{m_0}
   //    = \frac{ \sigma h^{2} }{ 2 pi e \lambda m_0 }

   //const double pi = 3.14159265369;

   double ksq; ksq = (kx * kx) + (ky * ky);

   double f_ej; f_ej = f_ej_single_scatterer( a, b, c, d, ksq);
   double trig_operand; trig_operand = 2 * PI * (kx * x + ky * y);
   sigmaVz_re += lambda * gamma * ab_inv * f_ej * cos( trig_operand ); 
   sigmaVz_im += lambda * gamma * ab_inv * f_ej * sin( trig_operand ); 

   return EXIT_SUCCESS;
}

double TEM_NS::f_ej_single_scatterer(
      const double* const a, const double* const b,
      const double* const c, const double* const d,
      const double& ksq
      )
{
   // f_{ej} = \sum_{i=1}^{3} \frac{a_{i}}{q^{2} + b_{i}} 
   //          + sum_{i=1}^{3} c_{i} e^{-d_{i} q^{2}}
   return 
      a[0]/(ksq + b[0])
      + a[1]/(ksq + b[1])
      + a[2]/(ksq + b[2])
      + c[0] * exp( - d[0] * ksq )
      + c[1] * exp( - d[1] * ksq )
      + c[2] * exp( - d[2] * ksq );
}
#endif
