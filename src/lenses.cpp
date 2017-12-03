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

#ifndef LENSES_CPP
#define LENSES_CPP

#include "lenses.hpp" 

void TEM_NS::objective_lens_transfer_function_coher(
          const double& ksqr, 
          const double& Cs3, 
          const double& defocus,
          const double& alpha_max_sqr,
          const double& lambda,
          const double& lambdasqr,
          double& H_0_re,
          double& H_0_im
          )
{
   // Kirkland (2010) eq. 5.27
   // H_{0}(\vec{k}) = \exp[ -i \chi(\vec{k})] A(\vec{k})
   // \chi(|k|) 
   //    = \pi \lambda k^{2} ( 0.5 C_{s} \lambda^{2} |k|^{2} - \Delta f)
   // where \Delta f is defocus, C_{s} is the coefficient of spherical 
   // aberration, and A(\vec{k}) is the aperture function:
   // A(\vec{k}) = 1 if \lambda k = \alpha < \alpha_{max}, 0 otherwise, 
   // where \alpha_{max} is the maximum semiangle allowed by the objective
   // aperture.
   //
   // NOTE: this version using eq 5.27 from Kirkland (2010) might cause
   //    a reflection about the origin to appear in the image.


   if ( lambdasqr * ksqr < alpha_max_sqr ) 
   {
      double chi = aberration_function_uncorrected(
            ksqr, Cs3, defocus, lambda, lambdasqr
            );
         //PI * lambda * ksqr * ( 0.5 * Cs * lambdasqr * ksqr - defocus);
      H_0_re = cos( chi );
      H_0_im = -1 * sin( chi );
   }
   else
   {
      H_0_re = 0.0;
      H_0_im = 0.0;
   }
   
   return ;
}

void TEM_NS::objective_lens_transfer_function_wp(
          const double& ksqr, 
          const double& Cs3, 
          const double& defocus,
          const double& alpha_max_sqr,
          const double& lambda,
          const double& lambdasqr,
          double& H_0_re,
          double& H_0_im
          )
{
   // Kirkland (2010) eq. 5.27
   // H_{0}(\vec{k}) = \exp[ -i \chi(\vec{k})] A(\vec{k})
   // \chi(|k|) 
   //    = \pi \lambda k^{2} ( 0.5 C_{s} \lambda^{2} |k|^{2} - \Delta f)
   // where \Delta f is defocus, C_{s} is the coefficient of spherical 
   // aberration, and A(\vec{k}) is the aperture function:
   // A(\vec{k}) = 1 if \lambda k = \alpha < \alpha_{max}, 0 otherwise, 
   // where \alpha_{max} is the maximum semiangle allowed by the objective
   // aperture.

   // This version of the function follows
   // Kirkland (2016) eq 12 which uses a linear approximation instead :
   // H_{WP}(\vec{k}) = FT[h_{WP}(\vec{x})]  = sin[ \Chi(\vec{k}) ]
   //
   // This approximation eliminates the imaginary part and thus the 
   //  centrosymmetry of the final image which was inappropriate anyway.

   if ( lambdasqr * ksqr < alpha_max_sqr ) 
   {
      //double chi = aberration_function_correctedtoCs5(
      //      ksqr, Cs3, Cs5, defocus, lambda, lambdasqr
      //      );
      double chi = aberration_function_uncorrected(
            ksqr, Cs3, defocus, lambda, lambdasqr
            );
      H_0_re = sin( chi );
      H_0_im = 0.0;
   }
   else
   {
      H_0_re = 0.0;
      H_0_im = 0.0;
   }
   
   return ;
}

void TEM_NS::objective_lens_transfer_function_wp_correctedtoCs5(
          const double& ksqr, 
          const double& Cs3, 
          const double& Cs5, 
          const double& defocus,
          const double& alpha_max_sqr,
          const double& lambda,
          const double& lambdasqr,
          double& H_0_re,
          double& H_0_im
          )
{
   // Kirkland (2010) eq. 5.27
   // H_{0}(\vec{k}) = \exp[ -i \chi(\vec{k})] A(\vec{k})
   // \chi(|k|) 
   //    = \pi \lambda k^{2} ( 0.5 C_{s} \lambda^{2} |k|^{2} - \Delta f)
   // where \Delta f is defocus, C_{s} is the coefficient of spherical 
   // aberration, and A(\vec{k}) is the aperture function:
   // A(\vec{k}) = 1 if \lambda k = \alpha < \alpha_{max}, 0 otherwise, 
   // where \alpha_{max} is the maximum semiangle allowed by the objective
   // aperture.
   //
   // Kirkland (2016) eq 12 uses a linear approximation instead :
   // H_{WP}(\vec{k}) = FT[h_{WP}(\vec{x})]  = sin[ \Chi(\vec{k}) ]
   //
   // This approximation eliminates the imaginary part and thus the 
   //  centrosymmetry of the final image which was inappropriate anyway.


   if ( lambdasqr * ksqr < alpha_max_sqr ) 
   {
      //double chi = aberration_function_correctedtoCs5(
      //      ksqr, Cs3, Cs5, defocus, lambda, lambdasqr
      //      );
      double chi = aberration_function_correctedtoCs5(
            ksqr, Cs3, Cs5, defocus, lambda, lambdasqr
            );
      H_0_re = sin( chi );
      H_0_im = 0.0;
   }
   else
   {
      H_0_re = 0.0;
      H_0_im = 0.0;
   }
   
   return ;
}

void TEM_NS::objective_lens_transfer_function_wp_correctedtoCs5_spread(
          const double& ksqr, 
          const double& Cs3, 
          const double& Cs5, 
          const double& defocus,
          const double& defocus_spread,   // \Delta_{0}
          const double& condenser_illumination_angle, // \lambda k_{s}
          const double& alpha_max_sqr,
          const double& lambda,
          const double& lambdasqr,
          double& H_0_re,
          double& H_0_im
          )
{
   // Kirkland (2016) eq. 14
   // H_{WP}(k) = \frac{1}{sqrt(1 + \epsilon k^{2})}
   //       \sin \left{
   //       \left[
   //          \frac{1}{3} C_{S5} 
   //          (1 - 2 \epsilon k^{2}) \lambda^{4} k^{4} 
   //          + 0.5 C_{S3} (1 - \epsilon k^{2}) \lambda^{2} k^{2} 
   //          - \Delta f
   //       \right]
   //       \right}
   //       \exp 
   //       \left( 
   //          - \left{
   //             [ 
   //                \pi \lambda k_{s} k ( C_{S5} \lambda^{4} k^{4}
   //                + C_{S3} \lambda^{2} k^{2} - \Delta f 
   //             ]^{2}
   //             + 0.25 ( \pi \lambda \Delta_{0} k^{2} )^{2}
   //          \right}
   //          / (1 + \epsilon k^{2})
   //       \right)
   // where \Delta f is defocus, C_{s} is the coefficient of spherical 
   // aberration, and A(\vec{k}) is the aperture function:
   // A(\vec{k}) = 1 if \lambda k = \alpha < \alpha_{max}, 0 otherwise, 
   // where \alpha_{max} is the maximum semiangle allowed by the objective
   // aperture.

   double epsilon = pow( PI * condenser_illumination_angle 
                           * defocus_spread, 2);

   const double onethird = 1.0/3.0;
   double epsilonksqr = epsilon * ksqr;
   double oneplusepsksqr = 1 + epsilonksqr;
   if ( lambdasqr * ksqr < alpha_max_sqr ) 
   {
      H_0_re = (1/(sqrt( oneplusepsksqr )))
            * sin( (PI * lambda * ksqr / oneplusepsksqr )
                     * (
                     onethird * Cs5 * (1 - 2 * epsilonksqr) 
                     * pow( lambdasqr * ksqr, 2)
                     + 0.5 * Cs3 * (1 - epsilonksqr) * lambdasqr * ksqr
                     - defocus
                     )
                 )
                  * exp(
                     -(
                        ksqr 
                        * pow(
                           PI * condenser_illumination_angle 
                           * ( Cs5 * pow( lambdasqr * ksqr, 2) 
                              + Cs3 * lambdasqr * ksqr - defocus
                             )
                         ,2)  
                        + 0.25 * pow( PI * lambda * defocus_spread * ksqr,2)
                     ) / oneplusepsksqr
                     );
   }
   else H_0_re = 0.0;

   H_0_im = 0.0;
   
   return ;
}

double TEM_NS::aberration_function_uncorrected( 
      const double& ksqr,
      const double& Cs3,
      const double& defocus,
      const double& lambda,
      const double& lambdasqr
      )
{
   // Kirkland (2010) eq 5.27 :
   // \Chi( k ) = \pi \lambda k^{2}(0.5 C_{s3} \lambda^{2} k^{2} - \delta f)
   return PI * lambda * ksqr
      * ( 0.5 * Cs3 * lambdasqr * ksqr - defocus);
}

double TEM_NS::aberration_function_correctedtoCs5( 
      const double& ksqr,
      const double& Cs3,
      const double& Cs5,
      const double& defocus,
      const double& lambda,
      const double& lambdasqr
      )
{
   // \Chi = \frac{2 \pi}{\lambda}
   //       \left( 
   //             \frac{1}{2} C_{1} \alpha^{2} 
   //             + \frac{1}{4} C_{S3} \alpha^{4} 
   //             + \frac{1}{6} C_{S5} \alpha^{6}  + ...
   //       \right)
   // where C_{1} = - \Delta f and \alpha = \lambda k , resulting in
   // \Chi = \frac{2 \pi}{\lambda}
   //       \left( 
   //             \frac{1}{2} C_{1} \alpha^{2} 
   //             + \frac{1}{4} C_{S3} \alpha^{4} 
   //             + \frac{1}{6} C_{S5} \alpha^{6}  + ...
   //       \right)
   double ksqrlambdasqr = ksqr * lambdasqr; 
   return PI * lambda * ksqr
      * ( - defocus + 0.5 * Cs3 * ksqrlambdasqr
         + (1.0/3.0) * Cs5 * pow( ksqrlambdasqr, 2)
         );
}

#endif
