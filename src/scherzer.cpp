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
// File: scherzer.cpp
// Purpose:

#ifndef SCHERZER_CPP
#define SCHERZER_CPP

#include "scherzer.hpp"


double TEM_NS::scherzer_defocus_bfctem_uncorrected(
      const double& lambda,
      const double& Cs3
      )
{
   // Kirkland (2016) Table 1
   // defocus = (1.5 C_{S3} \lambda) ^{1/2}
   return sqrt(1.5 * Cs3 * lambda);
}

double TEM_NS::scherzer_alphamax_bfctem_uncorrected(
      const double& lambda,
      const double& Cs3
      )
{
   // Kirkland (2016) Table 1
   // alphamax = (6 \lambda / C_{S3}) ^{1/4}
   return pow( 6.0 * lambda / Cs3, 0.25);
}

double TEM_NS::scherzer_approxresolution_bfctem_uncorrected(
      const double& lambda,
      const double& Cs3
      )
{
   // Kirkland (2016) Table 1
   // dmin = 0.67 (C_{S3} \lambda^{3}) ^{1/4}
   return 0.67 * pow( Cs3 * pow( lambda, 3.0), 0.25);
}

double TEM_NS::scherzer_defocus_adfstem_uncorrected(
      const double& lambda,
      const double& Cs3
      )
{
   // Kirkland (2016) Table 1
   // defocus = 0.87 (C_{S3} \lambda)^{1/2}
   return 0.87 * sqrt( Cs3 * lambda );
}

double TEM_NS::scherzer_alphamax_adfstem_uncorrected(
      const double& lambda,
      const double& Cs3
      )
{
   // Kirkland (2016) Table 1
   // alphamax = 1.34 ( \lambda / C_{S3}) ^{1/4}
   return 1.34 * pow( lambda / Cs3, 0.25);
}

double TEM_NS::scherzer_approxresolution_adfstem_uncorrected(
      const double& lambda,
      const double& Cs3
      )
{
   // Kirkland (2016) Table 1
   // dmin = 0.67 (C_{S3} \lambda^{3}) ^{1/4}
   return 0.43 * pow( Cs3 * pow( lambda, 3.0), 0.25);
}

double TEM_NS::scherzer_Cs3_bfctem_correctedtoCs5(
      const double& lambda,
      const double& Cs5
      )
{
   // Kirkland (2016) Table 2
   // C_{S3} = -3.2 ( \lambda C_{S5}^{2} )^{1/3}
   return -3.2 * pow( lambda * pow(Cs5, 2), 1.0/3.0);
}

double TEM_NS::scherzer_defocus_bfctem_correctedtoCs5(
      const double& lambda,
      const double& Cs5
      )
{
   // Kirkland (2016) Table 2
   // defocus =  -2 ( \lambda^{2} C_{S5} )^{1/3}
   return -2.0 * pow( pow( lambda, 2) * Cs5, 1.0/3.0);
}

double TEM_NS::scherzer_alphamax_bfctem_correctedtoCs5(
      const double& lambda,
      const double& Cs5
      )
{
   // Kirkland (2016) Table 2
   // alphamax = \frac{7}{4} ( \lambda / C_{S5} )^{1/6}
   return (7.0/4.0) * pow( lambda / Cs5, 1.0/6.0);
}

double TEM_NS::scherzer_approxresolution_bfctem_correctedtoCs5(
      const double& lambda,
      const double& Cs5
      )
{
   // Kirkland (2016) Table 2
   // dmin = \frac{4}{7} ( C_{S5} \lambda^{5} )^{1/6}
   return (4.0/7.0) * pow( Cs5 * pow( lambda, 5), 1.0/6.0);
}

double TEM_NS::scherzer_Cs3_adfstem_correctedtoCs5(
      const double& lambda,
      const double& Cs5
      )
{
   // Kirkland (2016) Table 2
   // C_{S3} = -2.289 ( \lambda C_{S5}^{2} )^{1/3}
   return -2.289 * pow( lambda * pow(Cs5, 2), 1.0/3.0);
}

double TEM_NS::scherzer_defocus_adfstem_correctedtoCs5(
      const double& lambda,
      const double& Cs5
      )
{
   // defocus =  -0.983 ( \lambda^{2} C_{S5} )^{1/3}
   return -0.983 * pow( pow( lambda, 2) * Cs5, 1.0/3.0);
}

double TEM_NS::scherzer_alphamax_adfstem_correctedtoCs5(
      const double& lambda,
      const double& Cs5
      )
{
   // Kirkland (2016) Table 2
   // alphamax = 1.513 ( \lambda / C_{S5} )^{1/6}
   return 1.513 * pow( lambda / Cs5, 1.0/6.0);
}

double TEM_NS::scherzer_approxresolution_adfstem_correctedtoCs5(
      const double& lambda,
      const double& Cs5
      )
{
   // Kirkland (2016) Table 2
   // dmin = 0.403 ( C_{S5} \lambda^{5} )^{1/6}
   return 0.403 * pow( Cs5 * pow( lambda, 5), 1.0/6.0);
}

#endif
