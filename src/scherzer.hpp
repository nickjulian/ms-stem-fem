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
// File: scherzer.hpp
// Purpose:

#ifndef SCHERZER_HPP
#define SCHERZER_HPP

#include <cstdlib>
#include <cmath>

#ifndef PI
#define PI 3.14159265369
#endif

namespace TEM_NS
{

   double scherzer_defocus_bfctem_uncorrected(
         const double& lambda,
         const double& Cs3
         );

   double scherzer_alphamax_bfctem_uncorrected(
         const double& lambda,
         const double& Cs3
         );

   double scherzer_approxresolution_bfctem_uncorrected(
         const double& lambda,
         const double& Cs3
         );

   double scherzer_defocus_adfstem_uncorrected(
         const double& lambda,
         const double& Cs3
         );

   double scherzer_alphamax_adfstem_uncorrected(
         const double& lambda,
         const double& Cs3
         );

   double scherzer_approxresolution_adfstem_uncorrected(
         const double& lambda,
         const double& Cs3
         );

   double scherzer_Cs3_bfctem_correctedtoCs5(
         const double& lambda,
         const double& Cs5
         );

   double scherzer_defocus_bfctem_correctedtoCs5(
         const double& lambda,
         const double& Cs5
         );

   double scherzer_alphamax_bfctem_correctedtoCs5(
         const double& lambda,
         const double& Cs5
         );

   double scherzer_approxresolution_bfctem_correctedtoCs5(
         const double& lambda,
         const double& Cs5
         );

   double scherzer_Cs3_adfstem_correctedtoCs5(
         const double& lambda,
         const double& Cs5
         );

   double scherzer_defocus_adfstem_correctedtoCs5(
         const double& lambda,
         const double& Cs5
         );

   double scherzer_alphamax_adfstem_correctedtoCs5(
         const double& lambda,
         const double& Cs5
         );

   double scherzer_approxresolution_adfstem_correctedtoCs5(
         const double& lambda,
         const double& Cs5
         );
}
#endif
