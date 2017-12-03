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
// File: scatterer.cpp
// Purpose:
//
// Use: 
//    1. instantiate a scatterer_param_LUT
//    2. instantiate your scatterers using the scatterer_param_LUT from 1, 
//       its position q[3], and its species (Z) 
// or:
//    1. instantiate a scatterer_param_LUT
//    2. For each given scattering position q[3] and corresponding species 
//       (Z), increment the scattering factor by fetching parameter 
//       addresses from the scatterer_param_LUT.get_param( Z )->a ,b,c,d ...

#ifndef SCATTERER_CPP
#define SCATTERER_CPP

#include "scatterer.hpp"



#endif
