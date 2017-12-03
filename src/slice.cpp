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
// File: slice.cpp
// Purpose:
//    - implement functions used to slice atomic systems and perform 
//       multislice TEM simulations

#ifndef SLICE_CPP
#define SLICE_CPP


#include <iostream>
#include <algorithm>    // std::sort

#include <iostream>  // cout
using std::cout;
using std::endl;

#include "slice.hpp"

#ifndef PI
#define PI 3.14159265369
#endif

int TEM_NS::delete_slices_from_list( list< slice* > sliceList )
{
   // Preconditions: 
   //    -  No memory location is pointed to by more than one member of 
   //        the slices contained in sliceList; each slice has unique pap,
   //        exp_i_sigma_v (transmission function), propagator, and 
   //        scatterers.
   //
   // This exists because calling delete on a list< slice* > won't delete
   // the slices pointed to by the listed pointers, but will only delete 
   // the pointers. 
   //
   // This function should be called before delete is called on 
   // list< slice* >.
   //
   // TODO: sort slices by memory location of their scatterers members and
   //       compare consecutive slices to decide if you're about to delete
   //       something twice.
   for ( list< slice* >::iterator sliceList_itr = sliceList.begin();
         sliceList_itr != sliceList.end();
         ++sliceList_itr )
   {
      delete[] (*sliceList_itr)->propagator_x_re;
      delete[] (*sliceList_itr)->propagator_x_im;
      delete[] (*sliceList_itr)->propagator_y_re;
      delete[] (*sliceList_itr)->propagator_y_im;

      fftw_free( (*sliceList_itr)->exp_i_sigma_v );
   }

   return EXIT_SUCCESS;
}


//int TEM_NS::delete_slices_from_list_ctem( list< slice* > sliceList )
//{
//   // Preconditions: 
//   //    -  No memory location is pointed to by more than one member of the
//   //        slices contained in sliceList; each slice has unique pap, 
//   //        exp_i_sigma_v (transmission function), propagator, and 
//   //        scatterers.
//   //
//   // This exists because calling delete on a list< slice* > won't delete
//   // the slices pointed to by the listed pointers, but will only delete the
//   // pointers. 
//   //
//   // This function should be called before delete is called on 
//   // list< slice* >.
//   //
//   // TODO: sort slices by memory location of their scatterers members and 
//   //       compare consequtive slices to decide if you're about to delete 
//   //       something twice.
//   for ( list< slice* >::iterator sliceList_itr = sliceList.begin();
//         sliceList_itr != sliceList.end();
//         ++sliceList_itr )
//   {
//      delete_members_of_slice_ctem( *sliceList_itr );
//      delete (*sliceList_itr); // (*sliceList_itr) is a pointer to a slice
//   }
//
//   return EXIT_SUCCESS;
//}

//int TEM_NS::delete_slices_from_list_stem( list< slice* > sliceList )
//{
//   // Preconditions: 
//   //    -  No memory location is pointed to by more than one member of the
//   //        slices contained in sliceList; each slice has unique pap, 
//   //        exp_i_sigma_v (transmission function), propagator, and 
//   //        scatterers.
//   //
//   // This exists because calling delete on a list< slice* > won't delete
//   // the slices pointed to by the listed pointers, but will only delete the
//   // pointers. 
//   //
//   // This function should be called before delete is called on 
//   // list< slice* >.
//   //
//   // TODO: sort slices by memory location of their scatterers members and 
//   //       compare consequtive slices to decide if you're about to delete 
//   //       something twice.
//   for ( list< slice* >::iterator sliceList_itr = sliceList.begin();
//         sliceList_itr != sliceList.end();
//         ++sliceList_itr )
//   {
//      delete_members_of_slice_stem( *sliceList_itr );
//      delete (*sliceList_itr); // (*sliceList_itr) is a pointer to a slice
//   }
//
//   return EXIT_SUCCESS;
//}


//int TEM_NS::delete_members_of_slice_ctem( slice* mySlice )
//{
//   // This function exists outside of the slice class because the members
//   //  of a slice are pointers to memory which is managed outside of the 
//   //  slice class.
//   
//   //TODO: check if the members to be deleted are allocated in memory
//
//   fftw_free( mySlice->pap ); // delete projected atomic potential
//
//   //fftw_free( mySlice->propagator );
//   //delete mySlice->propagator_x_re; // delete propagator
//   //delete mySlice->propagator_x_im;
//   //delete mySlice->propagator_y_re;
//   //delete mySlice->propagator_y_im;
//
//   //fftw_free( mySlice->exp_i_sigma_v );   // delete transmission function
//   delete mySlice->scatterers;// scatterers points to a vector<scatterer* >
//
//   return EXIT_SUCCESS;
//}
//
//int TEM_NS::delete_members_of_slice_stem( slice* mySlice )
//{
//   // This function exists outside of the slice class because the members
//   //  of a slice are pointers to memory which is managed outside of the 
//   //  slice class.
//   
//   //TODO: check if the members to be deleted are allocated in memory
//
//   fftw_free( mySlice->pap ); // delete projected atomic potential
//
//   //fftw_free( mySlice->propagator );
//   delete mySlice->propagator_x_re; // delete propagator
//   delete mySlice->propagator_x_im;
//   delete mySlice->propagator_y_re;
//   delete mySlice->propagator_y_im;
//
//   fftw_free( mySlice->exp_i_sigma_v );   // delete transmission function
//   delete mySlice->scatterers;// scatterers points to a vector<scatterer* >
//
//   return EXIT_SUCCESS;
//}


int TEM_NS::assign_atoms_to_slices_z(
      const vector< scatterer >& allScatterers,
      const vector< double >& slice_locations_z ,
      const double& xmin, const double& ymin, const double& zmin,
      const double& xperiod, const double& yperiod,
      const double& zperiod,
      list< slice* >& sliceList,
      const int& mynode,
      const int& rootnode,
      MPI_Comm comm
      )
{
   // Things to be aware of:
   // - slice_locations_z should be sorted elsewhere
   // - This function allocates memory of slice objects which need to
   //       be deleted externally by delete_slices_from_list() 
   ///////////////////////////////////////////////////////////////////////
   // Create a new slice for each valid slice_location and upper boundary
   ///////////////////////////////////////////////////////////////////////
   vector< double > slice_locations_copy;
   vector< double >::const_iterator sliceLocation_itr; 
   double zmax = zmin + zperiod;
   // copy and sort the locations of slice boundaries
   cout << "slice_locations_z : " << endl;//debug
   for ( sliceLocation_itr = slice_locations_z.begin(); 
         sliceLocation_itr != slice_locations_z.end(); 
         ++sliceLocation_itr )
   {
      if (
            (*sliceLocation_itr) >= zmin 
            && (*sliceLocation_itr) < zmax
         )
      {
         cout << (*sliceLocation_itr) << ", " ;//debug
         slice_locations_copy.push_back( *sliceLocation_itr );
      }
   }

   if ( slice_locations_copy.empty() )
   {  // require at least one slice 
      cout << "slice_locations_copy.size() : " //debug
         << slice_locations_copy.size() << endl; //debug
      slice_locations_copy.push_back( zmin );
   }

   //std::sort( slice_locations_copy.begin(), 
   //            slice_locations_copy.end() 
   //            );
   //debug
   cout << endl << "slice locations : " ;
   for ( sliceLocation_itr = slice_locations_copy.begin();
         sliceLocation_itr != slice_locations_copy.end();
         ++sliceLocation_itr)
   {
      cout << (*sliceLocation_itr) << ", ";
   }
   cout << endl << endl;
   //enddebug

   // instantiate the slices for each valid slice_location
   double thickness; 
   list< slice* >::iterator sliceList_itr;

   //  TODO: if first slice doesn't reach zmin, create a new slice having
   //          lower_bound = zmin 
   
   for ( sliceLocation_itr = slice_locations_copy.begin(); 
         sliceLocation_itr != slice_locations_copy.end(); // true if empty
         ++sliceLocation_itr )
   {
      // create a new slice, set its lower_bound and initiate its vector
      //  of scatterer pointers
      sliceList.push_back( new slice );
      (sliceList.back())->lower_bound = *sliceLocation_itr;
      (sliceList.back())->scatterers = new vector< const scatterer* >;
   }


   // Evaluate slice thickness
   if ( sliceList.front() == sliceList.back() ) 
   { // case in which only one slice exists
      (*(sliceList.begin()))->thickness = zperiod;
   }
   else
   {
      //list< slice* >::iterator sliceList_itr = sliceList.begin();
      sliceList_itr = sliceList.begin();

      double tmp_double;

      while ( *sliceList_itr != sliceList.back() )
      {
         tmp_double = (*sliceList_itr)->lower_bound;
         (*sliceList_itr)->thickness 
            = (*(++sliceList_itr))->lower_bound - tmp_double;
             //Note: questionable ++
            //+ (*(sliceList_itr + 1))->lower_bound;// does this work? No
         //++sliceList_itr;
      }
      // for the last slice assume periodicity in calculating thickness 
      (*sliceList_itr)->thickness // sliceList_itr == sliceList.back()
         = zperiod 
            - ((*sliceList_itr)->lower_bound) 
            + (*(sliceList.begin()))->lower_bound;
   }

   // Filling of sliceList should now be complete, and the slices it
   //  contains should all have true thickness and lower_bound member
   //  values.

   ///////////////////////////////////////////////////////////////////////
   // Create local list or vector of scatterers to be sorted and binned
   ///////////////////////////////////////////////////////////////////////
   // TODO: parallelize sorting and binning of the scatterers
   // Note: why is allScatterers a vector instead of a list? for sorting?
   //const vector< scatterer >& allScatterers,
   vector< const scatterer* > localScatterers; 
   //vector< scatterer* >::iterator scattererVec_itr; 
   
   for ( vector< scatterer >::const_iterator allScattererVec_itr 
            =  allScatterers.begin(); 
        allScattererVec_itr != allScatterers.end(); 
         ++allScattererVec_itr)
   {
      localScatterers.push_back( &(*allScattererVec_itr) ); 
   }
   
   ///////////////////////////////////////////////////////////////////////
   // Sort all scatterers using scatterer_sorting_lt_z() as "<" operator.
   ///////////////////////////////////////////////////////////////////////
   std::sort( localScatterers.begin(),
               localScatterers.end(),
               scatterer_sorting_lt_z );
   //debug: check that all local scatterers are sorted by z value
   //double tmp_z = ( *(localScatterers.begin() ) )->q[2];
   //cout << "z positions of sorted scatterers : " << endl;;
   //for (vector< const scatterer* >::iterator 
   //         vec_itr = localScatterers.begin(); 
   //      vec_itr != localScatterers.end(); ++vec_itr)
   //{
   //   if( tmp_z > (*vec_itr)->q[2] )
   //   {
   //      cout << "localScatterers not sorted,  tmp_z > (*vec_itr)->q[2] : "
   //         << tmp_z << " > " << (*vec_itr)->q[2] << endl;
   //      return EXIT_FAILURE;
   //   }
   //      cout << (*vec_itr)->q[2] << endl;
   //}
   //end debug

   // TODO: consider sorting slice_locations_z to ensure they're ordered

   ///////////////////////////////////////////////////////////////////////
   // Iterate through the scatterers and assign them to slices
   ///////////////////////////////////////////////////////////////////////
   
   // Access and assignment examples:
   // (*sliceList_itr) // pointer to a slice
   // (*sliceList_itr)->scatterers //pointer to vector of scatterer pointers
   //
   // *((*sliceList_itr)->scatterers) // vector of scatterer pointers
   // (*((*sliceList_itr)->scatterers)).push_back( somepointertoascatterer);
   //
   // ((*sliceList_itr)->scatterers)->push_back( somepointertoascatterer );
   //
   // someSlice = *sliceList_itr;
   // (someSlice.scatterers)->push_back( somePointerToAScatterer );
   // scattererPtr_itr = (someSlice.scatterers)->begin();
   // scattererPtr_itr = ((*sliceList_itr).scatterers)->begin();

   //sliceLocation_itr = slice_locations_z.begin();
   //vector< scatterer* > scattererVec_itr = localScatterers.begin();

   // list< vector< scatterer* >* > slicesOfScatterers
   // vector< scatterer* > localScatterers; 

   


   // If the first slice lower_bound is not equal to zmin, insert the 
   //  scatterers found below that first lower_bound into the last slice,
   //  not the first slice, as if they've travelled beyond a periodic 
   //  boundary.
   vector< const scatterer* >::iterator scattererVec_itr; 
   scattererVec_itr = localScatterers.begin();
   sliceList_itr = --(sliceList.end());
   for ( ;// if a scatterer has z < the first slice lower_bound, then
         (*scattererVec_itr)->q[2] < (sliceList.front())->lower_bound;
         ++scattererVec_itr )
   {// push scatterers onto the back slice
      ((*sliceList_itr)->scatterers)->push_back( *scattererVec_itr );
   }// TODO: consider replacing this for() loop with while().


   // For each remaining slice, iterate the scatterer iterator until they 
   //  exceed that slices upper bound.
   double slice_upper_bound;
   sliceList_itr = sliceList.begin();
   if ( sliceList_itr == --(sliceList.end()) )// TODO: is there a better way
   {
      slice_upper_bound = zmax;
   }
   else
   {
      ++sliceList_itr;
      slice_upper_bound = (*sliceList_itr)->lower_bound;
      --sliceList_itr;
   }
   // I don't have to worry about sliceList being empty because line 90
   //  instantiates a slice if the list of slice locations is empty.

   for ( ; //scattererVec_itr = localScatterers.begin();
         scattererVec_itr != localScatterers.end();
         ++scattererVec_itr )
   {
      // this assumes that all scatterers are ordered by q[2] position
      // TODO: generalize this to an arbitrary ordering orientation
      //       instead of along q[2]
      while ( (*scattererVec_itr)->q[2] >= slice_upper_bound )
      {
         ++sliceList_itr; // move on and begin assigning to the next slice
         // update slice_upper_bound
         if ( sliceList_itr == sliceList.end() )
         {
            slice_upper_bound = zmax;
            --sliceList_itr;  
            break;
         }
         else
         if ( sliceList_itr == --(sliceList.end()) )
         {
            slice_upper_bound = zmax;
         }
         else 
         {
            ++sliceList_itr;
            slice_upper_bound = (*sliceList_itr)->lower_bound;
            --sliceList_itr;
         }
         //if ( slice_upper_bound > zmax ) slice_upper_bound = zmax;
      }
      // add the current scatterer pointer to the current slice
      ((*sliceList_itr)->scatterers)->push_back( *scattererVec_itr );
   }

   // TODO: use a different scatterer_sorting_lt function to allow the 
   //       scatterers to be sorted along an arbitrary direction
   
   //debug
   //size_t slice_number = 0;
   //for ( sliceList_itr = sliceList.begin(); 
   //      sliceList_itr != sliceList.end(); 
   //      ++sliceList_itr )
   //{
   //   ++slice_number;
   //   cout << endl << "node " << mynode << ", ";
   //   cout << "slice #       : " << slice_number << endl;
   //   cout << "node " << mynode << ", ";
   //   cout << " lower_bound  : " << (*sliceList_itr)->lower_bound << endl;
   //   cout << "node " << mynode << ", ";
   //   cout << " thickness    : " << (*sliceList_itr)->thickness<< endl;
   //   cout << "node " << mynode << ", ";
   //   cout << " scatterer zs : " ;
   //   for ( scattererVec_itr = ((*sliceList_itr)->scatterers)->begin();
   //          scattererVec_itr != ((*sliceList_itr)->scatterers)->end();
   //          ++ scattererVec_itr )
   //   {
   //      cout << (*scattererVec_itr)->q[2] << ", " ;
   //   }
   //   cout << endl << endl;
   //}
   //end debug

   return EXIT_SUCCESS;
}

int TEM_NS::assign_atoms_to_slices_z_auto(
      const vector< scatterer >& allScatterers,
      const double& xmin, const double& ymin, const double& zmin,
      const double& xperiod, const double& yperiod,   //Note: these aren't
      const double& zperiod,                          // used unless there
      const double& min_thickness,                    // are no scatterers
      list< slice* >& sliceList,
      const unsigned int& input_flag_debug,
      const int& mynode,
      const int& rootnode,
      MPI_Comm comm
      )
{
   // Things to be aware of:
   // - This function allocates memory of slice objects which need to
   //       be deleted externally by delete_slices_from_list() 
   // -  This algorithm reduces the total sample thickness
   //       from zmax - zmin to a value dependent upon scatterer 
   //       locations.


   ///////////////////////////////////////////////////////////////////////
   // Create local list or vector of scatterers to be sorted and binned
   ///////////////////////////////////////////////////////////////////////
   // TODO: parallelize sorting and binning of the scatterers
   // Note: why is allScatterers a vector instead of a list? for sorting?
   //const vector< scatterer >& allScatterers,
   vector< const scatterer* > localScatterers; 
   //vector< scatterer* >::iterator scattererVec_itr; 
   
   for ( vector< scatterer >::const_iterator allScattererVec_itr 
            =  allScatterers.begin(); 
        allScattererVec_itr != allScatterers.end(); 
         ++allScattererVec_itr)
   {
      localScatterers.push_back( &(*allScattererVec_itr) ); 
   }
   
   ///////////////////////////////////////////////////////////////////////
   // Sort all scatterers using scatterer_sorting_lt_z() as "<" operator.
   ///////////////////////////////////////////////////////////////////////
   // TODO: use a different scatterer_sorting_lt function to allow the 
   //       scatterers to be sorted along an arbitrary direction. 
   //       Also replace every use of z coordinate with a function which 
   //       returns the scatterer position relative to an arbitrary 
   //       direction.
   
   std::sort( localScatterers.begin(),
               localScatterers.end(),
               scatterer_sorting_lt_z );
   //debug: check that all local scatterers are sorted by z value
   //double tmp_z = ( *(localScatterers.begin() ) )->q[2];
   //cout << "node " << mynode << ", ";
   //cout << "z positions of sorted scatterers : " << endl;;
   //for (vector< const scatterer* >::iterator 
   //         vec_itr = localScatterers.begin(); 
   //      vec_itr != localScatterers.end(); ++vec_itr)
   //{
   //   if( tmp_z > (*vec_itr)->q[2] )
   //   {
   //      cout << "node " << mynode << ", ";
   //      cout << "localScatterers not sorted,  tmp_z > (*vec_itr)->q[2] : "
   //         << tmp_z << " > " << (*vec_itr)->q[2] << endl;
   //      return EXIT_FAILURE;
   //   }
   //      cout << (*vec_itr)->q[2] << endl;
   //}
   //end debug

   ///////////////////////////////////////////////////////////////////////
   // Determine the locations to place slice boundaries
   ///////////////////////////////////////////////////////////////////////
   // When a large region of empty space exists between atoms, should a
   //  series of empty slices be made, or could it be accounted for in 
   //  the preceeding propagator?
   //
   // \psi_{n+1}(x,y) = FT^{-1}[
   //    P_{n}(k_{x}, k_{y}, \Delta z_{n})
   //    FT[ t_{n}(x,y) \psi_{n}(x,y) ]
   // ]
   // If the n'th slice is empty, then 
   //    t_{n}(x,y) = e^{i \sigma v(x,y)} = e^{i \sigma 0 } = 1
   // and 
   // \psi_{n+1}(x,y) 
   // = FT^{-1}[
   //       P_{n}(k_{x}, k_{y}, \Delta z_{n})
   //       FT[ 1 * \psi_{n}(x,y) ]
   //    ]
   // = FT^{-1}[
   //       P_{n}(k_{x}, k_{y}, \Delta z_{n})
   //       FT[ 
   //          FT^{-1}[
   //             P_{n-1}(k_{x}, k_{y}, \Delta z_{n-1})
   //             FT[ t_{n-1}(x,y) \psi_{n-1}(x,y) ]
   //          ]
   //       ]
   //    ]
   // = FT^{-1}[
   //       P_{n}(k_{x}, k_{y}, \Delta z_{n})
   //       P_{n-1}(k_{x}, k_{y}, \Delta z_{n-1})
   //       FT[ t_{n-1}(x,y) \psi_{n-1}(x,y) ]
   //    ]
   // = FT^{-1}[
   //       \exp[ -i \pi \lambda (k^{2}_{x} + k^{2}_{y}) \Delta z_{n}]
   //       \exp[ -i \pi \lambda (k^{2}_{x} + k^{2}_{y}) \Delta z_{n-1}]
   //       FT[ t_{n-1}(x,y) \psi_{n-1}(x,y) ]
   //    ] 
   // = FT^{-1}[
   //       \exp[ 
   //          -i \pi \lambda 
   //          (k^{2}_{x} + k^{2}_{y}) 
   //          (\Delta z_{n} + \Delta z_{n-1})
   //       ]
   //       FT[ t_{n-1}(x,y) \psi_{n-1}(x,y) ]
   //    ] 
   // \implies inserting empty slices is equivalent to using a 
   //          propagator whose thickness is the sum of the empty slice
   //          thicknesses.
   // The same equation shows that as long as the 3-D potential of a 
   //  slice at a certain position is approximately zero, then it has
   //  approximately no influence on the preceeding slice, which may
   //  thus be arbitrarily thin.
   //  

   // x Iterate through the scatterers tracking slice thickness. When 
   // x  thickness reaches a critical value (~1\AA), push a new cut 
   // x  location onto slice_locations_z. 
   // x  But, if some other atoms have very similar but lower z, move the
   // x  cut location to their position ??? No.
   // x  Also consider distance between consequtive scatterers, making 
   // x   sure not to put them on separate slices if they're z values are
   // x   very close together ??? No.
   // x  Ensure that the final slice thickness is at least min_thickness.
   //
   // Minimum slice thickness: ~ 1\AA.
   // x Maximum slice thickness: ??? if there are no more scatterers, then
   // x  Inf, or the upper bound on z should then be lowered...
   //
   // x Split z into increments of 1 \AA. Group scatterers by these 
   // x divisions. Assign each group of scatterers to a slice located at
   // x the average z position of the scatterers. Calculate slice 
   // x thicknesses.
   // x ^this could result in slice thicknesses < 1\AA.
   //
   // Iterate through the scatterers, adding them to a slice. If a
   // scatterrers position is > min_thickness away from the lower bound of
   // the group, then begin a new slice with it.
   // After creating the slices, calculate their average z and 
   // assign that to their cut location (lower_bound).
   // After identifying all cut locations, calculate slice thicknesses.
   // Ensure that the final slice thickness is at least min_thickness.
   //
   // Why shouldn't a slice be < 1\AA thick? Kirkland (2010):
   // "v_{\Delta z}(x,y) = \int^{z + \Delta z}_{z} V(x,y,z) dz
   // \approx
   // \int^{+Inf}_{-Inf} V(x,y,z) dz = v_{z}(x,y)
   // As long as the slice thickness is large compared to the effective
   // range of the atomic potential this approximation is valid."
   
   ///////////////////////////////////////////////////////////////////////
   // Handle the case of an empty scatterer list
   ///////////////////////////////////////////////////////////////////////

   if ( allScatterers.size() == 0 ) 
   {// add an empty slice occupying the entire zperiod
      sliceList.push_back( new slice );
      (sliceList.back())->lower_bound = zmin;
      (sliceList.back())->thickness = zperiod;//TODO: generalize direction
      (sliceList.back())->scatterers = new vector< const scatterer* >;
      return EXIT_SUCCESS;
   }//TODO: make the rest of the code perform equivalently and remove this

   ///////////////////////////////////////////////////////////////////////
   // Iterate through the scatterers, creating slices and pushing 
   //  scatterers onto them
   ///////////////////////////////////////////////////////////////////////

   sliceList.push_back( new slice );// ensure at least one slice exists

   list< slice* >::iterator sliceList_itr;
   sliceList_itr = sliceList.begin();
   (*sliceList_itr)->scatterers = new vector< const scatterer* >;
   (*sliceList_itr)->lower_bound = 0.0;
   //(sliceList.back())->scatterers = new vector< const scatterer* >;

   double tmp_lower_bound; 
   tmp_lower_bound = (localScatterers.front())->q[2];

   for ( vector< const scatterer* >::iterator 
            scattererVec_itr = localScatterers.begin();
         scattererVec_itr != localScatterers.end(); 
         ++scattererVec_itr )
   {
      if ( ((*scattererVec_itr)->q[2]) - tmp_lower_bound <= min_thickness)
      {// if 'current scatterer z' - 'first scatterers z' <= min_thickness
         // append the current scatterer to the current slice
         ((*sliceList_itr)->scatterers)->push_back( *scattererVec_itr );

         // accumulate z positions to be averaged later
         (*sliceList_itr)->lower_bound += (*scattererVec_itr)->q[2];
         // NOTE: (*sliceList_itr)->lower_bound will be divided by the 
         //       the number of scatterers later to set the cut locations
         //       to the average scatterer position in each slice.
      }
      else
      { // create a new slice containing the current scatterer
         sliceList.push_back( new slice );
         ++sliceList_itr;
         (*sliceList_itr)->scatterers = new vector< const scatterer* >;
         ((*sliceList_itr)->scatterers)->push_back( *scattererVec_itr );
         (*sliceList_itr)->lower_bound = (*scattererVec_itr)->q[2];
         tmp_lower_bound = (*scattererVec_itr)->q[2];
      }
   }
   
   ///////////////////////////////////////////////////////////////////////
   // Iterate through the slices, calculating the average z location
   //  of scatterers and assigning that to the slice lower_bound
   ///////////////////////////////////////////////////////////////////////

   tmp_lower_bound = 0.0;
   for ( sliceList_itr = sliceList.begin();
         sliceList_itr != sliceList.end();
         ++sliceList_itr )
   {  // For each slice, divide the sum of scatterer locations by the 
      // number of scatterers.
      (*sliceList_itr)->lower_bound = 
         (*sliceList_itr)->lower_bound
         /
         ((*sliceList_itr)->scatterers)->size(); 
   }

   ///////////////////////////////////////////////////////////////////////
   // Iterate through the slices again, calculating and assigning slice 
   //  thicknesses
   ///////////////////////////////////////////////////////////////////////
   //tmp_lower_bound = (sliceList.front())->lower_bound;
   list< slice* >::iterator sliceList_itr2;
   sliceList_itr2 = sliceList.begin();
   //TODO: would it be faster to use only one iterator instead of two?
   for ( sliceList_itr = sliceList.begin();
         sliceList_itr != sliceList.end();
         ++sliceList_itr )
   { // TODO: This algorithm reduces the total sample thickness
     //        from zmax - zmin to a value dependent upon scatterer 
     //        locations. Consider changing to preserve input thickness.
      ++sliceList_itr2;// itr2 is at the slice following the first itr
      if ( sliceList_itr2 == sliceList.end() )
      {// The final slice is made to extend min_thickness beyond the 
       //  position of it highest scatterer.
         (*sliceList_itr)->thickness
            = (((*sliceList_itr)->scatterers)->back())->q[2] 
               + min_thickness
               - (*sliceList_itr)->lower_bound;
      }
      else
      {
         (*sliceList_itr)->thickness 
            = (*sliceList_itr2)->lower_bound 
               - (*sliceList_itr)->lower_bound;
      }
   }

//   /////////////////////////////////////////////////////////////////////
//   // Create a new slice for each valid slice_location and upper 
//   //  boundary
//   /////////////////////////////////////////////////////////////////////
//   vector< double >::const_iterator sliceLocation_itr; 
//   double zmax = zmin + zperiod;
//   // copy and sort the locations of slice boundaries
//   cout << "slice_locations_z : " << endl;//debug
//   //for ( sliceLocation_itr = slice_locations_z.begin(); 
//   //      sliceLocation_itr != slice_locations_z.end(); 
//   //      ++sliceLocation_itr )
//   //{
//   //   if (
//   //         (*sliceLocation_itr) >= zmin 
//   //         && (*sliceLocation_itr) < zmax
//   //      )
//   //   {
//   //      cout << (*sliceLocation_itr) << ", " ;//debug
//   //      slice_locations_copy.push_back( *sliceLocation_itr );
//   //   }
//   //}
//
//   if ( slice_locations_z.empty() )
//   {  // require at least one slice 
//      cout << "slice_locations_z.size() : " //debug
//         << slice_locations_z.size() << endl; //debug
//      slice_locations_z.push_back( zmin );
//   }
//
//   //std::sort( slice_locations_copy.begin(), 
//   //            slice_locations_copy.end() 
//   //            );
//   //debug
//   cout << endl << "slice locations : " ;
//   for ( sliceLocation_itr = slice_locations_z.begin();
//         sliceLocation_itr != slice_locations_z.end();
//         ++sliceLocation_itr)
//   {
//      cout << (*sliceLocation_itr) << ", ";
//   }
//   cout << endl << endl;
//   //enddebug
//
//   // instantiate the slices for each valid slice_location
//   double thickness; 
//   list< slice* >::iterator sliceList_itr;
//
//   // TODO: if first slice doesn't reach zmin, create a new slice having
//   //          lower_bound = zmin 
//   
//   for ( sliceLocation_itr = slice_locations_z.begin(); 
//         sliceLocation_itr != slice_locations_z.end(); // true if empty
//         ++sliceLocation_itr )
//   {
//      // create a new slice, set its lower_bound and initiate its vector
//      //  of scatterer pointers
//      sliceList.push_back( new slice );
//      (sliceList.back())->lower_bound = *sliceLocation_itr;
//      (sliceList.back())->scatterers = new vector< const scatterer* >;
//   }
//
//
//   // Evaluate slice thickness
//   if ( sliceList.front() == sliceList.back() ) 
//   { // case in which only one slice exists
//      (*(sliceList.begin()))->thickness = zperiod;
//   }
//   else
//   {
//      //list< slice* >::iterator sliceList_itr = sliceList.begin();
//      sliceList_itr = sliceList.begin();
//
//      double tmp_double;
//
//      while ( *sliceList_itr != sliceList.back() )
//      {
//         tmp_double = (*sliceList_itr)->lower_bound;
//         (*sliceList_itr)->thickness 
//            = (*(++sliceList_itr))->lower_bound - tmp_double;
//             //TODO: questionable ++
//            //+ (*(sliceList_itr + 1))->lower_bound;// does this work? No
//         //++sliceList_itr;
//      }
//      // for the last slice assume periodicity in calculating thickness 
//      (*sliceList_itr)->thickness // sliceList_itr == sliceList.back()
//         = zperiod 
//            - ((*sliceList_itr)->lower_bound) 
//            + (*(sliceList.begin()))->lower_bound;
//   }
//
//   // Filling of sliceList should now be complete, and the slices it
//   //  contains should all have true thickness and lower_bound member
//   //  values.
//
//   // TODO: consider sorting slice_locations_z to ensure they're ordered
//
//   ///////////////////////////////////////////////////////////////////////
//   // Iterate through the scatterers and assign them to slices
//   ///////////////////////////////////////////////////////////////////////
//   
//   // Access and assignment examples:
//   // (*sliceList_itr) // pointer to a slice
//   // (*sliceList_itr)->scatterers //pointer to vector of scatterer pointers
//   //
//   // *((*sliceList_itr)->scatterers) // vector of scatterer pointers
//   // (*((*sliceList_itr)->scatterers)).push_back( somepointertoascatterer);
//   //
//   // ((*sliceList_itr)->scatterers)->push_back( somepointertoascatterer );
//   //
//   // someSlice = *sliceList_itr;
//   // (someSlice.scatterers)->push_back( somePointerToAScatterer );
//   // scattererPtr_itr = (someSlice.scatterers)->begin();
//   // scattererPtr_itr = ((*sliceList_itr).scatterers)->begin();
//
//   //sliceLocation_itr = slice_locations_z.begin();
//   //vector< scatterer* > scattererVec_itr = localScatterers.begin();
//
//   // list< vector< scatterer* >* > slicesOfScatterers
//   // vector< scatterer* > localScatterers; 
//
//   
//
//
//   // TODO: correct these assignments of scatterer pointers to slices
//   
//   // If the first slice lower_bound is not equal to zmin, insert the 
//   //  scatterers found below that first lower_bound into the last slice,
//   //  not the first slice, as if they've travelled beyond a periodic 
//   //  boundary.
//   vector< const scatterer* >::iterator scattererVec_itr; 
//   scattererVec_itr = localScatterers.begin();
//   sliceList_itr = --(sliceList.end());
//   for ( ;// if a scatterer has z < the first slice lower_bound, then
//         (*scattererVec_itr)->q[2] < (sliceList.front())->lower_bound;
//         ++scattererVec_itr )
//   {// push scatterers onto the back slice
//      ((*sliceList_itr)->scatterers)->push_back( *scattererVec_itr );
//   }// TODO: I don't think this is a good thing to have done
//
//
//   // For each remaining slice, iterate the scatterer iterator until they 
//   //  exceed that slices upper bound.
//   double slice_upper_bound;
//   sliceList_itr = sliceList.begin();
//   if ( sliceList_itr == --(sliceList.end()) )// TODO: is there a better way
//   {
//      slice_upper_bound = zmax;
//   }
//   else
//   {
//      ++sliceList_itr;
//      slice_upper_bound = (*sliceList_itr)->lower_bound;
//      --sliceList_itr;
//   }
//   // I don't have to worry about sliceList being empty because line 90
//   //  instantiates a slice if the list of slice locations is empty.
//
//   for ( ; //scattererVec_itr = localScatterers.begin();
//         scattererVec_itr != localScatterers.end();
//         ++scattererVec_itr )
//   {
//      // this assumes that all scatterers are ordered by q[2] position
//      // TODO: generalize this to an arbitrary ordering orientation
//      //       instead of along q[2]
//      while ( (*scattererVec_itr)->q[2] >= slice_upper_bound )
//      {
//         ++sliceList_itr; // move on and begin assigning to the next slice
//         // update slice_upper_bound
//         if ( sliceList_itr == sliceList.end() )
//         {
//            slice_upper_bound = zmax;
//            --sliceList_itr;  
//            break;
//         }
//         else
//         if ( sliceList_itr == --(sliceList.end()) )
//         {
//            slice_upper_bound = zmax;
//         }
//         else 
//         {
//            ++sliceList_itr;
//            slice_upper_bound = (*sliceList_itr)->lower_bound;
//            --sliceList_itr;
//         }
//         //if ( slice_upper_bound > zmax ) slice_upper_bound = zmax;
//      }
//      // add the current scatterer pointer to the current slice
//      ((*sliceList_itr)->scatterers)->push_back( *scattererVec_itr );
//   }

   //debug
   //size_t slice_number = 0;
   //for ( sliceList_itr = sliceList.begin(); 
   //      sliceList_itr != sliceList.end(); 
   //      ++sliceList_itr )
   //{
   //   ++slice_number;
   //   cout << "node " << mynode << ", ";
   //   cout << " slice #       : " << slice_number << endl;
   //   cout << "node " << mynode << ", ";
   //   cout << " lower_bound  : " << (*sliceList_itr)->lower_bound << endl;
   //   cout << "node " << mynode << ", ";
   //   cout << " thickness    : " << (*sliceList_itr)->thickness<< endl;
   //   cout << "node " << mynode << ", ";
   //   cout << " scatterer zs : " ;
   //   for ( vector< const scatterer* >::iterator 
   //            scattererVec_itr = ((*sliceList_itr)->scatterers)->begin();
   //          scattererVec_itr != ((*sliceList_itr)->scatterers)->end();
   //          ++scattererVec_itr )
   //   {
   //      cout << (*scattererVec_itr)->q[2] << ", " ;
   //   }
   //   cout << endl << endl;
   //}
   //end debug

   // debug
   if ( mynode == rootnode && input_flag_debug )
   {
      cout << "Created " << sliceList.size() << " slices:" << endl;
      size_t slice_number = 0;
      for ( sliceList_itr = sliceList.begin(); 
            sliceList_itr != sliceList.end(); 
            ++sliceList_itr )
      {
         ++slice_number;
         //cout << "node " << mynode << ", ";
         cout << " slice #       : " << slice_number << endl;
         //cout << "node " << mynode << ", ";
         cout << " lower_bound  : " << (*sliceList_itr)->lower_bound << endl;
         //cout << "node " << mynode << ", ";
         cout << " thickness    : " << (*sliceList_itr)->thickness<< endl;
         //cout << "node " << mynode << ", ";
         cout << " scatterer zs : " ;
         for ( vector< const scatterer* >::iterator 
                  scattererVec_itr = ((*sliceList_itr)->scatterers)->begin();
                scattererVec_itr != ((*sliceList_itr)->scatterers)->end();
                ++scattererVec_itr )
         {
            cout << (*scattererVec_itr)->q[2] << ", " ;
         }
         cout << endl << endl;
      }
   }
   // end debug

   return EXIT_SUCCESS;
}

//int TEM_NS::evaluate_propagator_ms(
//      const double& lambda,
//      const double* const kxdomain,
//      const ptrdiff_t& Nx,
//      const double* const kydomain,
//      const ptrdiff_t& Ny,
//      const double& delta_z,
//      fftw_complex* propagator_x,
//      fftw_complex* propagator_y
//      )
//{
//   // P_{n}(k_{x},k_{y},\Delta z_{n})
//   //    = \exp[-i \pi \lambda (k^{2}_{x} + k^{2}_{y}) \Delta z_{n}]
//   //    = \exp[-i \pi \lambda k^{2}_{x} \Delta z_{n}]
//   //      * \exp[-i \pi \lambda k^{2}_{y} \Delta z_{n}]
//   double trig_operand_x;
//   double trig_operand_y;
//   for ( ptrdiff_t i=0; i<Nx; ++i)
//      {
//         trig_operand_x = PI * lambda * delta_z *  pow(kxdomain[i], 2);
//         propagator_x[ i ][0] = cos( trig_operand );
//         propagator_x[ i ][1] = -sin( trig_operand );
//      }
//
//   for ( ptrdiff_t j=0; j<Ny; ++j)
//      {
//         trig_operand_y = PI * lambda * delta_z *  pow(kydomain[j], 2);
//         propagator_y[ j ][0] = cos( trig_operand );
//         propagator_y[ j ][1] = -sin( trig_operand );
//      }
//   return EXIT_SUCCESS;
//}

//int TEM_NS::propagate_ms(
//      const double& lambda,
//      const double* const kxdomain,
//      const ptrdiff_t& Nx,
//      const double* const kydomain,
//      const ptrdiff_t& Ny,
//      //const double& sqrtNxNy,
//      const double& bwcutoff_t,
//      const double& delta_z,
//      fftw_complex* ft_t_psi  // FT[t_{n}(x,y) \psi_{n}(x,y)]
//      )
//{
//   // Evaluate 
//   //    P_{n}(k_{x},k_{y},\Delta z_{n}) * FT[t_{n}(x,y) \psi_{n}(x,y)] 
//   // where
//   // P_{n}(k_{x},k_{y},\Delta z_{n})
//   //    = \exp[-i \pi \lambda (k^{2}_{x} + k^{2}_{y}) \Delta z_{n}]
//   double k_sqr;
//   double trig_operand;
//   double propagator_re;//, *propagator_x_re, *propagator_y_re;
//   double propagator_im;//, *propagator_x_im, *propagator_y_im;
//   double ft_t_psi_re;
//   double ft_t_psi_im;
//   const double bwcutoff_t_sqr = pow(bwcutoff_t, 2);
//   // Kirkland (2010) page 146:
//   // "A convenient way of limiting both functions is to set both
//   // the speciment transmission function and the propagator function
//   // to zero outside of (2/3)k_{max}." k_max = min{ kx_max/2, ky_max/2 }
//   // 
//   //  Condition for eliminating aliasing:
//   //  bwcutoff_t + bwcutoff_propagator + bwcutoff_psi < 2 * k_max,
//   //
//   //propagator_x_re = new double[Nx];
//   //propagator_x_im = new double[Nx];
//   //propagator_y_re = new double[Ny];
//   //propagator_y_im = new double[Ny];
//
//   //for ( ptrdiff_t i=0; i<Nx; ++i)  
//   //{
//   //   trig_operand = PI * lambda * delta_z * pow( kxdomain[i], 2) ;
//   //   propagator_x_re[i] = cos( trig_operand );
//   //   propagator_x_im[i] = -sin( trig_operand );
//   //}
//
//   //for ( ptrdiff_t j=0; j<Ny; ++j)
//   //{
//   //   trig_operand = PI * lambda * delta_z * pow( kydomain[j], 2) ;
//   //   propagator_y_re[j] = cos( trig_operand );
//   //   propagator_y_im[j] = -sin( trig_operand );
//   //}
//
//   //for ( ptrdiff_t i=0; i<Nx; ++i)  
//   //   for ( ptrdiff_t j=0; j<Ny; ++j)
//   //   {
//   //      k_sqr =  pow(kxdomain[i], 2) + pow(kydomain[j], 2);
//   //      if ( k_sqr < bwcutoff_t_sqr ) 
//   //      {
//   //         // complex multiplication:
//   //         //  (a + ib)(c + id) = ac - bd + i(ad + bc)
//   //         // P_{n} = P_{nx} * P_{ny} 
//   //         propagator_re 
//   //            = propagator_x_re[i] * propagator_y_re[j]
//   //               - propagator_x_im[i] * propagator_y_im[j];
//   //         propagator_im = -sin( trig_operand );
//
//   //         propagator_im
//   //            = propagator_x_re[i] * propagator_y_im[j]
//   //               + propagator_x_im[i] * propagator_y_re[j];
//
//   //         ft_t_psi_re = ft_t_psi[j + i * Ny][0];
//   //         ft_t_psi_im = ft_t_psi[j + i * Ny][1];
//
//   //         ft_t_psi[j + i * Ny][0] 
//   //            = ( propagator_re * ( ft_t_psi_re ) 
//   //                - propagator_im * ( ft_t_psi_im ) );
//   //                //- propagator_im * ( ft_t_psi_im ) ) / sqrtNxNy;
//
//   //         ft_t_psi[j + i * Ny][1] 
//   //            = ( propagator_re * ft_t_psi_im  
//   //                + propagator_im * ft_t_psi_re );
//   //                //+ propagator_im * ft_t_psi_re ) / sqrtNxNy;
//   //      }
//   //      else
//   //      {
//   //         ft_t_psi[j + i * Ny][0] = 0.0;
//   //         ft_t_psi[j + i * Ny][1] = 0.0;
//   //      }
//   //   }
//   // delete[] propagator_x_re; 
//   // delete[] propagator_x_im;
//   // delete[] propagator_y_re;
//   // delete[] propagator_y_im;
//
//
//   // The following version works, but calculates sin() & cos() Nx*Ny times,
//   //  whereas the above version calculates sin() & cos() Nx+Ny times at the
//   //  expense of Nx*Ny more complex multiplications.
//   // However, the above version spends time allocating memory which should
//   //  probably be allocated by the calling function and reused.
//   for ( ptrdiff_t i=0; i<Nx; ++i)  
//      for ( ptrdiff_t j=0; j<Ny; ++j)
//      {
//         k_sqr =  pow(kxdomain[i], 2) + pow(kydomain[j], 2);
//         if ( k_sqr < bwcutoff_t_sqr ) 
//         {
//            trig_operand
//               = PI * lambda * delta_z * ( k_sqr );
//            // TODO: would it be more efficient to simply multiply a double
//            //       by itself than to call pow( , 2) ?
//            propagator_re = cos( trig_operand );
//            propagator_im = -sin( trig_operand );
//
//            // complex multiplication:
//            //  (a + ib)(c + id) = ac - bd + i(ad + bc)
//            ft_t_psi_re = ft_t_psi[j + i * Ny][0];
//            ft_t_psi_im = ft_t_psi[j + i * Ny][1];
//
//            ft_t_psi[j + i * Ny][0] 
//               = ( propagator_re * ( ft_t_psi_re ) 
//                   - propagator_im * ( ft_t_psi_im ) );
//                   //- propagator_im * ( ft_t_psi_im ) ) / sqrtNxNy;
//
//            ft_t_psi[j + i * Ny][1] 
//               = ( propagator_re * ft_t_psi_im  
//                   + propagator_im * ft_t_psi_re );
//                   //+ propagator_im * ft_t_psi_re ) / sqrtNxNy;
//         }
//         else
//         {
//            ft_t_psi[j + i * Ny][0] = 0.0;
//            ft_t_psi[j + i * Ny][1] = 0.0;
//         }
//      }
//
//   return EXIT_SUCCESS;
//}

//int TEM_NS::propagate_ms(
//      const double& lambda,
//      const double* const kxdomain,
//      const ptrdiff_t& Nx,
//      const double* const kydomain,
//      const ptrdiff_t& Ny,
//      //const double& sqrtNxNy,
//      const double& bwcutoff_t,
//      const double& delta_z,
//      double* propagator_x_re,
//      double* propagator_x_im,
//      double* propagator_y_re,
//      double* propagator_y_im,
//      fftw_complex* ft_t_psi  // FT[t_{n}(x,y) \psi_{n}(x,y)]
//      )
//{
//   // Evaluate 
//   //    P_{n}(k_{x},k_{y},\Delta z_{n}) * FT[t_{n}(x,y) \psi_{n}(x,y)] 
//   // where
//   // P_{n}(k_{x},k_{y},\Delta z_{n})
//   //    = \exp[-i \pi \lambda (k^{2}_{x} + k^{2}_{y}) \Delta z_{n}]
//   double k_sqr;
//   double trig_operand;
//   double propagator_re;//, *propagator_x_re, *propagator_y_re;
//   double propagator_im;//, *propagator_x_im, *propagator_y_im;
//   double ft_t_psi_re;
//   double ft_t_psi_im;
//   const double bwcutoff_t_sqr = pow(bwcutoff_t, 2);
//   // Kirkland (2010) page 146:
//   // "A convenient way of limiting both functions is to set both
//   // the speciment transmission function and the propagator function
//   // to zero outside of (2/3)k_{max}." k_max = min{ kx_max/2, ky_max/2 }
//   // 
//   //  Condition for eliminating aliasing:
//   //  bwcutoff_t + bwcutoff_propagator + bwcutoff_psi < 2 * k_max,
//   //
//   //propagator_x_re = new double[Nx];
//   //propagator_x_im = new double[Nx];
//   //propagator_y_re = new double[Ny];
//   //propagator_y_im = new double[Ny];
//
//   for ( ptrdiff_t i=0; i<Nx; ++i)  
//   {
//      trig_operand = PI * lambda * delta_z * pow( kxdomain[i], 2) ;
//      propagator_x_re[i] = cos( trig_operand );
//      propagator_x_im[i] = -sin( trig_operand );
//   }
//
//   for ( ptrdiff_t j=0; j<Ny; ++j)
//   {
//      trig_operand = PI * lambda * delta_z * pow( kydomain[j], 2) ;
//      propagator_y_re[j] = cos( trig_operand );
//      propagator_y_im[j] = -sin( trig_operand );
//   }
//
//   for ( ptrdiff_t i=0; i<Nx; ++i)  
//      for ( ptrdiff_t j=0; j<Ny; ++j)
//      {
//         k_sqr =  pow(kxdomain[i], 2) + pow(kydomain[j], 2);
//         if ( k_sqr < bwcutoff_t_sqr ) 
//         {
//            // complex multiplication:
//            //  (a + ib)(c + id) = ac - bd + i(ad + bc)
//            // P_{n} = P_{nx} * P_{ny} 
//            propagator_re 
//               = propagator_x_re[i] * propagator_y_re[j]
//                  - propagator_x_im[i] * propagator_y_im[j];
//            propagator_im = -sin( trig_operand );
//
//            propagator_im
//               = propagator_x_re[i] * propagator_y_im[j]
//                  + propagator_x_im[i] * propagator_y_re[j];
//
//            ft_t_psi_re = ft_t_psi[j + i * Ny][0];
//            ft_t_psi_im = ft_t_psi[j + i * Ny][1];
//
//            ft_t_psi[j + i * Ny][0] 
//               = ( propagator_re * ( ft_t_psi_re ) 
//                   - propagator_im * ( ft_t_psi_im ) );
//                   //- propagator_im * ( ft_t_psi_im ) ) / sqrtNxNy;
//
//            ft_t_psi[j + i * Ny][1] 
//               = ( propagator_re * ft_t_psi_im  
//                   + propagator_im * ft_t_psi_re );
//                   //+ propagator_im * ft_t_psi_re ) / sqrtNxNy;
//         }
//         else
//         {
//            ft_t_psi[j + i * Ny][0] = 0.0;
//            ft_t_psi[j + i * Ny][1] = 0.0;
//         }
//      }
//   // delete[] propagator_x_re; 
//   // delete[] propagator_x_im;
//   // delete[] propagator_y_re;
//   // delete[] propagator_y_im;
//
//
//   // The following version works, but calculates sin() & cos() Nx*Ny times,
//   //  whereas the above version calculates sin() & cos() Nx+Ny times at the
//   //  expense of Nx*Ny more complex multiplications.
//   // However, the above version spends time allocating memory which should
//   //  probably be allocated by the calling function and reused.
//   //for ( ptrdiff_t i=0; i<Nx; ++i)  
//   //   for ( ptrdiff_t j=0; j<Ny; ++j)
//   //   {
//   //      k_sqr =  pow(kxdomain[i], 2) + pow(kydomain[j], 2);
//   //      if ( k_sqr < bwcutoff_t_sqr ) 
//   //      {
//   //         trig_operand
//   //            = PI * lambda * delta_z * ( k_sqr );
//   //         // TODO: would it be more efficient to simply multiply a double
//   //         //       by itself than to call pow( , 2) ?
//   //         propagator_re = cos( trig_operand );
//   //         propagator_im = -sin( trig_operand );
//
//   //         // complex multiplication:
//   //         //  (a + ib)(c + id) = ac - bd + i(ad + bc)
//   //         ft_t_psi_re = ft_t_psi[j + i * Ny][0];
//   //         ft_t_psi_im = ft_t_psi[j + i * Ny][1];
//
//   //         ft_t_psi[j + i * Ny][0] 
//   //            = ( propagator_re * ( ft_t_psi_re ) 
//   //                - propagator_im * ( ft_t_psi_im ) );
//   //                //- propagator_im * ( ft_t_psi_im ) ) / sqrtNxNy;
//
//   //         ft_t_psi[j + i * Ny][1] 
//   //            = ( propagator_re * ft_t_psi_im  
//   //                + propagator_im * ft_t_psi_re );
//   //                //+ propagator_im * ft_t_psi_re ) / sqrtNxNy;
//   //      }
//   //      else
//   //      {
//   //         ft_t_psi[j + i * Ny][0] = 0.0;
//   //         ft_t_psi[j + i * Ny][1] = 0.0;
//   //      }
//   //   }
//
//   return EXIT_SUCCESS;
//}

int TEM_NS::slice::update_propagator(
      const double& lambda,
      const double* const kxdomain,
      const ptrdiff_t& Nx,
      const double* const kydomain,
      const ptrdiff_t& Ny,
      const double& bwcutoff_t //,
      //const int& mynode,
      //const int& rootnode,
      //MPI_Comm comm
      )
{
   // Propagation :
   //    P_{n}(k_{x},k_{y},\Delta z_{n}) * FT[t_{n}(x,y) \psi_{n}(x,y)] 
   // where
   // P_{n}(k_{x},k_{y},\Delta z_{n})
   //    = \exp[-i \pi \lambda (k^{2}_{x} + k^{2}_{y}) \Delta z_{n}]
   double k_sqr;
   double trig_operand;
   const double bwcutoff_t_sqr = bwcutoff_t * bwcutoff_t;
   // Kirkland (2010) page 146:
   // "A convenient way of limiting both functions is to set
   //  both the specimen transmission function and the 
   //  propagator function to zero outside of (2/3)k_{max}."
   // k_max = min{ kx_max/2, ky_max/2 }
   // 
   //  Condition for eliminating aliasing:
   //  bwcutoff_t + bwcutoff_propagator + bwcutoff_psi 
   //    < 2 * k_max

   // TODO: bandwidth limit the propagator

   //      k_sqr =  pow(kxdomain[i], 2) + pow(kydomain[j], 2);
   //      if ( k_sqr < bwcutoff_t_sqr ) 
   
   for ( ptrdiff_t i=0; i<Nx; ++i)  
   {  
      if ( kxdomain[i] < bwcutoff_t ) 
      // computing over this domain is split across nodes
      {
         trig_operand 
            = PI * lambda * thickness * pow(kxdomain[i], 2);

         propagator_x_re[i] = cos( trig_operand );
         propagator_x_im[i] = -sin( trig_operand );
      }
      else
      {
         propagator_x_re[i] = 0.0;
         propagator_x_im[i] = 0.0;
      }
   }

   for ( ptrdiff_t j=0; j<Ny; ++j)
   {  
      // NOTE: computing this domain is duplicated accross nodes
      // TODO: maybe split this computation across nodes and broadcast it?
      if ( kydomain[j] < bwcutoff_t ) 
      {
         trig_operand = PI * lambda * thickness * pow( kydomain[j], 2) ;
         propagator_y_re[j] = cos( trig_operand );
         propagator_y_im[j] = -sin( trig_operand );
      }
      else
      {
         propagator_y_re[j] = 0.0;
         propagator_y_im[j] = 0.0;
      }
   }

   return EXIT_SUCCESS;
}


int TEM_NS::slice::propagate( 
      // This function requires separate evaluation of the 
      //  propagator.
      const double& lambda,
      const double* const kxdomain,
      const ptrdiff_t& Nx,
      const double* const kydomain,
      const ptrdiff_t& Ny,
      const double& bwcutoff_t, // transmission function bw
      fftw_complex* ft_t_psi  // FT[t_{n}(x,y) \psi_{n}(x,y)]
      )
{
   double k_sqr;
   double bwcutoff_t_sqr = bwcutoff_t * bwcutoff_t;
   double propagator_re;
   double propagator_im;
   double ft_t_psi_re;
   double ft_t_psi_im;

   for ( ptrdiff_t i=0; i<Nx; ++i)  
      for ( ptrdiff_t j=0; j<Ny; ++j)
      {
         k_sqr =  pow(kxdomain[i], 2) + pow(kydomain[j], 2);
         if ( k_sqr < bwcutoff_t_sqr ) 
         {
            // complex multiplication:
            //  (a + ib)(c + id) = ac - bd + i(ad + bc)
            // P_{n} = P_{nx} * P_{ny} 
            propagator_re 
               = propagator_x_re[i] * propagator_y_re[j]
                  - propagator_x_im[i] * propagator_y_im[j];

            propagator_im
               = propagator_x_re[i] * propagator_y_im[j]
                  + propagator_x_im[i] * propagator_y_re[j];

            ft_t_psi_re = ft_t_psi[j + i * Ny][0];
            ft_t_psi_im = ft_t_psi[j + i * Ny][1];

            ft_t_psi[j + i * Ny][0] 
               = ( propagator_re * ft_t_psi_re
                   - propagator_im * ft_t_psi_im );

            ft_t_psi[j + i * Ny][1] 
               = ( propagator_re * ft_t_psi_im  
                   + propagator_im * ft_t_psi_re );
         }
         else
         {
            ft_t_psi[j + i * Ny][0] = 0.0;
            ft_t_psi[j + i * Ny][1] = 0.0;
         }
      }

   return EXIT_SUCCESS;
}

//int TEM_NS::transmission_function_ms(
//      const ptrdiff_t& Nx,
//      const ptrdiff_t& Ny,
//      const double& sqrtNxNy,
//      const fftw_complex* const pap,
//      fftw_complex* exp_i_sigma_v
//      )
//{
//   double exp_neg_pap_im;
//
//   for ( ptrdiff_t i=0; i<Nx; ++i)
//      for ( ptrdiff_t j=0; j<Ny; ++j)
//      {
//         exp_neg_pap_im
//            = exp(-(pap[j + i * Ny][1]) / sqrtNxNy);
//
//         exp_i_sigma_v[j + i * Ny][0]
//            = exp_neg_pap_im 
//               * cos( (pap[j + i * Ny][0]) / sqrtNxNy );
//
//         exp_i_sigma_v[j + i * Ny][1]
//            = exp_neg_pap_im 
//               * sin( (pap[j + i * Ny][0]) / sqrtNxNy );
//      }
//
//   return EXIT_SUCCESS;
//}

int TEM_NS::slice::update_transmission_function(
            const double& lambda,
            const double& gamma,
            const double& ab_inv,
            //const double& bwcutoff_pap,
            const double& bwcutoff_t,
            //const scatterer_pap_LUT& myScattererPapLUT,
            const double* const kx_split,
            const ptrdiff_t& Nx_split,
            const double* const kx_joined,
            const double* const xx_joined,
            const ptrdiff_t& Nx_joined,
            const double* const ky, // reciprocal space y-domain
            const double* const yy,  // real space y-domain
            const ptrdiff_t& Ny,
            const double& sqrtNxNy,
            const unsigned int& input_flag_pap_tif,
            const string& outFileName,//debug
            const size_t& sliceNumber,//debug
            const ptrdiff_t& local_alloc_size_fftw,
            const ptrdiff_t& local_0_start_fftw,
            const int& mynode,
            const int& rootnode,
            MPI_Comm comm
      )
{
   // Summary : evaluate the transmission function of this slice given the
   //  list of scatterers by translating their pap_to_translate_{re,im} 
   //  members.

   // total projected atomic potential : exp_i_sigma_v[{0,1}]

   /////////////////////////////////////////////////////////////
   // Create forward and reverse fftw plans to enable bandwidth limiting
   /////////////////////////////////////////////////////////////
   fftw_plan pf_c2c_psi, pb_c2c_psi;
   // c2c in-place forward fft
   pf_c2c_psi = fftw_mpi_plan_dft_2d( // c2c in-place fft,
                           Nx_joined, Ny, 
                           exp_i_sigma_v, exp_i_sigma_v,
                           comm, FFTW_FORWARD, FFTW_MEASURE );


   // c2c in-place reverse fft
   pb_c2c_psi = fftw_mpi_plan_dft_2d( 
                           Nx_joined, Ny, 
                           exp_i_sigma_v, exp_i_sigma_v,
                           comm, FFTW_BACKWARD, FFTW_MEASURE );


   /////////////////////////////////////////////////////////////
   // Initialize the slice's transmission function to 0.0
   /////////////////////////////////////////////////////////////
   for ( size_t i=0; i < local_alloc_size_fftw; ++i)
   {
      exp_i_sigma_v[i][0] = 0.0;
      exp_i_sigma_v[i][1] = 0.0;
   }

   // debug
   //cout << "node " << mynode << ", initial exp_i_sigma_v[50][0,1] : " 
   //   << exp_i_sigma_v[50][0] << " + i" << exp_i_sigma_v[50][1]
   //   << endl;
   // end debug

   /////////////////////////////////////////////////////////////
   // For each scatterer, translate and add their projected atomic 
   //  potential to the total.
   /////////////////////////////////////////////////////////////

   // Identify the point in the joined domain at which to begin copying
   //  to the split domain (translate the sliced pap by selecting a 
   //  sub-region of the whole).
   
   size_t idx_x, idx_y;
   size_t idx_shift_x, idx_shift_y;
   size_t idx = 0;

   // debug
   //  Determine the value of idx_local_start_x such that
   //   kx_joined[i + idx_local_start_x] == kx_split[i]
   //   and thus
   //   pap_joined_x[i + idx_local_start_x] == pap_split_x[i]
   //size_t idx_local_start_x; 
   //while ( idx < Nx_joined )
   //{
   //   // NOTE: comparing reciprocal space domains to determine real space
   //   //  splitting location ...
   //   if ( kx_joined[idx] == kx_split[0] )
   //   {
   //      idx_local_start_x = idx;
   //      break;
   //   }
   //   else 
   //      ++idx;
   //}
   //if ( idx == Nx_joined )
   //{
   //   cout << "Error : update_transmission_function() failed;" 
   //      << " could not identify appropriate idx_local_start" << endl;
   //   return EXIT_FAILURE;
   //}
   //if ( idx_local_start_x != local_0_start_fftw )
   //{
   //   cout << "Error : idx_local_start_x != local_0_start_fftw"
   //      << endl << " idx_local_start_x: " << idx_local_start_x << endl
   //      << endl << " local_0_start_fftw: " << local_0_start_fftw << endl
   //      << " Shifting of cached projected atomic potentials and probes"
   //      << " will be erroneous." << endl;
   //   return( EXIT_FAILURE );
   //}
   // end debug

   //cout << "node " << mynode << ", idx_local_start_x : "  // debug
   //   << idx_local_start_x << endl; // debug

   double half_delta_x = 0.5 * (xx_joined[1] - xx_joined[0]);
   double half_delta_y = 0.5 * (yy[1] - yy[0]);

   // iterate over the scatterers
   for ( std::vector< const scatterer* >::const_iterator
          scatterer_ptr_itr = scatterers->begin();
          scatterer_ptr_itr != scatterers->end();
          ++scatterer_ptr_itr )
   {
      // TODO : figure out how to shift by (*scatterer_ptr_itr)->q[0]
      //         and (*scatterer_ptr_itr)->q[1]

      // Identify apropriate idx_shift_x idx_shift_y using kx_joined and 
      //  ky.

      // pap_shifted_y[i] = pap_template_y[i - idx_shift_y] with pbc
      // pap_shifted_split_x[i] 
      //    = pap_template_x[i + idx_local_start_x - idx_shift_x] with pbc
      
      // NOTE: the following assumes scatterer position is within 
      //       xx_joined & y domains
      idx_shift_x = 0;
      idx_shift_y = 0;
      for ( size_t i=0; i < Nx_joined; ++i)
      {
         if ( xx_joined[i] - (*scatterer_ptr_itr)->q[0] >= half_delta_x )
         {
            idx_shift_x = i-1;
            break;
         }
      }

      for ( size_t j=0; j < Ny; ++j)
      {
         if ( yy[j] - (*scatterer_ptr_itr)->q[1] >= half_delta_y )
         {
            idx_shift_y = j-1;
            break;
         }
      }

      //cout << "node " << mynode  // debug
      //   << ", scatterer Z : " << (*scatterer_ptr_itr)->Z // debug
      //   << ", scatterer pap_re : "  // debug
      //   << (*scatterer_ptr_itr)->pap_to_translate_re // debug
      //   << ", scatterer pap_re[0] : "  // debug
      //   << (*scatterer_ptr_itr)->pap_to_translate_re[0] // debug
      //   << ", scatterer pap_re[Nx_joined -1] : "  // debug
      //   << (*scatterer_ptr_itr)->pap_to_translate_re[Nx_joined -1] // debug
      //   << ", scatterer pap_im : "  // debug
      //   << (*scatterer_ptr_itr)->pap_to_translate_im // debug
      //   << ", scatterer pap_im[0] : "  // debug
      //   << (*scatterer_ptr_itr)->pap_to_translate_im[0] // debug
      //   << ", scatterer pap_im[Ny-1] : "  // debug
      //   << (*scatterer_ptr_itr)->pap_to_translate_im[Ny-1] // debug
      //   << ", idx_shift_x, idx_shift_y : " // debug
      //   << idx_shift_x << ", " << idx_shift_y << endl; // debug


      // iterate over the domains
      for ( size_t i=0; i < Nx_split; ++i)
      {
         //   pap_joined_x[i + idx_local_start_x] == pap_split_x[i]
         
         //if ( i + idx_local_start_x >= Nx_joined ) 
         //   idx_local_start_x = -i;

         // enforce periodic boundary condition, x direction
         //if ( i + idx_local_start_x < idx_shift_x )
         if ( i + local_0_start_fftw < idx_shift_x )
            idx_x = i + Nx_joined + local_0_start_fftw - idx_shift_x;
            //idx_x = i + Nx_joined + idx_local_start_x - idx_shift_x;
         else
            idx_x = i + local_0_start_fftw - idx_shift_x;
            //idx_x = i + idx_local_start_x - idx_shift_x;

         for ( size_t j=0; j < Ny; ++j)
         {
            // enforce periodic boundary condition, y direction 
            if ( j < idx_shift_y ) 
               idx_y = j + Ny - idx_shift_y;
            else
               idx_y = j - idx_shift_y;

            // debug
            if (
                  idx_x < 0 || idx_x > Nx_joined
                  || idx_y < 0 || idx_y > Ny
                  || idx_y + idx_x * Ny > Nx_joined * Ny
               )
               cout << "Error, (idx_x, idx_y, idx_y + idx_x*Ny, local_alloc_size_fftw): ("
                 << idx_x << ", " << idx_y 
                 << ", " << idx_y + idx_x * Ny << ", " 
                 << local_alloc_size_fftw
                 << ") exceeds bounds (Nx_joined,Ny)"
                 << Nx_joined << ", " << Ny << ")" << endl;
            // end debug

            // debug
            //   cout << "node " << mynode 
            //      << "; copying pap(" << idx_x << ", " << idx_y
            //      << ") :" 
            //      << (*scatterer_ptr_itr)
            //         ->pap_to_translate_re[ idx_y + idx_x * Ny ]
            //      << " + i" 
            //      << (*scatterer_ptr_itr)
            //         ->pap_to_translate_im[ idx_y + idx_x * Ny ]
            //      << endl;
            // end debug
            // accumulate the projected atomic potential onto the 
            //  transmission function
            exp_i_sigma_v[j + i * Ny][0] += 
               (*scatterer_ptr_itr)
                  ->pap_to_translate_re[ idx_y + idx_x * Ny ];

            exp_i_sigma_v[j + i * Ny][1] += 
               (*scatterer_ptr_itr)
                ->pap_to_translate_im[ idx_y + idx_x * Ny ];
         }
      }
      // debug
      //cout << "node " << mynode << ", exp_i_sigma_v[50][0,1] : " 
      //   << exp_i_sigma_v[50][0] << " + i" << exp_i_sigma_v[50][1]
      //   << endl;
      // end debug
   }

   /////////////////////////////////////////////////////////////
   // Write the projected atomic potential to tif if requested
   /////////////////////////////////////////////////////////////
   if ( input_flag_pap_tif )
   {
      double xResolution, yResolution; 
      if ( mynode == rootnode )
      {
         xResolution = 
            Nx_joined / (
                  (xx_joined[Nx_joined - 1] 
                   + (xx_joined[1] - xx_joined[0])
                   - xx_joined[0] ) * 1e-8
                 );
         yResolution = 
            Ny / ( (yy[Ny - 1] + (yy[1] - yy[0]) - yy[0] ) * 1e-8);
      }
      MPI_Bcast(&xResolution, 1, MPI_DOUBLE, rootnode, comm);
      MPI_Bcast(&yResolution, 1, MPI_DOUBLE, rootnode, comm);

      debug_output_complex_fftw_operand(
            exp_i_sigma_v, 
            1, 
            local_alloc_size_fftw,
            Nx_split, Nx_joined, Ny,
            3, // resolutionUnit = 3 indicates centimeters
            xResolution, yResolution,
            outFileName + "_pap_slice" 
               + TEM_NS::to_string(sliceNumber),
            mynode, rootnode, comm
            );
   }


   /////////////////////////////////////////////////////////////
   // Evaluate the transmission function in realspace given the 
   //  projected atomic potential 
   /////////////////////////////////////////////////////////////
   double exp_neg_pap_im;

   for ( ptrdiff_t i=0; i<Nx_split; ++i)
      for ( ptrdiff_t j=0; j<Ny; ++j)
      {
         exp_neg_pap_im
            = exp( -(exp_i_sigma_v[j + i * Ny][1])  / sqrtNxNy);

         exp_i_sigma_v[j + i * Ny][0]
            = exp_neg_pap_im 
               * cos( (exp_i_sigma_v[j + i * Ny][0]) / sqrtNxNy);

         exp_i_sigma_v[j + i * Ny][1]
            = exp_neg_pap_im 
               * sin( (exp_i_sigma_v[j + i * Ny][0]) / sqrtNxNy );
      }

   /////////////////////////////////////////////////////////////
   // Bandwidth limit the transmission function in reciprocal space
   /////////////////////////////////////////////////////////////
   fftw_execute( pf_c2c_psi );

   bw_limit(
         exp_i_sigma_v,
         Nx_split, kx_split,
         Ny, ky,
         bwcutoff_t
         );

   fftw_execute( pb_c2c_psi );

   // Remove factor of sqrt(Nx*Ny) caused by fftw
   for ( ptrdiff_t i=0; i<local_alloc_size_fftw; ++i)
   {
      exp_i_sigma_v[i][0] = exp_i_sigma_v[i][0] / sqrtNxNy;
      exp_i_sigma_v[i][1] = exp_i_sigma_v[i][1] / sqrtNxNy;
   }

   /////////////////////////////////////////////////////////////
   // clean up
   /////////////////////////////////////////////////////////////
   fftw_destroy_plan( pf_c2c_psi );
   fftw_destroy_plan( pb_c2c_psi );

   return EXIT_SUCCESS;
}
#endif
