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
// File: tem_io.cpp
// Purpose:
//

#ifndef TEM_IO_CPP
#define TEM_IO_CPP

#include "tem_io.hpp"

using std::getline;
using std::istringstream;
using std::setw;
using std::setprecision;

int TEM_NS::read_position_lammps_file( 
      const string& filename, 
      // If the upper and lower boundaries in the file are equal to
      // the periodic bounds, then Nx, Ny, Nz are not needed. 
      // However, if the boundaries are the greatest and least 
      // values contained in the discretized domain, then Nx, Ny, Nz
      // might be needed so that an extra delta_x = (xhi - xlo)/Nx 
      // may be added.
      double*& qq_contig,
      unsigned int*& Z_contig,
      unsigned int& total_population,
      double& xlo, double& ylo, double& zlo, 
      double& xperiod, double& yperiod, double& zperiod,
      const unsigned int& input_flag_debug,
      const int& mynode,
      const int& rootnode,
      MPI_Comm comm
      )
{
   // Precondition: 
   // - input is expected to be the same format as a lammps position
   //  input file, with an Atoms section and an enumerated list of 
   //  atoms types and coordinates
   // - qq_contig and Z_contig have not been allocated.
   
   // Postcondition:
   // - qq_contig and Z_contig have been allocated, and its memory 
   //    must be deleted elsewhere
   // - qq_contig[] and Z_contig[] all reside in contiguous regions 
   //    of memory for MPI
   
   if ( mynode == rootnode )
   {
      // declared_population 'atoms'
      size_t declared_population = 0;

      ifstream data_file( filename.c_str() );
      if ( data_file.is_open() )
      {

         //char data_line[256];// 256 is an arbitrary line length limit
         string data_line;
         getline(data_file, data_line); // first line
         
         // First line must begin with a space separated list of atomic
         //    numbers.
         size_t Z; 
         std::vector<size_t> Zcategories;
         if ( data_file.good() ) 
         {
            istringstream data_line_stream( data_line );
            while( data_line_stream >> Z ) Zcategories.push_back( Z );
         }

         if( getline( data_file, data_line) && data_file.good() )
         {

            // skip lines containing only whitespace
            size_t first = data_line.find_first_not_of(" \t" );
            while( first == std::string::npos ) 
                                       // npos : max size of a string
            {
               getline( data_file, data_line);
               first = data_line.find_first_not_of(" \t" );
            }

            istringstream data_line_stream( data_line );
            data_line_stream >> declared_population;

            string data_descriptor; 
            data_line_stream >> data_descriptor;
            // TODO: why can't I use std::tolower, but I can use tolower?
            transform( data_descriptor.begin(), data_descriptor.end(), 
                       data_descriptor.begin(), (int(*)(int))tolower );
            // TODO: this will only work with ASCII input, not UTF. 
            //       For Unicode support, use toLower from the ICU 
            //       library.
            
            if ( (! declared_population) 
                  || ( data_descriptor.compare("atoms") ) )
            {
               cerr << "node " << mynode << ", "
                  << "Error reading position data file; " 
                  //<< "declared_population, data_descriptor : "//debug
                  //<< declared_population  // debug
                  //<< ", " << data_descriptor  // debug
                  << endl;
               cerr << " current line should be : "<< endl;
               cerr << "<size_t> atoms " << endl;
               return EXIT_FAILURE;
            }
            if ( input_flag_debug )
               cout << "node " << mynode << ", "
                  << "declared_population : "
                  << declared_population << endl;
         }
         else
         {
            cerr << "node " << mynode << ", "
               << "Error reading position data file; data_line : "
               << data_line << endl;
            cerr << "data_file.good() : " << data_file.good();
            return EXIT_FAILURE;
         }

         // numberofspecies 'atom types'
         size_t numberofspecies = 0;
         if( getline( data_file, data_line) && data_file.good() )
         {
            // skip lines containing only whitespace
            size_t first = data_line.find_first_not_of(" \t" );
            while( first == std::string::npos ) // npos : max size of a string
            {
               getline( data_file, data_line);
               first = data_line.find_first_not_of(" \t" );
            }

            istringstream data_line_stream( data_line );
            data_line_stream >> numberofspecies;

            string data_descriptor1, data_descriptor2; 
            data_line_stream >> data_descriptor1 >> data_descriptor2;
            if ( data_descriptor1.empty() || data_descriptor2.empty() )
            {
               cerr << "node " << mynode << ", "
                  << "Error reading position data file; "
                  << "text following number of species must be 'atom types'" 
                  << endl;
            }
            // transform characters to lower case
            // TODO: why can't I use std::tolower, but I can use 
            //       tolower?
            // TODO: the following will only work with ASCII input, not
            //       UTF. 
            //       For Unicode support, use toLower from the ICU 
            //       library.
            transform( data_descriptor1.begin(), 
                        data_descriptor1.end(), 
                           data_descriptor1.begin(), 
                           (int(*)(int))tolower );
            transform( data_descriptor2.begin(), 
                        data_descriptor2.end(), 
                           data_descriptor2.begin(), 
                           (int(*)(int))tolower );
            
            if (  (! numberofspecies )
                    || ( data_descriptor1.compare("atom") ) 
                    || ( data_descriptor2.compare("types") ) 
               )
            {
               cerr << "node " << mynode << ", "
                  << "Error reading position data file; " 
                  << "numberofspecies, data_descriptor1, " // debug
                  << " data_descriptor2 : " // debug
                  << numberofspecies  // debug
                  << ", " << data_descriptor1  // debug
                  << ", " << data_descriptor2  // debug
                  << endl;
               cerr << " current line should be : "<< endl;
               cerr << "<double> atom types"
                  << endl;
               return EXIT_FAILURE;
            }
            if ( input_flag_debug )
               cout << "node " << mynode << ", "
                  << "numberofspecies : "
                  << numberofspecies << endl;

            if ( numberofspecies != Zcategories.size() )
            {
               cerr << "node " << mynode << ", "   
                  << "The first line of the input file must begin with "
                  << "a space separated list of atomic numbers having order "
                 << " corresponding with the lammps species enumeration"
                 << endl;
               return EXIT_FAILURE;
            }
         }
         else
         {
            cerr << "node " << mynode << ", "   
               << "Error reading position data file; data_line : "
               << data_line << endl;
            cerr << "data_file.good() : " << data_file.good();
            return EXIT_FAILURE;
         }

         // double double 'xlo xhi'
         double xlower, xupper;
         if( getline( data_file, data_line) && data_file.good() )
         {
            // skip lines containing only whitespace
            size_t first = data_line.find_first_not_of(" \t" );
            while( first == std::string::npos ) // npos : max size of a string
            {
               getline( data_file, data_line);
               first = data_line.find_first_not_of(" \t" );
            }

            string data_descriptor1, data_descriptor2; 
            istringstream data_line_stream( data_line );
            data_line_stream >> xlower >> xupper 
               >> data_descriptor1 >> data_descriptor2;

            if ( data_descriptor1.empty() || data_descriptor2.empty() )
            {
               cerr << "Error reading position data file;" 
                  << " current line should be : "<< endl;
               cerr << "<double> <double> xlo xhi"
                  << endl;
            }
            // transform characters to lower case
            // TODO: why can't I use std::tolower, but I can use tolower?
            // TODO: the following will only work with ASCII input, not UTF. 
            //       For Unicode support, use toLower from the ICU library.
            transform( data_descriptor1.begin(), data_descriptor1.end(), 
                           data_descriptor1.begin(), (int(*)(int))tolower );
                           //data_descriptor1.begin(), tolower );
            transform( data_descriptor2.begin(), data_descriptor2.end(), 
                           data_descriptor2.begin(), (int(*)(int))tolower );
                           //data_descriptor2.begin(), tolower );
            
            if (
                  ( data_descriptor1.compare("xlo") ) 
                  || ( data_descriptor2.compare("xhi") ) 
               )
            {
               cerr << "node " << mynode << ", "   // debug
                  << "Error reading position data file; " // debug
                  << "xlower, xupper, data_descriptor1, " // debug
                  << " data_descriptor2 : " // debug
                  << xlower << ", " << xupper // debug
                  << ", " << data_descriptor1 // debug
                  << ", " << data_descriptor2 // debug
                  << endl;
               return EXIT_FAILURE;
            }
            if ( input_flag_debug )
               cout << "node " << mynode << ", "
                  << "xlower, xupper : "
                  << xlower << ", " << xupper << endl;
         }
         else
         {
            cerr << "node " << mynode << ", "
               << "Error reading position data file; data_line : "
               << data_line << endl;
            cerr << "data_file.good() : " << data_file.good();
            return EXIT_FAILURE;
         }

         // double double 'ylo yhi'
         double ylower, yupper;
         if( getline( data_file, data_line) && data_file.good() )
         {
            // skip lines containing only whitespace
            size_t first = data_line.find_first_not_of(" \t" );
            while( first == std::string::npos ) // npos : max size of a string
            {
               getline( data_file, data_line);
               first = data_line.find_first_not_of(" \t" );
            }

            string data_descriptor1, data_descriptor2; 
            istringstream data_line_stream( data_line );
            data_line_stream >> ylower >> yupper 
               >> data_descriptor1 >> data_descriptor2;

            if ( data_descriptor1.empty() || data_descriptor2.empty() )
            {
               cerr << "node " << mynode << ", "
                  << "Error reading position data file;" 
                  << " current line should be : "<< endl;
               cerr << "<double> <double> ylo yhi"
                  << endl;
               return EXIT_FAILURE;
            }
            // transform characters to lower case
            // TODO: why can't I use std::tolower, but I can use tolower?
            // TODO: the following will only work with ASCII input, not UTF
            //       For Unicode support, use toLower from the ICU library.
            transform( data_descriptor1.begin(), data_descriptor1.end(), 
                           data_descriptor1.begin(), (int(*)(int))tolower);
                           //data_descriptor1.begin(), tolower );
            transform( data_descriptor2.begin(), data_descriptor2.end(), 
                           data_descriptor2.begin(), (int(*)(int))tolower);
                           //data_descriptor2.begin(), tolower );
            
            if (  
                    ( data_descriptor1.compare("ylo") ) 
                    || ( data_descriptor2.compare("yhi") ) 
               )
            {
               cerr << "node " << mynode << ", "
                  << "Error reading position data file; " 
                  << "ylower, yupper, data_descriptor1, "// debug
                  << " data_descriptor2 : "// debug
                  << ylower << ", " << yupper // debug
                  << ", " << data_descriptor1 // debug
                  << ", " << data_descriptor2 // debug
                  << endl;
               return EXIT_FAILURE;
            }
            if ( input_flag_debug )
               cout << "node " << mynode << ", "
                  << "ylower, yupper : "
                  << ylower << ", " << yupper << endl;
         }
         else
         {
            cerr << "node " << mynode << ", "
               << "Error reading position data file; data_line : "
               << data_line << endl;
            cerr << "data_file.good() : " << data_file.good();
            return EXIT_FAILURE;
         }

         // double double 'zlo zhi'
         double zlower, zupper;
         if( getline( data_file, data_line) && data_file.good() )
         {
            // skip lines containing only whitespace
            size_t first = data_line.find_first_not_of(" \t" );
            while( first == std::string::npos ) // npos : max size of a string
            {
               getline( data_file, data_line);
               first = data_line.find_first_not_of(" \t" );
            }

            string data_descriptor1, data_descriptor2; 
            istringstream data_line_stream( data_line );
            data_line_stream >> zlower >> zupper 
               >> data_descriptor1 >> data_descriptor2;

            if ( data_descriptor1.empty() || data_descriptor2.empty() )
            {
               cerr << "node " << mynode << ", "
                  << "Error reading position data file;" 
                  << " current line should be : "<< endl;
               cerr << "<double> <double> zlo zhi"
                  << endl;
            }
            // transform characters to lower case
            // TODO: why can't I use std::tolower, but I can use tolower?
            // TODO: the following will only work with ASCII input, not UTF. 
            //       For Unicode support, use toLower from the ICU library.
            transform( data_descriptor1.begin(), data_descriptor1.end(), 
                           data_descriptor1.begin(), (int(*)(int))tolower );
                           //data_descriptor1.begin(), tolower );
            transform( data_descriptor2.begin(), data_descriptor2.end(), 
                           data_descriptor2.begin(), (int(*)(int))tolower );
                           //data_descriptor2.begin(), tolower );
            
            if (  
                    ( data_descriptor1.compare("zlo") ) 
                    || ( data_descriptor2.compare("zhi") ) 
               )
            {
               cerr << "node " << mynode << ", "
                  << "Error reading position data file; " 
                  << "zlower, zupper, data_descriptor1, " // debug
                  << " data_descriptor2 : " // debug
                  << zlower << ", " << zupper  // debug
                  << ", " << data_descriptor1  // debug
                  << ", " << data_descriptor2  // debug
                  << endl; 
               return EXIT_FAILURE;
            }
            if ( input_flag_debug )
               cout << "node " << mynode << ", "
                  << "zlower, zupper : "
                  << zlower << ", " << zupper << endl;
         }
         else
         {
            cerr << "node " << mynode << ", "
               << "Error reading position data file; data_line : "
               << data_line << endl;
            cerr << "data_file.good() : " << data_file.good();
            return EXIT_FAILURE;
         }

         // double double double xy xz yz // tilt factors for triclinic system
         double tiltxy, tiltxz, tiltyz;
         if( getline( data_file, data_line) && data_file.good() )
         {
            // skip lines containing only whitespace
            size_t first = data_line.find_first_not_of(" \t" );
            while( first == std::string::npos ) // npos : max size of a string
            {
               getline( data_file, data_line);
               first = data_line.find_first_not_of(" \t" );
            }

            string data_descriptor1, data_descriptor2, data_descriptor3; 
            istringstream data_line_stream( data_line );
            data_line_stream >> tiltxy >> tiltxz >> tiltyz
               >> data_descriptor1 >> data_descriptor2 >> data_descriptor3;

            if ( data_descriptor1.empty() 
                  || data_descriptor2.empty() 
                  || data_descriptor3.empty() 
                  )
            {
               cerr << "node " << mynode << ", "
                  << "Error reading position data file;" 
                  << " current line should be : "<< endl
                  << "<double> <double> <double> xy xz yz"
                  << endl;
            }
            // transform characters to lower case
            // TODO: why can't I use std::tolower, but I can use tolower?
            // TODO: the following will only work with ASCII input, not UTF. 
            //       For Unicode support, use toLower from the ICU library.
            transform( data_descriptor1.begin(), data_descriptor1.end(), 
                           data_descriptor1.begin(), (int(*)(int))tolower );
                           //data_descriptor1.begin(), tolower );
            transform( data_descriptor2.begin(), data_descriptor2.end(), 
                           data_descriptor2.begin(), (int(*)(int))tolower );
                           //data_descriptor2.begin(), tolower );
            transform( data_descriptor3.begin(), data_descriptor3.end(), 
                           data_descriptor3.begin(), (int(*)(int))tolower );
                           //data_descriptor3.begin(), tolower );
            
            if (  
                    ( data_descriptor1.compare("xy") ) 
                    || ( data_descriptor2.compare("xz") ) 
                    || ( data_descriptor3.compare("yz") ) 
               )
            {
               cerr << "node " << mynode << ", "
                     << "Error reading position data file; " 
                     << "tiltxy, tiltxz, tiltyz, data_descriptor1, " // debug
                     << " data_descriptor2, data_descriptor3 : " // debug
                  << tiltxy << ", " << tiltxz << ", " << tiltyz << ", "//debug
                     << data_descriptor1 // debug
                     << ", " << data_descriptor2 // debug
                     << ", " << data_descriptor3 // debug
                     << endl
                     << " current line should be : "<< endl
                     << "<double> <double> <double> xy xz yz"
                     << endl
                     << "data_descriptor1.compar('xy') : " 
                     <<  data_descriptor1.compare("xy")  
                     << endl
                     << "data_descriptor2.compar('xz') : " 
                     <<  data_descriptor2.compare("xz")  
                     << endl
                     << "data_descriptor3.compar('yz') : " 
                     <<  data_descriptor3.compare("yz")  
                     << endl;
               return EXIT_FAILURE;
            }
            if ( input_flag_debug )
               cout << "node " << mynode << ", "
                     << "tiltxy, tiltxz, tiltyz: "
                     << tiltxy
                     << ", " << tiltxz
                     << ", " << tiltyz << endl;
         }
         else
         {
            cerr << "node " << mynode << ", "
               << "Error reading position data file; data_line : "
               << data_line << endl;
            cerr << "data_file.good() : " << data_file.good();
            return EXIT_FAILURE;
         }

         // check that boundary values are valid
         if ( xupper <= xlower ) 
         {
            cerr << "node " << mynode << ", "
               << "Error reading position data file; xhi <= xlo" << endl;
            return EXIT_FAILURE;
         }
         if ( yupper <= ylower ) 
         {
            cerr << "node " << mynode << ", "
               << "Error reading position data file; yhi <= ylo" << endl;
            return EXIT_FAILURE;
         }
         if ( zupper <= zlower ) 
         {
            cerr << "node " << mynode << ", "
               << "Error reading position data file; zhi <= zlo" << endl;
            return EXIT_FAILURE;
         }
         if ( tiltxy != 0.0 || tiltxz != 0.0 || tiltyz != 0.0 )
         {
            cerr << "node " << mynode << ", "
               << "Error reading position data file; xy, xz, yz, != 0.0" 
               << endl
               << "Cannot handle triclinic boundaries at the moment"
               << endl;
            return EXIT_FAILURE;
         }

         // Assign boundaries to output variables
         xperiod = xupper - xlower;
         yperiod = yupper - ylower;
         zperiod = zupper - zlower;
         xlo = xlower;
         ylo = ylower;
         zlo = zlower;
         
         if ( mynode == rootnode && input_flag_debug )
            cout << " From input file: (xperiod, xlower, xupper): ("
               << xlower << ", " << xupper << ")"
               << "(yperiod, ylower, yupper): ("
               << ylower << ", " << yupper << ")"
               << endl;

         // 'Atoms'
         
         if( getline( data_file, data_line) && data_file.good() )
         {
            // skip lines containing only whitespace
            size_t first = data_line.find_first_not_of(" \t" );
            // std::string::npos : max size of a string
            while( first == std::string::npos ) 
            {
               getline( data_file, data_line);
               first = data_line.find_first_not_of(" \t" );
            }

            string section_name; 
            istringstream data_line_stream( data_line );
            data_line_stream >> section_name;

            if ( section_name.empty() )
            {
               cerr << "node " << mynode << ", "
                  << "Error reading position data file;" 
                  << " current line should be : "<< endl
                  << "Atoms"
                  << endl;
            }
            // transform characters to lower case
            // TODO: why can't I use std::tolower, but I can use tolower?
            // TODO: the following will only work with ASCII input, not UTF
            //       For Unicode support, use toLower from the ICU library.
            transform( section_name.begin(), section_name.end(), 
                           section_name.begin(), (int(*)(int))tolower );
            if ( section_name.compare("atoms") )
            {
               cerr << "node " << mynode << ", "
                  << "Error reading position data file; " 
                  << " section_name : " // debug
                  << section_name //debug
                  << endl
                  << " current line should be : "<< endl
                  << "Atoms"
                  << endl;
               return EXIT_FAILURE;
            }
         }
         else
         {
            cerr << "node " << mynode << ", "
               << "Error reading position data file; data_line : "
               << data_line << endl
               << "data_file.good() : " << data_file.good();
            return EXIT_FAILURE;
         }


         // read types and positions of individual atoms

         size_t atom_number, atom_type; 
         size_t number_of_read_atoms = 0;
         double* qq;
         std::list<size_t> atom_number_list;
         std::list<size_t> atom_type_list;
         //std::vector<size_t> Zlist;
         std::vector<double*> qlist;

         // The following assignments will be overwritten since positions
         //  are ensured to be lower than (xupper, yupper, zupper).
         double xmin, ymin, zmin;
         xmin = xupper;
         ymin = yupper;
         zmin = zupper;

         while ( number_of_read_atoms < declared_population )
         {
            if( getline( data_file, data_line) && data_file.good() )
            {
               // skip lines containing only whitespace
               size_t first = data_line.find_first_not_of(" \t" );
               while( first == std::string::npos)//npos : max size of a string
               {
                  getline( data_file, data_line);
                  first = data_line.find_first_not_of(" \t" );
               } 

               qq = new double[3];
               // grab data from the current line
               istringstream data_line_stream( data_line );
               data_line_stream >> atom_number >> atom_type;
               data_line_stream >> qq[0] >> qq[1] >> qq[2] ;

               // check that the data is valid
               if ( atom_number > declared_population
                     || atom_type > numberofspecies
                     || atom_type <= 0
                     || qq[0] > xupper || qq[0] < xlower
                     || qq[1] > yupper || qq[1] < ylower
                     || qq[2] > zupper || qq[2] < zlower
                  )
               {
                  cerr << "node " << mynode << ", "
                     << "Error reading position data file;" 
                     << " on line with atom number : " << atom_number << endl
                     << " Check that atom number, type, and position are "
                     << "within appropriate boundaries." << endl;
                  return EXIT_FAILURE;
                  // TODO: clean allocated lists before returning failure
               }
               
               // Determine minimum positions so that all they may be 
               //  shifted to be all positive.
               if ( qq[0] < xmin ) xmin = qq[0];
               if ( qq[1] < ymin ) ymin = qq[1];
               if ( qq[2] < zmin ) zmin = qq[2];

               // TODO: assign data to qlist, atom_number_list
               qlist.push_back(qq);
               atom_number_list.push_back(atom_number);
               atom_type_list.push_back(atom_type );
               
               number_of_read_atoms++;
            }
            else
            {
               cerr << "node " << mynode << ", "
                  << "Error reading position data file;" 
                  << " atoms read : " 
                  << number_of_read_atoms  
                  << ", declared population : "
                  << declared_population
                  << endl
                  << "data_file.eof() : " << data_file.eof() << endl
                  << "data_file.good() : " << data_file.good() << endl
                  << "data_file.fail() : " << data_file.fail() << endl
                  << "data_file.bad() : " << data_file.bad() << endl;
               return EXIT_FAILURE;
            }
         }
         
         // Ensure that atom types (1 2 ...) match numberofspecies, and 
         //  that all atoms have unique indices (atom_number_list).

         std::list<size_t> atom_type_list_uniqued( atom_type_list );
         atom_type_list_uniqued.sort();
         atom_type_list_uniqued.unique();

         std::list<size_t> atom_number_list_uniqued( atom_number_list );
         atom_number_list_uniqued.sort();
         atom_number_list_uniqued.unique();

         ////////////////////////////////////////////////////////////////
         // shift all atoms so that all coordinates are positive
         ////////////////////////////////////////////////////////////////
         for ( size_t i=0; i < qlist.size() ; ++i )
         {
            qlist[i][0] = qlist[i][0] - xmin;
            qlist[i][1] = qlist[i][1] - ymin;
            qlist[i][2] = qlist[i][2] - zmin;
         }

         // Sequential pointers in qlist[] are contiguous, but the qq[] 
         //  arrays they point to are not contiguous to each other.
         // Here we transfer values from the qq[] arrays to a contiguous 
         //  array so that MPI may properly broadcast them between nodes 
         //  with ease.
         
         //std::vector<double*> qlist_contig;   // not needed
         size_t common_size = qlist.size();


         if ( 
               declared_population == common_size 
               && atom_type_list_uniqued.size() == numberofspecies 
               && atom_number_list_uniqued.size() == common_size 
            )
         {
            total_population = common_size;

            std::list<size_t>::iterator 
               atom_type_list_iterator = atom_type_list.begin();

            qq_contig = new double[ 3 * common_size ]; // 3-D
            Z_contig = new unsigned int[ total_population ];

            for (unsigned int i=0; i < total_population; i++)
            {
               qq_contig[ 3 * i ]       = qlist[i][0];
               qq_contig[ 3 * i + 1 ]   = qlist[i][1];
               qq_contig[ 3 * i + 2 ]   = qlist[i][2];

               delete[] qlist[i];// TODO: is this appropriate?
               qlist[i] = NULL;  // TODO: is this necessary?

               Z_contig[ i ] = Zcategories[ (*atom_type_list_iterator) -1];
               atom_type_list_iterator++;
            }
         }
         else
         {
            cerr << "node " << mynode << ", "
               << "Error reading position file: " << filename << endl
               << " declared_population : " << declared_population << endl
               << " positions read : " << qlist.size() << endl
               << " unique atom ID numbers : " 
               << atom_number_list_uniqued.size() << endl;
            for ( size_t i=0; i<qlist.size(); i++) 
            {
               delete[] qlist[i];
               qlist.pop_back();
            }
            return EXIT_FAILURE;
         }
//         ////////////////////////////////////////////////////////////////
//         // debug
//         for ( size_t i=0 ; i < Zcategories.size(); i++ )
//            cout << "Zcategories[" << i << "] : " << Zcategories[i] << endl;
//
//         cout << "declared_population : " << declared_population <<  endl;
//         cout << "atom types : " << numberofspecies << endl;
//         cout << "xlo xhi : " << xlower << " " << xupper << endl;
//         cout << "ylo yhi : " << ylower << " " << yupper << endl;
//         cout << "zlo zhi : " << zlower << " " << zupper << endl;
//         cout << "tiltxy tiltxz tiltyz : " 
//            << tiltxy << " " << tiltxz << " " << tiltyz << endl;
//         cout << "xperiod, yperiod, zperiod : " 
//            << xperiod << ", " << yperiod << ", " << zperiod << endl;
//         cout << "Atoms" << endl;
//
//
//         std::list<size_t>::iterator 
//            atom_number_list_iterator = atom_number_list.begin();
//         std::list<size_t>::iterator 
//            atom_type_list_iterator = atom_type_list.begin();
//         //std::vector<scatterer>::iterator 
//         //   myScatterers_iterator = myScatterers.begin();
//         cout << "Atom number, atom type, atomic number, position q" << endl;
//         int precision = 7;
//         int width = precision + 2;
//         for( size_t i=0; i< total_population; i++)
//         {
//            cout 
//               << std::setw(width) << std::setprecision(precision) 
//               << *atom_number_list_iterator 
//               << std::setw(width) << std::setprecision(precision) 
//               << *atom_type_list_iterator 
//               << std::setw(width) << std::setprecision(precision) 
//               << Z_contig[ i ]
//               << std::setw(width) << std::setprecision(precision) 
//               << qq_contig[ 3 * i ] 
//               << std::setw(width) << std::setprecision(precision) 
//               << qq_contig[ 3 * i + 1]
//               << std::setw(width) << std::setprecision(precision) 
//               << qq_contig[ 3 * i + 2]
//               << endl;
//               atom_type_list_iterator++;
//               atom_number_list_iterator++;
//         }
         // end debug
         //////////////////////////////////////////////////////////////////


         //////////////////////////////////////////////////////////////////
         // clean up
         for ( size_t i = 0; i < qlist.size(); i++)
         {
            //delete[] qlist.back();   // taken care of in previous loops
            //                      
            qlist.pop_back();
            // qq to be dealt with above after copying to qq_contig
            // qlist_contig.pop_back();
            // qq_contig to be deleted by calling function
         }
         for ( size_t i = 0; i < atom_type_list.size(); i++)
            atom_type_list.pop_back();
         for ( size_t i = 0; i < atom_number_list.size(); i++) 
            atom_number_list.pop_back(); 
         // end clean up
         //////////////////////////////////////////////////////////////////
      }
      else 
      {
         cerr << "node " << mynode << ", "
            << "Error opening position file: " << filename << endl;
         data_file.close();
         return EXIT_FAILURE;
      }

      data_file.close();
   } // end of the rootnode block
 
   // Broadcast results to the remaining nodes
   MPI_Bcast( &total_population, 1, MPI_UNSIGNED, rootnode, comm);

   if ( mynode != rootnode )
   {
      qq_contig = new double[ 3 * total_population ];
      Z_contig = new unsigned int[ total_population ];
   }

   MPI_Bcast( qq_contig, 3 * total_population, MPI_DOUBLE, rootnode, comm);

   MPI_Bcast( Z_contig, total_population, MPI_UNSIGNED, rootnode, comm);

   return EXIT_SUCCESS;
}





int TEM_NS::create_position_lammps_file_011diamond(
      const string& filename,
      const double& xmin, const double& ymin, const double& zmin,
      const size_t& unit_cells_x, 
      const size_t& unit_cells_y, 
      const size_t& unit_cells_z, 
      const int& Nx, const int& Ny, const int& Nz,
      const size_t& Z,
      const double& lattice_constant
      )
{
   std::ofstream data_file; 
   data_file.open( filename.c_str(), 
         std::ofstream::out  | std::ofstream::app );
   if ( ! data_file.good() )
   {
      cerr << "Error, create_position_lammps_file_011diamond : "
         << "could not open file '" << filename << "' for writing" << endl;
      data_file.close();
      return EXIT_FAILURE;
   }

   if ( data_file.is_open() )
   {
      if ( Nx <= 0 || Ny <= 0 || Nz <= 0
            || unit_cells_x <= 0 || unit_cells_y <= 0 || unit_cells_z <= 0
            || Z <= 0
            || lattice_constant <= 0
            )
      {
         cerr << "Error, create_position_lammps_file_011diamond : "
            << "N{x,y,z}, unit_cells_{x,y,z}, Z, and lattice_constant must"
            << " be greater than zero" << endl;
         return EXIT_FAILURE;
      }
      //size_t atoms_per_unit_cell = 8;// diamond viewed along 100
      size_t atoms_per_unit_cell = 4;// diamond viewed along 011

      const double oneoversqrttwo = 1.0/sqrt(2.0);
      double lattice_const_x = lattice_constant;
      double lattice_const_y = lattice_constant * oneoversqrttwo;
      double lattice_const_z = lattice_constant * oneoversqrttwo;

      double xperiod = unit_cells_x * lattice_const_x;
      double yperiod = unit_cells_y * lattice_const_y;
      double zperiod = unit_cells_z * lattice_const_z;

      unsigned int total_population = atoms_per_unit_cell
                              * unit_cells_x * unit_cells_y * unit_cells_z;

      double* q; q = new double[ 3 * total_population];
      //size_t * Z; Z = new size_t[ total_population ];

      double* x; x = new double[ Nx ];
      double* y; y = new double[ Ny ];
      double* z; z = new double[ Nz ];

      domain_3D_periodic( Nx, Ny, Nz,
                           xperiod, xmin,
                           yperiod, ymin,
                           zperiod, zmin,
                           x, y, z);

      positions_3D_lattice_diamond_011(
            unit_cells_x, unit_cells_y, unit_cells_z,
            x, Nx, y, Ny, z, Nz,
            q );

      ///////////////////////////////////////////////////////////////////
      // write to the file
      ///////////////////////////////////////////////////////////////////
      int precision = 7;
      int width = precision + 2;
      data_file 
         <<  setw(width) << setprecision(precision)
         << Z << " Z list" << endl << endl;

      data_file 
         <<  setw(width) << setprecision(precision)
         << total_population 
         <<  setw(width) << setprecision(precision)
         << "atoms" << endl << endl;
      data_file 
         <<  setw(width) << setprecision(precision)
         << "1" 
         <<  setw(width) << setprecision(precision)
         << " atom types" << endl << endl;
      data_file 
         <<  setw(width) << setprecision(precision)
         << xmin 
         <<  setw(width) << setprecision(precision)
         << xperiod + xmin 
         <<  setw(width) << setprecision(precision)
         << " xlo xhi" << endl;
      data_file 
         <<  setw(width) << setprecision(precision)
         << ymin 
         <<  setw(width) << setprecision(precision)
         << yperiod + ymin 
         <<  setw(width) << setprecision(precision)
         << " ylo yhi" << endl;
      data_file 
         <<  setw(width) << setprecision(precision)
         << zmin 
         <<  setw(width) << setprecision(precision)
         << zperiod + zmin 
         <<  setw(width) << setprecision(precision)
         << " zlo zhi" << endl;
      data_file 
         <<  setw(width) << setprecision(precision)
         << 0.0 
         <<  setw(width) << setprecision(precision)
         << 0.0 
         <<  setw(width) << setprecision(precision)
         << 0.0 
         <<  setw(width) << setprecision(precision)
         << " xy xz yz" << endl << endl;
      data_file << " Atoms" << endl << endl;

      for ( unsigned int i=0; i < total_population; i++)
      {
         data_file 
            <<  setw(width) << setprecision(precision)
            << i+1 
            <<  setw(width) << setprecision(precision)
            << 1 
            <<  setw(width) << setprecision(precision)
            << q[3 * i] 
            <<  setw(width) << setprecision(precision)
            << q[3 * i + 1] 
            <<  setw(width) << setprecision(precision)
            << q[3 * i + 2] << endl;
      }

      ///////////////////////////////////////////////////////////////////
      // clean up
      ///////////////////////////////////////////////////////////////////

      delete[] q;
      //delete[] Z;
      delete[] x;
      delete[] y;
      delete[] z;
   }
   else 
   {
      cerr << "Error opening position file: " << filename << endl;
      data_file.close();
      return EXIT_FAILURE;
   }

   data_file.close();

   return EXIT_SUCCESS;
}

int TEM_NS::create_position_lammps_file_001diamond(
      const string& filename,
      const double& xmin, const double& ymin, const double& zmin,
      const size_t& unit_cells_x, 
      const size_t& unit_cells_y, 
      const size_t& unit_cells_z, 
      const int& Nx, const int& Ny, const int& Nz,
      const size_t& Z,
      const double& lattice_constant
      )
{
   std::ofstream data_file; 
   data_file.open( filename.c_str(), 
         std::ofstream::out  | std::ofstream::app );
   if ( ! data_file.good() )
   {
      cerr << "Error, create_position_lammps_file_011diamond : "
         << "could not open file '" << filename << "' for writing" << endl;
      data_file.close();
      return EXIT_FAILURE;
   }

   if ( data_file.is_open() )
   {
      if ( Nx <= 0 || Ny <= 0 || Nz <= 0
            || unit_cells_x <= 0 || unit_cells_y <= 0 || unit_cells_z <= 0
            || Z <= 0
            || lattice_constant <= 0
            )
      {
         cerr << "Error, create_position_lammps_file_100diamond : "
            << "N{x,y,z}, unit_cells_{x,y,z}, Z, and lattice_constant must"
            << " be greater than zero" << endl;
         return EXIT_FAILURE;
      }
      size_t atoms_per_unit_cell = 8;// diamond viewed along 001
      //size_t atoms_per_unit_cell = 4;// diamond viewed along 011

      //const double oneoversqrttwo = 1.0/sqrt(2.0);
      //double lattice_const_x = lattice_constant;
      //double lattice_const_y = lattice_constant;
      //double lattice_const_z = lattice_constant;// * oneoversqrttwo;

      //double xperiod = unit_cells_x * lattice_const_x;
      //double yperiod = unit_cells_y * lattice_const_y;
      //double zperiod = unit_cells_z * lattice_const_z;
      double xperiod = unit_cells_x * lattice_constant;
      double yperiod = unit_cells_y * lattice_constant;
      double zperiod = unit_cells_z * lattice_constant;

      unsigned int total_population = atoms_per_unit_cell
                              * unit_cells_x * unit_cells_y * unit_cells_z;

      double* q; q = new double[ 3 * total_population];
      //size_t * Z; Z = new size_t[ total_population ];

      double* x; x = new double[ Nx ];
      double* y; y = new double[ Ny ];
      double* z; z = new double[ Nz ];

      domain_3D_periodic( Nx, Ny, Nz,
                           xperiod, xmin,
                           yperiod, ymin,
                           zperiod, zmin,
                           x, y, z);

      positions_3D_lattice_diamond_001(
            unit_cells_x, unit_cells_y, unit_cells_z,
            x, Nx, y, Ny, z, Nz,
            q );

      ///////////////////////////////////////////////////////////////////
      // write to the file
      ///////////////////////////////////////////////////////////////////
      int precision = 7;
      int width = precision + 2;
      data_file 
         <<  setw(width) << setprecision(precision)
         << Z << " Z list" << endl << endl;

      data_file 
         <<  setw(width) << setprecision(precision)
         << total_population 
         <<  setw(width) << setprecision(precision)
         << "atoms" << endl << endl;
      data_file 
         <<  setw(width) << setprecision(precision)
         << "1" 
         <<  setw(width) << setprecision(precision)
         << " atom types" << endl << endl;
      data_file 
         <<  setw(width) << setprecision(precision)
         << xmin 
         <<  setw(width) << setprecision(precision)
         << xperiod + xmin 
         <<  setw(width) << setprecision(precision)
         << " xlo xhi" << endl;
      data_file 
         <<  setw(width) << setprecision(precision)
         << ymin 
         <<  setw(width) << setprecision(precision)
         << yperiod + ymin 
         <<  setw(width) << setprecision(precision)
         << " ylo yhi" << endl;
      data_file 
         <<  setw(width) << setprecision(precision)
         << zmin 
         <<  setw(width) << setprecision(precision)
         << zperiod + zmin 
         <<  setw(width) << setprecision(precision)
         << " zlo zhi" << endl;
      data_file 
         <<  setw(width) << setprecision(precision)
         << 0.0 
         <<  setw(width) << setprecision(precision)
         << 0.0 
         <<  setw(width) << setprecision(precision)
         << 0.0 
         <<  setw(width) << setprecision(precision)
         << " xy xz yz" << endl << endl;
      data_file << " Atoms" << endl << endl;

      for ( unsigned int i=0; i < total_population; i++)
      {
         data_file 
            <<  setw(width) << setprecision(precision)
            << i+1 
            <<  setw(width) << setprecision(precision)
            << 1 
            <<  setw(width) << setprecision(precision)
            << q[ 3 * i] 
            <<  setw(width) << setprecision(precision)
            << q[ 3 * i + 1] 
            <<  setw(width) << setprecision(precision)
            << q[ 3 * i + 2] << endl;
      }

      ///////////////////////////////////////////////////////////////////
      // clean up
      ///////////////////////////////////////////////////////////////////

      delete[] q;
      //delete[] Z;
      delete[] x;
      delete[] y;
      delete[] z;
   }
   else 
   {
      cerr << "Error opening position file: " << filename << endl;
      data_file.close();
      return EXIT_FAILURE;
   }

   data_file.close();

   return EXIT_SUCCESS;
}

int TEM_NS::create_position_lammps_file_111diamond(
      const string& filename,
      const double& xmin, const double& ymin, const double& zmin,
      const size_t& unit_cells_x, 
      const size_t& unit_cells_y, 
      const size_t& unit_cells_z, 
      const int& Nx, const int& Ny, const int& Nz,
      const size_t& Z,
      const double& lattice_constant
      )
{
   // Initially treat it as though it's a 001 oriented crystal, using
   //  the positions_3D_lattice_diamond_111() function which will create
   //  the positions as though it were 001 oriented, but will also rotate
   //  the positions into a 111 orientation while enforcing periodic 
   //  boundary conditions.
   std::ofstream data_file; 
   data_file.open( filename.c_str(), 
         std::ofstream::out  | std::ofstream::app );
   if ( ! data_file.good() )
   {
      cerr << "Error, create_position_lammps_file_011diamond : "
         << "could not open file '" << filename << "' for writing" << endl;
      data_file.close();
      return EXIT_FAILURE;
   }

   if ( data_file.is_open() )
   {
      if ( Nx <= 0 || Ny <= 0 || Nz <= 0
            || unit_cells_x <= 0 || unit_cells_y <= 0 || unit_cells_z <= 0
            || Z <= 0
            || lattice_constant <= 0
            )
      {
         cerr << "Error, create_position_lammps_file_100diamond : "
            << "N{x,y,z}, unit_cells_{x,y,z}, Z, and lattice_constant must"
            << " be greater than zero" << endl;
         return EXIT_FAILURE;
      }
      size_t atoms_per_unit_cell = 8;// diamond viewed along 001
      //size_t atoms_per_unit_cell = 4;// diamond viewed along 011

      //const double oneoversqrttwo = 1.0/sqrt(2.0);
      //double lattice_const_x = lattice_constant;
      //double lattice_const_y = lattice_constant;
      //double lattice_const_z = lattice_constant;// * oneoversqrttwo;

      //double xperiod = unit_cells_x * lattice_const_x;
      //double yperiod = unit_cells_y * lattice_const_y;
      //double zperiod = unit_cells_z * lattice_const_z;
      //double xperiod = unit_cells_x * lattice_constant;
      //double yperiod = unit_cells_y * lattice_constant;
      //double zperiod = unit_cells_z * lattice_constant;

      //size_t total_population = atoms_per_unit_cell
      //                        * unit_cells_x * unit_cells_y * unit_cells_z;

      //double* q; q = new double[ 3 * total_population];
      std::list<double> q;
      std::list<double>::iterator q_itr;

      //size_t * Z; Z = new size_t[ total_population ];

      //double* x; x = new double[ Nx ];
      //double* y; y = new double[ Ny ];
      //double* z; z = new double[ Nz ];

      //domain_3D_periodic( Nx, Ny, Nz,
      //                     xperiod, xmin,
      //                     yperiod, ymin,
      //                     zperiod, zmin,
      //                     x, y, z);

      positions_3D_lattice_diamond_111(
            unit_cells_x, unit_cells_y, unit_cells_z,
            lattice_constant,
            //x, Nx, y, Ny, z, Nz,
            q );
      unsigned int total_population = q.size() / 3;

      // Determine periodic boundaries
      // x lattice constant is 6* the third atom's x-coordinate
      // y lattice constant is 4* the third atom's y-coordinate
      // z lattice constant is 3* the third atom's z-coordinate
      q_itr = q.begin();
      ++q_itr ; ++q_itr; ++q_itr;// first atom
      ++q_itr ; ++q_itr; ++q_itr;// 2nd atom
      // x component of the 3rd atom
      double xperiod = 6 * (*q_itr) * unit_cells_x;
      ++q_itr;// y component of the 3rd atom
      double yperiod = 4 * (*q_itr) * unit_cells_y;
      ++q_itr;// z component of the 3rd atom
      double zperiod = 3 * (*q_itr) * unit_cells_z;

      //double xperiod = unit_cells_x * lattice_constant_wrt111;
      //double yperiod = unit_cells_y * lattice_constant_wrt111;
      //double zperiod = unit_cells_z * lattice_constant_wrt111;

      // Shift atom positions to move [0,0,0] to [xmin,ymin,zmin]
      for ( q_itr = q.begin(); q_itr != q.end(); ++q_itr)
      {
         (*q_itr) = (*q_itr) + xmin;
         ++q_itr;
         (*q_itr) = (*q_itr) + ymin;
         ++q_itr;
         (*q_itr) = (*q_itr) + zmin;
         ++q_itr;
      }
      ///////////////////////////////////////////////////////////////////
      // write to the file
      ///////////////////////////////////////////////////////////////////
      int precision = 7;
      int width = precision + 6;
      data_file 
         <<  setw(width) << setprecision(precision)
         << Z << " Z list" << endl << endl;

      data_file 
         <<  setw(width) << setprecision(precision)
         << total_population 
         <<  setw(width) << setprecision(precision)
         << "atoms" << endl << endl;
      data_file 
         <<  setw(width) << setprecision(precision)
         << "1" 
         <<  setw(width) << setprecision(precision)
         << " atom types" << endl << endl;
      data_file 
         <<  setw(width) << setprecision(precision)
         << xmin 
         <<  setw(width) << setprecision(precision)
         << xperiod + xmin 
         <<  setw(width) << setprecision(precision)
         << " xlo xhi" << endl;
      data_file 
         <<  setw(width) << setprecision(precision)
         << ymin 
         <<  setw(width) << setprecision(precision)
         << yperiod + ymin 
         <<  setw(width) << setprecision(precision)
         << " ylo yhi" << endl;
      data_file 
         <<  setw(width) << setprecision(precision)
         << zmin 
         <<  setw(width) << setprecision(precision)
         << zperiod + zmin 
         <<  setw(width) << setprecision(precision)
         << " zlo zhi" << endl;
      data_file 
         <<  setw(width) << setprecision(precision)
         << 0.0 
         <<  setw(width) << setprecision(precision)
         << 0.0 
         <<  setw(width) << setprecision(precision)
         << 0.0 
         <<  setw(width) << setprecision(precision)
         << " xy xz yz" << endl << endl;
      data_file << " Atoms" << endl << endl;

      size_t atom_number=1;
      for ( q_itr = q.begin(); q_itr != q.end(); ++q_itr)
      {
         data_file 
            << setw(width) << setprecision(precision)
            << atom_number
            << setw(width) << setprecision(precision)
            << 1 // atom type
            << setw(width) << setprecision(precision)
            << *q_itr;
         ++q_itr;
         data_file 
            << setw(width) << setprecision(precision)
            << *q_itr;
         ++q_itr;
         data_file 
            << setw(width) << setprecision(precision)
            << *q_itr
            << endl;
         ++atom_number;
      }

      ///////////////////////////////////////////////////////////////////
      // clean up
      ///////////////////////////////////////////////////////////////////

      //delete[] q;
      //delete[] Z;
      //delete[] x;
      //delete[] y;
      //delete[] z;
   }
   else 
   {
      cerr << "Error opening position file: " << filename << endl;
      data_file.close();
      return EXIT_FAILURE;
   }

   data_file.close();

   return EXIT_SUCCESS;
}


int TEM_NS::create_position_lammps_file_111zincblende(
      const string& filename,
      const double& xmin, const double& ymin, const double& zmin,
      const size_t& unit_cells_x, 
      const size_t& unit_cells_y, 
      const size_t& unit_cells_z, 
      const int& Nx, const int& Ny, const int& Nz,
      const size_t& Z1,
      const size_t& Z2,
      const double& lattice_constant
      )
{
   // Initially treat it as though it's a 001 oriented crystal, using
   //  the positions_3D_lattice_diamond_111() function which will create
   //  the positions as though it were 001 oriented, but will also rotate
   //  the positions into a 111 orientation while enforcing periodic 
   //  boundary conditions.
   std::ofstream data_file; 
   data_file.open( filename.c_str(), 
         std::ofstream::out  | std::ofstream::app );
   if ( ! data_file.good() )
   {
      cerr << "Error, create_position_lammps_file_011diamond : "
         << "could not open file '" << filename << "' for writing" << endl;
      data_file.close();
      return EXIT_FAILURE;
   }

   if ( data_file.is_open() )
   {
      if ( Nx <= 0 || Ny <= 0 || Nz <= 0
            || unit_cells_x <= 0 || unit_cells_y <= 0 || unit_cells_z <= 0
            || Z1 <= 0
            || Z2 <= 0
            || lattice_constant <= 0
            )
      {
         cerr << "Error, create_position_lammps_file_100diamond : "
            << "N{x,y,z}, unit_cells_{x,y,z}, Z, and lattice_constant must"
            << " be greater than zero" << endl;
         return EXIT_FAILURE;
      }
      size_t atoms_per_unit_cell = 8;// diamond viewed along 001
      //size_t atoms_per_unit_cell = 4;// diamond viewed along 011

      //const double oneoversqrttwo = 1.0/sqrt(2.0);
      //double lattice_const_x = lattice_constant;
      //double lattice_const_y = lattice_constant;
      //double lattice_const_z = lattice_constant;// * oneoversqrttwo;

      //double xperiod = unit_cells_x * lattice_const_x;
      //double yperiod = unit_cells_y * lattice_const_y;
      //double zperiod = unit_cells_z * lattice_const_z;
      //double xperiod = unit_cells_x * lattice_constant;
      //double yperiod = unit_cells_y * lattice_constant;
      //double zperiod = unit_cells_z * lattice_constant;

      //size_t total_population = atoms_per_unit_cell
      //                        * unit_cells_x * unit_cells_y * unit_cells_z;

      //double* q; q = new double[ 3 * total_population];
      std::list<double> q;
      std::list<double>::iterator q_itr;
      std::list<size_t> Z;
      std::list<size_t>::iterator Z_itr;

      //size_t * Z; Z = new size_t[ total_population ];

      //double* x; x = new double[ Nx ];
      //double* y; y = new double[ Ny ];
      //double* z; z = new double[ Nz ];

      //domain_3D_periodic( Nx, Ny, Nz,
      //                     xperiod, xmin,
      //                     yperiod, ymin,
      //                     zperiod, zmin,
      //                     x, y, z);

      positions_3D_lattice_zincblende_111(
            unit_cells_x, unit_cells_y, unit_cells_z,
            lattice_constant,
            //x, Nx, y, Ny, z, Nz,
            Z1, Z2,
            Z,
            q );
      unsigned int total_population = q.size() / 3;

      // Determine periodic boundaries
      // x lattice constant is 6* the third atom's x-coordinate
      // y lattice constant is 4* the third atom's y-coordinate
      // z lattice constant is 3* the third atom's z-coordinate
      q_itr = q.begin();
      ++q_itr ; ++q_itr; ++q_itr;// first atom
      ++q_itr ; ++q_itr; ++q_itr;// 2nd atom
      // x component of the 3rd atom
      double xperiod = 6 * (*q_itr) * unit_cells_x;
      ++q_itr;// y component of the 3rd atom
      double yperiod = 4 * (*q_itr) * unit_cells_y;
      ++q_itr;// z component of the 3rd atom
      double zperiod = 3 * (*q_itr) * unit_cells_z;

      //double xperiod = unit_cells_x * lattice_constant_wrt111;
      //double yperiod = unit_cells_y * lattice_constant_wrt111;
      //double zperiod = unit_cells_z * lattice_constant_wrt111;

      // Shift atom positions to move [0,0,0] to [xmin,ymin,zmin]
      for ( q_itr = q.begin(); q_itr != q.end(); ++q_itr)
      {
         (*q_itr) = (*q_itr) + xmin;
         ++q_itr;
         (*q_itr) = (*q_itr) + ymin;
         ++q_itr;
         (*q_itr) = (*q_itr) + zmin;
         ++q_itr;
      }
      ///////////////////////////////////////////////////////////////////
      // write to the file
      ///////////////////////////////////////////////////////////////////
      int precision = 7;
      int width = precision + 6;
      data_file 
         <<  setw(width) << setprecision(precision)
         << Z1 << " " << Z2 << " Z list" << endl << endl;

      data_file 
         <<  setw(width) << setprecision(precision)
         << total_population 
         <<  setw(width) << setprecision(precision)
         << "atoms" << endl << endl;
      data_file 
         <<  setw(width) << setprecision(precision)
         << "2" 
         <<  setw(width) << setprecision(precision)
         << " atom types" << endl << endl;
      data_file 
         <<  setw(width) << setprecision(precision)
         << xmin 
         <<  setw(width) << setprecision(precision)
         << xperiod + xmin 
         <<  setw(width) << setprecision(precision)
         << " xlo xhi" << endl;
      data_file 
         <<  setw(width) << setprecision(precision)
         << ymin 
         <<  setw(width) << setprecision(precision)
         << yperiod + ymin 
         <<  setw(width) << setprecision(precision)
         << " ylo yhi" << endl;
      data_file 
         <<  setw(width) << setprecision(precision)
         << zmin 
         <<  setw(width) << setprecision(precision)
         << zperiod + zmin 
         <<  setw(width) << setprecision(precision)
         << " zlo zhi" << endl;
      data_file 
         <<  setw(width) << setprecision(precision)
         << 0.0 
         <<  setw(width) << setprecision(precision)
         << 0.0 
         <<  setw(width) << setprecision(precision)
         << 0.0 
         <<  setw(width) << setprecision(precision)
         << " xy xz yz" << endl << endl;
      data_file << " Atoms" << endl << endl;

      size_t atom_number=1;
      Z_itr = Z.begin();
      for ( q_itr = q.begin(); q_itr != q.end(); ++q_itr)
      {
         data_file 
            << setw(width) << setprecision(precision)
            << atom_number
            << setw(width) << setprecision(precision)
            << *Z_itr // atom type
            << setw(width) << setprecision(precision)
            << *q_itr;
         ++q_itr;
         data_file 
            << setw(width) << setprecision(precision)
            << *q_itr;
         ++q_itr;
         data_file 
            << setw(width) << setprecision(precision)
            << *q_itr
            << endl;
         ++atom_number;
         ++Z_itr;
      }

      ///////////////////////////////////////////////////////////////////
      // clean up
      ///////////////////////////////////////////////////////////////////

      //delete[] q;
      //delete[] Z;
      //delete[] x;
      //delete[] y;
      //delete[] z;
   }
   else 
   {
      cerr << "Error opening position file: " << filename << endl;
      data_file.close();
      return EXIT_FAILURE;
   }

   data_file.close();

   return EXIT_SUCCESS;
}

int TEM_NS::create_position_lammps_file_011zincblende(
      const string& filename,
      const double& xmin, const double& ymin, const double& zmin,
      const size_t& unit_cells_x, 
      const size_t& unit_cells_y, 
      const size_t& unit_cells_z, 
      const int& Nx, const int& Ny, const int& Nz,
      const size_t& Z1,
      const size_t& Z2,
      const double& lattice_constant
      )
{
   std::ofstream data_file; 
   data_file.open( filename.c_str(), 
         std::ofstream::out  | std::ofstream::app );
   if ( ! data_file.good() )
   {
      cerr << "Error, create_position_lammps_file_011diamond : "
         << "could not open file '" << filename << "' for writing" << endl;
      data_file.close();
      return EXIT_FAILURE;
   }

   if ( data_file.is_open() )
   {
      if ( Nx <= 0 || Ny <= 0 || Nz <= 0
            || unit_cells_x <= 0 || unit_cells_y <= 0 || unit_cells_z <= 0
            || Z1 <= 0
            || Z2 <= 0
            || lattice_constant <= 0
            )
      {
         cerr << "Error, create_position_lammps_file_100diamond : "
            << "N{x,y,z}, unit_cells_{x,y,z}, Z, and lattice_constant must"
            << " be greater than zero" << endl;
         return EXIT_FAILURE;
      }
      size_t atoms_per_unit_cell = 8;// diamond viewed along 001
      //size_t atoms_per_unit_cell = 4;// diamond viewed along 011

      //const double oneoversqrttwo = 1.0/sqrt(2.0);
      //double lattice_const_x = lattice_constant;
      //double lattice_const_y = lattice_constant;
      //double lattice_const_z = lattice_constant;// * oneoversqrttwo;

      //double xperiod = unit_cells_x * lattice_const_x;
      //double yperiod = unit_cells_y * lattice_const_y;
      //double zperiod = unit_cells_z * lattice_const_z;
      //double xperiod = unit_cells_x * lattice_constant;
      //double yperiod = unit_cells_y * lattice_constant;
      //double zperiod = unit_cells_z * lattice_constant;

      //size_t total_population = atoms_per_unit_cell
      //                        * unit_cells_x * unit_cells_y * unit_cells_z;

      //double* q; q = new double[ 3 * total_population];
      std::list<double> q;
      std::list<double>::iterator q_itr;
      std::list<size_t> Z;
      std::list<size_t>::iterator Z_itr;

      //size_t * Z; Z = new size_t[ total_population ];

      //double* x; x = new double[ Nx ];
      //double* y; y = new double[ Ny ];
      //double* z; z = new double[ Nz ];

      //domain_3D_periodic( Nx, Ny, Nz,
      //                     xperiod, xmin,
      //                     yperiod, ymin,
      //                     zperiod, zmin,
      //                     x, y, z);

      positions_3D_lattice_zincblende_011(
            unit_cells_x, unit_cells_y, unit_cells_z,
            lattice_constant,
            //x, Nx, y, Ny, z, Nz,
            Z1, Z2,
            Z,
            q );
      unsigned int total_population = q.size() / 3;

      // Determine periodic boundaries
      // x lattice constant is 2* the 2nd atom's x-coordinate
      // y lattice constant is 2* the 2nd atom's y-coordinate
      // z lattice constant is 2* the 3rd atom's z-coordinate
      q_itr = q.begin();
      ++q_itr ; ++q_itr; ++q_itr;// first atom
      // x component of the 2nd atom
      double xperiod = 2 * (*q_itr) * unit_cells_x;
      ++q_itr;// y component of the 2nd atom
      double yperiod = 2 * (*q_itr) * unit_cells_y;
      ++q_itr; ++q_itr ; ++q_itr; ++q_itr; // z component of the 3rd atom
      double zperiod = 2 * (*q_itr) * unit_cells_z;

      //double xperiod = unit_cells_x * lattice_constant_wrt111;
      //double yperiod = unit_cells_y * lattice_constant_wrt111;
      //double zperiod = unit_cells_z * lattice_constant_wrt111;

      // Shift atom positions to move [0,0,0] to [xmin,ymin,zmin]
      for ( q_itr = q.begin(); q_itr != q.end(); ++q_itr)
      {
         (*q_itr) = (*q_itr) + xmin;
         ++q_itr;
         (*q_itr) = (*q_itr) + ymin;
         ++q_itr;
         (*q_itr) = (*q_itr) + zmin;
         ++q_itr;
      }
      ///////////////////////////////////////////////////////////////////
      // write to the file
      ///////////////////////////////////////////////////////////////////
      int precision = 7;
      int width = precision + 6;
      data_file 
         <<  setw(width) << setprecision(precision)
         << Z1 << " " << Z2 << " Z list" << endl << endl;

      data_file 
         <<  setw(width) << setprecision(precision)
         << total_population 
         <<  setw(width) << setprecision(precision)
         << "atoms" << endl << endl;
      data_file 
         <<  setw(width) << setprecision(precision)
         << "2" 
         <<  setw(width) << setprecision(precision)
         << " atom types" << endl << endl;
      data_file 
         <<  setw(width) << setprecision(precision)
         << xmin 
         <<  setw(width) << setprecision(precision)
         << xperiod + xmin 
         <<  setw(width) << setprecision(precision)
         << " xlo xhi" << endl;
      data_file 
         <<  setw(width) << setprecision(precision)
         << ymin 
         <<  setw(width) << setprecision(precision)
         << yperiod + ymin 
         <<  setw(width) << setprecision(precision)
         << " ylo yhi" << endl;
      data_file 
         <<  setw(width) << setprecision(precision)
         << zmin 
         <<  setw(width) << setprecision(precision)
         << zperiod + zmin 
         <<  setw(width) << setprecision(precision)
         << " zlo zhi" << endl;
      data_file 
         <<  setw(width) << setprecision(precision)
         << 0.0 
         <<  setw(width) << setprecision(precision)
         << 0.0 
         <<  setw(width) << setprecision(precision)
         << 0.0 
         <<  setw(width) << setprecision(precision)
         << " xy xz yz" << endl << endl;
      data_file << " Atoms" << endl << endl;

      size_t atom_number=1;
      Z_itr = Z.begin();
      for ( q_itr = q.begin(); q_itr != q.end(); ++q_itr)
      {
         data_file 
            << setw(width) << setprecision(precision)
            << atom_number
            << setw(width) << setprecision(precision)
            << *Z_itr // atom type
            << setw(width) << setprecision(precision)
            << *q_itr;
         ++q_itr;
         data_file 
            << setw(width) << setprecision(precision)
            << *q_itr;
         ++q_itr;
         data_file 
            << setw(width) << setprecision(precision)
            << *q_itr
            << endl;
         ++atom_number;
         ++Z_itr;
      }

      ///////////////////////////////////////////////////////////////////
      // clean up
      ///////////////////////////////////////////////////////////////////

      //delete[] q;
      //delete[] Z;
      //delete[] x;
      //delete[] y;
      //delete[] z;
   }
   else 
   {
      cerr << "Error opening position file: " << filename << endl;
      data_file.close();
      return EXIT_FAILURE;
   }

   data_file.close();

   return EXIT_SUCCESS;
}


int TEM_NS::create_position_lammps_file_001zincblende(
      const string& filename,
      const double& xmin, const double& ymin, const double& zmin,
      const size_t& unit_cells_x, 
      const size_t& unit_cells_y, 
      const size_t& unit_cells_z, 
      const int& Nx, const int& Ny, const int& Nz,
      const size_t& Z1,
      const size_t& Z2,
      const double& lattice_constant
      )
{
   // Initially treat it as though it's a 001 oriented crystal, using
   //  the positions_3D_lattice_diamond_111() function which will create
   //  the positions as though it were 001 oriented, but will also rotate
   //  the positions into a 111 orientation while enforcing periodic 
   //  boundary conditions.
   std::ofstream data_file; 
   data_file.open( filename.c_str(), 
         std::ofstream::out  | std::ofstream::app );
   if ( ! data_file.good() )
   {
      cerr << "Error, create_position_lammps_file_011diamond : "
         << "could not open file '" << filename << "' for writing" << endl;
      data_file.close();
      return EXIT_FAILURE;
   }

   if ( data_file.is_open() )
   {
      if ( Nx <= 0 || Ny <= 0 || Nz <= 0
            || unit_cells_x <= 0 || unit_cells_y <= 0 || unit_cells_z <= 0
            || Z1 <= 0
            || Z2 <= 0
            || lattice_constant <= 0
            )
      {
         cerr << "Error, create_position_lammps_file_100diamond : "
            << "N{x,y,z}, unit_cells_{x,y,z}, Z, and lattice_constant must"
            << " be greater than zero" << endl;
         return EXIT_FAILURE;
      }
      size_t atoms_per_unit_cell = 8;// diamond viewed along 001
      //size_t atoms_per_unit_cell = 4;// diamond viewed along 011

      //const double oneoversqrttwo = 1.0/sqrt(2.0);
      //double lattice_const_x = lattice_constant;
      //double lattice_const_y = lattice_constant;
      //double lattice_const_z = lattice_constant;// * oneoversqrttwo;

      //double xperiod = unit_cells_x * lattice_const_x;
      //double yperiod = unit_cells_y * lattice_const_y;
      //double zperiod = unit_cells_z * lattice_const_z;
      //double xperiod = unit_cells_x * lattice_constant;
      //double yperiod = unit_cells_y * lattice_constant;
      //double zperiod = unit_cells_z * lattice_constant;

      //size_t total_population = atoms_per_unit_cell
      //                        * unit_cells_x * unit_cells_y * unit_cells_z;

      //double* q; q = new double[ 3 * total_population];
      std::list<double> q;
      std::list<double>::iterator q_itr;
      std::list<size_t> Z;
      std::list<size_t>::iterator Z_itr;

      //size_t * Z; Z = new size_t[ total_population ];

      //double* x; x = new double[ Nx ];
      //double* y; y = new double[ Ny ];
      //double* z; z = new double[ Nz ];

      //domain_3D_periodic( Nx, Ny, Nz,
      //                     xperiod, xmin,
      //                     yperiod, ymin,
      //                     zperiod, zmin,
      //                     x, y, z);

      positions_3D_lattice_zincblende_001(
            unit_cells_x, unit_cells_y, unit_cells_z,
            lattice_constant,
            //x, Nx, y, Ny, z, Nz,
            Z1, Z2,
            Z,
            q );
      unsigned int total_population = q.size() / 3;

      // Determine periodic boundaries
      double xperiod = unit_cells_x * lattice_constant;
      double yperiod = unit_cells_y * lattice_constant;
      double zperiod = unit_cells_z * lattice_constant;

      //double xperiod = unit_cells_x * lattice_constant_wrt111;
      //double yperiod = unit_cells_y * lattice_constant_wrt111;
      //double zperiod = unit_cells_z * lattice_constant_wrt111;

      // Shift atom positions to move [0,0,0] to [xmin,ymin,zmin]
      for ( q_itr = q.begin(); q_itr != q.end(); ++q_itr)
      {
         (*q_itr) = (*q_itr) + xmin;
         ++q_itr;
         (*q_itr) = (*q_itr) + ymin;
         ++q_itr;
         (*q_itr) = (*q_itr) + zmin;
         ++q_itr;
      }
      ///////////////////////////////////////////////////////////////////
      // write to the file
      ///////////////////////////////////////////////////////////////////
      int precision = 7;
      int width = precision + 6;
      data_file 
         <<  setw(width) << setprecision(precision)
         << Z1 << " " << Z2 << " Z list" << endl << endl;

      data_file 
         <<  setw(width) << setprecision(precision)
         << total_population 
         <<  setw(width) << setprecision(precision)
         << "atoms" << endl << endl;
      data_file 
         <<  setw(width) << setprecision(precision)
         << "2" 
         <<  setw(width) << setprecision(precision)
         << " atom types" << endl << endl;
      data_file 
         <<  setw(width) << setprecision(precision)
         << xmin 
         <<  setw(width) << setprecision(precision)
         << xperiod + xmin 
         <<  setw(width) << setprecision(precision)
         << " xlo xhi" << endl;
      data_file 
         <<  setw(width) << setprecision(precision)
         << ymin 
         <<  setw(width) << setprecision(precision)
         << yperiod + ymin 
         <<  setw(width) << setprecision(precision)
         << " ylo yhi" << endl;
      data_file 
         <<  setw(width) << setprecision(precision)
         << zmin 
         <<  setw(width) << setprecision(precision)
         << zperiod + zmin 
         <<  setw(width) << setprecision(precision)
         << " zlo zhi" << endl;
      data_file 
         <<  setw(width) << setprecision(precision)
         << 0.0 
         <<  setw(width) << setprecision(precision)
         << 0.0 
         <<  setw(width) << setprecision(precision)
         << 0.0 
         <<  setw(width) << setprecision(precision)
         << " xy xz yz" << endl << endl;
      data_file << " Atoms" << endl << endl;

      size_t atom_number=1;
      Z_itr = Z.begin();
      for ( q_itr = q.begin(); q_itr != q.end(); ++q_itr)
      {
         data_file 
            << setw(width) << setprecision(precision)
            << atom_number
            << setw(width) << setprecision(precision)
            << *Z_itr // atom type
            << setw(width) << setprecision(precision)
            << *q_itr;
         ++q_itr;
         data_file 
            << setw(width) << setprecision(precision)
            << *q_itr;
         ++q_itr;
         data_file 
            << setw(width) << setprecision(precision)
            << *q_itr
            << endl;
         ++atom_number;
         ++Z_itr;
      }

      ///////////////////////////////////////////////////////////////////
      // clean up
      ///////////////////////////////////////////////////////////////////

      //delete[] q;
      //delete[] Z;
      //delete[] x;
      //delete[] y;
      //delete[] z;
   }
   else 
   {
      cerr << "Error opening position file: " << filename << endl;
      data_file.close();
      return EXIT_FAILURE;
   }

   data_file.close();

   return EXIT_SUCCESS;
}


int TEM_NS::read_position_xyz_file(
      const string& filename,
      double*& qq_contig,
      unsigned int*& Z_contig,
      unsigned int& total_population,
      double& xlo, double& ylo, double& zlo,
      double& xperiod, double& yperiod, double& zperiod,
      const unsigned int& input_flag_debug,
      const int& mynode,
      const int& rootnode,
      MPI_Comm comm
      )

{
   // Currently elements 1 through 103 are modeled.

   if ( mynode == rootnode ) 
   {
      ifstream data_file( filename.c_str() );
      if ( data_file.is_open() )
      {
         string data_line;
         //getline(data_file, data_line);

         size_t declared_population;
         

         ////////////////////////////////////////////////////////////////
         // First line of an xyz file must begin with the number of atoms
         ////////////////////////////////////////////////////////////////
         size_t decleared_population = 0;
         if ( getline( data_file, data_line) && data_file.good() )
         {
            istringstream data_line_stream( data_line );
            data_line_stream >> declared_population;
            
            if ( !  declared_population )
            {
               cerr << "node " << mynode << ", "
                  << "Error reading position data file; "
                  << "declared_population : "
                  << declared_population << endl;
               cerr << " current line should be : " << endl;
               cerr << "<number of atoms>  " << endl;
               return EXIT_FAILURE;
            }
            if ( input_flag_debug )
               cout << "node " << mynode << ", "
                  << "declared_population : "
                  << declared_population << endl;
         }
         else
         {
            cerr << "node " << mynode << ", "
               << "Error reading position data file; data_line : "
               << data_line << endl;
            cerr << "data_file.good() : " << data_file.good();
            return EXIT_FAILURE;
         }

         ////////////////////////////////////////////////////////////////
         // The second line should begin with : 
         //    <xperiod> <yperiod> <zperiod>
         ////////////////////////////////////////////////////////////////
         if ( getline( data_file, data_line) && data_file.good() )
         {
            double xperiod_local, yperiod_local, zperiod_local;
            istringstream data_line_stream( data_line );
            data_line_stream >> xperiod_local 
               >> yperiod_local >> zperiod_local;

            if ( xperiod_local > 0.0 
                  && yperiod_local > 0.0 
                  && zperiod_local > 0.0 )
            {
               xperiod = xperiod_local;
               yperiod = yperiod_local;
               zperiod = zperiod_local;
            }
            else
            {
               cerr << "node " << mynode << ", "
                  << "Error reading xyz position data file;"
                  << " second line should begin with : " << endl;
               cerr << "<xperiod> <yperiod> <zperiod>" << endl;
               cerr << "data_line : " << endl << data_line << endl;
               cerr << "xperiod_local, yperiod_local, zperiod_local : " 
                  << xperiod_local << ", "
                  << yperiod_local << ", "
                  << zperiod_local << endl;
               return EXIT_FAILURE;
            }

         }
         else
         {
            cerr << "node " << mynode << ", "
               << "Error reading position data file; data_line : "
               << data_line << endl;
            cerr << "data_file.good() : " << data_file.good();
            return EXIT_FAILURE;
         }

         ////////////////////////////////////////////////////////////////
         // Read declared_population lines of the file to populate 
         //    positions and Zs
         ////////////////////////////////////////////////////////////////
         size_t atom_type;
         size_t number_of_read_atoms = 0;
         double* qq;
         std::vector<size_t> atom_type_list;
         std::vector<double*> qlist;

         double xmin = 0.0; double xmax = 0.0;
         double ymin = 0.0; double ymax = 0.0;
         double zmin = 0.0; double zmax = 0.0;

         if ( input_flag_debug )
            cout << "declared population : " 
               << declared_population << endl;

         while ( number_of_read_atoms < declared_population )
         {
            if ( getline( data_file, data_line) && data_file.good() )
            {
               // TODO : create a function to match element abbreviation with Z
               // TODO : populate qq_contig (atom positions) and Z_contig
               qq = new double[3];
               string element_name;
               istringstream data_line_stream( data_line );

               // grab element abbreviation from the current line
               data_line_stream >> element_name;

               // transform the read element name to lower case
               transform( element_name.begin(), element_name.end(),
                     element_name.begin(), (int(*)(int))tolower );

               // Translate the element name to atomic number Z
               atom_type = atom_element_abbrev_to_Z( element_name );
               if ( atom_type < 1 || atom_type > 103 )
               {
                  cerr << "node " << mynode << ", "
                     << "Error reading position data file; "
                     << "atom " << number_of_read_atoms << ", "
                     << "element " << element_name << " invalid; "
                     << "data_line : " << data_line << endl;
                   // clean up before exiting
                  for ( size_t i=0; i < qlist.size(); i++) 
                     qlist.pop_back();
                  for ( size_t i=0; i < atom_type_list.size(); i++)
                     atom_type_list.pop_back();
                  return EXIT_FAILURE;
               }

               data_line_stream >> qq[0] >> qq[1] >> qq[2] ;

               // keep track of atomic position maximums and minimums for 
               //  comparison to xperiod, yperiod, and zperiod
               if ( qq[0] < xmin ) xmin = qq[0];
               if ( qq[0] > xmax ) xmax = qq[0];
               if ( qq[1] < ymin ) ymin = qq[1];
               if ( qq[1] > ymax ) ymax = qq[1];
               if ( qq[2] < zmin ) zmin = qq[2];
               if ( qq[2] > zmax ) zmax = qq[2];

               qlist.push_back( qq );
               atom_type_list.push_back( atom_type );
               
               number_of_read_atoms++;

            }
            else
            {
               cerr << "node " << mynode << ", "
                  << "Error reading position data file; "
                  << "atoms read : " << number_of_read_atoms
                  << ", declared population : " << declared_population << endl
                  << ", data_file.good() : " << data_file.good()
                  << ", data_file.eof() : " << data_file.eof()
                  << ", data_file.fail() : " << data_file.fail()
                  << ", data_file.bad() : " << data_file.bad();
                // clean up before exiting
               for ( size_t i=0; i < qlist.size(); i++) 
                  qlist.pop_back();
               for ( size_t i=0; i < atom_type_list.size(); i++)
                  atom_type_list.pop_back();
               return EXIT_FAILURE;
            }

         }
         if ( input_flag_debug )
            cout << "number of atoms read: "
               << number_of_read_atoms << endl;
         
         // TODO : check that xmax - xmin <= xperiod, and similary for y & z
         if ( xmax - xmin > xperiod 
               || ymax - ymin > yperiod 
               || zmax - zmin > zperiod )
         {
            cerr << "node " << mynode << ", "
               << "Error reading position data file;"
               << " positions exceed period "
               << endl 
               << "xmax, xmin, xperiod : " 
               << xmax << ", " << xmin << ", " << xperiod << endl
               << "ymax, ymin, yperiod : " 
               << ymax << ", " << ymin << ", " << yperiod << endl
               << "zmax, zmin, zperiod : " 
               << zmax << ", " << zmin << ", " << zperiod << endl;

            // clean up before exiting
            for ( size_t i=0; i < qlist.size(); i++) 
               qlist.pop_back();
            for ( size_t i=0; i < atom_type_list.size(); i++)
               atom_type_list.pop_back();
            return EXIT_FAILURE;
         }

         // Calculate xlo, ylo, zlo from xmin, ymin, zmin and xperiod,
         //  yperiod, zperiod
         //xlo = xmin - ( xperiod - (xmax - xmin) )/2.0;
         //ylo = ymin - ( yperiod - (ymax - ymin) )/2.0;
         //zlo = zmin - ( zperiod - (zmax - zmin) )/2.0;
          
         // Note: The above causes the stem probe to be misplaced.
         // TODO: shift everything so that the lower left corner is (0,0,0)

         xlo = 0.0;
         ylo = 0.0;
         zlo = 0.0;
         //if ( input_flag_debug )
         //   cout << "node " << mynode << ", "
         //      << "ymin, yperiod : "
         //      << ylo  << ", " << yperiod << endl
         //      << "     " << mynode << ", "
         //      << "xmin, xperiod : "
         //      << xlo  << ", " << xperiod << endl;

         ////////////////////////////////////////////////////////////////
         // shift all atoms so that all coordinates are positive
         ////////////////////////////////////////////////////////////////
         for ( size_t i=0; i < qlist.size() ; ++i )
         {
            qlist[i][0] = qlist[i][0] - xmin;
            qlist[i][1] = qlist[i][1] - ymin;
            qlist[i][2] = qlist[i][2] - zmin;
         }

         ////////////////////////////////////////////////////////////////
         // copy atom positions and atomic number to qq_contig & Z_contig
         ////////////////////////////////////////////////////////////////
         if ( 
               qlist.size() != declared_population 
               ||
               qlist.size() != atom_type_list.size()
               )
         {
            cerr << "node " << mynode << ", "
               << "Error reading position data file; "
               << "number_of_read_atoms " << number_of_read_atoms 
              << ",  declared population " << declared_population 
              << ",  qlist.size() " << qlist.size()
              << ",  atom_type_list.size() " << atom_type_list.size()
              << endl;
            for ( size_t i=0; i < qlist.size(); i++) 
               qlist.pop_back();
            for ( size_t i=0; i < atom_type_list.size(); i++)
               atom_type_list.pop_back();
            return EXIT_FAILURE;
         }

         total_population = atom_type_list.size();

         qq_contig = new double[ 3 * qlist.size()]; // 3-D
         Z_contig = new unsigned int[ qlist.size() ];
         for ( size_t i=0; i < qlist.size() ; ++i )
         {
            qq_contig[ 3 * i ]      = qlist[i][0];
            qq_contig[ 3 * i + 1 ]  = qlist[i][1];
            qq_contig[ 3 * i + 2 ]  = qlist[i][2];

            delete[] qlist[i];// TODO : is this appropriate?
            qlist[i] = NULL;  // TODO : is this necessary?

            Z_contig[ i ] = atom_type_list[i];
         }

         ////////////////////////////////////////////////////////////////
         // clean up
         ////////////////////////////////////////////////////////////////
         for ( size_t i=0; i < qlist.size(); i++) 
            qlist.pop_back();
         for ( size_t i=0; i < atom_type_list.size(); i++)
            atom_type_list.pop_back();
         ////////////////////////////////////////////////////////////////
      }
      else
      {
         cerr << "node " << mynode << ", "
            << "Error opening position file: " << filename << endl;
         data_file.close();
         return EXIT_FAILURE;
      }

      data_file.close();
   }  // end of the rootnode block

   // Broadcast results to the remaining nodes
   MPI_Bcast( &total_population, 1, MPI_UNSIGNED, rootnode, comm);

   if ( mynode != rootnode )
   {
      qq_contig = new double[ 3 * total_population ];
      Z_contig = new unsigned int[ total_population ];
   }

   MPI_Bcast( qq_contig, 3 * total_population, MPI_DOUBLE, rootnode, comm);

   MPI_Bcast( Z_contig, total_population, MPI_UNSIGNED, rootnode, comm);

   return EXIT_SUCCESS;
}




size_t TEM_NS::atom_element_abbrev_to_Z( const string& element_name)
{
   // Precondition : 
   //    - element_name is a lower case abbreviation of an element
   //       having atomic number between 1 and 103
   if ( element_name ==  "h" ) return 1;
   if ( element_name ==  "1" ) return 1;
   if ( element_name ==  "he" ) return 2;
   if ( element_name ==  "2" ) return 2;
   if ( element_name ==  "li" ) return 3;
   if ( element_name ==  "3" ) return 3;
   if ( element_name ==  "be" ) return 4;
   if ( element_name ==  "4" ) return 4;
   if ( element_name ==  "b" ) return 5;
   if ( element_name ==  "5" ) return 5;
   if ( element_name ==  "c" ) return 6;
   if ( element_name ==  "6" ) return 6;
   if ( element_name ==  "n" ) return 7;
   if ( element_name ==  "7" ) return 7;
   if ( element_name ==  "o" ) return 8;
   if ( element_name ==  "8" ) return 8;
   if ( element_name ==  "f" ) return 9;
   if ( element_name ==  "9" ) return 9;
   if ( element_name ==  "ne" ) return 10;
   if ( element_name ==  "10" ) return 10;
   if ( element_name ==  "na" ) return 11;
   if ( element_name ==  "11" ) return 11;
   if ( element_name ==  "mg" ) return 12;
   if ( element_name ==  "12" ) return 12;
   if ( element_name ==  "al" ) return 13;
   if ( element_name ==  "13" ) return 13;
   if ( element_name ==  "si" ) return 14;
   if ( element_name ==  "14" ) return 14;
   if ( element_name ==  "p" ) return 15;
   if ( element_name ==  "15" ) return 15;
   if ( element_name ==  "s" ) return 16;
   if ( element_name ==  "16" ) return 16;
   if ( element_name ==  "cl" ) return 17;
   if ( element_name ==  "17" ) return 17;
   if ( element_name ==  "ar" ) return 18;
   if ( element_name ==  "18" ) return 18;
   if ( element_name ==  "k" ) return 19;
   if ( element_name ==  "19" ) return 19;
   if ( element_name ==  "ca" ) return 20;
   if ( element_name ==  "20" ) return 20;
   if ( element_name ==  "sc" ) return 21;
   if ( element_name ==  "21" ) return 21;
   if ( element_name ==  "ti" ) return 22;
   if ( element_name ==  "22" ) return 22;
   if ( element_name ==  "v" ) return 23;
   if ( element_name ==  "23" ) return 23;
   if ( element_name ==  "cr" ) return 24;
   if ( element_name ==  "24" ) return 24;
   if ( element_name ==  "mn" ) return 25;
   if ( element_name ==  "25" ) return 25;
   if ( element_name ==  "fe" ) return 26;
   if ( element_name ==  "26" ) return 26;
   if ( element_name ==  "co" ) return 27;
   if ( element_name ==  "27" ) return 27;
   if ( element_name ==  "ni" ) return 28;
   if ( element_name ==  "28" ) return 28;
   if ( element_name ==  "cu" ) return 29;
   if ( element_name ==  "29" ) return 29;
   if ( element_name ==  "zn" ) return 30;
   if ( element_name ==  "30" ) return 30;
   if ( element_name ==  "ga" ) return 31;
   if ( element_name ==  "31" ) return 31;
   if ( element_name ==  "ge" ) return 32;
   if ( element_name ==  "32" ) return 32;
   if ( element_name ==  "as" ) return 33;
   if ( element_name ==  "33" ) return 33;
   if ( element_name ==  "se" ) return 34;
   if ( element_name ==  "34" ) return 34;
   if ( element_name ==  "br" ) return 35;
   if ( element_name ==  "35" ) return 35;
   if ( element_name ==  "kr" ) return 36;
   if ( element_name ==  "36" ) return 36;
   if ( element_name ==  "rb" ) return 37;
   if ( element_name ==  "37" ) return 37;
   if ( element_name ==  "sr" ) return 38;
   if ( element_name ==  "38" ) return 38;
   if ( element_name ==  "y" ) return 39;
   if ( element_name ==  "39" ) return 39;
   if ( element_name ==  "zr" ) return 40;
   if ( element_name ==  "40" ) return 40;
   if ( element_name ==  "nb" ) return 41;
   if ( element_name ==  "41" ) return 41;
   if ( element_name ==  "mo" ) return 42;
   if ( element_name ==  "42" ) return 42;
   if ( element_name ==  "tc" ) return 43;
   if ( element_name ==  "43" ) return 43;
   if ( element_name ==  "ru" ) return 44;
   if ( element_name ==  "44" ) return 44;
   if ( element_name ==  "rh" ) return 45;
   if ( element_name ==  "45" ) return 45;
   if ( element_name ==  "pd" ) return 46;
   if ( element_name ==  "46" ) return 46;
   if ( element_name ==  "ag" ) return 47;
   if ( element_name ==  "47" ) return 47;
   if ( element_name ==  "cd" ) return 48;
   if ( element_name ==  "48" ) return 48;
   if ( element_name ==  "in" ) return 49;
   if ( element_name ==  "49" ) return 49;
   if ( element_name ==  "sn" ) return 50;
   if ( element_name ==  "50" ) return 50;
   if ( element_name ==  "sb" ) return 51;
   if ( element_name ==  "51" ) return 51;
   if ( element_name ==  "te" ) return 52;
   if ( element_name ==  "52" ) return 52;
   if ( element_name ==  "i" ) return 53;
   if ( element_name ==  "53" ) return 53;
   if ( element_name ==  "xe" ) return 54;
   if ( element_name ==  "54" ) return 54;
   if ( element_name ==  "cs" ) return 55;
   if ( element_name ==  "55" ) return 55;
   if ( element_name ==  "ba" ) return 56;
   if ( element_name ==  "56" ) return 56;
   if ( element_name ==  "la" ) return 57;
   if ( element_name ==  "57" ) return 57;
   if ( element_name ==  "ce" ) return 58;
   if ( element_name ==  "58" ) return 58;
   if ( element_name ==  "pr" ) return 59;
   if ( element_name ==  "59" ) return 59;
   if ( element_name ==  "nd" ) return 60;
   if ( element_name ==  "60" ) return 60;
   if ( element_name ==  "pm" ) return 61;
   if ( element_name ==  "61" ) return 61;
   if ( element_name ==  "sm" ) return 62;
   if ( element_name ==  "62" ) return 62;
   if ( element_name ==  "eu" ) return 63;
   if ( element_name ==  "63" ) return 63;
   if ( element_name ==  "gd" ) return 64;
   if ( element_name ==  "64" ) return 64;
   if ( element_name ==  "tb" ) return 65;
   if ( element_name ==  "65" ) return 65;
   if ( element_name ==  "dy" ) return 66;
   if ( element_name ==  "66" ) return 66;
   if ( element_name ==  "ho" ) return 67;
   if ( element_name ==  "67" ) return 67;
   if ( element_name ==  "er" ) return 68;
   if ( element_name ==  "68" ) return 68;
   if ( element_name ==  "tm" ) return 69;
   if ( element_name ==  "69" ) return 69;
   if ( element_name ==  "yb" ) return 70;
   if ( element_name ==  "70" ) return 70;
   if ( element_name ==  "lu" ) return 71;
   if ( element_name ==  "71" ) return 71;
   if ( element_name ==  "hf" ) return 72;
   if ( element_name ==  "72" ) return 72;
   if ( element_name ==  "ta" ) return 73;
   if ( element_name ==  "73" ) return 73;
   if ( element_name ==  "w" ) return 74;
   if ( element_name ==  "74" ) return 74;
   if ( element_name ==  "re" ) return 75;
   if ( element_name ==  "75" ) return 75;
   if ( element_name ==  "os" ) return 76;
   if ( element_name ==  "76" ) return 76;
   if ( element_name ==  "ir" ) return 77;
   if ( element_name ==  "77" ) return 77;
   if ( element_name ==  "pt" ) return 78;
   if ( element_name ==  "78" ) return 78;
   if ( element_name ==  "au" ) return 79;
   if ( element_name ==  "79" ) return 79;
   if ( element_name ==  "hg" ) return 80;
   if ( element_name ==  "80" ) return 80;
   if ( element_name ==  "tl" ) return 81;
   if ( element_name ==  "81" ) return 81;
   if ( element_name ==  "pb" ) return 82;
   if ( element_name ==  "82" ) return 82;
   if ( element_name ==  "bi" ) return 83;
   if ( element_name ==  "83" ) return 83;
   if ( element_name ==  "po" ) return 84;
   if ( element_name ==  "84" ) return 84;
   if ( element_name ==  "at" ) return 85;
   if ( element_name ==  "85" ) return 85;
   if ( element_name ==  "rn" ) return 86;
   if ( element_name ==  "86" ) return 86;
   if ( element_name ==  "fr" ) return 87;
   if ( element_name ==  "87" ) return 87;
   if ( element_name ==  "ra" ) return 88;
   if ( element_name ==  "88" ) return 88;
   if ( element_name ==  "ac" ) return 89;
   if ( element_name ==  "89" ) return 89;
   if ( element_name ==  "th" ) return 90;
   if ( element_name ==  "90" ) return 90;
   if ( element_name ==  "pa" ) return 91;
   if ( element_name ==  "91" ) return 91;
   if ( element_name ==  "u" ) return 92;
   if ( element_name ==  "92" ) return 92;
   if ( element_name ==  "np" ) return 93;
   if ( element_name ==  "93" ) return 93;
   if ( element_name ==  "pu" ) return 94;
   if ( element_name ==  "94" ) return 94;
   if ( element_name ==  "am" ) return 95;
   if ( element_name ==  "95" ) return 95;
   if ( element_name ==  "cm" ) return 96;
   if ( element_name ==  "96" ) return 96;
   if ( element_name ==  "bk" ) return 97;
   if ( element_name ==  "97" ) return 97;
   if ( element_name ==  "cf" ) return 98;
   if ( element_name ==  "98" ) return 98;
   if ( element_name ==  "es" ) return 99;
   if ( element_name ==  "99" ) return 99;
   if ( element_name ==  "fm" ) return 100;
   if ( element_name ==  "100" ) return 100;
   if ( element_name ==  "md" ) return 101;
   if ( element_name ==  "101" ) return 101;
   if ( element_name ==  "no" ) return 102;
   if ( element_name ==  "102" ) return 102;
   if ( element_name ==  "lr" ) return 103;
   if ( element_name ==  "103" ) return 103;
   
   // If none of the modeled elements matched, return the invalid
   //  atomic number 0 .
   return 0;
}

#endif
