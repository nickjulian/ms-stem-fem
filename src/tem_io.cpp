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
      const string& filename,  // only valid on root node
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
                     << "tiltxy, tiltxz, tiltyz, data_descriptor1, "//debug
                     << " data_descriptor2, data_descriptor3 : "// debug
                  << tiltxy << ", " << tiltxz << ", " // debug
                  << tiltyz << ", "// debug
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
               << "Error reading position file" << endl
               << "   declared_population : " << declared_population 
               << endl
               << "   positions read : " << qlist.size() << endl
               << "   unique atom ID numbers : " 
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
         ///////////////////////////////////////////////////////////////


         ///////////////////////////////////////////////////////////////
         // clean up
         for ( size_t i = 0; i < qlist.size(); i++)
         {
            //delete[] qlist.back(); // taken care of in previous loops
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
         ///////////////////////////////////////////////////////////////
      }
      else 
      {
         cerr << "node " << mynode << ", "
            << "Error opening position file" << endl;
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


int TEM_NS::read_position_lammps_file_nonmpi( 
      const string& filename,  // only valid on root node
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
      const unsigned int& input_flag_debug
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
   
   if ( 1 )
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
               cerr //<< "node " << mynode << ", "
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
               cout //<< "node " << mynode << ", "
                  << "declared_population : "
                  << declared_population << endl;
         }
         else
         {
            cerr //<< "node " << mynode << ", "
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
               cerr //<< "node " << mynode << ", "
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
               cerr //<< "node " << mynode << ", "
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
               cout //<< "node " << mynode << ", "
                  << "numberofspecies : "
                  << numberofspecies << endl;

            if ( numberofspecies != Zcategories.size() )
            {
               cerr //<< "node " << mynode << ", "   
                  << "The first line of the input file must begin with "
                  << "a space separated list of atomic numbers having order "
                 << " corresponding with the lammps species enumeration"
                 << endl;
               return EXIT_FAILURE;
            }
         }
         else
         {
            cerr //<< "node " << mynode << ", "   
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
               cerr //<< "node " << mynode << ", "   // debug
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
               cout //<< "node " << mynode << ", "
                  << "xlower, xupper : "
                  << xlower << ", " << xupper << endl;
         }
         else
         {
            cerr //<< "node " << mynode << ", "
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
               cerr //<< "node " << mynode << ", "
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
               cerr //<< "node " << mynode << ", "
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
               cout //<< "node " << mynode << ", "
                  << "ylower, yupper : "
                  << ylower << ", " << yupper << endl;
         }
         else
         {
            cerr //<< "node " << mynode << ", "
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
               cerr //<< "node " << mynode << ", "
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
               cerr //<< "node " << mynode << ", "
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
               cout //<< "node " << mynode << ", "
                  << "zlower, zupper : "
                  << zlower << ", " << zupper << endl;
         }
         else
         {
            cerr //<< "node " << mynode << ", "
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
               cerr //<< "node " << mynode << ", "
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
               cerr //<< "node " << mynode << ", "
                     << "Error reading position data file; " 
                     << "tiltxy, tiltxz, tiltyz, data_descriptor1, "//debug
                     << " data_descriptor2, data_descriptor3 : "// debug
                  << tiltxy << ", " << tiltxz << ", " // debug
                  << tiltyz << ", "// debug
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
               cout //<< "node " << mynode << ", "
                     << "tiltxy, tiltxz, tiltyz: "
                     << tiltxy
                     << ", " << tiltxz
                     << ", " << tiltyz << endl;
         }
         else
         {
            cerr //<< "node " << mynode << ", "
               << "Error reading position data file; data_line : "
               << data_line << endl;
            cerr << "data_file.good() : " << data_file.good();
            return EXIT_FAILURE;
         }

         // check that boundary values are valid
         if ( xupper <= xlower ) 
         {
            cerr //<< "node " << mynode << ", "
               << "Error reading position data file; xhi <= xlo" << endl;
            return EXIT_FAILURE;
         }
         if ( yupper <= ylower ) 
         {
            cerr //<< "node " << mynode << ", "
               << "Error reading position data file; yhi <= ylo" << endl;
            return EXIT_FAILURE;
         }
         if ( zupper <= zlower ) 
         {
            cerr //<< "node " << mynode << ", "
               << "Error reading position data file; zhi <= zlo" << endl;
            return EXIT_FAILURE;
         }
         if ( tiltxy != 0.0 || tiltxz != 0.0 || tiltyz != 0.0 )
         {
            cerr //<< "node " << mynode << ", "
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
         
         if ( //mynode == rootnode && 
               input_flag_debug )
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
               cerr //<< "node " << mynode << ", "
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
               cerr //<< "node " << mynode << ", "
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
            cerr //<< "node " << mynode << ", "
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
                  cerr //<< "node " << mynode << ", "
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
               cerr //<< "node " << mynode << ", "
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
            cerr //<< "node " << mynode << ", "
               << "Error reading position file" << endl
               << "   declared_population : " << declared_population 
               << endl
               << "   positions read : " << qlist.size() << endl
               << "   unique atom ID numbers : " 
               << atom_number_list_uniqued.size() << endl;
            for ( size_t i=0; i<qlist.size(); i++) 
            {
               delete[] qlist[i];
               qlist.pop_back();
            }
            return EXIT_FAILURE;
         }

         ///////////////////////////////////////////////////////////////
         // clean up
         for ( size_t i = 0; i < qlist.size(); i++)
         {
            //delete[] qlist.back(); // taken care of in previous loops
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
         ///////////////////////////////////////////////////////////////
      }
      else 
      {
         cerr //<< "node " << mynode << ", "
            << "Error opening position file" << endl;
         data_file.close();
         return EXIT_FAILURE;
      }

      data_file.close();
   } // end of the rootnode block
 
   // Broadcast results to the remaining nodes
   //MPI_Bcast( &total_population, 1, MPI_UNSIGNED, rootnode, comm);

   //if ( mynode != rootnode )
   //{
   //   qq_contig = new double[ 3 * total_population ];
   //   Z_contig = new unsigned int[ total_population ];
   //}

   //MPI_Bcast( qq_contig, 3 * total_population, MPI_DOUBLE, rootnode, comm);

   //MPI_Bcast( Z_contig, total_population, MPI_UNSIGNED, rootnode, comm);

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

int TEM_NS::read_parameter_file(
         const string& parameter_filename, 
         string& model_file,
         input_flags& flags,
         string& output_prefix,
         ptrdiff_t& Nx,
         ptrdiff_t& Ny,
         double& VV,
         double& defocus,
         double& alpha_max,
         double& defocus_spread,
         double& condenser_illumination_angle,
         double& Cs3,
         double& Cs5,
         double& detector_inner_angle,
         double& detector_outer_angle,
         double& raster_spacing,
         double& azimuthal_binning_size_factor,
         double& minSliceThickness,
         const int& mynode,
         const int& rootnode,
         MPI_Comm comm
      )
{
   // sufficiency and conflicts of the input should be checked by caller
   int Nx_int;
   int Ny_int;
   //if ( mynode == rootnode )
   //{  // NOTE:
   //    If only root node is allowed to access the file, then
   //    root node will have to find a way to Bcast the strings it 
   //    reads for output_prefix and model_file. Unfortunately, I've
   //    not yet found a way to do Bcast strings without assuming that
   //    they're ASCII. If they happen to be UTF, then I don't know
   //    how to measure their length and thus can't Bcast them yet.
   ifstream data_file( parameter_filename.c_str() );
   if ( data_file.is_open() )
   {
      //cout << " reading " << parameter_filename.c_str() << endl; // debug
      string data_line;
      while( getline( data_file, data_line) && data_file.good() )
      {
         // skip lines containing only whitespace
         size_t first = data_line.find_first_not_of(" \t" );
         while( first == std::string::npos)//npos: max size of a string
         {
            getline( data_file, data_line);
            first = data_line.find_first_not_of(" \t" );
         }

         istringstream data_line_stream( data_line );
         string data_descriptor; 
         data_line_stream >> data_descriptor;
         // Transform ASCII characters to lower case.
         //    TODO: find another solution for transforming case
         //    Supposedly this won't work for UTF. ?!
         //    For Unicode support, I've heard that I should use 
         //    toLower from the ICU library.
         //    Should they even be transformed to lower case?
         transform( data_descriptor.begin(), data_descriptor.end(), 
                    data_descriptor.begin(), (int(*)(int))tolower );

         if ( ! data_descriptor.compare("debug") )
         {
            flags.debug = 1;
         }
         else if ( ! data_descriptor.compare("output_prefix") )
         {
            data_line_stream >> output_prefix;
            flags.o = 1;
         }
         else if ( ! data_descriptor.compare("samples_x") )
         {
            data_line_stream >> Nx_int;
            Nx = (ptrdiff_t) Nx_int;
            flags.nx = 1;
         }
         else if ( ! data_descriptor.compare("samples_y") )
         {
            data_line_stream >> Ny_int;
            Ny = (ptrdiff_t) Ny_int;
            flags.ny = 1;
         }
         else if ( ! data_descriptor.compare("model_file") )
         {
            data_line_stream >> model_file;
            flags.a = 1;
         }
         else if ( ! data_descriptor.compare("microscope_voltage") )
         {
            data_line_stream >> VV;
            flags.microscope_voltage = 1;
         }
         else if ( ! data_descriptor.compare("scherzer_defocus") )
         {
            flags.scherzer_defocus = 1;
         }
         else if ( ! data_descriptor.compare("scherzer_alphamax") )
         {
            flags.scherzer_alphamax = 1;
         }
         else if ( ! data_descriptor.compare("scherzer_cs3") )
         {
            flags.scherzer_cs3 = 1;
         }
         else if ( ! data_descriptor.compare("defocus") )
         {
            data_line_stream >> defocus;
            flags.defocus= 1;
         }
         else if ( ! data_descriptor.compare("alphamax") )
         {
            data_line_stream >> alpha_max;
            flags.alpha_max = 1;
         }
         else if ( ! data_descriptor.compare("spread") )
         {
            data_line_stream >> defocus_spread;
            data_line_stream >> condenser_illumination_angle;
            flags.spread= 1;
         }
         else if ( ! data_descriptor.compare("cs3") )
         {
            data_line_stream >> Cs3;
            flags.cs3 = 1;
         }
         else if ( ! data_descriptor.compare("cs5") )
         {
            data_line_stream >> Cs5;
            flags.cs5 = 1;
         }
         else if ( ! data_descriptor.compare("detectorangles") )
         {
            data_line_stream >> detector_inner_angle;
            data_line_stream >> detector_outer_angle;

            if ( detector_inner_angle > detector_outer_angle )
            {
               double tmp_detector_angle;
               tmp_detector_angle = detector_outer_angle;
               detector_outer_angle = detector_inner_angle;
               detector_inner_angle = tmp_detector_angle;
            }
            flags.adfstem_detector_angles = 1;
         }
         else if ( ! data_descriptor.compare("raster_spacing") )
         {
            data_line_stream >> raster_spacing;
            flags.raster_spacing = 1;
         }
         else if ( ! data_descriptor.compare("adfstemuncorrfem") )
         {
            flags.adfstem_uncorrected = 1;
            flags.fem = 1;
         }
         else if ( ! data_descriptor.compare("adfstemcorrfem") )
         {
            flags.adfstem_corrected = 1;
            flags.fem = 1;
         }
         else if ( ! data_descriptor.compare("vgt17") )
         {
            flags.gt17 = 1;
         }
         else if ( ! data_descriptor.compare("vomega") )
         {
            flags.d1 = 1;
         }
         else if ( ! data_descriptor.compare("vr") )
         {
            flags.d2 = 1;
         }
         else if ( ! data_descriptor.compare("vre") )
         {
            flags.d3 = 1;
         }
         else if ( ! data_descriptor.compare("omegavimage") )
         {
            flags.d4 = 1;
         }
         else if ( ! data_descriptor.compare("vphi") )
         {
            flags.rva = 1;
         }
         else if ( ! data_descriptor.compare("adfstemcorr") )
         {
            flags.adfstem_corrected = 1;
         }
         else if ( ! data_descriptor.compare("adfstemuncorr") )
         {
            flags.adfstem_uncorrected = 1;
         }
         else if ( ! data_descriptor.compare("bfctemcorr") )
         {
            flags.bfctem_corrected = 1;
         }
         else if ( ! data_descriptor.compare("bfctemuncorr") )
         {
            flags.bfctem_uncorrected = 1;
         }
         else if ( ! data_descriptor.compare("paptif") )
         {
            flags.pap_tif = 1;
         }
         else if ( ! data_descriptor.compare("dupe") )
         {
            flags.dupe = 1;
         }
         else if ( ! data_descriptor.compare("dr") )
         {
            data_line_stream >> azimuthal_binning_size_factor;
         }
         else if ( ! data_descriptor.compare("minslice") )
         {
            data_line_stream >> minSliceThickness; // unexpected camel
         }
         else if ( ! data_descriptor.compare("images") )
         {
            flags.image_output = 1;
         }
         else if ( ! data_descriptor.compare("diffractionimages") )
         {
            flags.diffraction_output = 1;
         }
         else if ( ! data_descriptor.compare("netcdfimages") )
         {
            flags.netcdf_images = 1;
         }
         else if ( ! data_descriptor.compare("netcdfvariance") )
         {
            flags.netcdf_variance = 1;
         }
         else if ( ! data_descriptor.compare("correlograph") )
         {
            flags.correlograph = 1;
         }
         else if ( ! data_descriptor.compare("correlograph_variance") )
         {
            flags.correlograph_variance = 1;
            flags.correlograph = 1;
         }
         else if ( ! data_descriptor.compare("correlograph_everyimage") )
         {
            flags.correlograph_everyimage = 1;
            flags.correlograph = 1;
         }
         else if ( ! data_descriptor.compare("correlograph_everytxt") )
         {
            flags.correlograph_everytxt = 1;
            flags.correlograph = 1;
         }
         //else if ( ! data_descriptor.compare("correlograph_everynetcdf") )
         //{
         //   flags.correlograph_everynetcdf = 1;
         //   flags.correlograph = 1;
         //}
         else
         {
            if ( mynode == rootnode )
            {
               cerr << "Error, unexpected argument in parameter file: "
                  << data_line << endl;
               cout << "Acceptable phrases in the parameter file:" << endl
                  << "  samples_x <integer> " << endl
                  << "  samples_y <integer> " << endl
                  << "    WARNING: it's recommended that samples_x and "
                  << "    samples_y be identical and even, otherwise "
                  << "    bad data will creep in from the edges of "
                  << "    reciprocal space" << endl
                  << "  model_file <path/to/model/file> " << endl
                  << "  microscope_voltage <voltage value> " << endl
                  << "  scherzer_defocus" << endl
                  << "  scherzer_alphamax" << endl
                  << "  scherzer_cs3" << endl
                  << "  defocus <defocus value>" << endl
                  << "  alphamax <maximum alpha>" << endl
                  << "  spread <defocus spread>" 
                  << " <condenser illumination angle>" << endl
                  << "  cs3 <cs3 value>" << endl
                  << "  cs5 <cs5 value>" << endl
                  << "  detectorangles <inner_angle>"
                  <<    " <outer_angle [radians]>" << endl
                  << "  raster_spacing <distance [A] between probe"
                  << "  positions>" << endl
                  << "  adfstemcorrfem" << endl
                  << "  adfstemuncorrfem" << endl
                  << "  Vgt17" << endl
                  << "  Vomega" << endl
                  << "  Vr" << endl
                  << "  Vre" << endl
                  << "  Omegavimage" << endl
                  << "  Vphi" << endl
                  << "  correlograph" << endl
                  << "  correlograph_variance" << endl
                  << "  correlograph_everyimage" << endl
                  << "  correlograph_everytxt" << endl
                  //<< "  correlograph_everynetcdf" << endl
                  << "  adfstemcorr" << endl
                  << "  adfstemuncorr" << endl
                  //<< "  bfctemcorr" << endl
                  //<< "  bfctemuncorr" << endl
                  << "  paptif" << endl
                  << "  dupe" << endl
                  << "  dr <azimuthal binning size factor>" << endl
                  << "  minslice <minimum slice thickness>" << endl
                  << "  images" << endl
                  << "  diffractionimages" << endl
                  << "  netcdfimages" << endl
                  << "  netcdfvariance" << endl
                  << "  debug" << endl;
            }
            data_file.close();
            return EXIT_FAILURE;
         }
      }
   }
   else 
   {
      cerr << "Error opening parameter file: " << parameter_filename 
         << endl;
      data_file.close();
      return EXIT_FAILURE;
   }
   data_file.close();
   //} // rootnode

   if ( flags.fem && (!( flags.d1  || flags.d2 || flags.d3 || flags.d4 || flags.gt17 || flags.rva || flags.correlograph)) )
   {
      flags.d1 = 1; // default fem mode
   }

   //// Broadcast results to the remaining nodes
   //
   //// pack the boolean values into an array and then broadcast it
   //unsigned int* flags_to_send ;
   //flags_to_send = new unsigned int[ 30 ];
   //if ( mynode == rootnode ) 
   //{
   //   flags_to_send[0] = flags.o;
   //   flags_to_send[1] = flags.nx;
   //   flags_to_send[2] = flags.ny;
   //   flags_to_send[3] = flags.a;
   //   flags_to_send[4] = flags.microscope_voltage;
   //   flags_to_send[5] = flags.scherzer_defocus;
   //   flags_to_send[6] = flags.scherzer_alphamax;
   //   flags_to_send[7] = flags.scherzer_cs3;
   //   flags_to_send[8] = flags.defocus;
   //   flags_to_send[9] = flags.alpha_max;
   //   flags_to_send[10] = flags.spread;
   //   flags_to_send[11] = flags.cs3;
   //   flags_to_send[12] = flags.cs5;
   //   flags_to_send[13] = flags.raster_spacing;
   //   flags_to_send[14] = flags.fem;
   //   flags_to_send[15] = flags.adfstem_uncorrected;
   //   flags_to_send[16] = flags.adfstem_corrected;
   //   flags_to_send[17] = flags.gt17;
   //   flags_to_send[18] = flags.d1;
   //   flags_to_send[19] = flags.d2;
   //   flags_to_send[20] = flags.d3;
   //   flags_to_send[21] = flags.d4;
   //   flags_to_send[22] = flags.bfctem_corrected;
   //   flags_to_send[23] = flags.bfctem_uncorrected;
   //   flags_to_send[24] = flags.pap_tif;
   //   flags_to_send[25] = flags.dupe;
   //   flags_to_send[26] = flags.image_output;
   //   flags_to_send[27] = flags.netcdf_images;
   //   flags_to_send[28] = flags.netcdf_variance;
   //   flags_to_send[29] = flags.debug;
   //}
   //MPI_Bcast( flags_to_send , 30, MPI_UNSIGNED, rootnode, comm);

   //int* N_int;
   //N_int = int[2];
   //if ( mynode == rootnode )
   //{
   //   N_int[0] = Nx_int;
   //   N_int[1] = Ny_int;
   //}
   //MPI_Bcast( N_int, 2, MPI_INT, rootnode, comm);

   //// TODO: broadcast output_prefix string
   //// TODO: broadcast model_file string

   //// pack the double values into an array and then broadcast it
   //double* doubles_to_send;
   //doubles_to_send = new double[ 10 ];
   //if ( mynode == rootnode ) 
   //{
   //   doubles_to_send[0] = VV;
   //   doubles_to_send[1] = defocus;
   //   doubles_to_send[2] = alpha_max;
   //   doubles_to_send[3] = defocus_spread;
   //   doubles_to_send[4] = condenser_illumination_angle;
   //   doubles_to_send[5] = Cs3;
   //   doubles_to_send[6] = Cs5;
   //   doubles_to_send[7] = raster_spacing;
   //   doubles_to_send[8] = azimuthal_binning_size_factor;
   //   doubles_to_send[9] = minSliceThickness;
   //}
   //MPI_Bcast( doubles_to_send, 10, MPI_DOUBLE, rootnode, comm);

   //// unpack the broadcast values
   //if ( mynode != rootnode ) 
   //{
   //   VV = doubles_to_send[0];
   //   defocus = doubles_to_send[1];
   //   alpha_max = doubles_to_send[2];
   //   defocus_spread = doubles_to_send[3];
   //   condenser_illumination_angle = doubles_to_send[4];
   //   Cs3 = doubles_to_send[5];
   //   Cs5 = doubles_to_send[6];
   //   raster_spacing = doubles_to_send[7];
   //   azimuthal_binning_size_factor = doubles_to_send[8];
   //   minSliceThickness = doubles_to_send[9];

   //   Nx = (ptrdiff_t) N_int[0];
   //   Ny = (ptrdiff_t) N_int[1];

   //   flags.o = flags_to_send[0];
   //   flags.nx = flags_to_send[1];
   //   flags.ny = flags_to_send[2];
   //   flags.a = flags_to_send[3];
   //   flags.microscope_voltage = flags_to_send[4];
   //   flags.scherzer_defocus = flags_to_send[5];
   //   flags.scherzer_alphamax = flags_to_send[6];
   //   flags.scherzer_cs3 = flags_to_send[7];
   //   flags.defocus = flags_to_send[8];
   //   flags.alpha_max = flags_to_send[9];
   //   flags.spread = flags_to_send[10];
   //   flags.cs3 = flags_to_send[11];
   //   flags.cs5 = flags_to_send[12];
   //   flags.raster_spacing = flags_to_send[13];
   //   flags.fem = flags_to_send[14];
   //   flags.adfstem_uncorrected = flags_to_send[15];
   //   flags.adfstem_corrected = flags_to_send[16];
   //   flags.gt17 = flags_to_send[17];
   //   flags.d1 = flags_to_send[18];
   //   flags.d2 = flags_to_send[19];
   //   flags.d3 = flags_to_send[20];
   //   flags.d4 = flags_to_send[21];
   //   flags.bfctem_corrected = flags_to_send[22];
   //   flags.bfctem_uncorrected = flags_to_send[23];
   //   flags.pap_tif = flags_to_send[24];
   //   flags.dupe = flags_to_send[25];
   //   flags.image_output = flags_to_send[26];
   //   flags.netcdf_images = flags_to_send[27];
   //   flags.netcdf_variance = flags_to_send[28];
   //   flags.debug = flags_to_send[29];
   //}

   //delete[] flags_to_send;
   //delete[] doubles_to_send;
   //delete[] N_int;

   return EXIT_SUCCESS;
}

int TEM_NS::read_position_xyz_file(
      const string& filename, // only valid on root node
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


int TEM_NS::read_position_xyz_file_nonmpi(
      const string& filename, // only valid on root node
      double*& qq_contig,
      unsigned int*& Z_contig,
      unsigned int& total_population,
      double& xlo, double& ylo, double& zlo,
      double& xperiod, double& yperiod, double& zperiod,
      const unsigned int& input_flag_debug
      )
{
   // Currently elements 1 through 103 are modeled.

   if ( 1 )//mynode == rootnode ) 
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
               cerr //<< "node " << mynode << ", "
                  << "Error reading position data file; "
                  << "declared_population : "
                  << declared_population << endl;
               cerr << " current line should be : " << endl;
               cerr << "<number of atoms>  " << endl;
               return EXIT_FAILURE;
            }
            if ( input_flag_debug )
               cout //<< "node " << mynode << ", "
                  << "declared_population : "
                  << declared_population << endl;
         }
         else
         {
            cerr //<< "node " << mynode << ", "
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
               cerr //<< "node " << mynode << ", "
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
            cerr //<< "node " << mynode << ", "
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
                  cerr //<< "node " << mynode << ", "
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
               cerr //<< "node " << mynode << ", "
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
            cerr //<< "node " << mynode << ", "
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
            cerr //<< "node " << mynode << ", "
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
         cerr //<< "node " << mynode << ", "
            << "Error opening position file: " << filename << endl;
         data_file.close();
         return EXIT_FAILURE;
      }

      data_file.close();
   }  // end of the rootnode block

   // Broadcast results to the remaining nodes
   //MPI_Bcast( &total_population, 1, MPI_UNSIGNED, rootnode, comm);

   //if ( mynode != rootnode )
   //{
   //   qq_contig = new double[ 3 * total_population ];
   //   Z_contig = new unsigned int[ total_population ];
   //}

   //MPI_Bcast( qq_contig, 3 * total_population, MPI_DOUBLE, rootnode, comm);

   //MPI_Bcast( Z_contig, total_population, MPI_UNSIGNED, rootnode, comm);

   return EXIT_SUCCESS;
}

unsigned int TEM_NS::read_cmdline_options( 
   const std::vector<string>& args,
   string& model_file_name,
   input_flags& flags,
   string& output_prefix,
   ptrdiff_t& Nx,
   ptrdiff_t& Ny,
   double& VV,
   double& defocus,
   double& alpha_max,
   double& defocus_spread,
   double& condenser_illumination_angle,
   double& Cs3,
   double& Cs5,
   double& detector_inner_angle,
   double& detector_outer_angle,
   double& raster_spacing,
   double& azimuthal_binning_size_factor,
   double& minSliceThickness,
   unsigned int& dupe_x,
   unsigned int& dupe_y,
   unsigned int& dupe_z,
   const int& mynode,
   const int& rootnode,
   MPI_Comm comm
)
{
   unsigned int failflag;
   int Nx_int;
   int Ny_int;
   failflag = 0;
   for ( size_t idx=1; idx < args.size(); idx++)
   {
      if (args[idx] == "--parameter_file")// microscope parameter file name
      {
         string parameter_file_name;
         if (idx + 1 < args.size()) 
            parameter_file_name = string(args[idx + 1]);
            //istringstream( args[idx + 1] ) >> parameter_file_name;
         if ( 
               read_parameter_file(
                  parameter_file_name,
                  model_file_name,
                  flags,
                  output_prefix,
                  Nx,
                  Ny,
                  VV,
                  defocus,
                  alpha_max,
                  defocus_spread,
                  condenser_illumination_angle,
                  Cs3,
                  Cs5,
                  detector_inner_angle,
                  detector_outer_angle,
                  raster_spacing,
                  azimuthal_binning_size_factor,
                  minSliceThickness,
                  mynode,
                  rootnode,
                  MPI_COMM_WORLD
               ) == EXIT_FAILURE
         )
         {
            if ( mynode == rootnode )
            {
               cout << "Error, could not read parameter file : " 
                  << parameter_file_name << endl;
            }
            failflag = 1 ;
         }

         if ( flags.microscope_voltage && flags.nx && flags.ny ) 
            flags.m = 1;
         else
            flags.m = 0;

         flags.pf = 1;
         idx += 1;

      }
      else if ( args[idx] == "--debug" )
      {
         flags.debug  = 1;
      }
      else if ( args[idx] == "--output_prefix" )
      {
         if (idx + 1 < args.size()) 
            output_prefix = string(args[idx + 1]);//=string(argv[2]);
         flags.o = 1;
         idx += 1;
      }
      else if ( args[idx] == "-m" )
      {
         if (idx + 1 < args.size()) 
            istringstream(args[idx + 1]) >> Nx_int;
         if (idx + 2 < args.size()) 
            istringstream(args[idx + 2]) >> Ny_int;
         //if (idx + 3 < args.size()) 
         // istringstream(args[idx + 3]) >> Nz_int;
         if (idx + 3 < args.size()) istringstream(args[idx + 3]) >> VV;

         Nx = (ptrdiff_t) Nx_int;
         Ny = (ptrdiff_t) Ny_int;
         // Nz = (ptrdiff_t) Nz_int;

         flags.m = 1;
         idx += 3;
      }
      else if ( args[idx] == "-a" )
      {
         //if ( mynode == rootnode )
         //{
         //istringstream( args[idx + 1] ) >> model_file_name;
         if (idx + 1 < args.size()) 
            model_file_name = string( args[ idx + 1] );

         flags.a = 1;
         idx += 1;
      }
      else if ( args[idx] == "--scherzer_defocus")
      {
         flags.scherzer_defocus = 1;
      }
      else if ( args[idx] == "--scherzer_alphamax")
      {
         flags.scherzer_alphamax = 1;
      }
      else if ( args[idx] == "--scherzer_cs3")
      {
         flags.scherzer_cs3 = 1;
      }
      else if ( args[idx] == "--defocus" )
      {
         if (idx + 1 < args.size()) 
            istringstream( args[idx + 1] ) >> defocus;
         flags.defocus= 1;
         idx += 1;
      }
      else if ( args[idx] == "--alphamax" )
      {
         if (idx + 1 < args.size()) 
            istringstream( args[idx + 1] ) >> alpha_max;
         flags.alpha_max = 1;
         idx += 1;
      }
      else if ( args[idx] == "--spread" )
      {
         if (idx + 1 < args.size()) 
            istringstream( args[idx + 1] ) >> defocus_spread;
         if (idx + 2 < args.size()) 
            istringstream( args[idx + 2] ) >> condenser_illumination_angle;
         flags.spread = 1;
         idx += 2;
      }
      else if ( args[idx] == "--cs3" )
      {
         if (idx + 1 < args.size()) 
            istringstream( args[idx + 1] ) >> Cs3;
         flags.cs3 = 1;
         idx += 1;
      }
      else if ( args[idx] == "--cs5" )
      {
         if (idx + 1 < args.size()) 
            istringstream( args[idx + 1] ) >> Cs5;
         flags.cs5 = 1;
         idx += 1;
      }
      else if ( args[idx] == "--detectorangles" )
      {
         if (idx + 1 < args.size()) 
            istringstream( args[idx + 1] ) >> detector_inner_angle;
         if (idx + 2 < args.size()) 
            istringstream( args[idx + 1] ) >> detector_outer_angle;

         if ( detector_inner_angle > detector_outer_angle )
         {
            double tmp_detector_angle;
            tmp_detector_angle = detector_outer_angle;
            detector_outer_angle = detector_inner_angle;
            detector_inner_angle = tmp_detector_angle;
         }
         flags.adfstem_detector_angles = 1;
         idx += 1;
      }
      else if ( args[idx] == "--rasterspacing" )
      {
         if ( idx + 1 < args.size())
            istringstream( args[idx + 1] ) >> raster_spacing;
         flags.raster_spacing = 1;
         idx += 1;
      }
      else if ( args[idx] == "--adfstemuncorrfem" )
      {
         flags.adfstem_uncorrected = 1;
         flags.fem = 1;
      }
      else if ( args[idx] == "--adfstemcorrfem" )
      {
         flags.adfstem_corrected = 1;
         flags.fem = 1;
      }
      else if ( args[idx] == "--Vgt17" )
      {
         flags.gt17 = 1;
      }
      else if ( args[idx] == "--Vomega" )
      {
         flags.d1 = 1;
      }
      else if ( args[idx] == "--Vr" )
      {
         flags.d2 = 1;
      }
      else if ( args[idx] == "--Vre" )
      {
         flags.d3 = 1;
      }
      else if ( args[idx] == "--Omegavimage" )
      {
         flags.d4 = 1;
      }
      else if ( args[idx] == "--Vphi" )
      {
         flags.rva = 1;
      }
      else if ( args[idx] == "--adfstemcorr" )
      {
         flags.adfstem_corrected = 1;
      }
      else if ( args[idx] == "--adfstemuncorr" )
      {
         flags.adfstem_uncorrected = 1;
      }
      else if ( args[idx] == "--bfctemcorr" )
      {
         flags.bfctem_corrected = 1;
      }
      else if ( args[idx] == "--bfctemuncorr" )
      {
         flags.bfctem_uncorrected = 1;
      }
      else if ( args[idx] == "--paptif" )
      {
          flags.pap_tif = 1;
      }
      else if ( args[idx] == "--dupe" )
      {
         if (idx + 3 < args.size()) 
            istringstream( args[idx + 1] ) >> dupe_x;
            istringstream( args[idx + 2] ) >> dupe_y;
            istringstream( args[idx + 3] ) >> dupe_z;
         flags.dupe = 1;
         idx += 3;
      }
      else if ( args[idx] == "--dr" )
      {
         if ( idx + 1 < args.size())
            istringstream( args[idx + 1] ) 
               >> azimuthal_binning_size_factor;
         idx += 1;
      }
      else if ( args[idx] == "--minslice" )
      {
         if ( idx + 1 < args.size())
            istringstream( args[idx + 1] ) 
               >> minSliceThickness;
         idx += 1;
      }
      else if ( args[idx] == "--images" )
      {
         flags.image_output = 1;
      }
      else if ( args[idx] == "--netcdfimages" )
      {
         flags.netcdf_images = 1;
      }
      else if ( args[idx] == "--netcdfvariance" )
      {
         flags.netcdf_variance = 1;
      }
      else if ( args[idx] == "--correlograph" )
      {
         flags.correlograph = 1;
      }
      else if ( args[idx] == "--correlograph_variance" )
      {
         flags.correlograph_variance = 1;
         flags.correlograph = 1;
      }
      else if ( args[idx] == "--correlograph_everyimage" )
      {
         flags.correlograph_everyimage = 1;
         flags.correlograph = 1;
      }
      else if ( args[idx] == "--correlograph_everytxt" )
      {
         flags.correlograph_everytxt = 1;
         flags.correlograph = 1;
      }
      //else if ( args[idx] == "--correlograph_everynetcdf" )
      //{
      //   flags.correlograph_everynetcdf = 1;
      //   flags.correlograph = 1;
      //}
      else
      {
         if ( mynode == rootnode )
         {
            cerr << "Error, unexpected argument : " << args[idx] << endl;
         }
         failflag = 1;
      }
   } // iteration over command line arguments

   if ( (! flags.gt17 ) 
         && (! flags.d1 ) && (! flags.d2 ) 
         && (! flags.d3 ) && (! flags.d4 )
         && (! flags.rva) && (! flags.correlograph) )
   {
      flags.d1 = 1; // default variance calculation mode
   }

   if ( failflag ==1 )
   {
      return EXIT_FAILURE;
   }

   return EXIT_SUCCESS;
}

unsigned int TEM_NS::check_runtime_flags(
   const input_flags& flags,
   const string& args0,
   const int& mynode,
   const int& rootnode
)
{
   unsigned int failflag;
   failflag = 0;
   // debug
   if( mynode == rootnode  && flags.debug )
   {
      cout << 
      "flags.m " << 
      flags.m << endl
      << "flags.pf " << 
      flags.pf << endl <<
      "flags.o " << 
      flags.o << endl <<
      "flags.a " << 
      flags.a << endl <<
      "flags.defocus " << 
      flags.defocus << endl <<
      "flags.spread " << 
      flags.spread << endl <<
      "flags.dupe " << 
      flags.dupe << endl <<
      "flags.image_output " << 
      flags.image_output << endl <<
      "flags.diffraction_output " << 
      flags.diffraction_output << endl <<
      "flags.adfstem_corrected " << 
      flags.adfstem_corrected << endl <<
      "flags.adfstem_uncorrected  " << 
      flags.adfstem_uncorrected  << endl <<
      "flags.bfctem_corrected  " << 
      flags.bfctem_corrected  << endl <<
      "flags.bfctem_uncorrected  " << 
      flags.bfctem_uncorrected  << endl <<
      "flags.fem  " << 
      flags.fem  << endl <<
      "flags.gt17 " << 
      flags.gt17  << endl <<
      "flags.d1 " << 
      flags.d1  << endl <<
      "flags.d2 " << 
      flags.d2  << endl <<
      "flags.d3 " << 
      flags.d3  << endl <<
      "flags.d4 " << 
      flags.d4  << endl <<
      "flags.rva " << 
      flags.rva  << endl <<
      "flags.scherzer_defocus " << 
      flags.scherzer_defocus << endl <<
      "flags.scherzer_alphamax " << 
      flags.scherzer_alphamax << endl <<
      "flags.scherzer_cs3 " << 
      flags.scherzer_cs3 << endl <<
      "flags.cs3 " << 
      flags.cs3 << endl <<
      "flags.cs5 " << 
      flags.cs5 << endl <<
      "flags.alpha_max " << 
      flags.alpha_max << endl <<
      "flags.aberration_correction " << 
      flags.aberration_correction << endl <<
      "flags.raster_spacing " << 
      flags.raster_spacing << endl <<
      "flags.correlograph " << 
      flags.correlograph << endl <<
      "flags.correlograph_variance" << 
      flags.correlograph_variance << endl <<
      "flags.correlograph_everyimage " << 
      flags.correlograph_everyimage << endl <<
      "flags.correlograph_everytxt" << 
      flags.correlograph_everytxt << endl <<
      //"flags.correlograph_everynetcdf" << 
      //flags.correlograph_everynetcdf << endl <<
      "flags.debug " << 
      flags.debug 
      << endl;
   }
   // end debug

   if (
         !(
            flags.o   // output file name prefix
            &&
            flags.a   // atom position and species input file
            &&
            flags.m
            //( // command line xor file input of microscope parameters
            //   // exclusive or:
            //   (! flags.m) != (! flags.pf) 
            //)
            &&
            // if aberration correction is used with bfctem, require 
            //    flags.spread, flags.cs3, flags.cs5
            !(
               (
                  flags.adfstem_corrected 
                  || 
                  flags.bfctem_corrected
               )
               &&
               (!
                  (
                     //flags.spread
                     //&&
                                                // defocus_spread is 
                                                // currently only used
                                                // in bfctem
                     flags.cs5
                     &&
                     (
                        flags.scherzer_cs3
                        ||
                        flags.cs3
                     )
                  )
               )
            )
            &&
            !(
               (// require specifying Cs3 if using uncorrected TEM
                  flags.adfstem_uncorrected
                  ||
                  flags.bfctem_uncorrected
               )
               && 
               (! flags.cs3)
            )
            &&
            (
               flags.scherzer_defocus 
               || 
               flags.defocus
            )
            &&
            (
               flags.scherzer_alphamax
               || 
               flags.alpha_max
            )
         )
      )
   {
      if ( mynode == rootnode )
      {
         cout << args0 << " : lacking required parameters" << endl;
         // TODO: the following two if(){} statements duplicate tests above
         //       Perhaps implement these error messages inside the above
         //       tests or vice versa.
         if (
              (flags.adfstem_corrected || flags.bfctem_corrected)
               && 
               ! (
                  //flags.spread
                  //&&
                                          // defocus_spread is 
                                          // currently only used 
                                          // in bfctem
                  (flags.cs3 || flags.scherzer_cs3)
                  && 
                  flags.cs5
               )
                
            )
         {
            cout << "Use of aberraction correction requires"
               << " specifying --defocus_spread and --cs5. If not using"
               << " --scherzer_cs3, then --cs3 is also required." << endl;
         }
         if (
               (// require specifying Cs3 if using uncorrected TEM
                  flags.adfstem_uncorrected
                  ||
                  flags.bfctem_uncorrected
               )
               && 
               (! flags.cs3)
            )
         {
            cout << "Use of uncorrected aberration"
             << " requires specifying --cs3." << endl;
         }

         if (
               !(
                  flags.scherzer_defocus 
                  || 
                  flags.defocus
               )
            )
         {
            cout << "You must specify either --defocus <value [angstrom]>"
               << " or --scherzer_defocus" << endl;
         }
         if (
               !(
                  flags.scherzer_alphamax
                  || 
                  flags.alpha_max
               )
            )
         {
            cout << "You must specify either --alphamax <value [radians]>"
               << " or --scherzer_alphamax" << endl;
         }
      }
      failflag = 1;
   }

   if (
         failflag == 0
         &&
         ! flags.adfstem_uncorrected
         &&
         ! flags.adfstem_corrected
         //&&
         //! flags.bfctem_corrected
         //&&
         //! flags.bfctem_uncorrected
      )
   {
      if ( mynode == rootnode )
      {
         cout << args0 << " : you must specify at least one of the following" << endl
            << "   [--adfstemcorrfem] simulate fluctuation microscopy using aberration corrected adfstem mode" 
            << endl 
            << "   [--adfstemuncorrfem] simulate fluctuation microscopy using adfstem mode without aberration correction" 
            << endl 
            << "   [--adfstemcorr] simulate aberration corrected adfstem" 
            << endl 
            << "   [--adfstemuncorr] simulate adfstem mode without aberration correction" 
            //<< endl 
            //<< "   [--bfctemcorr] simulate bright field TEM with aberration correction" 
            //<< endl 
            //<< "   [--bfctemuncorr] simulate bright field TEM without aberration correction" 
            << endl;
      }
      failflag = 1;
   }
   if ( flags.dupe && ! flags.fem && mynode == rootnode )
   {
      cout << "Note: --dupe flag only duplicates the sample for FTEM."
         << " The sample will not be duplicated this time." 
         << endl;
   }
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

string TEM_NS::atom_element_Z_to_abbrev( const size_t& element_Z)
{
   // Precondition : 
   //    - element_Z is an atomic number between 1 and 103
   if ( element_Z == 1 ) return "H";
   if ( element_Z == 2 ) return "He";
   if ( element_Z == 3 ) return "Li";
   if ( element_Z == 4 ) return "Be";
   if ( element_Z == 5 ) return "B";
   if ( element_Z == 6 ) return "C";
   if ( element_Z == 7 ) return "N";
   if ( element_Z == 8 ) return "O";
   if ( element_Z == 9 ) return "F";
   if ( element_Z == 10 ) return "Ne";
   if ( element_Z == 11 ) return "Na";
   if ( element_Z == 12 ) return "Mg";
   if ( element_Z == 13 ) return "Al";
   if ( element_Z == 14 ) return "Si";
   if ( element_Z == 15 ) return "P";
   if ( element_Z == 16 ) return "S";
   if ( element_Z == 17 ) return "Cl";
   if ( element_Z == 18 ) return "Ar";
   if ( element_Z == 19 ) return "K";
   if ( element_Z == 20 ) return "Ca";
   if ( element_Z == 21 ) return "Sc";
   if ( element_Z == 22 ) return "Ti";
   if ( element_Z == 23 ) return "V";
   if ( element_Z == 24 ) return "Cr";
   if ( element_Z == 25 ) return "Mn";
   if ( element_Z == 26 ) return "Fe";
   if ( element_Z == 27 ) return "Co";
   if ( element_Z == 28 ) return "Ni";
   if ( element_Z == 29 ) return "Cu";
   if ( element_Z == 30 ) return "Zn";
   if ( element_Z == 31 ) return "Ga";
   if ( element_Z == 32 ) return "Ge";
   if ( element_Z == 33 ) return "As";
   if ( element_Z == 34 ) return "Se";
   if ( element_Z == 35 ) return "Br";
   if ( element_Z == 36 ) return "Kr";
   if ( element_Z == 37 ) return "Rb";
   if ( element_Z == 38 ) return "Sr";
   if ( element_Z == 39 ) return "Y";
   if ( element_Z == 40 ) return "Zr";
   if ( element_Z == 41 ) return "Nb";
   if ( element_Z == 42 ) return "Mo";
   if ( element_Z == 43 ) return "Tc";
   if ( element_Z == 44 ) return "Ru";
   if ( element_Z == 45 ) return "Rh";
   if ( element_Z == 46 ) return "Pd";
   if ( element_Z == 47 ) return "Zg";
   if ( element_Z == 48 ) return "Cd";
   if ( element_Z == 49 ) return "In";
   if ( element_Z == 50 ) return "Sn";
   if ( element_Z == 51 ) return "Sb";
   if ( element_Z == 52 ) return "Te";
   if ( element_Z == 53 ) return "I";
   if ( element_Z == 54 ) return "Xe";
   if ( element_Z == 55 ) return "Cs";
   if ( element_Z == 56 ) return "Ba";
   if ( element_Z == 57 ) return "La";
   if ( element_Z == 58 ) return "Ce";
   if ( element_Z == 59 ) return "Pr";
   if ( element_Z == 60 ) return "Nd";
   if ( element_Z == 61 ) return "Pm";
   if ( element_Z == 62 ) return "Sm";
   if ( element_Z == 63 ) return "Eu";
   if ( element_Z == 64 ) return "Gd";
   if ( element_Z == 65 ) return "Tb";
   if ( element_Z == 66 ) return "Dy";
   if ( element_Z == 67 ) return "Ho";
   if ( element_Z == 68 ) return "Er";
   if ( element_Z == 69 ) return "Tm";
   if ( element_Z == 70 ) return "Yb";
   if ( element_Z == 71 ) return "Lu";
   if ( element_Z == 72 ) return "Hf";
   if ( element_Z == 73 ) return "Ta";
   if ( element_Z == 74 ) return "W";
   if ( element_Z == 75 ) return "Re";
   if ( element_Z == 76 ) return "Os";
   if ( element_Z == 77 ) return "Ir";
   if ( element_Z == 78 ) return "Pt";
   if ( element_Z == 79 ) return "Au";
   if ( element_Z == 80 ) return "Hg";
   if ( element_Z == 81 ) return "Tl";
   if ( element_Z == 82 ) return "Pb";
   if ( element_Z == 83 ) return "Bi";
   if ( element_Z == 84 ) return "Po";
   if ( element_Z == 85 ) return "At";
   if ( element_Z == 86 ) return "Rn";
   if ( element_Z == 87 ) return "Fr";
   if ( element_Z == 88 ) return "Ra";
   if ( element_Z == 89 ) return "Ac";
   if ( element_Z == 90 ) return "Th";
   if ( element_Z == 91 ) return "Pa";
   if ( element_Z == 92 ) return "U";
   if ( element_Z == 93 ) return "Np";
   if ( element_Z == 94 ) return "Pu";
   if ( element_Z == 95 ) return "Am";
   if ( element_Z == 96 ) return "Cm";
   if ( element_Z == 97 ) return "Bk";
   if ( element_Z == 98 ) return "Cf";
   if ( element_Z == 99 ) return "Es";
   if ( element_Z == 100 ) return "Fm";
   if ( element_Z == 101 ) return "Md";
   if ( element_Z == 102 ) return "No";
   if ( element_Z == 103 ) return "Lr";
   
   // If none of the modeled elements matched, return the invalid
   //  atomic number 0 .
   return 0;
}


#endif
