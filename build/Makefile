ROOT = ..

#############################
# Hoffman2 configuration
#############################
#CXX		= mpiicpc
#INCLUDE_NETCDF = -I$(NETCDF_HOME)/include 
#INCLUDE_FFTW = -I/u/local/apps/fftw3/3.3.6-pl2_intel_17_mpi/include 
#CXXFLAGS	=  -I$(ROOT)/include $(INCLUDE_NETCDF) $(INCLUDE_FFTW) -fp-model precise -fimf-arch-consistency=true
#LIBS_NETCDF = -L$(NETCDF_HOME)/lib -lnetcdf_c++ -lnetcdf -lnetcdf_c++ # -lm
#LIBS_FFTW = -L/u/local/apps/fftw3/3.3.6-pl2_intel_17_mpi/lib -Wl,-rpath,/u/local/apps/fftw3/3.3.6-pl2_intel_17_mpi/lib -lfftw3_mpi -lfftw3 
##LIBS		=  -L$(ROOT)/lib -L/usr/lib64/ $(LIBS_NETCDF) $(LIBS_FFTW) -ltiff -lm 
#LIBS		=  -L$(ROOT)/lib -L/usr/lib64/ $(LIBS_FFTW) -ltiff -lm 

#############################
# gcc + openmpi configuration
#############################
CXX = mpicxx
CXXFLAGS = -I$(ROOT)/include
#LIBS_NETCDF	= -lnetcdf_c++ -lnetcdf -lnetcdf_c++
LIBS_FFTW = -lfftw3_mpi -lfftw3
#LIBS_LMP = -L$(ROOT)/<put path to lammps library here> -llammps
LIBS_LMP = -llammps # to use this, copy all lammps header files to include/
#LIBS = -L$(ROOT)/lib -L/usr/lib64 $(LIBS_NETCDF) $(LIBS_FFTW) -ltiff -lm 
LIBS = -L$(ROOT)/lib -L/usr/lib64 $(LIBS_FFTW) $(LIBS_LMP) -ltiff -lm 
#############################

OBJS = func.o projected_atomic_potential.o write_tif.o scherzer.o bw_limit.o lenses.o tem_io.o probe_wavefunction.o adfstem.o slice.o variance.o io_txt.o 
#io_netcdf.o 
SRCS = $(ROOT)/src/func.cpp $(ROOT)/src/projected_atomic_potential.cpp $(ROOT)/src/write_tif.cpp $(ROOT)/src/scherzer.cpp $(ROOT)/src/bw_limit.cpp  $(ROOT)/src/lenses.cpp $(ROOT)/src/tem_io.cpp $(ROOT)/src/slice.cpp $(ROOT)/src/probe_wavefunction.cpp $(ROOT)/src/adfstem.cpp $(ROOT)/src/variance.cpp $(ROOT)/src/io_txt.cpp 
#$(ROOT)/src/io_netcdf.cpp 

MSSTEMFEM = $(ROOT)/bin/ms-stem-fem.x

MSSTEMFEMMD = $(ROOT)/bin/ms-stem-fem-md.x

$(MSSTEMFEM):	$(OBJS)
	$(CXX) -g $(ROOT)/src/tem.cpp $(OBJS) $(CXXFLAGS) $(LIBS) -o $(MSSTEMFEM)


$(MSSTEMFEMMD):	$(OBJS)
	$(CXX) -g $(ROOT)/src/ftem_lammps.cpp $(OBJS) $(CXXFLAGS) $(LIBS) -o $(MSSTEMFEM)

#$(CXX) -g $(ROOT)/src/tem.cpp $(OBJS) $(CXXFLAGS) $(LIBS) -o $(MSSTEMFEM)

$(OBJS):			
	$(CXX) -g -c $(CXXFLAGS) $(SRCS)

ms-stem-fem: $(MSSTEMFEM)

ms-stem-fem-md: $(MSSTEMFEMMD)

clean: 
	rm -f $(OBJS) $(MSSTEMFEM) $(MSSTEMFEMMD)
