# It is possible to generate a MATLAB or an Octave interface thanks to
# the Octave MEX file compatibility. Comment/uncomment the lines below
# depending on whether you want to generate the MATLAB or the Octave
# interface

# Location of the MATLAB mex compiler; change it based on the MATLAB version on your machine.
# For Apple silicon Macs, currently R2022a_Beta is the only compatible version.
# ( the use of -largeArrayDims is necessary to work with sparse
# matrices since R2006b)
MEX = /Applications/MATLAB_R2022a_Beta.app/bin/mex -g -largeArrayDims

# To generate the Octave interface uncomment the following line
# MEX         = mkoctfile -g --mex

# Change to the directory where you compiled MUMPS
MUMPS_DIR = /Users/chsu/install/compile/MUMPS_5.5.0

# Orderings (see main Makefile.inc file from MUMPS)
# Uncomment the two lines below if you installed METIS to the default /usr/local folder
#LMETISDIR = /usr/local/lib
#LMETIS     = -L$(LMETISDIR) -lmetis
LPORDDIR   = $(MUMPS_DIR)/PORD/lib
LPORD      = -L$(LPORDDIR) -lpord
LORDERINGS = $(LPORD) $(LMETIS)

# Fortran runtime library
# Please find out the path and name of your
# Fortran runtime, examples below:
# g95:
# LIBFORT = /usr/lib/libf95.a /usr/lib/libgcc.a 
# Intel:
# LIBFORT = /opt/intel80/lib/libifcore.a /opt/intel80/lib/libifport.a /opt/intel80/lib/libirc.a
# PGI:
# LIBFORT = -L/usr/local/pgi/linux86/5.2/lib -llapack -lblas -lpgf90 -lpgc -lpgf90rtl -lpgftnrtl -lpgf902 -lpgf90_rpm1 -lpghpf2
# SGI 32-bit
# LIBFORT = -L/usr/lib32 -lblas -L/usr/lib32/mips4 -lfortran
# Sun
# LIBFORT = -L/opt2/SUNWspro7/lib -lsunperf -lfminvai -lfai2 -lfsu -lfmaxvai -lfmaxlai -lfai -lfsumai -lLIBFORT = /usr/local/lib/libgfortran.a

# gfortran
# Use the following if you have an Apple silicon Mac and have gfortran installed through "brew install gcc"; change to your gcc version.
LIBFORT = -L/opt/homebrew/lib/gcc/12 -lgfortran

# Use the following if you have an Intel Mac and have gfortran installed through "brew install gcc"; change to your gcc version.
#LIBFORT = -L/usr/local/opt/gcc/lib/gcc/12 -lgfortran

# BLAS library:
# Use the following to link to BLAS and LAPACK in Apple's Accelerate framework.
# Note that mex does not support "-framework", so instead of "-framework Accelerate", here we use "-lBLAS -lLAPACK", which will search for BLAS and LAPACK (and the one in in Apple's Accelerate framework will be found).
LIBBLAS = -lBLAS -lLAPACK

# Use the following if you have an Apple silicon Mac, have openblas installed through "brew install openblas", and want to use that.
#LIBBLAS = -L/opt/homebrew/opt/openblas/lib -lLAPACK

# Use the following if you have an Intel Mac, have openblas installed through "brew install openblas", and want to use that.
#LIBBLAS = -L/usr/local/opt/openblas/lib -lLAPACK

# extra options passed via mex command
# Add -DINTSIZE64 if MUMPS was compiled with 64-bit integers (BLAS
# library should then have 64-bit integers as well)
OPTC    = -g

