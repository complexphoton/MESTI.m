# It is possible to generate a MATLAB or an Octave interface thanks to
# the Octave MEX file compatibility. Comment/uncomment the lines below
# depending on whether you want to generate the MATLAB or the Octave
# interface

# To generate the MATLAB interface uncomment the following line
# ( the use of -largeArrayDims is necessary to work with sparse
# matrices since R2006b)
#MEX         = /opt/matlab/matlab/bin/mex -g -largeArrayDims
MEX         = /Applications/MATLAB_R2018a.app/bin/mex -g -largeArrayDims

# To generate the Octave interface uncomment the following line
# MEX         = mkoctfile -g --mex

# Main MUMPS_DIR
#MUMPS_DIR = $(HOME)/MUMPS_5.1.2
MUMPS_DIR = /Users/chsu/install/compile/MUMPS_5.1.2

# Orderings (see main Makefile.inc file from MUMPS)
#LMETISDIR = ${HOME}/parmetis-4.0.3/build/Linux-x86_64/libmetis
LMETISDIR = /Users/chsu/install/compile/metis-5.1.0/build/Darwin-x86_64/libmetis
LMETIS     = -L$(LMETISDIR) -lmetis
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
#  LIBFORT = /usr/lib/gcc/x86_64-linux-gnu/4.7/libgfortran.so
#LIBFORT = /usr/local/opt/gcc/lib/gcc/8/libgfortran.a
#LIBFORT = -L/usr/local/opt/gcc/lib/gcc/8
#LIBFORT = -L/usr/local/opt/gcc/lib
LIBFORT = /usr/local/Cellar/gcc/8.2.0/lib/gcc/8/libgfortran.dylib

# BLAS library:
# LIBBLAS = -L/usr/lib/atlas -lblas
# LIBBLAS = -lsunperf -lf77compat
#LIBBLAS = -lblas
# LIBBLAS = /home/jylexcel/libs_courge/libgoto_opteronp-r1.26.a
# LIBBLAS = /usr/local/opt/openblas/lib/libopenblas.a
#LIBBLAS = -L/usr/local/opt/openblas/lib
LIBBLAS = -L/usr/local/opt/openblas/lib -llapack
# LIBBLAS = /home/jylexcel/libs_courge/libblas.a  # -fPIC missing
# LIBBLAS = /usr/lib/libblas.so

#LAPACK = -llapack
#LAPACK = -L/usr/local/opt/lapack/lib

# extra options passed via mex command
# Add -DINTSIZE64 if MUMPS was compiled with 64-bit integers (BLAS
# library should then have 64-bit integers as well)
OPTC    = -g

