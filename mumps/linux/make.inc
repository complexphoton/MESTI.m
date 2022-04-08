# It is possible to generate a MATLAB or an Octave interface thanks to
# the Octave MEX file compatibility. Comment/uncomment the lines below
# depending on whether you want to generate the MATLAB or the Octave
# interface

# To generate the MATLAB interface uncomment the following line
# ( the use of -largeArrayDims is necessary to work with sparse
# matrices since R2006b)
MEX         = mex -g -largeArrayDims

# To generate the Octave interface uncomment the following line
# MEX         = mkoctfile -g --mex

# Main MUMPS_DIR
MUMPS_DIR = /project/cwhsu_38/shared/software/mumps-5.4.1-seq_wo_metis

# Orderings (see main Makefile.inc file from MUMPS)
# If METIS has been installed, you can uncomment the following two lines to specify
# the directory and library of METIS.
#LMETISDIR = /project/cwhsu_38/shared/software/metis-5.1.0/build/Linux-x86_64/libmetis
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
LIBFORT = -lgfortran
#LIBFORT = /usr/lib/gcc/x86_64-linux-gnu/4.7/libgfortran.so

# BLAS library:
# LIBBLAS = -L/usr/lib/atlas -lblas
# LIBBLAS = -lsunperf -lf77compat
# LIBBLAS = -lblas
LIBBLAS = /spack/apps/linux-centos7-x86_64/intel-19.0.4/intel-mkl-2019.5.281-ftihihwas7hnpmex6tt3luvqzuy6fgb7/compilers_and_libraries_2019.5.281/linux/mkl/lib/intel64_lin/libmkl_intel_lp64.so
#-lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core
# LIBBLAS = /home/jylexcel/libs_courge/libblas.a  # -fPIC missing
# LIBBLAS = /usr/lib/libblas.so

# extra options passed via mex command
# Add -DINTSIZE64 if MUMPS was compiled with 64-bit integers (BLAS
# library should then have 64-bit integers as well)
OPTC    = -g
