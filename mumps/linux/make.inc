# It is possible to generate a MATLAB or an Octave interface thanks to
# the Octave MEX file compatibility. Comment/uncomment the lines below
# depending on whether you want to generate the MATLAB or the Octave
# interface

# Location of the MATLAB mex compiler; change it based on the MATLAB version on your 
# machine. For example, you may change the line to 
# Mex         = /spack/apps/matlab/2020b/bin/mex -g -largeArrayDims
# If you are using in Lmod module system and MATLAB is installed, you can use
# "module load matlab". You do not need to change the following line.
# ( the use of -largeArrayDims is necessary to work with sparse
# matrices since R2006b)
MEX         = mex -g -largeArrayDims

# To generate the Octave interface uncomment the following line
# MEX         = mkoctfile -g --mex

# Change to the directory where you compiled MUMPS
MUMPS_DIR = /project/cwhsu_38/shared/software/mumps-5.5.0-seq_wo_metis

# Orderings (see main Makefile.inc file from MUMPS)
# If METIS has been installed, you can uncomment the following two lines to specify
# the directory and library of METIS.
# The following is the path of the METIS on Hsu's group in USC Discovery cluster:
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
# Make this variable point to the path where the gfortran library is
LIBFORT = /spack/apps/gcc/8.3.0/lib64/libgfortran.so

# BLAS library:
# Make this variable point to the path where the Intel MKL library is
# installed. Default install directory for Intel MKL is following:
# MKLROOT=/opt/intel/mkl/lib/intel64
# However, in different Linux clusters, the Intel MKL library may be put in different
# path. If you are using in Lmod module system and Intel libraries are installed,
# you can use "module load intel" and "printenv PATH". One of the PATH should be the
# path of Intel directory, then typically the MKLROOT is under the path
# ".../mkl/lib/intel64", where "..." is the Intel directory.
# If Lmod module system is not used in your cluster, then you can also try to find the
# path like ".../mkl/lib/intel64", where "..." is the Intel directory.
# If you still cannot figure out the path of the libraries, you may need to contact
# the maintenance team of your cluster.
# The following is the path of the Intel MKL libraries in USC Discovery cluster:
MKLROOT=/spack/apps/linux-centos7-x86_64/intel-19.0.4/intel-mkl-2019.5.281-ftihihwas7hnpmex6tt3luvqzuy6fgb7/compilers_and_libraries_2019.5.281/linux/mkl/lib/intel64
LIBBLAS = -L$(MKLROOT) -lmkl_intel_lp64

# extra options passed via mex command
# Add -DINTSIZE64 if MUMPS was compiled with 64-bit integers (BLAS
# library should then have 64-bit integers as well)
OPTC    = -g
