#
#  This file is part of MUMPS 5.5.0, released
#  on Thu Apr 14 11:45:33 UTC 2022
#
################################################################################
#
#   Makefile.inc.generic.SEQ
#
#    Generic Makefile.inc for sequential (MPI free, Scalapack free) version
#
#
#   This defines some parameters dependent on your platform; you should
#   look for the appropriate file in the directory ./Make.inc/ and copy it
#   into a file called Makefile.inc. For example, from the MUMPS root
#   directory, use 
#   "cp Make.inc/Makefile.inc.generic.SEQ ./Makefile.inc"
#   (see the main README file for details)
#
#   If you do not find any suitable Makefile in Makefile.inc, use this file:
#   "cp Make.inc/Makefile.inc.generic ./Makefile.inc" and modify it according
#   to the comments given below. If you manage to build MUMPS on a new platform,
#   and think that this could be useful to others, you may want to send us
#   the corresponding Makefile.inc file.
#
################################################################################


########################################################################
#Begin orderings

# NOTE that PORD is distributed within MUMPS by default. It is recommended to
# install other orderings. For that, you need to obtain the corresponding package
# and modify the variables below accordingly.
# For example, to have Metis available within MUMPS:
#          1/ download Metis and compile it
#          2/ uncomment (suppress # in first column) lines
#             starting with LMETISDIR,  LMETIS
#          3/ add -Dmetis in line ORDERINGSF
#             ORDERINGSF  = -Dpord -Dmetis
#          4/ Compile and install MUMPS
#             make clean; make   (to clean up previous installation)
#
#          Metis/ParMetis and SCOTCH/PT-SCOTCH (ver 6.0 and later) orderings are recommended.
#

#SCOTCHDIR  = ${HOME}/scotch_6.0
#ISCOTCH    = -I$(SCOTCHDIR)/include
#
# You have to choose one among the following two lines depending on
# the type of analysis you want to perform. If you want to perform only
# sequential analysis choose the first (remember to add -Dscotch in the ORDERINGSF
# variable below); for both parallel and sequential analysis choose the second 
# line (remember to add -Dptscotch in the ORDERINGSF variable below)

#LSCOTCH    = -L$(SCOTCHDIR)/lib -lesmumps -lscotch -lscotcherr
#LSCOTCH    = -L$(SCOTCHDIR)/lib -lptesmumps -lptscotch -lptscotcherr


LPORDDIR = $(topdir)/PORD/lib/
IPORD    = -I$(topdir)/PORD/include/
LPORD    = -L$(LPORDDIR) -lpord

# If METIS has been installed, you can uncomment the following two lines to specify 
# your own path of METIS.
# The following is the path of the METIS on Hsu's group in USC Discovery cluster:
#LMETISDIR = /project/cwhsu_38/shared/software/metis-5.1.0/build/Linux-x86_64/libmetis
#IMETIS    = -I/project/cwhsu_38/shared/software/metis-5.1.0/include

# You have to choose one among the following two lines depending on
# the type of analysis you want to perform. If you want to perform only
# sequential analysis choose the first (remember to add -Dmetis in the ORDERINGSF
# variable below); for both parallel and sequential analysis choose the second 
# line (remember to add -Dparmetis in the ORDERINGSF variable below)

# If METIS has been installed, you can uncomment the following line to specify 
# the library of METIS.
#LMETIS    = -L$(LMETISDIR) -lmetis

# The following variables will be used in the compilation process.
# Please note that -Dptscotch and -Dparmetis imply -Dscotch and -Dmetis respectively.
# If you want to use Metis 4.X or an older version, you should use -Dmetis4 instead of -Dmetis
# or in addition with -Dparmetis (if you are using parmetis 3.X or older).
#ORDERINGSF = -Dscotch -Dmetis -Dpord -Dptscotch -Dparmetis

# Add -Dmetis to the following line if you installed METIS
ORDERINGSF  = -Dpord 
ORDERINGSC  = $(ORDERINGSF)

LORDERINGS = $(LMETIS) $(LPORD) $(LSCOTCH)
IORDERINGSF = $(ISCOTCH)
IORDERINGSC = $(IMETIS) $(IPORD) $(ISCOTCH)

#End orderings
########################################################################

########################################################################
# DEFINE HERE SOME COMMON COMMANDS, THE COMPILER NAMES, ETC...

# PLAT : use it to add a default suffix to the generated libraries
PLAT =
# Library extension, + C and Fortran "-o" option
# may be different under Windows
LIBEXT  = .a
LIBEXT_SHARED  = .so
FPIC_OPT = -fPIC
OUTC    = -o 
OUTF    = -o 
# RM : remove files
RM      = /bin/rm -f
# CC : C compiler
CC      = gcc
# FC : Fortran 90 compiler
FC = gfortran
# FL : Fortran linker
FL = gfortran
# AR : Archive object in a library
#      keep a space at the end if options have to be separated from lib name
AR = ar vr 
# RANLIB : generate index of an archive file
#   (optional use "RANLIB = echo" in case of problem)
RANLIB  = ranlib
#RANLIB  = echo

# DEFINE HERE YOUR LAPACK LIBRARY

# Make MKLROOT variable point to the path where the Intel MKL is
# installed. Default install directory for Intel MKL is following:
# MKLROOT=/opt/intel/mkl
# However, in different Linux clusters, the Intel MKL may be put in different
# path. If you are using in Lmod module system and Intel MKL are installed,
# you can use "module load intel", "module load intel-mkl", and "echo $MKLROOT".
# MKLROOT path should be printed out, if Intel MKL is available.
# The MKLROOT should be under the path ".../mkl", where "..." is the Intel path.
#
# If Lmod module system is not used in your cluster, then you can also try to
# find the path like ".../mkl", where "..." is the Intel directory path.
# You can specify the MKLROOT by "source .../mkl/bin/mklvars.sh intel64"
# and type "echo $MKLROOT". MKLROOT path should be printed out.
#
# If you still cannot figure out the path of the libraries, you may need to
# contact the maintenance team of your cluster.
#
# The following is the path of the Intel MKL path in USC Discovery cluster:
MKLROOT=/spack/apps/linux-centos7-x86_64/intel-19.0.4/intel-mkl-2019.5.281-ftihihwas7hnpmex6tt3luvqzuy6fgb7/compilers_and_libraries_2019.5.281/linux/mkl
# LAPACK library is under the folder of $(MKLROOT)/lib/intel64
LAPACK = -L$(MKLROOT)/lib/intel64 -lmkl_sequential -lmkl_intel_lp64 -lmkl_core

# The next two lines should not be modified. They concern
# the sequential library provided by MUMPS, to use instead
# of ScaLAPACK and MPI.
INCSEQ = -I$(topdir)/libseq
LIBSEQ  = $(LAPACK) -L$(topdir)/libseq -lmpiseq

# DEFINE HERE YOUR BLAS LIBRARY
# BLAS library is under the folder of $(MKLROOT)/lib/intel64
LIBBLAS = -L$(MKLROOT)/lib/intel64 -lmkl_intel_lp64

# DEFINE HERE YOUR PTHREAD LIBRARY
LIBOTHERS = -lpthread

# FORTRAN/C COMPATIBILITY:
#  Use:
#    -DAdd_ if your Fortran compiler adds an underscore at the end
#              of symbols,
#     -DAdd__ if your Fortran compiler adds 2 underscores,
#
#     -DUPPER if your Fortran compiler uses uppercase symbols
#
#     leave empty if your Fortran compiler does not change the symbols.
#

CDEFS   = -DAdd_

#COMPILER OPTIONS
# The -fallow-argument-mismatch is needed for gfortran 10.0 and later for OPTF.
# If you use gfortran 10.0 and later, uncomment -fallow-argument-mismatch.
OPTF    = -O -fPIC #-fallow-argument-mismatch 
OPTL    = -O -I. -fPIC
OPTC    = -O -fPIC 

#Sequential:

INCS = $(INCSEQ)
LIBS = $(LIBSEQ)
LIBSEQNEEDED = libseqneeded
