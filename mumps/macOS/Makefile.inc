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
#   look for the approriate file in the directory ./Make.inc/ and copy it
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
#
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

# Uncomment the two lines below if you installed METIS to the default /usr/local folder
#LMETISDIR = /usr/local/lib
#IMETIS    = -I/usr/local/include

# You have to choose one among the following two lines depending on
# the type of analysis you want to perform. If you want to perform only
# sequential analysis choose the first (remember to add -Dmetis in the ORDERINGSF
# variable below); for both parallel and sequential analysis choose the second 
# line (remember to add -Dparmetis in the ORDERINGSF variable below)

# Uncomment the following line if you installed METIS
#LMETIS    = -L$(LMETISDIR) -lmetis

# The following variables will be used in the compilation process.
# Please note that -Dptscotch and -Dparmetis imply -Dscotch and -Dmetis respectively.
# If you want to use Metis 4.X or an older version, you should use -Dmetis4 instead of -Dmetis
# or in addition with -Dparmetis (if you are using parmetis 3.X or older).
#ORDERINGSF = -Dscotch -Dmetis -Dpord -Dptscotch -Dparmetis

# Add -Dmetis to the following line if you installed METIS
ORDERINGSF = -Dpord
ORDERINGSC  = $(ORDERINGSF)

LORDERINGS = $(LMETIS) $(LPORD) $(LSCOTCH)
IORDERINGSF = $(ISCOTCH)
IORDERINGSC = $(IMETIS) $(IPORD) $(ISCOTCH)

#End orderings
########################################################################

########################################################################
# DEFINE HERE SOME COMMON COMMANDS, THE COMPILER NAMES, ETC...

# PLAT : use it to add a default suffix to the generated libraries
PLAT    = 
# Library extension, + C and Fortran "-o" option
# may be different under Windows
LIBEXT  = .a
OUTC    = -o 
OUTF    = -o 
# RM : remove files
RM      = /bin/rm -f
# CC : C compiler
CC      = clang
# FC : Fortran 90 compiler
FC      = gfortran
# FL : Fortran linker
FL      = gfortran
# AR : Archive object in a library
#      keep a space at the end if options have to be separated from lib name
AR      = ar vr 
# RANLIB : generate index of an archive file
#   (optionnal use "RANLIB = echo" in case of problem)
RANLIB  = ranlib
#RANLIB  = echo

# DEFINE HERE YOUR LAPACK LIBRARY

# Use LAPACK in Apple's Accelerate framework; alternately we can set "LAPACK = -lLAPACK".
LAPACK = -framework Accelerate

# The next two lines should not be modified. They concern
# the sequential library provided by MUMPS, to use instead
# of ScaLAPACK and MPI.
INCSEQ  = -I$(topdir)/libseq
LIBSEQ  = $(LAPACK) -L$(topdir)/libseq -lmpiseq

# DEFINE HERE YOUR BLAS LIBRARY

# Use the following to use BLAS in Apple's Accelerate framework; alternately we can set "LIBBLAS = -lBLAS".
LIBBLAS = -framework Accelerate

# Use the following if you have an Apple silicon Mac, have openblas installed through "brew install openblas", and want to use that.
#LIBBLAS = -L/opt/homebrew/opt/openblas/lib

# Use the following if you have an Intel Mac, have openblas installed through "brew install openblas", and want to use that.
#LIBBLAS = -L/usr/local/opt/openblas/lib

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

CDEFS = -DAdd_

#COMPILER OPTIONS
# The -fallow-argument-mismatch is needed for gfortran 10.0 and later; remove that if you use an earlier version of gfortran
OPTF    = -O -fPIC -fallow-argument-mismatch
OPTC    = -O -I. -fPIC
OPTL    = -O -fPIC

#Sequential:

INCS = $(INCSEQ)
LIBS = $(LIBSEQ)
LIBSEQNEEDED = libseqneeded
