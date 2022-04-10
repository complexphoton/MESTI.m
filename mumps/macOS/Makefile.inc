#
#  This file is part of MUMPS 5.2.1, released
#  on Fri Jun 14 14:46:05 UTC 2019
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
#SCOTCHDIR  = /Users/chsu/install/compile/scotch_6.0.6
#ISCOTCH    = -I$(SCOTCHDIR)/include
#ISCOTCH    = -I/usr/local/include

# You have to choose one among the following two lines depending on
# the type of analysis you want to perform. If you want to perform only
# sequential analysis choose the first (remember to add -Dscotch in the ORDERINGSF
# variable below); for both parallel and sequential analysis choose the second 
# line (remember to add -Dptscotch in the ORDERINGSF variable below)

#LSCOTCH    = -L$(SCOTCHDIR)/lib -lesmumps -lscotch -lscotcherr
#LSCOTCH    = -L$(SCOTCHDIR)/lib -lptesmumps -lptscotch -lptscotcherr
#LSCOTCH    = -L/usr/local/lib/libscotch.a -lesmumps -lscotch -lscotcherr
#LSCOTCH    = -L/usr/local/lib -lesmumps -lscotch -lscotcherr
#LSCOTCH    = -L/Users/chsu/install/compile/scotch_6.0.6/lib -lesmumps -lscotch -lscotcherr


LPORDDIR = $(topdir)/PORD/lib/
IPORD    = -I$(topdir)/PORD/include/
LPORD    = -L$(LPORDDIR) -lpord

#LMETISDIR = /opt/metis-5.1.0/build/Linux-x86_64/libmetis
#LMETISDIR = /Users/chsu/install/compile/metis-5.1.0/build/Darwin-x86_64/libmetis
#IMETIS    = /opt/metis-5.1.0/include
#LMETISDIR = /usr/local/lib/libmetis.a
#LMETISDIR = /usr/local/lib
#IMETIS    = -I/usr/local/include
#IMETIS    = -I/Users/chsu/install/compile/metis-5.1.0/include

# You have to choose one among the following two lines depending on
# the type of analysis you want to perform. If you want to perform only
# sequential analysis choose the first (remember to add -Dmetis in the ORDERINGSF
# variable below); for both parallel and sequential analysis choose the second 
# line (remember to add -Dparmetis in the ORDERINGSF variable below)

#LMETIS    = -L$(LMETISDIR) -lmetis
#LMETIS    = -L$(LMETISDIR) -lparmetis -lmetis

# The following variables will be used in the compilation process.
# Please note that -Dptscotch and -Dparmetis imply -Dscotch and -Dmetis respectively.
# If you want to use Metis 4.X or an older version, you should use -Dmetis4 instead of -Dmetis
# or in addition with -Dparmetis (if you are using parmetis 3.X or older).
#ORDERINGSF = -Dscotch -Dmetis -Dpord -Dptscotch -Dparmetis
#ORDERINGSF = -Dmetis -Dpord # add -Dmetis if you have Metis installed.
ORDERINGSF  = -Dpord 
ORDERINGSC  = $(ORDERINGSF)

LORDERINGS = $(LMETIS) $(LPORD)

#IORDERINGSF = $(ISCOTCH)
IORDERINGSF = $(IMETIS) $(IPORD)

IORDERINGSC = $(IMETIS) $(IPORD)

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
CC      = cc -fPIC
# FC : Fortran 90 compiler
FC      = gfortran -fallow-argument-mismatch -fPIC
# FL : Fortran linker
FL      = gfortran -fPIC
# AR : Archive object in a library
#      keep a space at the end if options have to be separated from lib name
AR      = ar vr 
# RANLIB : generate index of an archive file
#   (optionnal use "RANLIB = echo" in case of problem)
RANLIB  = ranlib
#RANLIB  = echo

# DEFINE HERE YOUR LAPACK LIBRARY

LAPACK = -lLAPACK
#LAPACK = -llapack
#LAPACK = -L/usr/local/opt/lapack/lib

# The next two lines should not be modified. They concern
# the sequential library provided by MUMPS, to use instead
# of ScaLAPACK and MPI.
INCSEQ  = -I$(topdir)/libseq
LIBSEQ  = $(LAPACK) -L$(topdir)/libseq -lmpiseq

# DEFINE HERE YOUR BLAS LIBRARY

LIBBLAS = -lBLAS 
#LIBBLAS = -L/usr/local/opt/openblas/lib # on intel-based systems, openblas is faster.

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
OPTF    = -O
# changed OPTS to below to disable pthread
#OPTC    = -O -I. -DWITHOUT_PTHREAD
OPTC    = -O -I.
OPTL    = -O 

#Sequential:

INCS = $(INCSEQ)
LIBS = $(LIBSEQ)
LIBSEQNEEDED = libseqneeded
