## Prerequisite 

Users will need to do the following prepares for compiling MUMPS and running MATLAB interface for MUMPS on Linux cluster.

### Compilers

The compilation of MUMPS requires both C and Fortran compilers. Both C and Fortran compliers should be included in GNU Compiler Collection (GCC), which should be installed on Linux cluster. In Lmod module system, it can be done by typing  

```shell
module load gcc
```

### BLAS and LAPACK

MUMPS requires both BLAS and LAPACK libraries, which are standard libraries on Linux cluster. These libraries are also included in many implementations, such as Intel MKL and OpenBLAS. For example, if Intel MKL is included on the Linux cluster, in Lmod module system, it can be done by typing 

```shell
module load intel
```

### Running MATLAB interface for MUMPS

In some cases, system cannot find the BLAS and LAPACK libraries by itself when users run MATLAB interface for MUMPS. To solve this issue, users can `export LD_PRELOAD` before running MATLAB. For example, if BLAS and LAPACK libraries are used through Intel MKL, users can typing,

```shell
export LD_PRELOAD=$MKLROOT/libmkl_intel_lp64.so:$MKLROOT/libmkl_sequential.so:$MKLROOT/libmkl_intel_thread.so:$MKLROOT/libmkl_core.so
```

Note that `MKLROOT` is the path for Intel MKL libraries and exact libraries to be preload may be different from Intel MKL versions and cluster systems.