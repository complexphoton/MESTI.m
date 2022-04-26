## Prerequisites for MUMPS installation on a Linux cluster

We need the following tools before compiling MUMPS and its MATLAB interface in the Linux cluster.

### MATLAB 

The MATLAB should be installed before compiling MUMPS-MATLAB interface. The MATLAB compiler `mex` is required when compile the MUMPS-MATLAB interface. 

If MATLAB is already installed in the cluster, in Lmod module system, you can type  
```shell
module load matlab
```

If MATLAB is already installed in the cluster, but it does not used Lmod module system, you may add the path of your MATLAB to `PATH` by
```shell
export PATH=".../matlab/bin/:$PATH"
```
where `...` is the path to MATLAB in the cluster.

### GNU compilers collection

The compilation of MUMPS requires both C and Fortran compilers. Both C and Fortran compliers should be included in GNU Compiler Collection (GCC), which should be installed in the cluster. 

It typically that users in clusters have default GCC, you can check it by 
```shell
gcc -v
```
Then you can see the details of your GCC.

If you do not have default GCC, in Lmod module system, it can be done by typing  
```shell
module load gcc
```
In non-Lmod module system, you may type,
```shell
.../gcc/setup.sh
```
where `...` is the path to GCC in the cluster.

### BLAS and LAPACK

MUMPS requires both BLAS and LAPACK libraries, which are standard libraries on Linux cluster. These libraries are also included in many implementations, such as Intel MKL and OpenBLAS. For example, if Intel MKL is included on the Linux cluster, in Lmod module system, it can be done by typing 

```shell
module load intel
```

If Lmod module system is not used in your cluster, then you can try to find the path like `.../mkl/lib/intel64`, where `...` is the Intel directory.

If you still cannot figure out the path of the libraries, you may need to contact the maintenance team of your cluster.

### Running MATLAB interface for MUMPS

In some cases, your cluster cannot find the BLAS and LAPACK libraries by itself when users run MATLAB interface for MUMPS. To solve this issue, users can `export LD_PRELOAD` before running MATLAB. For example, if BLAS and LAPACK libraries are used through Intel MKL, users can type,

```shell
export LD_PRELOAD=$MKLROOT/libmkl_intel_lp64.so:$MKLROOT/libmkl_sequential.so:$MKLROOT/libmkl_intel_thread.so:$MKLROOT/libmkl_core.so
```

Note that `MKLROOT` is the path for Intel MKL libraries, where you find your BLAS and LAPACK libraries.

 In different Linux clusters, the Intel MKL library may be put in different
path. If you are using in Lmod module system and Intel libraries are installed, you can use `module load intel` and `printenv PATH`. One of the PATH should be the path of Intel directory, then typically the MKLROOT is under the path `.../mkl/lib/intel64`, where `...` is the Intel directory. If Lmod module system is not used in your cluster, then you can also try to find the path like `.../mkl/lib/intel64`, where `...` is the Intel directory.