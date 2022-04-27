## Prerequisites for MUMPS installation on a Linux cluster

We need the following tools before compiling MUMPS and its MATLAB interface in the Linux cluster.

### MATLAB 

The MATLAB should be installed before compiling MUMPS-MATLAB interface. The MATLAB compiler `mex` is required when compiling the MUMPS-MATLAB interface. 

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

The compilation of MUMPS requires both C and Fortran compilers. Both C and Fortran compliers are  included in GNU Compiler Collection (GCC), which should be installed in the cluster. 

It is common that GCC is loaded by default on a cluster. You can type
```shell
gcc -v
```
to see if GCC is available.

If GCC is not loaded, you can load it on a Lmod module system by typing
```shell
module load gcc
```
In non-Lmod module system, you may type,
```shell
.../gcc/setup.sh
```
where `...` is the path to GCC in the cluster.

### BLAS and LAPACK

MUMPS requires both BLAS and LAPACK libraries, which are standard libraries on Linux cluster. These libraries are also included in many implementations, such as Intel MKL and OpenBLAS. 



In this example, we use Intel MKL. 



In different Linux clusters, the Intel MKL may be put in different path. Here we show how to find the path of `MKLROOT`, where is the path for Intel MKL, including BLAS and LAPACK libraries in `$(MKLROOT)/lib/intel64`. 



If you are using in Lmod module system and Intel MKL are installed, you can use 

```shell
module load intel
module load intel-mkl
echo $MKLROOT
```

`MKLROOT`  is assigned and should be printed out, if Intel MKL is available. The `MKLROOT` should be under the path `.../mkl`, where `...` is the Intel directory path.



If Lmod module system is not used in your cluster, then you can also try to find the path like `.../mkl`, where `...` is the Intel directory path. Then you can assign the `MKLROOT` by 

```shell
source .../mkl/bin/mklvars.sh intel64
echo $MKLROO
```

`MKLROOT`  is assigned and should be printed out.



In some cases, your cluster cannot find the BLAS and LAPACK libraries by itself when you run MATLAB interface for MUMPS. To solve this issue, you can `export LD_PRELOAD` and `LD_LIBRARY_PATH` before running MATLAB. For example, if BLAS and LAPACK libraries are used through Intel MKL, you can type,

```shell
export LD_PRELOAD=$MKLROOT/lib/intel64/libmkl_intel_lp64.so:$MKLROOT/lib/intel64/libmkl_sequential.so:$MKLROOT/lib/intel64/libmkl_intel_thread.so:$MKLROOT/lib/intel64/libmkl_core.so
export LD_LIBRARY_PATH=$MKLROOT/lib/intel64
```



If you still cannot figure out the path of the libraries, you may need to contact the maintenance team of your cluster.