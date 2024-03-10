## Prerequisites for MUMPS installation on a Mac

We need the following tools before compiling MUMPS and its MATLAB interface on macOS.

Much of the following (and the whole MUMPS compilation process) needs to be done in a command-line environment. You can either use Apple's built-in Terminal app (in the /Applications/Utilities/ folder) or a 3rd-party one like iTerm2.

### MATLAB for Apple silicon Macs

If you are not sure whether your Mac runs on an Intel processor or Apple silicon, click the Apple logo on the top left corner of the menu bar, and click About This Mac. A Mac with Intel processor will show an item **Processor** (*e.g.*, Intel Core i7); a Mac with Apple silicon will show an item **Chip** (*e.g.*, Apple M1) instead.

If your Mac runs on Apple silicon, the version of your MATLAB must be R2023b or later. Older versions of MATLAB can run on Apple silicon through Rosetta 2 but do not support building mex files as native ARM binaries, which we need to compile the MATLAB interface for MUMPS.

### Xcode Command Line Tools

Apple's Xcode Command Line Tools (CLT) include <code>make</code>, <code>ar</code>, and <code>ranlib</code>, Apple's C compiler <code>clang</code>, and Apple's implementation of BLAS and LAPACK, [vecLib](https://developer.apple.com/documentation/accelerate/veclib), within its Accelerate framework. We need those.

In the next step, we will install Homebrew, which will install CLT (if not already installed). So, nothing needs to be done here.

### Homebrew

Xcode and CLT do not include a Fortran compiler. Here, we use [Homebrew](https://brew.sh/) to install one; Homebrew can also be used for the optional installation of <code>openblas</code> and <code>cmake</code> below. If you already have Homebrew installed, run <code>brew update</code> to update it. If you don't have Homebrew installed, copy the following line and paste it in terminal.
```
/bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"
```
Follow instructions from the script to install Homebrew and then to add it to your PATH.

If CLT was not installed prior, Homebrew will install it as part of the script above. 

As described in the [Homebrew installation page](https://docs.brew.sh/Installation), this installs Homebrew to <code>/opt/homebrew</code> for an Apple Silicon Mac, <code>/usr/local</code> for an Intel Mac.

After installation, enter
```
brew doctor
```
to make sure there's no outstanding issues.

### GNU compiler collection

After installing Homebrew, enter
```
brew install gcc
```
in terminal. This will install the [GNU compiler collection (GCC)](https://gcc.gnu.org/), which includes the Fortran compiler <code>gfortran</code>. If you already installed GCC through Homebrew before, this will update GCC to the most current stable version.

If you update GCC in the future (through <code>brew upgrade</code> or by running <code>brew install gcc</code> again), and the first number of the GCC version changes (_e.g._ from gcc 11.x.x to gcc 12.x.x), then you'll need to update that version number in the <code>LIBFORT</code> path within the <code>make.inc</code> of the MATLAB interface of MUMPS, and recompile that MATLAB interface.

### (Optional) OpenBLAS

MUMPS uses BLAS extensively, so we need a BLAS library. One option is Apple's [vecLib](https://developer.apple.com/documentation/accelerate/veclib) within its Accelerate framework (which is already installed above with Xcode and CLT); another is [OpenBLAS](https://www.openblas.net/); another is [MKL](https://www.intel.com/content/www/us/en/developer/tools/oneapi/onemkl.html).

If you have an Apple silicon Mac, you should use vecLib for BLAS. On an Apple silicon Mac, vecLib utilizes the Apple Matrix coprocessor (AMX), so it is much faster than single-threaded OpenBLAS (*e.g.*, see [this benchmark](https://github.com/danielchalef/openblas-benchmark-m1)) and has a comparable speed as multithreaded OpenBLAS (*e.g.*, see [this benchmark](https://github.com/OpenMathLib/OpenBLAS/issues/2814#issuecomment-771505972)). In this case, you don't need to install OpenBLAS, and you can skip this section.

If you have an Intel Mac, we recommend installing OpenBLAS. We found on an Intel Macbook Pro that MESTI with the APF method is about 30% faster when MUMPS is compiled with single-threaded OpenBLAS compared to with vecLib. Most likely MKL will be even faster than OpenBLAS on an Intel Mac, but we did not test it; you can follow the instructions [here](https://www.intel.com/content/www/us/en/developer/tools/oneapi/onemkl.html) if you want to use MKL instead (but note the macOS `Makefile.inc` and `make.inc` files we provide only consider vecLib and OpenBLAS.) To install OpenBLAS, enter
```
brew install openblas
```
in the terminal. This will install OpenBLAS in <code>/usr/local/opt/openblas</code> for an Intel Mac (or <code>/opt/homebrew/opt/openblas</code> for an Apple Silicon Mac). You can find this information in the future with the <code>brew info openblas</code> command.

### (Optional) CMake

If you plan to install [METIS](http://glaros.dtc.umn.edu/gkhome/metis/metis/overview) for matrix ordering, you'll also need [CMake](https://cmake.org/). Enter
```
brew install cmake
```
in terminal to install CMake.

## Troubleshooting

During the MATLAB interface compilation, you may encounter a warning with Xcode `license has not been accepted` and then an error with `no supported compiler was found`. This issue can be resolved by entering
```
/usr/libexec/PlistBuddy -c 'Add :IDEXcodeVersionForAgreedToGMLicense string 10.0' ~/Library/Preferences/com.apple.dt.Xcode.plist
```
in terminal. Alternatively, you can install the full Xcode and accept the Xcode license agreement if this trick does not work.
