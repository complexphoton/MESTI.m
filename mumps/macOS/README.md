## Prerequisites for MUMPS installation on a Mac

We need the following tools before compiling MUMPS and its MATLAB interface on macOS.

Much of the following (and the whole MUMPS compilation process) needs to be done in a command-line environment. You can either use Apple's built-in Terminal app (in the /Applications/Utilities/ folder) or a 3rd-party one like iTerm2.

Note that now `Makefile.inc` and `make.inc` on macOS we provided do not activate the OpenMP feature.

### MATLAB for Apple silicon Macs

If you are not sure whether your Mac runs on an Intel processor or Apple silicon, click the Apple logo on the top left corner of the menu bar, and click About This Mac. A Mac with Intel processor will show an item **Processor** (*e.g.*, Intel Core i7); a Mac with Apple silicon will show an item **Chip** (*e.g.*, Apple M1) instead.

If your Mac runs on an Intel processor, you can skip this part :)

If your Mac runs on Apple silicon&mdash;congratulations, it's faster (we've found an M1 Macbook Pro to be about twice as fast as an Intel Macbook Pro when running MESTI with the APF method due to the improved performance of vecLib), but there's extra work for you. None of the official MATLAB releases to date provide native support for Apple silicon; they run on Apple silicon Macs through Rosetta 2. That means when we compile the MATLAB interface for MUMPS, the <code>mex</code> compiler will try to compile for an Intel architecture, but Apple silicon is an ARM architecture, resulting in an error like [this](https://www.mathworks.com/matlabcentral/answers/1696860-use-gsl-compiled-on-apple-silicon-with-mex-function-on-matlab-2021b).

MATLAB recently released a public beta with native Apple silicon support [here](https://www.mathworks.com/support/apple-silicon-r2022a-beta.html). Go to this page (you'll need to log in with a Mathworks account), fill out the required fields and Submit, and follow the instructions there. You'll need to install Azul Zulu OpenJDK 8 with the .dmg option, then download and install the MATLAB R2022a beta.

### Xcode

We need to install [Xcode](https://developer.apple.com/xcode/) because it is [required](https://www.mathworks.com/support/requirements/supported-compilers.html) by the MATLAB compiler <code>mex</code>. 

Supposedly, <code>mex</code> only needs the <code>clang</code> compiler of Xcode, which can be obtained through the Xcode Command Line Tools (CLT)&mdash;a much smaller installation compared to the full Xcode. However, if you only install CLT, you'll get errors like what's described in [this thread](https://www.mathworks.com/matlabcentral/answers/307362-mex-on-macosx-without-xcode) when using <code>mex</code>. If you don't have enough disk space for Xcode, you can skip this part and use the workarounds described in that thread when compiling the MATLAB interface for MUMPS.

To install Xcode, open the App Store app, search for Xcode, and install it. This is a large download (12.7GB for Xcode 13), and the resulting Xcode.app takes up 33 GB of disk space.

After installation, you'll need to accept the Xcode license agreement. You can open Xcode.app from the /Applications folder, upon which a prompt will ask you to accept the license. Alternatively, you can enter the following in the terminal
```
sudo xcodebuild -license accept
```

### Xcode Command Line Tools

Apple's Xcode Command Line Tools (CLT) include <code>make</code>, <code>ar</code>, and <code>ranlib</code>, Apple's C compiler <code>clang</code>, and Apple's implementation of BLAS and LAPACK, [vecLib](https://developer.apple.com/documentation/accelerate/veclib), within its Accelerate framework. We need those.

Technically, the same tools are already included in the full Xcode installation. However, there can be issues linking to vecLib when Xcode is installed but not CLT. Plus, Homebrew requires CLT.

In the next step, we will install Homebrew, which will install CLT (if not already installed). So, nothing needs to be done here.

### Homebrew

Xcode and CLT do not include a Fortran compiler. Here, we use [Homebrew](https://brew.sh/) to install one; Homebrew can also be used for the optional installation of <code>openblas</code> and <code>cmake</code> below. If you already have Homebrew installed, run <code>brew update</code> to update it. If you don't have Homebrew installed, copy the following line and paste it in terminal.
```
/bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"
```
Follow instructions from the script to install Homebrew and then to add it to your PATH.

If CLT was not installed prior, Homebrew will install it as part of the script above. Note that even though we already installed Xcode, Homebrew will still ask for CLT to be installed, for [various reasons](https://github.com/Homebrew/brew/issues/10714#issuecomment-786663987).

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

MUMPS uses BLAS extensively, so we need a BLAS library. One option is Apple's [vecLib](https://developer.apple.com/documentation/accelerate/veclib) (which is already installed above with Xcode and CLT); another is [OpenBLAS](https://www.openblas.net/); another is [MKL](https://www.intel.com/content/www/us/en/developer/tools/oneapi/onemkl.html).

If you have an Apple silicon Mac, you should use vecLib, which we found to be much faster than OpenBLAS on an M1 Macbook Pro, consistent with [what others have found](https://github.com/danielchalef/openblas-benchmark-m1). (In fact, it is the vecLib improvement that makes an M1 Mac faster than an Intel Mac; when OpenBLAS is used, we found there to be no performance difference between an M1 Macbook Pro and an Intel Mackbook Pro when running MESTI with APF.) There's no reason to use MKL since MKL is only optimized for Intel machines. So, skip this part.

If you have an Intel Mac, you can consider installing OpenBLAS or MKL. It's optional, since vecLib is still available with no extra work. But we did find on an Intel Macbook Pro that MESTI with the APF method is about 30% faster when MUMPS is compiled with OpenBLAS compared to with vecLib. Most likely MKL will be even faster than OpenBLAS on an Intel Mac, but we did not test it; you can follow the instructions [here](https://www.intel.com/content/www/us/en/developer/tools/oneapi/onemkl.html) if you want to use MKL instead (but note the macOS `Makefile.inc` and `make.inc` files we provide only consider vecLib and OpenBLAS.)

To install OpenBLAS, enter
```
brew install openblas
```
in the terminal. This will install OpenBLAS in <code>/opt/homebrew/opt/openblas</code> for an Apple Silicon Mac, <code>/usr/local/opt/openblas</code> for an Intel Mac. You can find this information in the future with the <code>brew info openblas</code> command.

### (Optional) CMake

If you plan to install [METIS](http://glaros.dtc.umn.edu/gkhome/metis/metis/overview) for matrix ordering, you'll also need [CMake](https://cmake.org/). Enter
```
brew install cmake
```
in terminal to install CMake.
