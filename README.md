# MESTI.m

**MESTI** (Maxwell's Equations Solver with Thousands of Inputs) is a finite-difference frequency-domain (FDFD) solver that implements the **Schur complement scattering analysis (SCSA)** method described in [arXiv2203.xxxxx](https://arxiv.org/abs/2203.xxxxx). It allows users to solve Maxwell's equations with thousands of inputs while using computing time and memory that is even less than what a typical direct solver uses to solve with just one input.

MESTI.m is written in MATLAB and considers transverse-magnetic waves in 2D. A 3D vectorial version of MESTI written in Julia is under development and will be released later.  

With MESTI.m, the user can specify arbitrary permittivity profiles of the system, arbitrary lists of input sources (user-specified or automatically built), and arbitrary lists of output projections (or no projection, in which case the complete field profiles are returned). MESTI.m implements all of the common boundary conditions, perfectly matched layer (PML) with both imaginary and real coordinate stretching, as well as exact outgoing boundaries in two-sided or one-sided geometries. In addition to SCSA, MESTI.m also implements several conventional direct solvers.

## Installation

No installation is required for MESTI.m itself; just download it and add the <code>MESTI.m/src</code> folder to the search path using the <code>addpath</code> command in MATLAB. The MATLAB version should be R2019b or later. (Using an earlier version is possible but requires minor edits.)

However, to use the SCSA method (which is the key distinguishing feature of MESTI), the user needs to install the MUMPS package and its MATLAB interface. Without MUMPS, MESTI.m will still run but will only use slower conventional solvers that are not based on the Schur complement. So, MUMPS installation is strongly recommended.  See this [MUMPS installation](./mumps/README.md) page for steps to install MUMPS.

## Available Functions 

The function <code>mesti(syst, B, C, D)</code> provides the most flexibility. The user can specify arbitrary permittivity profiles, arbitrary boundary conditions, PML on any or all sides, with any lists of sources and any lists of output projections (or no output projection). But the user needs to build the source and/or projection profiles beforehand.

The function <code>mesti2s(syst, in, out)</code> deals specifically with two-sided or one-sided geometries where the boundary condition in y is closed (e.g. periodic or PEC)  and the boundary condition in x is outgoing. It builds the lists of sources and lists of output projections for the user, so it is simpler to use when such geometries are of interest. <code>mesti2s()</code> also provides (1) exact outgoing boundaries implemented through the self-energy, which is not available in <code>mesti()</code>, and (2) the recursive Green's function method from the [RGF](https://github.com/chiaweihsu/RGF) repository which is efficient in 1D and for systems with small-to-medium transverse widths. 

For field-profile computations, simply omit the argument <code>C</code> or  <code>out</code>, or set it to <code>[]</code>.

Detailed usage of these functions are given in the documentation section (comment lines at the beginning) of the <code>.m</code> function files. One can, for example, type <code>help mesti</code> in MATLAB to see such documentation.

In addition, the user can use the function <code>mesti_build_channels()</code> to build the input and/or output matrices when using <code>mesti()</code>, or to determine which subset of the channels are of interest when using <code>mesti2s()</code>.

MESTI.m has two other general-purpose functions <code>mesti_build_fdfd_matrix()</code> and <code>mesti_matrix_solver()</code>, which are called by <code>mesti()</code>. A typical user shouldn't need to directly use these two functions, but they are available.

## Examples

Several examples are provided in the <code>examples</code> folder, which are useful for learning the functionalities and usage of MESTI.m. Each example has its own folder, with its <code>.m</code> script, a Markdown page that shows the output of the example, and auxiliary files specific to that example. Below are links to their Markdown pages:

1. [Fabry–Pérot etalon](./examples/1d_fabry_perot/fabry_perot.md) (1D)
2. [Distributed Bragg reflector](./examples/1d_distributed_bragg_reflector/distributed_bragg_reflector.md) (1D)
3. [Open channel in a scattering medium](./examples/2d_open_channel_through_disorder/open_channel_through_disorder.md) (2D)
4. [Reflection Matrix in Gaussian-Beam Basis](./examples/2d_reflection_matrix_Gaussian_beams/reflection_matrix_Gaussian_beams.md) (2D) 
5. [Meta-atom design for metasurfaces](./examples/2d_meta_atom_design/meta_atom_design.md) (2D)
6. [Millimeter-wide metalens](./examples/2d_metalens_full/metalens_full.md) (2D)

More examples will be added later.

## Gallery
Below are some images from the <code>examples</code> folder:

1. Propagation through a Fabry–Pérot etalon
<img src="./examples/1d_fabry_perot/fabry_perot_field_profile.gif" width="336" height="264"> 

2. Open channel propagating through a strongly scattering medium
<img src="./examples/2d_open_channel_through_disorder/disorder_open_channel.gif" width="530" height="398"> 

3. Field profile of Gaussian beams focused at different locations
<img src="./examples/2d_reflection_matrix_Gaussian_beams/reflection_matrix_Gaussian_beams.gif" width="596" height="343"> 

4. Full-wave simulation of a mm-wide hyperbolic metalens 
<img src="./examples/2d_metalens_full/metalens_intensity_profile_0_degree.png" width="504" height="378"> 

5. Angle dependence of the metalens above
<img src="./examples/2d_metalens_full/metalens_Strehl_ratio_and_transmission_eff.png" width="504" height="378"> 
