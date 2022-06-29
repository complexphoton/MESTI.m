function [S, channels, info] = mesti2s(syst, in, out, opts)
%MESTI2S Solves frequency-domain scattering problems in a two-sided geometry.
%   [field_profiles, channels, info] = MESTI2S(syst, in) returns the spatial
%   field profiles of Ez(x,y) for scattering problems of 2D transverse-magnetic
%   (TM) fields satisfying
%      [- (d/dx)^2 - (d/dy)^2 - (omega/c)^2*epsilon(x,y)] Ez(x,y) = 0,
%   or of Hz(x,y) for 2D transverse-electric (TE) fields satisfying
%      [- (d/dx)*(1/epsilon(x,y))_yy*(d/dx) - (d/dy)*(1/epsilon(x,y))_xx*(d/dy)
%       + (d/dy)*(1/epsilon(x,y))_xy*(d/dx) + (d/dx)*(1/epsilon(x,y))_yx*(d/dy)
%       - (omega/c)^2] Hz(x,y) = 0,
%   where Ez or Hz is each a superposition of incident and scattered waves.
%      The polarization (TM or TE), relative permittivity profile epsilon(x,y),
%   frequency omega, and boundary conditions are specified by structure 'syst';
%   epsilon(x,y) must have homogeneous spaces on the left (-x) and right (+x)
%   sides, with an outgoing boundary in x for the scattered waves and a closed
%   boundary in y.
%      The incident wavefronts from left and/or right are specified by variable
%   'in'.
%      The returned 'field_profiles' is a 3D array, with field_profiles(:,:,i)
%   being the total (incident + scattered) field profile of Ez or Hz given the
%   i-th incident wavefront. The returned 'channels' is a structure containing
%   properties of the propagating and evanescent channels in the homogeneous
%   spaces on the left and right. The information of the computation is returned
%   in structure 'info'.
%
%   [S, channels, info] = MESTI2S(syst, in, out) returns the scattering matrix
%   S, where 'in' and 'out' specify either the list of input/output channels or
%   the input/output wavefronts. When the MUMPS function zmumps() is available,
%   this is typically done by computing the Schur complement of an augmented
%   matrix K through a partial factorization.
%
%   [field_profiles, channels, info] = MESTI2S(syst, in, [], opts) and
%   [S, channels, info] = MESTI2S(syst, in, out, opts) allow detailed options to
%   be specified with structure 'opts'.
%
%   In mesti2s(), the boundary condition in y must be closed (e.g., periodic or
%   PEC). Given the closed boundary, the set of transverse modes forms a
%   complete and orthonormal basis of propagating and evanescent channels.
%   The inputs and outputs are specified in the basis of these propagating 
%   channels, with coefficients normalized with respect to the flux in the
%   longitudinal (x) direction. Properties of those channels are given by
%   mesti_build_channels().
%
%   When an open boundary in y is of interest, the appropriate input/output
%   basis depends on the problem, so the user should use the more general
%   function mesti() and will need to specify the input and output matrices B
%   and C.
%
%   MESTI considers nonmagnetic materials with isotropic (i.e., scalar)
%   permittivity. Even though epsilon(x,y) is originally a scalar in the
%   continuous system, subpixel smoothing (which is needed to ensure smooth
%   variation with respect to parameters and to reach second-order accuracy for
%   TE polarization) creates anisotropy where epsilon(x,y) becomes a symmetric
%   tensor. TM polarization uses the zz component of epsilon(x,y). TE
%   polarization uses the xx, yy, and xy (equals yx) components of epsilon(x,y).
%   The other components of epsilon(x,y) are zero.
%
%   This file builds the input and output channels using mesti_build_channels(),
%   builds the matrices B and C, and then calls mesti() to solve the scattering
%   problems. For scattering-matrix computations in TM polarization, mesti2s()
%   can also use the recursive Green's function method instead, which is
%   efficient in 1D and for 2D systems with a small width in y.
%
%   === Input Arguments ===
%   syst (scalar structure; required):
%      A structure that specifies the system, used to build the FDFD matrix A.
%      It contains the following fields:
%      syst.polarization (character vector; optional):
%         Polarization. Possible choices are:
%            'TM' - Transverse-magnetic field (Hx, Hy, Ez)
%            'TE' - Transverse-electric field (Ex, Ey, Hz)
%         TM field uses syst.epsilon, and TE field uses syst.inv_epsilon. If
%         only one of syst.epsilon and syst.inv_epsilon is given,
%         syst.polarization is optional and will be automatically picked based
%         on which one is given. If syst.epsilon and syst.inv_epsilon are both
%         given, then syst.polarization must be specified.
%      syst.epsilon (numeric matrix, real or complex; required for TM):
%         An ny_Ez-by-nx_Ez matrix discretizing the relative permittivity
%         profile epsilon(x,y) of the scattering region. Specifically,
%         syst.epsilon(m,n) is the scalar epsilon(x,y) averaged over a square
%         with area (syst.dx)^2 centered at the point (x_n, y_m) where Ez(x,y)
%         is located on the Yee lattice. It is the zz component of the
%         discretized epsilon(x,y) tensor from subpixel smoothing, used by TM
%         fields. We choose
%            (x_n, y_m) = (n-0.5, m-0.5)*syst.dx,
%            with n = 1, ..., nx_Ez, m = 1, ..., ny_Ez.
%         such that the lower corner of the first pixel syst.epsilon(m=1,n=1) is
%         at (x,y) = (0,0), and the back side of the last slice
%         syst.epsilon(m,n=nx_Ez) is at the edge of the scattering region, x = L
%         = nx_Ez*syst.dx;
%            Outside the scattering region (with x < 0 and x > L), epsilon(x,y)
%         is given by scalars syst.epsilon_L and syst.epsilon_R for the
%         homogeneous spaces on the two sides. Note that nx_Ez = 0 corresponds
%         to L = 0 (ie, no scattering region) and is allowed.
%            Note that y corresponds to the first index m, and x corresponds to
%         the second index n.
%            The positive imaginary part of syst.epsilon describes absorption,
%         and the negative imaginary part describes linear gain.
%            One can use syst.epsilon with ny_Ez = 1 and with a periodic or
%         Bloch periodic boundary in y to simulate 1D systems where the relative
%         permittivity profile is translationally invariant in y.
%      syst.inv_epsilon (cell array; required for TE):
%         The xx, yy, and xy components of the discretized inverse relative
%         permittivity 1/epsilon(x,y) tensor of the scattering region from
%         subpixel smoothing, used by TE fields. It has three elements
%            inv_epsilon{1}: (1/epsilon(x,y))_xx, size [ny_Ez, nx_Hz]
%            inv_epsilon{2}: (1/epsilon(x,y))_yy, size [ny_Hz, nx_Ez]
%            inv_epsilon{3}: (1/epsilon(x,y))_xy, size [ny_Ez, nx_Ez]
%         The third element, inv_epsilon{3}, is optional and is treated as zero
%         when syst.inv_epsilon only has two elements. The yx component is not
%         specified since we only consider symmetric 1/epsilon tensors where
%         (1/epsilon)_yx = (1/epsilon)_xy.
%            The different components are located at different points:
%            - Hz(x,y) at (x_{n-0.5}, y_{m-0.5}); size [ny_Hz, nx_Hz].
%            - (1/epsilon(x,y))_xx and Dx ~ dHz/dy at (x_{n+0.5}, y_m).
%            - (1/epsilon(x,y))_yy and Dy ~ dHz/dx at (x_n, y_{m+0.5}).
%            - (1/epsilon(x,y))_xy at (x_n, y_m), same as Ez.
%         Here, (x_n, y_m) is the location of Ez and syst.epsilon above.
%            inv_epsilon{1}, inv_epsilon{2}, and inv_epsilon{3} should each be
%         (1/epsilon(x,y))_xx, (1/epsilon(x,y))_yy, and (1/epsilon(x,y))_xy from
%         subpixel smoothing, averaged over a square with area (syst.dx)^2
%         centered at the points where each of them is located. The
%         1/epsilon(x,y) tensor from subpixel smoothing is given by Eq. (1) of
%         Farjadpour et al, Optics Letters 31, 2972 (2006).
%            In x direction for a two-sided geometry, nx_Hz = nx_Ez + 1 =
%         L/syst.dx + 1, and the sites on x_{n+0.5} start from x = 0 (i.e., half
%         a pixel before the first site of x_n) and end on x = L (i.e., half a
%         pixel after the last site of x_n). That means the first and last
%         slices of Hz sit exactly on the front and back surfaces of the
%         scattering region, and the first and last slices of inv_epsilon{1}
%         extend half a pixel outside the scattering region. For a one-sided
%         geometry (where syst.epsilon_R is not given), nx_Hz = nx_Ez =
%         L/syst.dx, and the last slice of Hz sits on x = L - syst.dx (since Hz
%         = 0 at x = L).
%            In y direction, where these points start and end depend on the
%         boundary condition.
%            For periodic, Bloch periodic, and PMCPEC boundary conditions in y,
%         ny_Hz = ny_Ez, and all of the sites on y_{n+0.5} are half a pixel
%         after the corresponding sites on y_n.
%            For PEC boundary condition in y, ny_Hz = ny_Ez + 1, and the sites
%         on y_{n+0.5} start from half a pixel before the first site of y_n and
%         end on half a pixel after the last site of y_n.
%            For PMC boundary condition in y, ny_Hz = ny_Ez - 1, and the sites
%         on y_{n+0.5} start from half a pixel after the first site of y_n and
%         end on half a pixel before the last site of y_n.
%            For PECPMC boundary condition in y, ny_Hz = ny_Ez, and all of the
%         sites on y_{n+0.5} are half a pixel before the corresponding sites on
%         y_n.
%            For a 1D problem with ny_Hz = 1 and periodic boundary condition in
%         y, the y derivative vanishes, so inv_epsilon{1} is not used and does
%         not need to be specified.
%      syst.epsilon_L (real scalar; required):
%         Relative permittivity of the homogeneous space on the left.
%      syst.epsilon_R (real scalar or []; optional):
%         Relative permittivity of the homogeneous space on the right. If
%         syst.epsilon_R is not given or is empty, the system will be one-sided,
%         terminated on the right with a PEC boundary with Ez = 0 at x = L +
%         0.5*syst.dx, or wth a PMC boundary with Hz = 0 at x = L.
%      syst.length_unit (anything; optional):
%         Length unit, such as micron, nm, or some reference wavelength. This
%         code only uses dimensionless quantities, so syst.length_unit is never
%         used. This syst.length_unit is meant to help the user interpret the
%         units of (x,y), dx, wavelength, ky_B, etc.
%      syst.wavelength (numeric scalar, real or complex; required):
%         Vacuum wavelength 2*pi*c/omega, in units of syst.length_unit.
%      syst.dx (positive scalar; required):
%         Discretization grid size, in units of syst.length_unit.
%      syst.yBC (character vector; required unless syst.ky_B is specified):
%         Boundary condition (BC) at the two ends in y direction, effectively
%         specifying Ez(m,n) or Hz(m,n) at m=0 and m=ny_Ez+1 or ny_Hz+1, one
%         pixel beyond the computation domain. Available choices are:
%           'Bloch'    - Ez(m+ny_Ez,n) = Ez(m,n)*exp(1i*syst.ky_B*ny_Ez*syst.dx)
%                        Hz(m+ny_Hz,n) = Hz(m,n)*exp(1i*syst.ky_B*ny_Hz*syst.dx)
%           'periodic' - equivalent to 'Bloch' with syst.ky_B = 0
%           'PEC'      - Ez(0,n) = Ez(ny_Ez+1,n) = 0
%                        Hz(0,n) = Hz(1,n); Hz(ny_Hz+1,n) = Hz(ny_Hz,n)
%           'PMC'      - Ez(0,n) = Ez(1,n); Ez(ny_Ez+1,n) = Ez(ny_Ez,n)
%                        Hz(0,n) = Hz(ny_Hz+1,n) = 0
%           'PECPMC'   - Ez(0,n) = 0; Ez(ny_Ez+1,n) = Ez(ny_Ez,n)
%                        Hz(0,n) = Hz(1,n); Hz(ny_Hz+1,n) = 0
%           'PMCPEC'   - Ez(0,n) = Ez(1,n); Ez(ny_Ez+1,n) = 0
%                        Hz(0,n) = 0; Hz(ny_Hz+1,n) = Hz(ny_Hz,n)
%         where PEC stands for perfect electric conductor (for which Ez = 0 and
%         Ex ~ dHz/dy = 0 at the boundary) and PMC stands for perfect magnetic
%         conductor (for which Hz = 0 and Hx ~ dEz/dy = 0 at the boundary).
%            Note that this yBC also defines a complete and orthonormal set of
%         transverse modes, upon which the input and output channels in input
%         arguments 'in' and 'out' are defined; mesti2s() does not support PML
%         in y direction because a closed boundary is necessary for defining
%         such a transverse basis.
%            Here, syst.yBC is required, with no default choice (except when
%         syst.ky_B is given, in which case syst.yBC = 'Bloch' is automatically
%         used).
%      syst.ky_B (real scalar; optional):
%         Bloch wave number in y direction, in units of 1/syst.length_unit.
%         syst.ky_B is only used when syst.yBC = 'Bloch'. It is allowed to
%         specify a complex-valued syst.ky_B, but a warning will be displayed.
%      syst.xBC (character vector or scalar structure or cell array; optional):
%         In mesti(), outgoing boundary condition is always used in x direction.
%         But there are different options for implementing such outgoing
%         boundary. By default, 
%            syst.xBC = 'outgoing';
%         in which case self-energy based on the retarded Green's function of a
%         semi-infinite discrete homogeneous space is used to enforce an exact
%         outgoing boundary for all propagating and evanescent waves on both
%         sides. Doing so doesn't increase the simulation domain; For TM, only a
%         single-pixel slice of syst.epsilon_L and a single-pixel slice of
%         syst.epsilon_R are added on the left and on the right to place the
%         self-energy and the input and output, with a total simulation domain
%         size of [ny_Ez, nx_Ez+2]. For TE, the self-energy and input/output are
%         directly placed on the surface where the first and last pixels of Hz
%         are, so the total simulation domain size remains [ny_Hz, nx_Hz].
%            A drawback of self-energy is that it introduces coupling between
%         every pair of pixels in y direction on that slice, making matrix A
%         less sparse. When the system is wide (roughly when ny > 800), it is
%         more efficient to implement the outgoing boundary using perfectly
%         matched layer (PML), which attenuates outgoing waves with minimal
%         reflection.
%            Note that in mesti2s(), the PML is placed in the homogeneous spaces
%         specified by syst.epsilon_L and syst.epsilon_R, outside of the
%         scattering region specified by syst.epsilon or syst.inv_epsilon. (This
%         is different from the more general function mesti(), where
%         syst.epsilon specifies the entire simulation domain, so PML is placed
%         inside syst.epsilon.)
%            To use the same PML on both sides or to use PML on the left of a
%         one-sided geometry, set syst.xBC to be a scalar structure with the
%         following fields:
%            npixels_PML (positive integer scalar; required): Number of PML
%               pixels. This number of pixels is added in addition to the
%               scattering region.
%            npixels_spacer (non-negative integer scalar; optional): Number of
%               homogeneous-space pixels to be added between the PML pixels and
%               syst.epsilon/syst.inv_epsilon, used to attenuate evanescent
%               waves. Defaults to 0.
%            power_sigma (non-negative scalar; optional): Power of the
%               polynomial grading for the conductivity sigma; defaults to 3.
%            sigma_max_over_omega (non-negative scalar; optional):
%               Conductivity at the end of the PML; defaults to
%                  0.8*(power_sigma+1)/((2*pi/wavelength)*dx*sqrt(epsilon_bg)).
%               where epsilon_bg is either syst.epsilon_L or syst.epsilon_R.
%               This is used to attenuate propagating waves.
%            power_kappa (non-negative scalar; optional): Power of the
%               polynomial grading for the real-coordinate-stretching factor
%               kappa; defaults to 3.
%            kappa_max (real scalar no smaller than 1; optional):
%               Real-coordinate-stretching factor at the end of the PML;
%               defaults to 15. This is used to accelerate the attenuation of
%               evanescent waves. kappa_max = 1 means no real-coordinate
%               stretching.
%            power_alpha (non-negative scalar; optional): Power of the
%               polynomial grading for the CFS alpha factor; defaults to 1.
%            alpha_max_over_omega (non-negative scalar; optional): Complex-
%               frequency-shifting (CFS) factor at the beginning of the PML.
%               This is typically used in time-domain simulations to suppress
%               late-time (low-frequency) reflections. We don't use it by
%               default (alpha_max_over_omega = 0) since we are in frequency
%               domain.
%         We use the following PML coordinate-stretching factor:
%            s(u) = kappa(u) + sigma(u)./(alpha(u) - i*omega)
%         with
%            sigma(u)/omega = sigma_max_over_omega*(u.^power_sigma),
%            kappa(u) = 1 + (kappa_max-1)*(u.^power_kappa),
%            alpha(u)/omega = alpha_max_over_omega*((1-u).^power_alpha),
%         where omega is frequency, and u goes linearly from 0 at the beginning
%         of the PML to 1 at the end of the PML. 
%            The syntax above (setting syst.xBC to 'outgoing' or a scalar
%         structure) uses the same parameters on the two sides. For a two-sided
%         geometry (i.e., when syst.epsilon_L and syst.epsilon_R are both
%         provided), one can also use different parameters on the two sides. To
%         do so, set syst.xBC to a two-element cell array, with the first
%         element specifying parameters on the left, the second element
%         specifying parameters on the right. For example,
%            syst.xBC = {'outgoing', struct('npixels_PML', 40)};
%         specifies exact outgoing boundary on the left, 40 pixels of PML on the
%         right.
%            With real-coordinate stretching, PML can attenuate evanescent waves
%         more efficiently than free space, so npixels_spacer defaults to 0.
%            The PML thickness should be chosen based on the acceptable level of
%         reflectivity given the discretization resolution and the range of wave
%         numbers (i.e., angles) involved; more PML pixels gives lower
%         reflectivity. Typically 10-40 pixels are sufficient.
%            PML is a layer, not technically a boundary condition. In mesti2s(),
%         the boundary condition behind the PML is always set to Dirichlet (ie,
%         PEC for TM, PMC for TE).
%            Note that PML cannot be used when opts.method = 'RGF'.
%      syst.use_continuous_dispersion (logical scalar; optional):
%         Whether to use the dispersion equation of the continuous wave equation
%         when building the input/output channels. Defaults to false, in which
%         case the finite-difference dispersion is used.
%      syst.m0 (real numeric scalar; optional, defaults to 0):
%         Center of the transverse mode profile with periodic or Bloch periodic
%         boundary condition, phi_{m,a} = exp(i*ky(a)*syst.dx*(m-m0))/sqrt(ny),
%         where ky(a) = syst.ky_B + a*(2*pi/ny*syst.dx) and ny = ny_Ez for TM,
%         ny_Hz for TE.
%   in (cell array or scalar structure; required):
%      The set of input channels or input wavefronts.
%         To specify all propagating channels on one side or on both sides, use
%            in = {'left'}, or 
%            in = {'right'}, or
%            in = {'left', 'right'}.
%      The character vectors 'left' and 'right' can be abbreviated as 'L' and
%      'R' respectively.
%         To specify a subset of the propagating channels, let 'in' be a scalar
%      structure with the following fields:
%            in.ind_L (integer vector): Vector containing the indices of
%               propagating channels incident on the left side.
%            in.ind_R (integer vector): Vector containing the indices of
%               propagating channels incident on the right side.
%      One can provide only in.ind_L or only in.ind_R or both of them.
%         The above generates flux-normalized single-channel inputs of the form
%            f_a^L(x,y) = 1/sqrt(nu_L(a))*phi_a(y)*exp(i*kx_L(a)*x)
%      for input channel 'a' from the left, and/or
%            f_a^R(x,y) = 1/sqrt(nu_R(a))*phi_a(y)*exp(-i*kx_R(a)*(x-L))
%      for input channel 'a' from the right, where nu normalizes flux in the x
%      direction, kx is the longitudinal wave number, and phi_a(y) is the
%      transverse profile of the a-th propagating channel; nu = sin(kx(a)*dx)
%      for TM polarization, nu = sin(kx(a)*dx)/epsilon_bg for TE polarization.
%         The user can first use mesti_build_channels() to get the indices, wave
%      numbers, and transverse profiles of the propagating channels; based on
%      that, the user can specify the list of channels of interest through
%      in.ind_L and in.ind_R above, or a list of customized wavefronts through
%      in.v_L and in.v_R below.
%         When the input(s) of interests are not a single channel but a
%      superposition of multiple propagating channels, such custom input
%      wavefronts can be specified by letting 'in' be a scalar structure with
%      the following fields:
%            in.v_L (numeric matrix): Matrix where each column specifies the
%               coefficients of propagating channels on the left for one input
%               wavefront from the left; the wavefront is a superposition of all
%               propagating channels of the form f_a^L(x,y) above, with the
%               superposition coefficients given by that column of in.v_L.
%               size(in.v_L, 1) must equal N_prop_L, the total number of
%               propagating channels on the left; size(in.v_L, 2) is the number
%               of input wavefronts.
%            in.v_R (numeric matrix): Analogous to to in.v_L, but specifying
%               input wavefronts from the right instead.
%      Note that the input wavefronts from the left and the input wavefronts
%      from the right are treated as separate inputs. In other words, each input
%      either comes from the left or comes from the right, but not from both
%      sides. If an input with incidence from both sides is of interest, the
%      user can manually superimpose results from the separate-side
%      computations.
%   out (cell array or scalar structure or []; optional):
%      The set of output channels or output wavefronts.
%         When out=[] or when out is omitted as in input argument (i.e., when
%      nargin <= 2), no output projection is used, and the spatial field
%      profiles Ez(x,y) or Hz(x,y) corresponding to the set of inputs are
%      returned.
%         When out is given, the scattering matrix is returned, with the output
%      basis of the scattering matrix specified by out. In this case, out
%      follows the same syntax as the input argument 'in'. Specifically, one can
%      specify all propagating channels on one side or on both sides with
%            out = {'left'}, or 
%            out = {'right'}, or
%            out = {'left', 'right'}.
%      The character vectors 'left' and 'right' can be abbreviated as 'L' and
%      'R' respectively.
%         One can alternatively specify a subset of propagating
%      channels with
%            out.ind_L (integer vector), and/or
%            out.ind_R (integer vector).
%      For example, if r_full is the full reflection matrix computed with in =
%      {'L'} and out = {'L'}, then the reflection matrix computed using in.ind_L
%      and out.ind_L is equivalent to r_full(out.ind_L, in.ind_L).
%         One can also let the output basis of the scattering matrix be a
%      superposition of multiple propagating channels, with such custom output
%      wavefronts specified by
%            out.v_L (numeric matrix), and/or
%            out.v_R (numeric matrix).
%      Here, each row of the scattering matrix corresponds to projection onto an
%      output wavefront specified by a column of out.v_L or of out.v_R. For
%      example, if r_full is the full reflection matrix computed with in = {'L'}
%      and out = {'L'}, then the reflection matrix computed using in.v_L and
%      out.v_L is equivalent to (out.v_L)'*r_full*in.v_L where ' is the
%      conjugate transpose.
%   opts (scalar structure; optional, defaults to an empty struct):
%      A structure that specifies the options of computation; defaults to an
%      empty structure. It can contain the following fields (all optional):
%      opts.verbal (logical scalar; optional, defaults to true):
%         Whether to print system information and timing to the standard output.
%      opts.nx_L (non-negative integer scalar; optional, defaults to 0):
%         Number of pixels of homogeneous space on the left (syst.epsilon_L) to
%         include when returning the spatial field profile; not used for
%         scattering matrix computations.
%      opts.nx_R (non-negative integer scalar; optional, defaults to 0):
%         Number of pixels of homogeneous space on the right (syst.epsilon_R) to
%         include when returning the spatial field profile; not used for
%         scattering matrix computations. Note that opts.nx_R can still be used
%         in one-sided geometries where syst.epsilon_R is not given; the field
%         profile on the right is simply zero in such case.
%      opts.solver (character vector; optional):
%         The software used for sparse matrix factorization. Available choices
%         are (case-insensitive):
%            'MUMPS'  - (default) Uses MUMPS. Its MATLAB interface zmumps.m must
%                       be in MATLAB's search path. This is much faster and uses
%                       less memory.
%            'MATLAB' - Uses the built-in lu() function in MATLAB, which uses
%                       UMFPACK with AMD ordering. This requires no installation
%                       but is much slower. This is be used by default if
%                       zmumps.m is not found in the search path.
%         opts.solver is not used when opts.method = 'RGF'.
%      opts.method (character vector; optional):
%         The solution method. Available choices are (case-insensitive):
%            'APF' - Augmented partial factorization. When opts.solver =
%                    'MUMPS', C*inv(A)*B is obtained through the Schur
%                    complement of an augmented matrix K = [A,B;C,0] using a
%                    partial factorization; this is the true APF. When
%                    opts.solver = 'MATLAB', C*inv(A)*B is obtained as
%                    C*inv(U)*inv(L)*B with optimized grouping, which is not the
%                    true APF but is slightly better than factorize_and_solve.
%                    Cannot be used for computing the full field profile
%                    inv(A)*B or with iterative refinement.
%            'FS'  - Factorize and solve. Factorize A=L*U, solve for inv(A)*B
%                    with forward and backward substitutions, and optionally
%                    project with C.
%            'RGF' - Recursive Green's function method. Cannot be used for
%                    computing the full field profile or with iterative
%                    refinement, and cannot be used with PML or with TE
%                    polarization.
%            'factorize_and_solve' - Same as 'FS'.
%         By default, if the input argument 'out' is not given or if
%         opts.iterative_refinement = true, then 'factorize_and_solve' is used.
%         Otherwise, 'APF' is used if ny is large or if TE polarization is used
%         or if PML has been specified; 'RGF' is used otherwise.
%      opts.symmetrize_K (logical scalar; optional):
%         Whether or not to pad input and/or output channels and perform
%         permutations to make matrix K = [A,B;C,0] symmetric when computing its
%         Schur complement, which lowers computing time and memory usage. Such
%         channel padding and permutation is reversed afterwards and does not
%         affect what mesti2s() returns.
%            opts.symmetrize_K can only be used when all of the following are
%         met: (1) input argument 'out' is given, (2) 'in' and 'out' are not
%         specified as wavefronts, (3) opts.solver = 'MUMPS', (4) opts.method =
%         'APF', and (5) the boundary condition in y is not Bloch periodic.
%         When all of these conditions are met, opts.symmetrize_K defaults to
%         true; otherwise it is not used.
%      opts.clear_syst (logical scalar; optional, defaults to false):
%         When opts.clear_syst = true and opts.method is not 'RGF', variable
%         'syst' will be cleared in the caller's workspace to reduce peak memory
%         usage. This can be used when syst.epsilon or syst.inv_epsilon takes up
%         significant memory and is not needed after calling mesti2s().
%      opts.clear_memory (logical scalar; optional, defaults to true):
%         Whether or not to clear variables inside mesti2s() to reduce peak
%         memory usage.
%      opts.verbal_solver (logical scalar; optional, defaults to false):
%         Whether to have the solver print detailed information to the standard
%         output. Note the behavior of output from MUMPS depends on compiler.
%      opts.use_METIS (logical scalar; optional, defaults to false):
%         Whether to use METIS (instead of the default AMD) to compute the
%         ordering in MUMPS. Using METIS can sometimes reduce memory usage
%         and/or factorization and solve time, but it typically takes longer at
%         the analysis (i.e., ordering) stage.
%      opts.nrhs (positive integer scalar; optional):
%         The number of right-hand sides (number of columns of the input matrix
%         B) to consider simultaneously, used only when opts.method =
%         'factorize_and_solve' and input argument 'out' is given. Defaults to 1
%         if opts.iterative_refinement = true, 10 if opts.solver = 'MUMPS' with
%         opts.iterative_refinement = false, 4 otherwise.
%      opts.store_ordering (logical scalar; optional, defaults to false):
%         Whether to store the ordering sequence (permutation) for matrix A or
%         matrix K; only possible when opts.solver = 'MUMPS'. If
%         opts.store_ordering = true, the ordering will be returned in
%         info.ordering.
%      opts.ordering (positive integer vector; optional):
%         A user-specified ordering sequence for matrix A or matrix K, used only
%         when opts.solver = 'MUMPS'. Using the ordering from a previous
%         computation can speed up the analysis stage, but the matrix size must
%         be the same.
%      opts.analysis_only (logical scalar; optional, defaults to false):
%         When opts.analysis_only = true, the factorization and solution steps
%         will be skipped, and S = [] will be returned. The user can use
%         opts.analysis_only = true with opts.store_ordering = true to return
%         the ordering for A or K; only possible when opts.solver = 'MUMPS'.
%      opts.nthreads_OMP (positive integer scalar; optional):
%         Number of OpenMP threads used in MUMPS; overwrites the OMP_NUM_THREADS
%         environment variable.
%      opts.iterative_refinement (logical scalar; optional, defaults to false):
%         Whether to use iterative refinement in MUMPS to lower round-off
%         errors. Iterative refinement can only be used when opts.solver =
%         'MUMPS' and opts.method = 'factorize_and_solve' and input argument
%         'out' is given, in which case opts.nrhs must equal 1. When iterative
%         refinement is used, the relevant information will be returned in
%         info.itr_ref_nsteps, info.itr_ref_omega_1, and info.itr_ref_omega_2.
%
%   === Output Arguments ===
%   field_profiles (3D array):
%      For field-profile computations (i.e., when 'out' is not given), the
%      returned field_profiles is a 3D array containing the field profiles, such
%      that field_profiles(:,:,a) is the total (incident + scattered) field of
%      Ez or Hz corresponding to the a-th input wavefront, with
%         size(field_profiles, 1) = ny
%         size(field_profiles, 2) = opts.nx_L + nx + opts.nx_R
%      where [ny, nx] = [ny_Ez, nx_Ez] for TM, [ny_Hz, nx_Hz] for TE. Recall
%      that nx_Hz = nx_Ez + 1 = L/syst.dx + 1 for a two-sided geometry; nx_Hz =
%      nx_Ez = L/syst.dx for a one-sided geometry.
%         By default, opts.nx_L = opts.nx_R = 0, and field_profiles only contain
%      Ez or Hz in the scattering region. When opts.nx_L and opts.nx_R are
%      specified, field_profiles will also contain opts.nx_L pixels of Ez or Hz
%      in the homogeneous space on the left, opts.nx_R pixels of Ez or Hz
%      in the homogeneous space on the right.
%   S (numeric matrix):
%      For scattering-matrix computations (i.e., when 'out' is given), S is the
%      scattering matrix, such that S(b,a) is the flux-normalized field
%      coefficient in the b-th propagating output channel (or the b-th output
%      wavefront) given incident wave in the a-th propagating input channel (or
%      the a-th input wavefront).
%         When all propagating channels on one side or on both sides are
%      requested, e.g. with in = {'left'} or out = {'left', 'right'}, matrix S
%      includes channels.L.N_prop propagating channels on the left and/or
%      channels.R.N_prop on the right, with transverse wave numbers given by
%      vectors channels.L.kydx_prop and channels.R.kydx_prop, longitudinal wave
%      numbers given by vectors channels.L.kxdx_prop and channels.R.kxdx_prop.
%      The transverse profiles of these channels are channels.fun_phi(kydx).
%      Matrix S is then one block of S_full = [[r_L; t_L], [t_R; r_R]] depending
%      on which side(s) are specified for the input and output, with r_L and r_R
%      being reflection matrices on the left and right sides, t_L and t_R being
%      transmission matrices with input from the left and right respectively.
%      Note that channels on the left always come before channels on the right
%      in matrix S, so the first channels.L.N_prop elements of the matrix are
%      always channels on the left (if requested), with the last
%      channels.R.N_prop elements being channels on the right (if requested);
%      the ordering of elements within the cell arrays 'in' and/or 'out' is not
%      considered.
%         The phases of the elements of the scattering matrix depend on the
%      reference plane. For channels on the left, the reference plane is at
%      x = 0. For channels on the right, the reference plane is at x = L =
%      nx_Ez*syst.dx.
%         When a subset of the propagating channels are requested, with
%      in.ind_L, in.ind_R, out.ind_L, and/or out.ind_R, matrix S includes such
%      subset of the propagating channels. For example, if r_full is the full
%      reflection matrix computed with in = {'L'} and out = {'L'}, then the
%      reflection matrix computed using in.ind_L and out.ind_L is equivalent to
%      r_full(out.ind_L, in.ind_L).
%         When the input wavefronts and/or output basis is specified by in.v_L,
%      in.v_R, out.v_L, and/or out.v_R, matrix S is the full scattering matrix
%      with the input and/or output basis changed. For example, if r_full is the
%      full reflection matrix computed with in = {'L'} and out = {'L'}, then the
%      reflection matrix computed using in.v_L and out.v_L is equivalent to
%      (out.v_L)'*r_full*in.v_L where ' is the conjugate transpose.
%   channels (structure):
%      A structure returned by function mesti_build_channels() that contains
%      properties of the propagating and evanescent channels of the homogeneous
%      spaces on the left and right. Type "help mesti_build_channels" for more
%      information.
%   info (scalar structure):
%      A structure that contains the following fields:
%      info.opts (scalar structure):
%         The final 'opts' used, excluding the user-specified matrix ordering.
%      info.timing (scalar structure):
%         A structure containing timing of the various stages, in seconds, in
%         fields 'total', 'init', 'build', 'analyze', 'factorize', 'solve'.
%      info.xPML (two-element cell array; optional);
%         PML parameters on the low and high sides of x direction, if used.
%      info.ordering_method (character vector; optional):
%         Ordering method used in MUMPS.
%      info.ordering (positive integer vector; optional):
%         Ordering sequence returned by MUMPS when opts.store_ordering = true.
%      info.itr_ref_nsteps (integer vector; optional):
%         Number of steps of iterative refinement for each input, if
%         opts.iterative_refinement = true; 0 means no iterative refinement.
%      info.itr_ref_omega_1 (real vector; optional):
%         Scaled residual omega_1 at the end of iterative refinement for each
%         input; see MUMPS user guide section 3.3.2 for definition.
%      info.itr_ref_omega_2 (real vector; optional):
%         Scaled residual omega_2 at the end of iterative refinement for each
%         input; see MUMPS user guide section 3.3.2 for definition.
%
%   See also: mesti_build_channels, mesti

%% Part 1.1: Check validity of syst, assign default values to its fields, and parse BC and PML specifications

t0 = clock;

if nargin < 2; error('Not enough input arguments.'); end
if ~(isstruct(syst) && isscalar(syst)); error('Input argument ''syst'' must be a scalar structure.'); end
if ~isfield(syst, 'epsilon_L');         error('Input argument ''syst'' must have field ''epsilon_L''.'); end
if ~isfield(syst, 'wavelength');        error('Input argument ''syst'' must have field ''wavelength''.'); end
if ~isfield(syst, 'dx');                error('Input argument ''syst'' must have field ''dx''.'); end
if ~(isreal(syst.epsilon_L)     && isscalar(syst.epsilon_L));  error('syst.epsilon_L must be a real scalar.'); end
if ~(isnumeric(syst.wavelength) && isscalar(syst.wavelength)); error('syst.wavelength must be a numeric scalar.'); end
if ~(isreal(syst.dx) && isscalar(syst.dx) && syst.dx > 0);     error('syst.dx must be a positive scalar.'); end

% Pick the polarization to use; assign use_TM
if isfield(syst, 'polarization')
    if strcmpi(syst.polarization, 'TM')
        use_TM = true;
    elseif strcmpi(syst.polarization, 'TE')
        use_TM = false;
    else
        error('syst.polarization, if given, must be ''TM'' or ''TE''.');
    end
else
    % When syst.polarization is not given, we automatically pick based on whether syst.epsilon or syst.inv_epsilon is given.
    if isfield(syst, 'epsilon') && ~isfield(syst, 'inv_epsilon')
        use_TM = true;
    elseif ~isfield(syst, 'epsilon') && isfield(syst, 'inv_epsilon')
        use_TM = false;
    elseif isfield(syst, 'epsilon') && isfield(syst, 'inv_epsilon')
        error('syst.polarization must be given when syst.epsilon and syst.inv_epsilon both exist.');
    else % neither syst.epsilon nor syst.inv_epsilon exists
        error('Input argument ''syst'' must have field ''epsilon'' or ''inv_epsilon''.');
    end
end

% Check syst.epsilon (for TM) and syst.inv_epsilon (for TE)
if use_TM
    if ~isfield(syst, 'epsilon')
        error('syst.epsilon must be given when syst.polarization = ''TM''.');
    elseif ~(isnumeric(syst.epsilon) && ismatrix(syst.epsilon))
        error('syst.epsilon must be a numeric matrix, if given.');
    end
    syst.polarization = 'TM';
    str_pol = 'Ez'; % for printing system info
else
    if ~isfield(syst, 'inv_epsilon')
        error('syst.inv_epsilon must be given when syst.polarization = ''TE''.');
    elseif ~iscell(syst.inv_epsilon) || (numel(syst.inv_epsilon) ~= 2 && numel(syst.inv_epsilon) ~= 3)
        error('syst.inv_epsilon must be a two-element or three-element cell array, if given.');
    elseif ~(ismatrix(syst.inv_epsilon{1}) && isnumeric(syst.inv_epsilon{1}))
        error('syst.inv_epsilon{1} must be a numeric matrix.');
    elseif ~(ismatrix(syst.inv_epsilon{2}) && isnumeric(syst.inv_epsilon{2}))
        error('syst.inv_epsilon{2} must be a numeric matrix.');
    elseif numel(syst.inv_epsilon) == 3 && ~(ismatrix(syst.inv_epsilon{3}) && isnumeric(syst.inv_epsilon{3}))
        error('syst.inv_epsilon{3} must be a numeric matrix, if given.');
    end
    syst.polarization = 'TE';
    str_pol = 'Hz'; % for printing system info
end

% Check that the user did not accidentally use options only in mesti()
if isfield(syst, 'PML') && ~isempty(syst.PML)
    error('syst.PML is not used in mesti2s(); use syst.xBC instead.');
elseif isfield(syst, 'self_energy') && ~isempty(syst.self_energy)
    error('syst.self_energy is not used in mesti2s(); use syst.xBC = ''outgoing'' instead.');
elseif isfield(syst, 'kx_B') && ~isempty(syst.kx_B)
    error('syst.kx_B is not supported in mesti2s(); use mesti() if Bloch periodic BC in x is needed.');
elseif isfield(syst, 'PML_type') && ~isempty(syst.PML_type)
    warning('syst.PML_type is not supported in mesti2s(); UPML will be used. Use mesti() if SC-PML is needed.')
    syst = rmfield(syst, 'PML_type');
end

% syst.epsilon_R is an optional argument; determines whether the system is one-sided
if ~isfield(syst, 'epsilon_R') || isempty(syst.epsilon_R)
    two_sided = false;
    syst.epsilon_R = [];
elseif ~(isreal(syst.epsilon_R) && isscalar(syst.epsilon_R))
    error('syst.epsilon_R must be a real scalar, if given.');
else
    two_sided = true;
end

% Number of sites in y and x
if use_TM
    [ny_Ez, nx_Ez] = size(syst.epsilon);
    if ny_Ez==0; error('ny_Ez = size(syst.epsilon,1) cannot be zero.'); end

    nx_s = nx_Ez;    % number of pixels of Ez in the scattering region
    dn = 0.5;     % the source/detection is half a pixel away from x=0 and x=L

    nx = nx_Ez;
    ny = ny_Ez;
else
    [ny_Hz, nx_Ez] = size(syst.inv_epsilon{2}); % (1/epsilon)_yy
    if ny_Hz==0; error('ny_Hz = size(syst.inv_epsilon{2},1) cannot be zero.'); end

    if two_sided
        % nx_Hz = nx_Ez + 1 = L/syst.dx + 1 for a two-sided geometry
        nx_Hz = nx_Ez + 1;
    else
        % nx_Hz = nx_Ez = L/syst.dx for a one-sided geometry
        nx_Hz = nx_Ez;
    end

    nx_s = nx_Ez - 1; % number of pixels of Hz in the scattering region, excluding the two surfaces
    dn = 0;        % the source/detection is already at x=0 and x=L

    nx = nx_Hz;
    ny = ny_Hz;
end

% Check boundary condition in y
if isfield(syst, 'ky_B') && ~isempty(syst.ky_B)
    if ~(isnumeric(syst.ky_B) && isscalar(syst.ky_B))
        error('syst.ky_B must be a numeric scalar, if given.');
    elseif (isfield(syst, 'yBC') && ~isempty(syst.yBC)) && (iscell(syst.yBC) || ~strcmpi(syst.yBC, 'Bloch'))
        error('When syst.ky_B is given, syst.yBC must be ''Bloch'' if specified.');
    end
    syst.yBC = 'Bloch';
    % mesti_build_channels() uses ky_B*periodicity as the input arguments yBC for Bloch BC
    yBC = (syst.ky_B)*(ny*syst.dx); % dimensionless
else
    if ~isfield(syst, 'yBC') || isempty(syst.yBC)
        error('Input argument ''syst'' must have non-empty field ''yBC'' when syst.ky_B is not given.');
    elseif ~((ischar(syst.yBC) && isrow(syst.yBC)) || (isstring(syst.yBC) && isscalar(syst.yBC)))
        error('syst.yBC must be a character vector or string, if given.');
    elseif ~ismember(lower(syst.yBC), lower({'Bloch', 'periodic', 'PEC', 'PMC', 'PECPMC', 'PMCPEC'}))
        error('syst.yBC = ''%s'' is not a supported option; type ''help mesti2s'' for supported options.', syst.yBC);
    elseif strcmpi(syst.yBC, 'Bloch')
        error('syst.yBC = ''Bloch'' but syst.ky_B is not given.');
    end
    yBC = syst.yBC;
end

if ~use_TM
    if strcmpi(syst.yBC, 'PMC')
        ny_Ez = ny_Hz + 1;
    elseif strcmpi(syst.yBC, 'PEC')
        ny_Ez = ny_Hz - 1;
    else
        ny_Ez = ny_Hz;
    end

    % Check that size(syst.inv_epsilon{1}) = [ny_Ez, nx_Hz] since we will append to it later
    if ny_Hz == 1 && (strcmpi(syst.yBC, 'periodic') || isequal(yBC, 0))
        syst.inv_epsilon{1} = sparse(ny_Ez, nx_Hz); % (1/epsilon)_xx is not used when d/dy = 0
    else
        if isequal(size(syst.inv_epsilon{1}), [0, 0]) && ~(ny_Ez==0 && nx_Hz==0)
            error('syst.inv_epsilon{1} must be given (unless ny_Hz = 1 and yBC = ''periodic'').');
        elseif size(syst.inv_epsilon{1}, 1) ~= ny_Ez
            error('size(syst.inv_epsilon{1}, 1) must equal ny_Ez = %d.', ny_Ez);
        elseif size(syst.inv_epsilon{1}, 2) ~= nx_Hz
            if two_sided
                error('size(syst.inv_epsilon{1}, 2) must equal nx_Hz = size(syst.inv_epsilon{2}, 2) + 1 = %d for a two-sided geometry.', nx_Hz);
            else
                error('size(syst.inv_epsilon{1}, 2) must equal nx_Hz = size(syst.inv_epsilon{2}, 2) = %d for a one-sided geometry.', nx_Hz);
            end
        end
    end

    % Check that size(syst.inv_epsilon{3}) = [ny_Ez, nx_Ez] since we will append to it later
    if numel(syst.inv_epsilon) == 3
        if size(syst.inv_epsilon{3}, 1) ~= ny_Ez
            error('size(syst.inv_epsilon{3}, 1) must equal ny_Ez = %d, if given.', ny_Ez);
        elseif size(syst.inv_epsilon{3}, 2) ~= nx_Ez
            error('size(syst.inv_epsilon{3}, 2) must equal nx_Ez = size(syst.inv_epsilon{2}, 2) = %d, if given.', nx_Ez);
        end
    end
end

% Defaults to using an exact outgoing boundary in x
if ~isfield(syst, 'xBC') || isempty(syst.xBC)
    syst.xBC = 'outgoing';
elseif ~((~iscell(syst.xBC) && strcmpi(syst.xBC, 'outgoing')) || (isstruct(syst.xBC) && isscalar(syst.xBC)) || (iscell(syst.xBC) && numel(syst.xBC)==2))
    error('syst.xBC must be ''outgoing'' or a scalar structure or a two-element cell array, if given.');
end

% Apply the same xBC on all sides; the second element will be ignored if two_sided = false
if ~(iscell(syst.xBC) && numel(syst.xBC)==2)
    syst.xBC = {syst.xBC, syst.xBC};
elseif ~two_sided
    error('For a one-sided geometry, boundary condition on the right is PEC for TM, PMC for TE; syst.xBC only specifies boundary condition on the left and must be ''outgoing'' or a scalar structure, if given.');
end

% Start with default setting: exact outgoing boundary using self-energy, with no PML and no spacer on all sides
% nx_extra is the number of homogeneous-space pixels to be added to syst.epsilon or syst.inv_epsilon{2} in x direction.
% The two elements of nx_extra and use_self_energy correspond to the left and right sides.
if two_sided
    n_sides = 2;
    nx_extra = [1, 1];
    use_self_energy = [true, true];
    str_xBC = {'outgoing', 'outgoing'}; % to be used for printing later
else
    n_sides = 1;
    nx_extra = [1, 0];
    use_self_energy = [true, false];
    if use_TM
        str_xBC = {'outgoing', 'PEC'};
    else
        str_xBC = {'outgoing', 'PMC'};
    end
end
syst.PML = {}; % to be used in mesti()
n_PML = 0; % number of sides where PML is used

% Loop over syst.xBC and handle PML parameters, if specified
str_sides = {'low', 'high'};
for ii = 1:n_sides
    if ~strcmpi(syst.xBC{ii}, 'outgoing')
        use_self_energy(ii) = false;
        n_PML = n_PML + 1;
        str_xBC{ii} = 'PML';
        PML_ii = syst.xBC{ii};
        if ~(isstruct(PML_ii) && isscalar(PML_ii))
            error('syst.xBC{%d} must be ''outgoing'' or a scalar structure.', ii)
        end

        % Check that the user did not accidentally use options only in mesti()
        if isfield(PML_ii, 'npixels') && ~isempty(PML_ii.npixels)
            warning('syst.xBC{%d}.npixels is ignored in mesti2s(); use field ''npixels_PML''.', ii);
        end
        if isfield(PML_ii, 'direction') && ~isempty(PML_ii.direction)
            warning('syst.xBC{%d}.direction is ignored in mesti2s().', ii);
        end
        if isfield(PML_ii, 'side') && ~isempty(PML_ii.side)
            warning('syst.xBC{%d}.side is ignored in mesti2s().', ii);
        end

        % Number of PML pixels must be given
        % Other fields are optional and will be checked in mesti_build_fdfd_matrix()
        if ~isfield(PML_ii, 'npixels_PML') || isempty(PML_ii.npixels_PML)
            error('When syst.xBC{%d} is a structure, it must contain field ''npixels_PML''.', ii);
        end
        npix_PML = PML_ii.npixels_PML;
        if ~(isreal(npix_PML) && isscalar(npix_PML) && npix_PML >= 0 && round(npix_PML) == npix_PML)
            error('syst.xBC{%d}.npixels_PML must be a non-negative integer scalar.', ii);
        end

        % Defaults to no spacer
        if ~isfield(PML_ii, 'npixels_spacer') || isempty(PML_ii.npixels_spacer)
            PML_ii.npixels_spacer = 0;
        end
        npix_spacer = PML_ii.npixels_spacer;
        if ~(isreal(npix_spacer) && isscalar(npix_spacer) && npix_spacer >= 0 && round(npix_spacer) == npix_spacer)
            error('syst.xBC{%d}.npixels_spacer must be a non-negative integer scalar, if given.', ii);
        end

        % number of pixels in x to be added outside of syst.epsilon or syst.inv_epsilon
        nx_extra(ii) = 1 + npix_PML + npix_spacer;

        PML_ii = rmfield(PML_ii, {'npixels_PML', 'npixels_spacer'}); % these won't be used in mesti()
        PML_ii.npixels = npix_PML;   % to be used in mesti()
        PML_ii.direction = 'x';      % to be used in mesti()
        PML_ii.side = str_sides{ii}; % to be used in mesti()

        syst.PML{n_PML} = PML_ii; % to be used in mesti()
    end
end

% Total number of pixels in x
nx_tot = nx_s + sum(nx_extra);

% syst.use_continuous_dispersion and syst.m0 will be initialized/checked in mesti_build_channels()
if ~isfield(syst, 'use_continuous_dispersion') || isempty(syst.use_continuous_dispersion)
    syst.use_continuous_dispersion = false; % Use finite-difference dispersion by default
end
if ~isfield(syst, 'm0'); syst.m0 = []; end

%% Part 1.2: Check the input argument 'opts' and assign default values

if ~((isstruct(in) && isscalar(in)) || (iscell(in) && numel(in) <= 2))
    error('Input argument ''in'' must be a scalar structure or a cell array with no more than two elements.');
end

% out is an optional argument
if nargin < 3
    out = [];
end
if ~((isstruct(out) && isscalar(out)) || (iscell(out) && numel(out) <= 2) || isempty(out))
    error('Input argument ''out'' must be a scalar structure or a cell array with no more than two elements or [], if given.');
end

% opts is an optional argument
if nargin < 4 || isempty(opts)
    opts = struct();
end
if ~(isstruct(opts) && isscalar(opts))
    error('Input argument ''opts'' must be a scalar structure or [], if given.');
end

% Check that the user did not accidentally use options only in mesti()
if isfield(opts, 'prefactor') && ~isempty(opts.prefactor)
    error('opts.prefactor is not used in mesti2s(); the -2i prefactor is automatically included.');
end

% We return the scattering matrix S=C*inv(A)*B if out is given from input argument; else we return the spatial field profiles
% This opts.return_field_profile is not specified by the user; it will be returned as info.opts.return_field_profile to help debugging
opts.return_field_profile = isempty(out);

% By default, for field-profile computation, we only return result within the scattering region, setting nx_L = nx_R = 0.
if opts.return_field_profile
    if ~isfield(opts, 'nx_L') || isempty(opts.nx_L)
        opts.nx_L = 0;
    elseif ~(isreal(opts.nx_L) && isscalar(opts.nx_L) && opts.nx_L >= 0 && round(opts.nx_L) == opts.nx_L)
        error('opts.nx_L must be a non-negative integer scalar, if given.');
    end
    if ~isfield(opts, 'nx_R') || isempty(opts.nx_R)
        opts.nx_R = 0;
    elseif ~(isreal(opts.nx_R) && isscalar(opts.nx_R) && opts.nx_R >= 0 && round(opts.nx_R) == opts.nx_R)
        error('opts.nx_R must be a non-negative integer scalar, if given.');
    end
else
    if isfield(opts, 'nx_L') && ~isempty(opts.nx_L)
        warning('opts.nx_L is not used for scattering matrix computation; will be ignored.');
        opts = rmfield(opts, 'nx_L');
    end
    if isfield(opts, 'nx_R') && ~isempty(opts.nx_R)
        warning('opts.nx_R is not used for scattering matrix computation; will be ignored.');
        opts = rmfield(opts, 'nx_R');
    end
end

% Turn on verbal output by default
if ~isfield(opts, 'verbal') || isempty(opts.verbal)
    opts.verbal = true;
elseif ~(islogical(opts.verbal) && isscalar(opts.verbal))
    error('opts.verbal must be a logical scalar, if given.');
end

% By default, we don't clear syst in the caller's workspace
if ~isfield(opts, 'clear_syst') || isempty(opts.clear_syst)
    opts.clear_syst = false;
elseif ~(islogical(opts.clear_syst) && isscalar(opts.clear_syst))
    error('opts.clear_syst must be a logical scalar, if given.');
end

% By default, we will clear internal variables
if ~isfield(opts, 'clear_memory') || isempty(opts.clear_memory)
    opts.clear_memory = true;
elseif ~(islogical(opts.clear_memory) && isscalar(opts.clear_memory))
    error('opts.clear_memory must be a logical scalar, if given.');
end

% Use MUMPS for opts.solver when it is available
MUMPS_available = exist('zmumps','file');
if ~isfield(opts, 'solver') || isempty(opts.solver)
    if MUMPS_available
        opts.solver = 'MUMPS';
    else
        opts.solver = 'MATLAB';
    end
elseif ~((ischar(opts.solver) && isrow(opts.solver)) || (isstring(opts.solver) && isscalar(opts.solver)))
    error('opts.solver must be a character vector or string, if given.');
elseif ~ismember(lower(opts.solver), lower({'MUMPS', 'MATLAB'}))
    error('opts.solver = ''%s'' is not a supported option; use ''MUMPS'' or ''MATLAB''.', opts.solver);
elseif strcmpi(opts.solver, 'MUMPS') && ~MUMPS_available
    error('opts.solver = ''%s'' but function zmumps() is not found.', opts.solver)
end

if ~isfield(opts, 'method') || isempty(opts.method)
    % By default, if the input argument 'out' is not given or if
    % opts.iterative_refinement = true, then 'factorize_and_solve' is used.
    % Otherwise, 'APF' is used if ny is large or if TE polarization is used
    % or if PML has been specified; 'RGF' is used otherwise.
    if opts.return_field_profile || (isfield(opts, 'iterative_refinement') && isequal(opts.iterative_refinement, true))
        opts.method = 'factorize_and_solve';
    elseif n_PML > 0 || ~use_TM
        opts.method = 'APF';
    elseif isfield(opts, 'store_ordering') && isequal(opts.store_ordering, true)
        opts.method = 'APF';
    elseif isfield(opts, 'analysis_only') && isequal(opts.analysis_only, true)
        opts.method = 'APF';
    else
        % APF is more efficient when ny is large. RGF is more efficient when ny is small.
        if strcmpi(opts.solver, 'MUMPS')
            if ny > 80
                opts.method = 'APF';
            else
                opts.method = 'RGF';
            end
        else
            if ny > 200
                opts.method = 'APF';
            else
                opts.method = 'RGF';
            end
        end
    end
elseif ~((ischar(opts.method) && isrow(opts.method)) || (isstring(opts.method) && isscalar(opts.method)))
    error('opts.method must be a character vector or string, if given.');
elseif ~ismember(lower(opts.method), lower({'APF', 'factorize_and_solve', 'FS', 'RGF'}))
    error('opts.method = ''%s'' is not a supported option; use ''APF'' or ''factorize_and_solve'' or ''RGF''.', opts.method);
elseif opts.return_field_profile && (strcmpi(opts.method, 'APF') || strcmpi(opts.method, 'RGF'))
    error('opts.method = ''%s'' cannot be used for field-profile computations where out = []; use opts.method = ''factorize_and_solve'' instead.', opts.method)
elseif (isfield(opts, 'iterative_refinement') && isequal(opts.iterative_refinement, true)) && (strcmpi(opts.method, 'APF') || strcmpi(opts.method, 'RGF'))
    error('opts.method = ''%s'' cannot be used when opts.iterative_refinement = true; use opts.method = ''factorize_and_solve'' instead.', opts.method)
elseif strcmpi(opts.method, 'FS')
    opts.method = 'factorize_and_solve';  % opts.method = 'FS' is short for opts.method = 'factorize_and_solve'
end

if strcmpi(opts.method, 'RGF')
    % check options not used in RGF
    opts = rmfield(opts, 'solver'); % not used in RGF
    opts = rmfield(opts, 'clear_syst'); % not used in RGF
    if n_PML > 0
        error('PML cannot be used when opts.method = ''RGF''; use syst.xBC = ''outgoing'' or change opts.method.');
    elseif ~use_TM
        error('Our RGF implementation only works with TM polarization; choose a different opts.method.');
    end
    if isfield(opts, 'nrhs') && ~isempty(opts.nrhs)
        warning('opts.nrhs is not used when opts.method = ''RGF''; will be ignored.');
        opts = rmfield(opts, 'nrhs');
    end
    if isfield(opts, 'use_METIS') && ~isempty(opts.use_METIS)
        warning('opts.use_METIS is not used when opts.method = ''RGF''; will be ignored.');
        opts = rmfield(opts, 'use_METIS');
    end
    if isfield(opts, 'store_ordering') && isequal(opts.store_ordering, true)
        error('opts.store_ordering = true cannot be used when opts.method = ''RGF''.');
    end
    if isfield(opts, 'ordering') && ~isempty(opts.ordering)
        warning('opts.ordering is not used when opts.method = ''RGF''; will be ignored.');
        opts = rmfield(opts, 'ordering');
    end
    if isfield(opts, 'analysis_only') && isequal(opts.analysis_only, true)
        error('opts.analysis_only = true cannot be used when opts.method = ''RGF''.');
    end
    opts.analysis_only = false;
    if isfield(opts, 'nthreads_OMP') && ~isempty(opts.nthreads_OMP)
        warning('opts.nthreads_OMP is not used when opts.method = ''RGF''; will be ignored.');
        opts = rmfield(opts, 'nthreads_OMP');
    end
    use_RGF = true;
else
    use_RGF = false;
end

% opts.symmetrize_K will be checked/initialized later

% The following fields of opts will be checked/initialized in mesti_matrix_solver():
%    opts.verbal_solver
%    opts.use_METIS
%    opts.nrhs
%    opts.store_ordering
%    opts.ordering
%    opts.analysis_only
%    opts.nthreads_OMP
%    opts.iterative_refinement

% Set up the homogeneous-space channels on the two sides
k0dx = (2*pi/syst.wavelength)*(syst.dx);
if ~two_sided
    % For convenience, we set syst.epsilon_R so it is not [], and so we can still use channels.L below
    syst.epsilon_R = NaN;
end
% we can also use channels = mesti_build_channels(syst);
channels = mesti_build_channels(ny, syst.polarization, yBC, k0dx, syst.epsilon_L, syst.epsilon_R, syst.use_continuous_dispersion, syst.m0);
N_prop_L = channels.L.N_prop;
if two_sided
    N_prop_R = channels.R.N_prop;
end

if opts.verbal
    % print basic system info
    fprintf('System size: ny = %d, nx = %d => %d', ny, nx, nx_tot);
    if two_sided
        fprintf('; N_prop= {%d, %d}\n', N_prop_L, N_prop_R);
    else
        fprintf('; one-sided; N_prop_L = %d\n', N_prop_L);
    end
    fprintf('xBC = {%s, %s}; yBC = %s', str_xBC{1}, str_xBC{2}, syst.yBC);
    if strcmpi(syst.yBC, 'Bloch'); fprintf(' (ky_B = %.4f)', syst.ky_B); end
    fprintf('; %s polarization\n', str_pol);
end

t1 = clock; timing_init = etime(t1,t0); % Initialization time

%% Part 2.1: Build the self-energy matrix, if needed

if use_self_energy(1) || use_self_energy(2)
    if syst.use_continuous_dispersion
        warning('The outgoing boundary condition based on self-energy is not exact when syst.use_continuous_dispersion = true; consider using PML or setting syst.use_continuous_dispersion = false.');
    end
    if opts.verbal; fprintf('Building G0 ... '); end

    % Build ny-by-ny unitary matrix phi where the a-th column is the a-th transverse mode
    % This includes all the propagating and evanescent channels
    % With periodic boundary condition, phi is the shifted DFT matrix
    phi = channels.fun_phi(channels.kydx_all);

    % Retarded Green's function G0 of a semi-infinite homogeneous space, evaluated at the surface (just before the space is terminated)
    % A few notes about TE polarization:
    % In homogeneous space, A0_TE = (A0_TM)/epsilon_bg, so G0_TE = G0_TM*epsilon_bg.
    % When we compute the self-energy = V_SF*inv(A_F)*V_FS below, we get an additional (1/epsilon_bg)^2 prefactor for TE polarization from V_SF and V_FS.
    % So, overall, (self-energy)_TE = (self-energy)_TM/epsilon_bg.
    % Here, we absorb the additional (1/epsilon_bg)^2 prefactor into G0_TE. So, for TE polarization, the G0_L and G0_R below is actually G0_L/syst.epsilon_L^2 and G0_R/syst.epsilon_R^2. This is fine since we won't use G0 anywhere else; RGF is currently not supported for TE.
    if use_self_energy(1)
        if use_TM
            G0_L = (phi.*exp(1i*channels.L.kxdx_all))*(phi'); % use implicit expansion; note that channels.L.kxdx_all is a 1-by-ny row vector
        else
            G0_L = (phi.*(exp(1i*channels.L.kxdx_all)/syst.epsilon_L))*(phi'); % use implicit expansion; note that channels.L.kxdx_all is a 1-by-ny row vector
        end
    end
    if all(use_self_energy) && syst.epsilon_R == syst.epsilon_L
        G0_R = G0_L;
    elseif use_self_energy(2)
        if use_TM
            G0_R = (phi.*exp(1i*channels.R.kxdx_all))*(phi'); % use implicit expansion; note that channels.L.kxdx_all is a 1-by-ny row vector
        else
            G0_R = (phi.*(exp(1i*channels.R.kxdx_all)/syst.epsilon_R))*(phi'); % use implicit expansion; note that channels.L.kxdx_all is a 1-by-ny row vector
        end
    end

    % phi_prop_L and phi_prop_R will be used later for building B and/or C
    phi_prop_L = phi(:, channels.L.ind_prop);
    has_phi_prop_L = true;
    if two_sided
        if syst.epsilon_R == syst.epsilon_L
            phi_prop_R = phi_prop_L;
        else
            phi_prop_R = phi(:, channels.R.ind_prop);
        end
        has_phi_prop_R = true;
    end
    if opts.clear_memory; clear phi; end

    % self-energy can be used instead of PML to implement exact outgoing boundary condition.
    % self-energy is the retarded Green's function of the surrounding space (with a Dirichlet boundary surrounding the scattering region) evaluated on the surface.
    % Specifically, self-energy = V_SF*inv(A_F)*V_FS where F denotes the surrounding infinite free space, S denotes the scattering region as given by syst.epsilon or syst.inv_epsilon, A_F is the differential operator of the free space with a Dirichlet boundary surrounding the scattering region S, and V_FS is the coupling matrix between F and S (which is nonzero only along the interface between F and S).
    % For the geometry here, self-energy equals the retarded Green's function of a semi-infinite homogeneous space.
    if ~use_RGF
        % syst.self_energy will be used in mesti()
        if all(use_self_energy)
            syst.self_energy = [[G0_L; sparse((nx_s+nx_extra(2))*ny, ny)], sparse(nx_tot*ny, nx_s*ny), [sparse((nx_extra(1)+nx_s)*ny, ny); G0_R]];
        elseif use_self_energy(1)
            syst.self_energy = [[G0_L; sparse((nx_s+nx_extra(2))*ny, ny)], sparse(nx_tot*ny, (nx_s+nx_extra(2))*ny)];
        elseif use_self_energy(2)
            syst.self_energy = [sparse(nx_tot*ny, (nx_extra(1)+nx_s)*ny), [sparse((nx_extra(1)+nx_s)*ny, ny); G0_R]];
        end
        if opts.clear_memory; clear G0_L G0_R; end
    end

    t2 = clock; timing_build_G0 = etime(t2,t1);
    if opts.verbal; fprintf('elapsed time: %7.3f secs\n', timing_build_G0); end
else
    has_phi_prop_L = false;
    has_phi_prop_R = false;
    t2 = t1;
    timing_build_G0 = 0;
end

%% Part 2.2: Parse the input argument 'in'

if opts.verbal; fprintf('Building B,C... '); end

% Indices of input channels, in row vectors
ind_in_L = zeros(1, 0);
if two_sided
    ind_in_R = zeros(1, 0);
end
use_ind_in = true; % whether to use indices of input channels (if not, wavefront coefficients will be used)

if iscell(in)
    % Take all input channels on the left and/or right
    for ii = 1:numel(in)
        element = in{ii};
        if ~((ischar(element) && isrow(element)) || (isstring(element) && isscalar(element)))
            error('in{%d} must be a character vector or string.', ii);
        elseif ~ismember(element, {'left', 'right', 'L', 'R'})
            error('in{%d} = ''%s'' is not a supported option; use ''left'' or ''right''.', ii, element);
        elseif ismember(element, {'left', 'L'})
            ind_in_L = 1:N_prop_L;
        else % element = 'right' or 'R'
            if two_sided
                ind_in_R = 1:N_prop_R;
            else
                error('in{%d} = ''%s'' cannot be used in a one-sided geometry.', ii, element);
            end
        end
    end
else % isstruct(in)
    field_names = setdiff(fieldnames(in), {'ind_L','ind_R', 'v_L', 'v_R'});
    if numel(field_names) > 0
        error('Input argument ''in'' contains unrecognized field %s.', field_names{1});
    end

    % Use the user-specified set of input channels on the left and/or right
    if isfield(in, 'ind_L')
        ind_in_L = reshape(in.ind_L, 1, []);
        if ~(isreal(ind_in_L) && isequal(ind_in_L, round(ind_in_L)) && min(ind_in_L)>0 && max(ind_in_L)<=N_prop_L)
            error('in.ind_L, when specified, must be an array of positive integers not exceeding N_prop_L = %d.', N_prop_L);
        end
    end
    if isfield(in, 'ind_R')
        if two_sided
            ind_in_R = reshape(in.ind_R, 1, []);
            if ~(isreal(ind_in_R) && isequal(ind_in_R, round(ind_in_R)) && min(ind_in_R)>0 && max(ind_in_R)<=N_prop_R)
                error('in.ind_R, when specified, must be an array of positive integers not exceeding N_prop_R = %d.', N_prop_R);
            end
        else
            error('in.ind_R cannot be used in a one-sided geometry.');
        end
    end

    % Use user-specified input wavefronts
    if isfield(in, 'v_L') || isfield(in, 'v_R')
        use_ind_in = false;
        % The user can specify either the set of input channels or the input wavefronts, but not both
        if ~(isempty(ind_in_L) && (~two_sided || isempty(ind_in_R)))
            error('in.v_L and in.v_R cannot be used together with in.ind_L and in.ind_R.');
        end
        v_in_L = sparse(N_prop_L, 0);
        v_in_R = sparse(N_prop_R, 0);
        if isfield(in, 'v_L')
            v_in_L = in.v_L;
            if ~(ismatrix(v_in_L) && size(v_in_L, 1) == N_prop_L && isnumeric(v_in_L))
                error('in.v_L, when specified, must be a numeric matrix with size [N_prop_L, M_in_L] where N_prop_L = %d.', N_prop_L);
            end
        end
        if isfield(in, 'v_R')
            if two_sided
                v_in_R = in.v_R;
                if ~(ismatrix(v_in_R) && size(v_in_R, 1) == N_prop_R && isnumeric(v_in_R))
                    error('in.v_R, when specified, must be a numeric matrix with size [N_prop_R, M_in_R] where N_prop_R = %d.', N_prop_R);
                end
            else
                error('in.v_R cannot be used in a one-sided geometry.');
            end
        end
    end
end

% Number of inputs
if use_ind_in
    M_in_L = numel(ind_in_L);
    if two_sided
        M_in_R = numel(ind_in_R);
    end
else
    M_in_L = size(v_in_L, 2);
    % if nnz(v_in_L) is large, we will need phi_prop_L when we build B_L later
    if (nnz(v_in_L) >= N_prop_L) && ~has_phi_prop_L
        phi_prop_L = channels.fun_phi(channels.L.kydx_prop);
        has_phi_prop_L = true;
    end
    if two_sided
        M_in_R = size(v_in_R, 2);
        % if nnz(v_in_R) is large, we will need phi_prop_R when we build B_R later
        if (nnz(v_in_R) >= N_prop_R) && ~has_phi_prop_R
            if (syst.epsilon_R == syst.epsilon_L) && has_phi_prop_L
                phi_prop_R = phi_prop_L;
            else
                phi_prop_R = channels.fun_phi(channels.R.kydx_prop);
            end
            has_phi_prop_R = true;
        end
    end
end

%% Part 2.3: Parse the input argument 'out'

if ~isempty(out)
    % Indices of output channels, in row vectors
    ind_out_L = zeros(1, 0);
    if two_sided
        ind_out_R = zeros(1, 0);
    end
    use_ind_out = true; % whether to use indices of output channels (if not, wavefront coefficients will be used)

    if iscell(out)
        % Take all output channels on the left and/or right
        for ii = 1:numel(out)
            element = out{ii};
            if ~((ischar(element) && isrow(element)) || (isstring(element) && isscalar(element)))
                error('out{%d} must be a character vector or string.', ii);
            elseif ~ismember(element, {'left', 'right', 'L', 'R'})
                error('out{%d} = ''%s'' is not a supported option; use ''left'' or ''right''.', ii, element);
            elseif ismember(element, {'left', 'L'})
                ind_out_L = 1:N_prop_L;
            else % element = 'right' or 'R'
                if two_sided
                    ind_out_R = 1:N_prop_R;
                else
                    error('out{%d} = ''%s'' cannot be used in a one-sided geometry.', ii, element);
                end
            end
        end
    else % isstruct(out)
        field_names = setdiff(fieldnames(out), {'ind_L','ind_R', 'v_L', 'v_R'});
        if numel(field_names) > 0
            error('Input argument ''out'' contains unrecognized field %s.', field_names{1});
        end

        % Use the user-specified set of output channels on the left and/or right
        if isfield(out, 'ind_L')
            ind_out_L = reshape(out.ind_L, 1, []);
            if ~(isreal(ind_out_L) && isequal(ind_out_L, round(ind_out_L)) && min(ind_out_L)>0 && max(ind_out_L)<=N_prop_L)
                error('out.ind_L, when specified, must be an array of positive integers not exceeding N_prop_L = %d.', N_prop_L);
            end
        end
        if isfield(out, 'ind_R')
            if two_sided
                ind_out_R = reshape(out.ind_R, 1, []);
                if ~(isreal(ind_out_R) && isequal(ind_out_R, round(ind_out_R)) && min(ind_out_R)>0 && max(ind_out_R)<=N_prop_R)
                    error('out.ind_R, when specified, must be an array of positive integers not exceeding N_prop_R = %d.', N_prop_R);
                end
            else
                error('out.ind_R cannot be used in a one-sided geometry.');
            end
        end

        % Use user-specified output wavefronts
        if isfield(out, 'v_L') || isfield(out, 'v_R')
            use_ind_out = false;
            % The user can specify either the set of output channels or the output wavefronts, but not both
            if ~(isempty(ind_out_L) && (~two_sided || isempty(ind_out_R)))
                error('out.v_L and out.v_R cannot be used together with out.ind_L and out.ind_R.');
            end
            v_out_L = sparse(N_prop_L, 0);
            v_out_R = sparse(N_prop_R, 0);
            if isfield(out, 'v_L')
                v_out_L = out.v_L;
                if ~(ismatrix(v_out_L) && size(v_out_L, 1) == N_prop_L && isnumeric(v_out_L))
                    error('out.v_L, when specified, must be a numeric matrix with size [N_prop_L, M_out_L] where N_prop_L = %d.', N_prop_L);
                end
            end
            if isfield(out, 'v_R')
                if two_sided
                    v_out_R = out.v_R;
                    if ~(ismatrix(v_out_R) && size(v_out_R, 1) == N_prop_R && isnumeric(v_out_R))
                        error('out.v_R, when specified, must be a numeric matrix with size [N_prop_R, M_out_R] where N_prop_R = %d.', N_prop_R);
                    end
                else
                    error('out.v_R cannot be used in a one-sided geometry.');
                end
            end
        end
    end

    % Number of outputs
    if use_ind_out
        M_out_L = numel(ind_out_L);
        if two_sided
            M_out_R = numel(ind_out_R);
        end
    else
        M_out_L = size(v_out_L, 2);
        % if nnz(v_out_L) is large, we will need phi_prop_L when we build C_L later
        if (nnz(v_out_L) >= N_prop_L) && ~has_phi_prop_L
            phi_prop_L = channels.fun_phi(channels.L.kydx_prop);
            has_phi_prop_L = true;
        end
        if two_sided
            M_out_R = size(v_out_R, 2);
            % if nnz(v_out_R) is large, we will need phi_prop_R when we build C_R later
            if (nnz(v_out_R) >= N_prop_R) && ~has_phi_prop_R
                if (syst.epsilon_R == syst.epsilon_L) && has_phi_prop_L
                    phi_prop_R = phi_prop_L;
                else
                    phi_prop_R = channels.fun_phi(channels.R.kydx_prop);
                end
                has_phi_prop_R = true;
            end
        end
    end
end

%% Part 2.4: Build the source B_L/B_R and the projection C_L/C_R on the two surfaces

% First, we symmetrize the set of input channels and output channels, if possible

% Bloch periodic boundary condition with ky_B*periodicity != 0 or pi breaks the symmetry of A
% mesti_build_channels() currently does not return channels.L.ind_prop_conj when yBC == pi (even though such permutation does exist), so we cannot symmetrize K when yBC == pi.
% When ny == 1, even though matrix A is symmetric, there is still no permutation that flips the sign of ky, so we still set is_symmetric_A to false.
if isnumeric(yBC) && yBC ~= 0 % && yBC ~= pi && ny > 1
    is_symmetric_A = false;
else
    is_symmetric_A = true;
end

% opts.symmetrize_K can only be used when all of the following are met: (1) input argument 'out' is given, (2) 'in' and 'out' are not specified as wavefronts, (3) opts.solver = 'MUMPS', (4) opts.method = 'APF', and (5) the boundary condition in y is not Bloch periodic.
if ~isempty(out) && use_ind_in && use_ind_out && strcmpi(opts.method, 'APF') && strcmpi(opts.solver, 'MUMPS') && is_symmetric_A
    % By default, opts.symmetrize_K = true if it's possible
    if ~isfield(opts, 'symmetrize_K') || isempty(opts.symmetrize_K)
        opts.symmetrize_K = true;
    elseif ~(islogical(opts.symmetrize_K) && isscalar(opts.symmetrize_K))
        error('opts.symmetrize_K must be a logical scalar, if given.');
    end
    use_transpose_B = opts.symmetrize_K;
    opts = rmfield(opts, 'symmetrize_K'); % opts.symmetrize_K is not used in mesti()
else
    % symmetrize K is not possible
    if isfield(opts, 'symmetrize_K')
        if isequal(opts.symmetrize_K, true)
            warning('opts.symmetrize_K = true is only available when isempty(out) = false, use_ind_in = true, use_ind_out = true, opts.solver = ''MUMPS'', opts.method = ''APF'', and syst.yBC is not ''Bloch''. Here isempty(out) = %d, use_ind_in = %d, use_ind_out = %d, opts.solver = ''%s'', opts.method = ''%s'', syst.yBC = ''%s''; opts.symmetrize_K will be ignored', isempty(out), use_ind_in, use_ind_out, opts.solver, opts.method, syst.yBC);
        end
        opts = rmfield(opts, 'symmetrize_K');
    end
    use_transpose_B = false;
end

% No need to build C if we symmetrize K
if opts.return_field_profile || use_transpose_B
    build_C = false;
else
    build_C = true;
end

% Here we build:
% (1) the blocks of input matrix B on the left and right surfaces: B_L and B_R
% (2) the blocks of output matrix C on the left and right surfaces: C_L and C_R
% B_L, C_L, B_R, C_R are all dense matrices with size(..., 1) = ny. They are inputs/outputs on Ez for TM, Hz for TE.
% For TM polarization, these inputs/outputs are placed at one pixel outside syst.epsilon (at n=0 and n=nx_Ez+1), which is at x=-0.5*dx and x=L+0.5*dx.
% For TE polarization, these inputs/outputs are placed at the first and last pixels of syst.inv_epsilon{1} (at n=0.5 and n=nx_Ez+0.5), which is at x=0 and x=L.
% We want the reference plane to be at x=0 and x=L. For TM, we need to shift the phase by back propagating half a pixel (dn = 0.5). For TE, dn = 0.
% A line source of -2i*sqrt(nu)*phi(m) at n=0 will generate an x-flux-normalized incident field of exp(i*kxdx*|n|)*phi(m)/sqrt(nu), where nu = sin(kxdx) for TM polarization, nu = sin(kxdx)/epsilon_bg for TE polarization.
% Including the back propagation, we want B_L(m,a) = -2i*sqrt(nu(a))*exp(-i*kxdx(a)*dn)*phi(m,a); we will multiple the -2i prefactor at the end.
% The flux-normalized output projection is sqrt(nu(b))*conj(phi(:,b)); it will be transposed in mesti() or before rgf(). We also need to shift the phase by back propagating dn pixel.
% Therefore, we want C_L(m,b) = sqrt(nu(b))*exp(-i*kxdx(b)*dn)*conj(phi(m,b)).
% Note that the complex conjugation only applies to phi; the sqrt(nu) prefactor is not conjugated. (At real-valued frequency, nu is real-valued, so this doesn't matter. But at complex-valued frequency, nu is complex-valued, and we should not conjugate it. Note that the transverse basis is complete and orthonormal even when the frequency is complex, so the output projection doesn't need to be modified when the frequency is complex.)
% When the input/output channels are specified by channel indices, we will multiply the sqrt(nu(a))*exp(-i*kxdx(a)*dn) and sqrt(nu(b))*exp(-i*kxdx(b)*dn) prefactors at the end, after C*inv(A)*B is computed.
% When the input/output wavefronts are specified by in.v_L, in.v_R, out.v_L, out.v_R, we take superpositions of the channels using the v coefficients, with the sqrt(nu)*exp(-i*kxdx*dn) prefactors included.
if use_transpose_B % when opts.symmetrize_K = true
    % Here, we pad channels and/or permutate them such that C = transpose(B); this makes matrix K = [A,B;C,0] symmetric.
    % To have C=transpose(B), the complex conjugate of the transverse field profiles of the list of output channels must equal the transverse field profiles of the list of input channels, and the list of input channels and the list of output channels must have the same prefactor nu.
    % So, we expand the list of input channels (ind_in_L) to include the conjugate pairs of the output channels (channels.L.ind_prop_conj(ind_out_L)). The conjugate pairs correspond to flipping the sign of ky, and they share the same nu (which only depends on ky^2).
    % We only build B_L and B_R here, from which matrix B will be built in mesti(). C_L and C_R are not needed since matrix C will not be used until in mesti_matrix_solver().

    % We only need to keep the unique channel indices, since ind_in_L and channels.L.ind_prop_conj(ind_out_L) are likely to contain the same indices.
    % ind_in_out_L satisfies ind_L(ind_in_out_L) = [ind_in_L, channels.L.ind_prop_conj(ind_out_L)]. The computations will be done in ind_L, so later we can use ind_in_out_L to retrieve the original lists of input channels (with the first half of ind_in_out_L) and the original list of output channels (with the second half of ind_in_out_L).
    [ind_L, ~, ind_in_out_L] = unique([ind_in_L, channels.L.ind_prop_conj(ind_out_L)]);
    if has_phi_prop_L
        B_L = phi_prop_L(:,ind_L);
    else
        B_L = channels.fun_phi(channels.L.kydx_prop(ind_L));
    end
    if two_sided
        [ind_R, ~, ind_in_out_R] = unique([ind_in_R, channels.R.ind_prop_conj(ind_out_R)]);
        if has_phi_prop_R
            B_R = phi_prop_R(:,ind_R);
        else
            B_R = channels.fun_phi(channels.R.kydx_prop(ind_R));
        end
    end
else % without opts.symmetrize_K
    % Build input matrices B_L and B_R
    if use_ind_in % input channels specified by ind_in_L and ind_in_R
        if has_phi_prop_L
            B_L = phi_prop_L(:,ind_in_L);
        else
            B_L = channels.fun_phi(channels.L.kydx_prop(ind_in_L));
        end
        if two_sided
            if has_phi_prop_R
                B_R = phi_prop_R(:,ind_in_R);
            else
                B_R = channels.fun_phi(channels.R.kydx_prop(ind_in_R));
            end
        end
    else % input wavefronts specified by v_in_L and v_in_R
        if has_phi_prop_L % must be so when nnz(v_in_L) >= N_prop_L
            B_L = phi_prop_L*(reshape(channels.L.sqrt_nu_prop.*exp((-1i*dn)*channels.L.kxdx_prop), N_prop_L, 1).*v_in_L); % use implicit expansion
        else
            % build B_L using only the channels used
            B_L = zeros(ny, M_in_L);
            for ii = 1:M_in_L
                ind = find(v_in_L(:,ii));
                B_L(:,ii) = channels.fun_phi(channels.L.kydx_prop(ind))*(reshape(channels.L.sqrt_nu_prop(ind).*exp((-1i*dn)*channels.L.kxdx_prop(ind)),[],1).*v_in_L(ind,ii));
            end
        end
        if two_sided
            if has_phi_prop_R % must be so when nnz(v_in_R) >= N_prop_R
                B_R = phi_prop_R*(reshape(channels.R.sqrt_nu_prop.*exp((-1i*dn)*channels.R.kxdx_prop), N_prop_R, 1).*v_in_R); % use implicit expansion
            else
                B_R = zeros(ny, M_in_R);
                for ii = 1:M_in_R
                    ind = find(v_in_R(:,ii));
                    B_R(:,ii) = channels.fun_phi(channels.R.kydx_prop(ind))*(reshape(channels.R.sqrt_nu_prop(ind).*exp((-1i*dn)*channels.R.kxdx_prop(ind)),[],1).*v_in_R(ind,ii));
                end
            end
        end
    end

    % Build output matrices C_L and C_R
    if build_C
        if use_ind_out % output channels specified by ind_out_L and ind_out_R
            if has_phi_prop_L
                C_L = conj(phi_prop_L(:,ind_out_L));
            else
                C_L = conj(channels.fun_phi(channels.L.kydx_prop(ind_out_L)));
            end
            if two_sided
                if has_phi_prop_R
                    C_R = conj(phi_prop_R(:,ind_out_R));
                else
                    C_R = conj(channels.fun_phi(channels.R.kydx_prop(ind_out_R)));
                end
            end
        else % output wavefronts specified by v_out_L and v_out_R
            if has_phi_prop_L % must be so when nnz(v_out_L) >= N_prop_L
                C_L = conj(phi_prop_L*(reshape(conj(channels.L.sqrt_nu_prop.*exp((-1i*dn)*channels.L.kxdx_prop)), N_prop_L, 1).*v_out_L)); % use implicit expansion
            else
                % build C_L using only the channels used
                C_L = zeros(ny, M_out_L);
                for ii = 1:M_out_L
                    ind = find(v_out_L(:,ii));
                    C_L(:,ii) = conj(channels.fun_phi(channels.L.kydx_prop(ind))*(reshape(conj(channels.L.sqrt_nu_prop(ind).*exp((-1i*dn)*channels.L.kxdx_prop(ind))),[],1).*v_out_L(ind,ii)));
                end
            end
            if two_sided
                if has_phi_prop_R % must be so when nnz(v_out_R) >= N_prop_R
                    C_R = conj(phi_prop_R*(reshape(conj(channels.R.sqrt_nu_prop.*exp((-1i*dn)*channels.R.kxdx_prop)), N_prop_R, 1).*v_out_R)); % use implicit expansion
                else
                    C_R = zeros(ny, M_out_R);
                    for ii = 1:M_out_R
                        ind = find(v_out_R(:,ii));
                        C_R(:,ii) = conj(channels.fun_phi(channels.R.kydx_prop(ind))*(reshape(conj(channels.R.sqrt_nu_prop(ind).*exp((-1i*dn)*channels.R.kxdx_prop(ind))),[],1).*v_out_R(ind,ii)));
                    end
                end
            end
        end
    end
end

if opts.clear_memory && (has_phi_prop_L || has_phi_prop_R); clear phi_prop_L phi_prop_R; end

t3 = clock; timing_build_BC = etime(t3,t2);
if opts.verbal; fprintf('elapsed time: %7.3f secs\n', timing_build_BC); end

%% Part 3: call rgf() or mesti() to do the computation

if use_RGF
    if ~two_sided
        % A Dirichlet boundary condition at n = nx_Ez+1 corresponds to G0_R = 0
        G0_R = sparse(ny, ny);
        B_R = sparse(ny, 0);
        C_R = sparse(ny, 0);
    end
    C_L = C_L.';
    C_R = C_R.';
    syst.yBC = yBC;
    [S, info] = rgf(syst, G0_L, G0_R, B_L, B_R, C_L, C_R, opts);
else
    % Add syst.epsilon_L and syst.epsilon_R, to be used in mesti()
    if use_TM
        if two_sided
            syst.epsilon = [syst.epsilon_L*ones(ny,nx_extra(1)), syst.epsilon, syst.epsilon_R*ones(ny,nx_extra(2))];
        else
            syst.epsilon = [syst.epsilon_L*ones(ny,nx_extra(1)), syst.epsilon];
        end
    else
        % Recall:
        %   size(inv_epsilon{1}) = [ny_Ez, nx_Hz]
        %   size(inv_epsilon{2}) = [ny_Hz, nx_Ez]
        %   size(inv_epsilon{3}) = [ny_Ez, nx_Ez]
        % In syst.inv_epsilon, nx_Hz extends half a pixel beyond the scattering region, with
        %   two-sided: nx_Hz = nx_Ez + 1
        %   one_sided: nx_Hz = nx_Ez
        % But in mesti(), we use xBC = PMC, so nx_Hz = nx_Ez - 1. So, here we pad one less pixel per side to syst.inv_epsilon{1}.
        if two_sided
            syst.inv_epsilon{1} = [(1/syst.epsilon_L)*ones(ny_Ez,nx_extra(1)-1), syst.inv_epsilon{1}, (1/syst.epsilon_R)*ones(ny_Ez,nx_extra(2)-1)];
            syst.inv_epsilon{2} = [(1/syst.epsilon_L)*ones(ny_Hz,nx_extra(1)  ), syst.inv_epsilon{2}, (1/syst.epsilon_R)*ones(ny_Hz,nx_extra(2)  )];
            if numel(syst.inv_epsilon) == 3
                syst.inv_epsilon{3} = [zeros(ny_Ez,nx_extra(1)), syst.inv_epsilon{3}, zeros(ny_Ez,nx_extra(2))];
            end
        else
            syst.inv_epsilon{1} = [(1/syst.epsilon_L)*ones(ny_Ez,nx_extra(1)-1), syst.inv_epsilon{1}];
            syst.inv_epsilon{2} = [(1/syst.epsilon_L)*ones(ny_Hz,nx_extra(1)  ), syst.inv_epsilon{2}];
            if numel(syst.inv_epsilon) == 3
                syst.inv_epsilon{3} = [zeros(ny_Ez,nx_extra(1)), syst.inv_epsilon{3}];
            end
        end
    end
    syst.epsilon_L = []; % mesti() will throw warning if syst.epsilon_L is given
    syst.epsilon_R = []; % mesti() will throw warning if syst.epsilon_R is given

    % The original syst.epsilon or syst.inv_epsilon is no longer needed but still exists in the caller's workspace; we may clear it to reduce memory usage
    if opts.clear_syst
        syst_name = inputname(1); % name of the variable we call syst in the caller's workspace; will be empty if there's no variable for it in the caller's workspace
        if ~isempty(syst_name)
            evalin('caller', ['clear ', syst_name]); % do 'clear syst' in caller's workspace
        end
    end

    % Whether we use self-energy or PML or have a one-sided geometry, the BC in x direction will be Dirichlet when building matrix A in mesti()
    if use_TM
        syst.xBC = 'PEC';
    else
        syst.xBC = 'PMC';
    end

    % location of inputs/outputs on the left surface (n=0 for TM; n=0.5 for TE)
    n_L = nx_extra(1);

    % location of inputs/outputs on the right surface (n=nx_Ez+1 for TM; n=nx_Ez+0.5 for TE)
    n_R = nx_extra(1) + nx_s + 1;

    % Specify inputs
    if two_sided; B = struct('pos', {[],[]}, 'data', {[],[]}); end % pre-allocate
    % inputs on the left surface
    B(1).pos  = [1, n_L, ny, 1];
    B(1).data = reshape(B_L, ny, 1, []); % we don't specify the number of inputs since it can be M_in_L or length(ind_L)
    % inputs on the right surface
    if two_sided
        B(2).pos  = [1, n_R, ny, 1];
        B(2).data = reshape(B_R, ny, 1, []); % we don't specify the number of inputs since it can be M_in_R or length(ind_R)
    end

    % Specify outputs
    if build_C
        if two_sided; C = struct('pos', {[],[]}, 'data', {[],[]}); end % pre-allocate
        % outputs on the left surface
        C(1).pos  = [1, n_L, ny, 1];
        C(1).data = reshape(C_L, ny, 1, M_out_L);
        % outputs on the right surface
        if two_sided
            C(2).pos  = [1, n_R, ny, 1];
            C(2).data = reshape(C_R, ny, 1, M_out_R);
        end
    elseif use_transpose_B
        C = 'transpose(B)';
    else
        C = [];
    end

    if opts.clear_memory
        clear B_L B_R C_L C_R
    end

    % D will be subtracted later
    D = [];

    % variables syst, B, and C can be cleared from mesti() since we don't need them anymore
    opts.clear_syst = opts.clear_memory;
    opts.clear_BC = opts.clear_memory;

    % Main computation happens here
    % Note that we should no longer use syst, B, and C beyond this point since they may be cleared inside mesti()
    [S, info] = mesti(syst, B, C, D, opts);
end

info.timing.init  = info.timing.init  + timing_init;
info.timing.build = info.timing.build + timing_build_G0 + timing_build_BC; % combine with build time for A and K

if info.opts.analysis_only
    t3 = clock;
    info.timing.total = etime(t3,t0);
    if opts.verbal; fprintf('          Total elapsed time: %7.3f secs\n', info.timing.total); end
    return
end

t1 = clock;

% Include the -2i prefactor that should have been in the input matrix B
S = (-2i)*S;

%% Part 4: wrap up

% Recover the original list of input and output channels if we symmetrized K = [A,B;C,0]
if use_transpose_B % when opts.symmetrize_K = true
    % Indices for the original list of input channels on the left
    ind_in = ind_in_out_L(1:M_in_L);

    % Indices for the original list of output channels on the left
    % There is no need to use ind_prop_conj again since the use of C = transpose(B) already compensates the previous use of ind_prop_conj.
    ind_out = ind_in_out_L(M_in_L+(1:M_out_L));

    if two_sided
        % Include channels on the right; note ind_in_out_L and ind_in_out_R are column vectors
        ind_in  = [ind_in ; length(ind_L) + ind_in_out_R(1:M_in_R)];
        ind_out = [ind_out; length(ind_L) + ind_in_out_R(M_in_R+(1:M_out_R))];
    end

    S = S(ind_out, ind_in);
end

% Include the sqrt(nu(a))*exp(-i*kxdx*dn) prefactors that should have been in the input matrix B
if use_ind_in % input channels specified by ind_in_L and ind_in_R
    prefactor = channels.L.sqrt_nu_prop(ind_in_L).*exp((-1i*dn)*channels.L.kxdx_prop(ind_in_L));
    if two_sided
        prefactor = [prefactor, channels.R.sqrt_nu_prop(ind_in_R).*exp((-1i*dn)*channels.R.kxdx_prop(ind_in_R))];
    end
    if opts.return_field_profile
        S = S.*reshape(prefactor, 1, 1, []); % use implicit expansion
    else
        S = S.*prefactor; % use implicit expansion
    end
end

if ~opts.return_field_profile

    % Include the sqrt(nu(b))*exp(-i*kxdx*dn) prefactors that should have been in the output matrix C
    if use_ind_out % output channels specified by ind_out_L and ind_out_R
        prefactor = channels.L.sqrt_nu_prop(ind_out_L).*exp((-1i*dn)*channels.L.kxdx_prop(ind_out_L));
        if two_sided
            prefactor = [prefactor, channels.R.sqrt_nu_prop(ind_out_R).*exp((-1i*dn)*channels.R.kxdx_prop(ind_out_R))];
        end
        S = reshape(prefactor, [], 1).*S; % use implicit expansion
    end

    % Subtract D = C*inv(A_0)*B - S_0 where A_0 is a reference system and S_0 is its scattering matrix
    % When use_ind_in = use_ind_out = true, D is basically the identity matrix.
    % But we need to shift the phase back by dn, which means D_ba = delta_ba*exp(-i*kxdx(a)*2*dn)
    % When user-specified input and output wavefronts are used, we have D_L = (v_out_L')*diag(exp(-i*kxdx*2*dn))*v_in_L
    phase_factor = reshape(exp((-1i*2*dn)*channels.L.kxdx_prop), [], 1);
    if use_ind_in
        D_L = spdiags(phase_factor, 0, N_prop_L, N_prop_L);
        D_L = D_L(:, ind_in_L);
    else
        D_L = phase_factor .* v_in_L; % use implicit expansion
    end
    if use_ind_out
        D_L = D_L(ind_out_L,:);
    else
        D_L = (v_out_L')*D_L;
    end
    if two_sided && (M_in_R ~= 0 || M_out_R ~= 0)
        phase_factor = reshape(exp((-1i*2*dn)*channels.R.kxdx_prop), [], 1);
        if use_ind_in
            D_R = spdiags(phase_factor, 0, N_prop_R, N_prop_R);
            D_R = D_R(:, ind_in_R);
        else
            D_R = phase_factor .* v_in_R; % use implicit expansion
        end
        if use_ind_out
            D_R = D_R(ind_out_R,:);
        else
            D_R = (v_out_R')*D_R;
        end
        S = S - [[D_L; sparse(M_out_R, M_in_L)], [sparse(M_out_L, M_in_R); D_R]];
    else
        S = S - D_L;
    end
else % when opts.return_field_profile = true
    % The S returned by mesti() has size [ny, nx_tot, M_in_L+M_in_R] where nx_tot = nx_s + sum(nx_extra);
    % Here, we remove the npixels_PML and npixels_spacer pixels
    if n_PML > 0
        nx_remove_L = nx_extra(1) - 1; % we keep the surface pixel
        if two_sided
            nx_remove_R = nx_extra(2) - 1;
        else
            nx_remove_R = 0;
        end
        ind_n = (1+nx_remove_L):(nx_tot-nx_remove_R);
        S = S(:, ind_n, :);
    end
    % At this point, S has this number of pixels in x:
    % TM, two-sided: 1+nx_Ez+1
    % TM, one-sided: 1+nx_Ez
    % TE, two-sided: nx_Hz = nx_Ez+1
    % TE, one-sided: nx_Hz = nx_Ez

    if use_TM && opts.nx_L <= 1 && (~two_sided || opts.nx_R <= 1) % nothing to add (except possibly zeros)
        if ~(opts.nx_L == 1 && ((two_sided && opts.nx_R == 1) || (~two_sided && opts.nx_R == 0)))
            % Remove pixels such that we have nx_L pixels of syst.epsilon_L and nx_R pixels of syst.epsilon_R
            if two_sided
                ind_n = (2-opts.nx_L):(1+nx_Ez+opts.nx_R);
                S = S(:, ind_n, :);
            else
                ind_n = (2-opts.nx_L):(1+nx_Ez);
                S = S(:, ind_n, :);
                if opts.nx_R > 0
                    % Add nx_R slices of zero on the right
                    S = cat(2, S, zeros(ny, opts.nx_R, M_in_L));
                end
            end
        end
    elseif ~use_TM && opts.nx_L == 0 && (~two_sided || opts.nx_R == 0) % nothing to add (except possibly zeros)
        if ~two_sided && opts.nx_R > 0
            % Add nx_R slices of zero on the right
            S = cat(2, S, zeros(ny, opts.nx_R, M_in_L));
        end
    else
        if opts.verbal; fprintf('            ... '); end

        % phi is a ny-by-ny unitary matrix where the a-th column is the a-th transverse mode; it includes all the propagating and evanescent channels
        phi = channels.fun_phi(channels.kydx_all);

        % We use phi' to project the field on the left or right surface onto the complete and orthonormal set of transverse modes.
        phi_prime = phi';

        % Ez already includes one extra pixel on the left
        if use_TM
            nx_L_extra = opts.nx_L - 1;
        else
            nx_L_extra = opts.nx_L;
        end

        if nx_L_extra == -1
            % Remove the single pixel of syst.epsilon_L on the left
            S = S(:, 2:end, :);
        elseif nx_L_extra > 0
            % Add pixels such that we have opts.nx_L pixels of syst.epsilon_L on the left
            S_L = zeros(ny, nx_L_extra, size(S,3));
            % The n below is n = n - n_L, where n_L is the location of the inputs/outputs on the left surface.
            n = (-nx_L_extra):1:-1;  % 1-by-nx_L_extra row vector
            kx_x = reshape(channels.L.kxdx_all, ny, 1).*n; % kx*x; ny-by-nx_L_extra matrix through implicit expansion
            exp_pikx = exp( 1i*kx_x); % exp(+i*kx*x)
            exp_mikx = exp(-1i*kx_x); % exp(-i*kx*x)
            u_in = zeros(ny, 1);
            n_L = 1; % index for the inputs/outputs on the left surface
            if M_in_L > 0
                prefactor = reshape(exp((-1i*dn)*channels.L.kxdx_prop)./channels.L.sqrt_nu_prop, N_prop_L, 1);
            end
            for ii = 1:M_in_L % input from left
                u = phi_prime*S(:, n_L, ii);  % u is a ny-by-1 column vector of transverse mode coefficients
                % u_in is the incident wavefront at n_L; note we need to back propagate dn pixel from x=0
                u_in(:) = 0;
                if use_ind_in
                    u_in(channels.L.ind_prop(ind_in_L(ii))) = prefactor(ind_in_L(ii));
                else
                    u_in(channels.L.ind_prop) = prefactor.*v_in_L(:,ii);
                end
                u_out = u - u_in;
                S_L(:, :, ii) = phi*(u_in.*exp_pikx + u_out.*exp_mikx);
            end
            if two_sided
                for ii = 1:M_in_R % input from right
                    u_out = phi_prime*S(:, n_L, M_in_L+ii); % u_out = u because there is no input on the left side
                    S_L(:, :, M_in_L+ii) = phi*(u_out.*exp_mikx);
                end
            end
            S = cat(2, S_L, S);
        end

        if two_sided

            % Ez already includes one extra pixel on the right
            if use_TM
                nx_R_extra = opts.nx_R - 1;
            else
                nx_R_extra = opts.nx_R;
            end

            if nx_R_extra == -1
                % Remove the single pixel of syst.epsilon_R on the right
                S = S(:, 1:(end-1), :);
            elseif nx_R_extra > 0
                % Add pixels such that we have opts.nx_R pixels of syst.epsilon_R on the right
                S_R = zeros(ny, nx_R_extra, size(S,3));
                % The n below is n = n - n_R, where n_R is the location of the inputs/outputs on the right surface.
                n = 1:nx_R_extra;  % 1-by-nx_R_extra row vector
                kx_x = reshape(channels.R.kxdx_all, ny, 1).*n; % kx*x; ny-by-nx_R_extra matrix through implicit expansion
                exp_pikx = exp( 1i*kx_x); % exp(+i*kx*x)
                exp_mikx = exp(-1i*kx_x); % exp(-i*kx*x)
                u_in = zeros(ny, 1);
                n_R = size(S, 2); % index for the inputs/outputs on the right surface
                for ii = 1:M_in_L % input from left
                    u_out = phi_prime*S(:, n_R, ii); % u_out = u because there is no input on the right side
                    S_R(:, :, ii) = phi*(u_out.*exp_pikx);
                end
                if M_in_R > 0
                    prefactor = reshape(exp((-1i*dn)*channels.R.kxdx_prop)./channels.R.sqrt_nu_prop, N_prop_R, 1);
                end
                for ii = 1:M_in_R % input from right
                    u = phi_prime*S(:, n_R, M_in_L+ii);  % u is a ny-by-1 column vector of transverse mode coefficients
                    % u_in is the incident wavefront at n_R; note we need to back propagate dn pixel from x=L
                    u_in(:) = 0;
                    if use_ind_in
                        u_in(channels.R.ind_prop(ind_in_R(ii))) = prefactor(ind_in_R(ii));
                    else
                        u_in(channels.R.ind_prop) = prefactor.*v_in_R(:,ii);
                    end
                    u_out = u - u_in;
                    S_R(:, :, M_in_L+ii) = phi*(u_in.*exp_mikx + u_out.*exp_pikx);
                end
                S = cat(2, S, S_R);
            end
        elseif opts.nx_R > 0 % one-sided
            % Add nx_R slices of zero on the right
            S = cat(2, S, zeros(ny, opts.nx_R, M_in_L));
        end

        t2 = clock;
        if opts.verbal; fprintf('elapsed time: %7.3f secs\n', etime(t2,t1)); end
    end
end

t3 = clock;
info.timing.solve = info.timing.solve + etime(t3,t1); % Add the post-processing time
info.timing.total = etime(t3,t0);

if opts.verbal; fprintf('          Total elapsed time: %7.3f secs\n', info.timing.total); end

end


% Compute scattering matrix using the recursive Green's function (RGF) method.
% RGF proceeds by incrementally adding slices of the system and evaluating the Green's function on the two surfaces.
% This implementation only works for TM polarization.
function [S, info] = rgf(syst, G0_L, G0_R, B_L, B_R, C_L, C_R, opts)

t1 = clock;

% Number of sites in y and x
[ny, nx] = size(syst.epsilon);

% Finite-difference wave operator for one slice, without the (k0dx)^2*epsilon term
A0 = mesti_build_fdfd_matrix(zeros(ny,1), 0, 'PEC', syst.yBC);
k0dx2 = ((2*pi/syst.wavelength)*(syst.dx))^2;

t2 = clock; info.timing.build = etime(t2,t1);  % the build time for B, C will be added in addition to this
t1 = clock;

M_in_L  = size(B_L, 2);
M_in_R  = size(B_R, 2);
M_out_L = size(C_L, 1);
M_out_R = size(C_R, 1);

% Iterate through by adding slices
% To avoid unnecessary steps, we pick the scheme based on which parts of the S-matrix are needed
if opts.verbal
    fprintf('< Method: RGF >\n');
    fprintf('Iterating   ... ');
end
if M_in_R==0 && M_out_R==0
    % input and output both on the left; loop from right to left
    G_LL = G0_R;
    if opts.clear_memory; clear G0_R; evalin('caller', 'clear G0_R'); end
    for n = nx:-1:1
        A_nn = A0 - spdiags(k0dx2*(syst.epsilon(:,n)), 0, ny, ny);
        G_LL = inv(A_nn - G_LL);
    end
    G_LL = (eye(ny) - G0_L * G_LL) \ G0_L;
    if opts.clear_memory; clear A0 A_nn G0_L; evalin('caller', 'clear G0_L'); end
    S = C_L * (G_LL * B_L);
elseif M_in_R==0
    % input from the left; loop from right to left
    G_LL = G0_R; G_RL = G0_R;
    if opts.clear_memory; clear G0_R; evalin('caller', 'clear G0_R'); end
    for n = nx:-1:1
        A_nn = A0 - spdiags(k0dx2*(syst.epsilon(:,n)), 0, ny, ny);
        G_LL = inv(A_nn - G_LL);
        G_RL = G_RL * G_LL;
    end
    G_LL = (eye(ny) - G0_L * G_LL) \ G0_L;
    if opts.clear_memory; clear A0 A_nn G0_L; evalin('caller', 'clear G0_L'); end
    G_RL = G_RL * G_LL;
    S = [(C_L * G_LL) * B_L; C_R * (G_RL * B_L)]; % also works when M_out_L==0
elseif M_in_L==0 && M_out_L==0
    % input and output both on the right; loop from left to right
    G_RR = G0_L;
    if opts.clear_memory; clear G0_L; evalin('caller', 'clear G0_L'); end
    for n = 1:nx
        A_nn = A0 - spdiags(k0dx2*(syst.epsilon(:,n)), 0, ny, ny);
        G_RR = inv(A_nn - G_RR);
    end
    G_RR = (eye(ny) - G0_R * G_RR) \ G0_R;
    if opts.clear_memory; clear A0 A_nn G0_R; evalin('caller', 'clear G0_R'); end
    S = C_R * (G_RR * B_R);
elseif M_in_L==0
    % input from the right; loop from left to right
    G_RR = G0_L; G_LR = G0_L;
    if opts.clear_memory; clear G0_L; evalin('caller', 'clear G0_L'); end
    for n = 1:nx
        A_nn = A0 - spdiags(k0dx2*(syst.epsilon(:,n)), 0, ny, ny);
        G_RR = inv(A_nn - G_RR);
        G_LR = G_LR * G_RR;
    end
    G_RR = (eye(ny) - G0_R * G_RR) \ G0_R;
    G_LR = G_LR * G_RR;
    if opts.clear_memory; clear A0 A_nn G0_R; evalin('caller', 'clear G0_R'); end
    S = [C_L * (G_LR * B_R); (C_R * G_RR) * B_R]; % also works when M_out_R==0
else
    % general case: input from both sides; loop from left to right
    % Initial step: retarded Green's function for an semi-infinite homogeneous space on the left
    G_RR = G0_L; G_RL = G0_L; G_LL = G0_L; G_LR = G0_L;
    if opts.clear_memory; clear G0_L; evalin('caller', 'clear G0_L'); end

    % Main loop of RGF: attach the scattering region slice by slice
    for n = 1:nx
        A_nn = A0 - spdiags(k0dx2*(syst.epsilon(:,n)), 0, ny, ny);
        G_RR = inv(A_nn - G_RR);
        G_RL = G_RR * G_RL;
        G_LL = G_LL + G_LR * G_RL;
        G_LR = G_LR * G_RR;
    end

    % Final step: attach homogeneous space on the right
    G_RR = (eye(ny) - G0_R * G_RR) \ G0_R;
    G_RL = G_RR * G_RL;
    G_LL = G_LL + G_LR * G_RL;
    G_LR = G_LR * G_RR;
    if opts.clear_memory; clear A0 A_nn G0_R; evalin('caller', 'clear G0_R'); end

    % Fisher-Lee relation for S = [[r;t],[tp,rp]]; the -2i factor and the Kronecker delta term will be included later
    % Multiply to the right first, since typically the number of input channels equals or is less than the number of output channels
    S = [[C_L * (G_LL * B_L); C_R * (G_RL * B_L)], [C_L * (G_LR * B_R); C_R * (G_RR * B_R)]];
end
if issparse(S); S = full(S); end
t2 = clock;
info.timing.init  = 0;
info.timing.factorize = 0;
info.timing.solve = etime(t2,t1);
if opts.verbal; fprintf('elapsed time: %7.3f secs\n', info.timing.solve); end
info.opts = opts; % Return the parameters used for user's reference
end
