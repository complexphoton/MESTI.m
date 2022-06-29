function channels = mesti_build_channels(ny, polarization, yBC, k0dx, epsilon_L, epsilon_R, use_continuous_dispersion, m0)
%MESTI_BUILD_CHANNELS Set up properties of channels in the homogeneous space.
%   MESTI_BUILD_CHANNELS(ny, polarization, yBC, k0dx, epsilon_bg) returns a
%   structure containing properties of the propagating and evanescent channels
%   in a discretized homogeneous space with ny pixels in the transverse (y)
%   direction, with the given polarization, boundary condition yBC along y,
%   background relative permittivity epsilon_bg, and dimensionless frequency
%   k0dx = (2*pi/vacuum_wavelength)*dx where dx is the discretization grid size.
%
%   MESTI_BUILD_CHANNELS(ny, polarization, yBC, k0dx, epsilon_L, epsilon_R)
%   returns a structure containing the properties of discretized homogeneous
%   spaces with relative permittivities epsilon_L on the left, epsilon_R on the
%   right.
%
%   MESTI_BUILD_CHANNELS(syst) does the same but with the parameters extracted
%   from structure syst; see 'help mesti2s' for the fields required for
%   structure 'syst'.
%
%   === Input Arguments ===
%   ny (positive integer scalar; required):
%      Number of sites in the transverse (y) direction; ny = ny_Ez for TM
%      polarization, ny = ny_Hz for TE polarization.
%   polarization (character vector; required):
%      Polarization. Possible choices are:
%         'TM' - Ez component of transverse-magnetic field (Hx, Hy, Ez)
%         'TE' - Hz component of transverse-electric field (Ex, Ey, Hz)
%   yBC (character vector or scalar number; required):
%      Boundary condition (BC) at the two ends in y direction, effectively
%      specifying Ez(m,n) or Hz(m,n) at m=0 and m=ny+1 which are one pixel
%      beyond the computation domain. For character inputs, the options are:
%         'periodic'         - f(m+ny,n) = f(m,n)
%         'Dirichlet'        - f(0,n) = f(ny+1,n) = 0
%         'Neumann'          - f(0,n) = f(1,n); f(ny+1,n) = f(ny,n)
%         'DirichletNeumann' - f(0,n) = 0; f(ny+1,n) = f(ny,n)
%         'NeumannDirichlet' - f(0,n) = f(1,n); f(ny+1,n) = 0
%      where f = Ez or Hz depending on the polarization.
%         One can also specify 'PEC', 'PMC', 'PECPMC', or 'PMCPEC'. For TM
%      polarization, PEC is equivalent to Dirichlet for which Ez = 0 at the
%      boundary, and PMC is equivalent to Neumann for which dEz/dx or dEz/dy = 0
%      at the boundary. For TE polarization, PMC is equivalent to Dirichlet for
%      which Hz = 0 at the boundary, and PEC is equivalent to Neumann for which
%      dHz/dx or dHz/dy = 0 at the boundary.
%         When yBC is a numeric scalar, the Bloch periodic boundary condition is
%      used with
%                              f(m+ny,n) = f(m,n)*exp(1i*yBC)
%      In other words, yBC = ky_B*p where ky_B is the Bloch wave number, and p =
%      ny*dx is the periodicity in y.
%   k0dx (numeric scalar, real or complex; required):
%      Dimensionless frequency, k0*dx = (2*pi/vacuum_wavelength)*dx.
%   epsilon_L (numeric scalar, real or complex; required):
%      Relative permittivity of the homogeneous space on the left.
%   epsilon_R (numeric scalar, real or complex, or []; optional):
%      Relative permittivity of the homogeneous space on the right. Only the
%      left side will be considered if epsilon_R is not given or is [].
%   use_continuous_dispersion (logical scalar; optional, defaults to false):
%      Whether to use the dispersion equation of the continuous wave equation
%      when building the input/output channels. Defaults to false, in which case
%      the finite-difference dispersion is used.
%   m0 (real numeric scalar; optional, defaults to 0):
%      Center of the transverse mode profile with periodic or Bloch periodic
%      boundary condition, phi_{m,a} = exp(i*ky(a)*dx*(m-m0))/sqrt(ny), where
%      ky(a) = ky_B + a*(2*pi/ny*dx).
%
%   === Output Arguments ===
%   channels (scalar structure):
%      channels.kydx_all (1-by-ny real row vector):
%         Dimensionless transverse wave number ky*dx for all ny channels,
%         including both propagating and evanescent ones. They are real-valued
%         and are ordered from small to large.
%      channels.fun_phi (function_handle):
%         A function that, given one element of kydx_all as the input, returns
%         its normalized transverse field profile as an ny-by-1 column vector;
%         when the input is a row vector, it returns a matrix where each column
%         is the respective transverse profile. The transverse modes form a
%         complete and orthonormal set, so the ny-by-ny matrix
%         channels.fun_phi(channels.kydx_all) is unitary.
%      channels.L (scalar structure):
%         When epsilon_L and epsilon_R are both given (i.e., epsilon_R is given
%         and is not []), the properties specific to the left and right sides
%         are returned in channels.L and channels.R; channels.L and channels.R
%         are both scalar structures, and their fields are described below.
%            When only epsilon_L = epsilon_bg is given; there will be no
%         channels.L and channels.R; instead, the fields that would have been
%         assigned to channels.L will be assigned to channels directly. For
%         example, channels.L.N_prop will be channels.N_prop instead.
%      channels.L.N_prop (integer scalar):
%         Number of propagating channels.
%      channels.L.kxdx_all (1-by-ny complex row vector):
%         Dimensionless longitudinal wave number kx*dx for all channels,
%         including both propagating and evanescent ones. Due to the
%         discretization, kxdx is equivalent to kxdx + 2*pi. Whenever kxdx is a
%         solution, -kxdx is also a solution. Here, we choose the sign of kxdx
%         such that:
%         When k0dx is real, we have
%            - Propagating channels: 0 < Re(kxdx) < pi, Im(kxdx) = 0.
%            - Evanescent channels: Im(kxdx) >= 0, mod(Re(kxdx),pi) = 0.
%         When k0dx is complex, we analytically continue the above choice onto
%         the complex-frequency plane. Specifically, we pick the kxdx that is
%         continuously connected to one with real k0dx through a vertical line
%         in the complex-(k0dx^2) plane.
%      channels.L.ind_prop (1-by-N_prop integer row vector):
%         Indices of the N_prop propagating channels among all ny channels.
%      channels.L.kxdx_prop (1-by-N_prop real row vector):
%         Dimensionless longitudinal wave number kx*dx for the N_prop propagating
%         channels, equal to channels.L.kxdx_all(channels.L.ind_prop).
%      channels.L.kydx_prop (1-by-N_prop real row vector):
%         Dimensionless transverse wave number ky*dx for the N_prop propagating
%         channels, equal to channels.kydx_all(channels.L.ind_prop).
%      channels.L.sqrt_nu_prop (1-by-N_prop row vector):
%         Normalization factor for the longitudinal flux of the propagating
%         transverse modes: sqrt_nu = sqrt(sin(kxdx)) for TM polarization,
%         sqrt_nu = sqrt(sin(kxdx)/epsilon_L) for TE polarization, such that
%         the x component of the Poynting vector of a transverse mode integrated
%         over y is proportional to sqrt_nu^2. The longitudinal group velocity
%         is v_g = (sin(kxdx)/k0dx)*(c/epsilon_L).
%      channels.L.ind_prop_conj (1-by-N_prop integer row vector; optional):
%         A permutation vector that switches one propagating channel with one
%         having a complex-conjugated transverse profile. In particular,
%            channels.fun_phi(kydx_prop(channels.L.ind_prop_conj))
%         equals
%            conj(channels.fun_phi(kydx_prop)).
%         For periodic boundary, this flips the sign of ky. For Dirichlet and
%         Neumann boundary, fun_phi is real, so there is no permutation.
%            Given a periodic boundary condition in y, when the set of output
%         channels is the same as the set of input channels and they have the
%         same order, the reflection matrix r is not symmetric because the
%         output projection involves a complex conjugation. But
%         r(channels.L.ind_prop_conj,:) is symmetric.
%            For Bloch periodic boundary with nonzero ky_B, complex conjugation
%         maps ky to -ky where ky_B is flipped, so such permutation does not
%         exist, and ind_prop_conj is not given.
%      channels.R (scalar structure; optional):
%         Structure containing properties specific to the right (R) side,
%         similar to channels.R; only provided when epsilon_R is given.
%
%   See also: mesti2s

if nargin == 1
    % Extract parameters from syst
    syst = ny;
    if ~(isstruct(syst) && isscalar(syst)); error('Input argument ''syst'' must be a scalar structure.'); end

    % Assign polarization
    % We don't check syst.epsilon and syst.inv_epsilon since those are checked in mesti() or mesti2s().
    if isfield(syst, 'polarization')
        polarization = syst.polarization;
        if strcmpi(polarization, 'TM')
            ny = size(syst.epsilon, 1); % ny_Ez
        elseif strcmpi(polarization, 'TE')
            ny = size(syst.inv_epsilon{2}, 1); % ny_Hz
        else
            error('syst.polarization must be either ''TM'' or ''TE''.');
        end
    else
        % When syst.polarization is not given, we automatically pick based on whether syst.epsilon or syst.inv_epsilon is given.
        if isfield(syst, 'epsilon') && ~isfield(syst, 'inv_epsilon')
            polarization = 'TM';
            ny = size(syst.epsilon, 1); % ny_Ez
        elseif ~isfield(syst, 'epsilon') && isfield(syst, 'inv_epsilon')
            polarization = 'TE';
            ny = size(syst.inv_epsilon{2}, 1); % ny_Hz
        elseif isfield(syst, 'epsilon') && isfield(syst, 'inv_epsilon')
            error('syst.polarization must be given when syst.epsilon and syst.inv_epsilon both exist.');
        else % neither syst.epsilon nor syst.inv_epsilon exists
            error('Input argument ''syst'' must have field ''epsilon'' or ''inv_epsilon''.');
        end
    end

    if ~isfield(syst, 'epsilon_L'); error('Input argument ''syst'' must have field ''epsilon_L''.'); end
    epsilon_L = syst.epsilon_L;

    if ~isfield(syst, 'wavelength'); error('Input argument ''syst'' must have field ''wavelength''.'); end
    if ~isfield(syst, 'dx'); error('Input argument ''syst'' must have field ''dx''.'); end
    if ~(isnumeric(syst.wavelength) && isscalar(syst.wavelength)); error('syst.wavelength must be a numeric scalar.'); end
    if ~(isreal(syst.dx) && isscalar(syst.dx) && syst.dx > 0); error('syst.dx must be a positive scalar.'); end
    k0dx = (2*pi/syst.wavelength)*(syst.dx);

    % Check boundary condition in y
    if isfield(syst, 'ky_B') && ~isempty(syst.ky_B)
        if ~(isnumeric(syst.ky_B) && isscalar(syst.ky_B))
            error('syst.ky_B must be a numeric scalar, if given.');
        elseif (isfield(syst, 'yBC') && ~isempty(syst.yBC)) && ~strcmpi(syst.yBC, 'Bloch')
            error('When syst.ky_B is given, syst.yBC must be ''Bloch'' if specified.');
        end
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

    if ~isfield(syst, 'epsilon_R')
        epsilon_R = [];
    else
        epsilon_R = syst.epsilon_R;
    end
    if ~isfield(syst, 'use_continuous_dispersion')
        use_continuous_dispersion = [];
    else
        use_continuous_dispersion = syst.use_continuous_dispersion;
    end
    if ~isfield(syst, 'm0')
        m0 = [];
    else
        m0 = syst.m0;
    end
else
    if nargin < 5
        error('Not enough input arguments.');
    end
    if nargin < 6
        epsilon_R = [];
    end
    if nargin < 7
        use_continuous_dispersion = [];
    end
    if nargin < 8
        m0 = [];
    end
end

% Check input arguments and set default values
if ~(isreal(ny) && isscalar(ny) && round(ny)==ny && ny>0)
    error('Input argument ny must be a positive integer scalar.');
elseif ~(strcmpi(polarization, 'TM') || strcmpi(polarization, 'TE'))
    error('Input argument polarization must be either ''TM'' or ''TE''.');
elseif ~((ischar(yBC) && isrow(yBC)) || ((isstring(yBC) || isnumeric(yBC)) && isscalar(yBC)))
    error('Input argument yBC must be a character vector or string, or numeric scalar (for Bloch periodic boundary).');
elseif ~(isscalar(k0dx) && isnumeric(k0dx))
    error('Input argument k0dx must be a numeric scalar.');
elseif ~(isscalar(epsilon_L) && isnumeric(epsilon_L))
    error('Input argument epsilon_L or syst.epsilon_L must be a numeric scalar.');
end
if isempty(epsilon_R)
    two_sided = false;
else
    two_sided = true;
    if ~(isscalar(epsilon_R) && isnumeric(epsilon_R)); error('Input argument epsilon_R or syst.epsilon_R must be a numeric scalar, if given.'); end
end
if isempty(use_continuous_dispersion)
    use_continuous_dispersion = false;
elseif ~(islogical(use_continuous_dispersion) && isscalar(use_continuous_dispersion))
    error('Input argument use_continuous_dispersion or syst.use_continuous_dispersion must be a logical scalar, if given.');
end
if isempty(m0)
    m0 = 0;
elseif ~(isreal(m0) && isscalar(m0))
    error('Input argument m0 or syst.m0 must be a real scalar, if given.');
end

use_TM = strcmpi(polarization, 'TM');

% Convert PEC/PMC to Dirichlet/Neumann based on whether TE or TM polarization is used
if use_TM
    % For TM field, we consider Ez(x,y), so PEC => Dirichlet, PMC => Neumann
    if strcmpi(yBC, 'PEC')
        yBC = 'Dirichlet';
    elseif strcmpi(yBC, 'PMC')
        yBC = 'Neumann';
    elseif strcmpi(yBC, 'PECPMC')
        yBC = 'DirichletNeumann';
    elseif strcmpi(yBC, 'PMCPEC')
        yBC = 'NeumannDirichlet';
    end
else
    % For TE field, we consider Hz(x,y), so PMC => Dirichlet, PEC => Neumann
    if strcmpi(yBC, 'PMC')
        yBC = 'Dirichlet';
    elseif strcmpi(yBC, 'PEC')
        yBC = 'Neumann';
    elseif strcmpi(yBC, 'PMCPEC')
        yBC = 'DirichletNeumann';
    elseif strcmpi(yBC, 'PECPMC')
        yBC = 'NeumannDirichlet';
    end
end

% These are used only for periodic and Bloch periodic boundary conditions; otherwise they stay empty
ka = []; % ky_B*periodicity
ind_zero_ky = [];

% Handle periodic and Bloch periodic boundary conditions
if strcmpi(yBC, 'Bloch')
    error('To use Bloch periodic boundary condition in mesti_build_channels(), set the second input argument yBC to ky_B*p where ky_B is the Bloch wave number and p is the periodicity.');
elseif isnumeric(yBC)
    ka = yBC;
    yBC = 'Bloch';
    % ka must be real for channels.fun_phi(channels.kydx_all) to be unitary
    if ~isreal(ka)
        warning('ky_B*a = %g + 1i*%g is a complex number; must be real for a complete orthonormal transverse basis.', real(ka), imag(ka));
    end
elseif strcmpi(yBC, 'periodic')
    ka = 0;
    yBC = 'Bloch';
end

% f = [f(1), ..., f(ny)].'; 
% For periodic and Bloch periodic boundary, we order channels.kydx_all such that it increases monotonically from negative to positive
% For other boundary conditions, ky >= 0, and we order channels.kydx_all such that it increases monotonically from smallest to largest

% Transverse modes (form a complete basis and are independent of epsilon_L/R)
if strcmpi(yBC, 'Bloch')
    % f(ny+1) = f(1)*exp(1i*ka); f(0) = f(ny)*exp(-1i*ka)
    % The transverse mode index where kydx = ka/ny
    if mod(ny,2) == 1
        ind_zero_ky = round((ny+1)/2);
    else
        ind_zero_ky = round(ny/2);
    end
    channels.kydx_all = (ka/ny) + ((1:ny)-ind_zero_ky)*(2*pi/ny);
    % Dimensionless transverse mode profile: phi_{m,a} = exp(i*(m-m0)*kydx(a))/sqrt(ny)
    channels.fun_phi = @(kydx) exp(((1:ny).'-m0)*(1i*kydx))/sqrt(ny);
elseif strcmpi(yBC, 'Dirichlet') % Dirichlet on both sides
    % f(0) = f(ny+1) = 0
    channels.kydx_all = (1:ny)*(pi/(ny+1));
    % Dimensionless transverse mode profile: phi_{m,a} = sin(m*kydx(a))*sqrt(2/(ny+1))
    channels.fun_phi = @(kydx) sin(((1:ny).')*kydx)*sqrt(2/(ny+1)); 
elseif strcmpi(yBC, 'Neumann') % Neumann on both sides
    % f(0) = f(1), f(ny+1) = f(ny)
    channels.kydx_all = ((1:ny)-1)*(pi/ny);
    % Dimensionless transverse mode profile:
    % When kydx == 0: phi_{m,a} = sqrt(1/ny)
    % When kydx != 0: phi_{m,a} = cos((m-0.5)*kydx(a))*sqrt(2/ny)
    % We subtract (~kydx)*(1-sqrt(1/2)) from the cos() which is nonzero only when kydx=0
    channels.fun_phi = @(kydx) (cos(((0.5:ny).')*kydx)-((~kydx)*(1-sqrt(1/2))))*sqrt(2/ny); 
elseif strcmpi(yBC, 'DirichletNeumann') % Dirichlet on the low side, Neumann on the high side
    % f(0) = 0, f(ny+1) = f(ny)
    channels.kydx_all = (0.5:ny)*(pi/(ny+0.5));
    % Dimensionless transverse mode profile: phi_{m,a} = sin(m*kydx(a))*sqrt(2/(ny+0.5))
    channels.fun_phi = @(kydx) sin(((1:ny).')*kydx)*sqrt(2/(ny+0.5)); 
elseif strcmpi(yBC, 'NeumannDirichlet') % Neumann on the low side, Dirichlet on the high side
    % f(0) = f(1), f(ny+1) = 0
    channels.kydx_all = (0.5:ny)*(pi/(ny+0.5));
    % Dimensionless transverse mode profile: phi_{m,a} = cos((m-0.5)*kydx(a))*sqrt(2/(ny+0.5))
    channels.fun_phi = @(kydx) cos(((0.5:ny).')*kydx)*sqrt(2/(ny+0.5)); 
else
    error('Input argument yBC = ''%s'' is not a supported option.', yBC);
end

% Properties for the homogeneous space on the left (kxdx, sqrt_nu_prop, number of propagating channels, etc; depends on epsilon_L/R)
side = setup_longitudinal(k0dx, epsilon_L, channels.kydx_all, use_TM, ka, ind_zero_ky, use_continuous_dispersion);

if two_sided
    channels.L = side;
    % Homogeneous space on the right
    if epsilon_R == epsilon_L
        channels.R = side;
    elseif ~isnan(epsilon_R)
        channels.R = setup_longitudinal(k0dx, epsilon_R, channels.kydx_all, use_TM, ka, ind_zero_ky, use_continuous_dispersion);
    end
else
    % add the fields of 'side' to 'channels'
    fnames = fieldnames(side);
    for ii = 1:length(fnames)
        channels.(fnames{ii}) = side.(fnames{ii});
    end
end

end


function side = setup_longitudinal(k0dx, epsilon_bg, kydx_all, use_TM, ka, ind_zero_ky, use_continuous_dispersion)
% Returns a structure 'side'. See comments at the beginning of this file for more info.

k0dx2_epsilon = (k0dx^2)*epsilon_bg;

if ~use_continuous_dispersion
    % use the finite-difference dispersion for homogeneous space: k0dx2_epsilon = 4*sin^2(kxdx/2) + 4*sin^2(kydx/2)

    % sin_kxdx_over_two_sq = sin^2(kxdx/2)
    sin_kxdx_over_two_sq = 0.25*k0dx2_epsilon - sin(kydx_all/2).^2;

    % Dimensionless longitudinal wave number
    % asin(sqrt(z)) has two branch points (at z=0 and z=1) and with the branch cuts going outward on the real-z axis; we will address the branch choice below
    % Note kxdx is only defined up to modulo 2*pi (ie, kxdx is equivalent to kxdx + 2*pi, kxdx - 2*pi, etc) because sin(kxdx) and exp(1i*kxdx) are both invariant under 2pi shifts.
    side.kxdx_all = 2*asin(sqrt(sin_kxdx_over_two_sq));

    % Indices of the propagating channels
    % When k0dx2_epsilon is real, these are indices of the channels with real-valued kxdx
    % When k0dx2_epsilon is complex, these are indices of the channels we consider "propagating-like"; they have complex kxdx with 0 < real(kxdx) < pi. When k0dx2_epsilon is tuned to a real number continuously, this set continuously becomes that at real k0dx2_epsilon.
    side.ind_prop = find((real(sin_kxdx_over_two_sq) > 0) & (real(sin_kxdx_over_two_sq) < 1));

    % Here we address the sign choice of kxdx, namely its branch
    % When k0dx2_epsilon is real, we choose the sign of kxdx such that:
    % 1) 0 < kxdx < pi for propagating channels (where kxdx is real)
    % 2) Im(kxdx) >= 0 for evanescent channels
    % Using the correct sign is important when we build the retarded Green's function of the semi-infinite homogeneous space.
    % The default branch choice of asin(sqrt(z)) returns the correct sign for the most part, except when z > 1. We need to flip the sign of those (which can only occur if k0dx2_epsilon > 4).
    % When k0dx2_epsilon is complex-valued, it is not always possible to unambiguously choose the sign that is "physical", because kxdx will be complex-valued, and the sign we "want" for real(kxdx) and the sign we want for imag(kxdx) may be incompatible.
    % What we do with complex k0dx2_epsilon is that we choose the sign for the (complex-valued) kxdx such that when k0dx2_epsilon is tuned to a real number continuously by fixing Re(k0dx2_epsilon) and varying Im(k0dx2_epsilon), the kxdx we choose continuously becomes the "correct" one at real k0dx2_epsilon without crossing any branch cut. To do so, we rotate the two branch cuts of asin(sqrt(z)) by 90 degrees to the lower part of the complex-z plane (ie, the lower part of the complex-k0dx2_epsilon plane), and we pick the branch with the correct sign when k0dx2_epsilon is real. This is implemented by flipping the sign of kxdx for the ones that require flipping.
    % Note that we will get a discontinuity whenever k0dx2_epsilon crosses one of those vertical-pointing branch cuts. That is unavoidable.
    % The following few lines implement the "flipping".
    if ~isreal(k0dx2_epsilon) || (isreal(k0dx2_epsilon) && k0dx2_epsilon > 4)
        % Note that when imag(sin_kxdx_over_two_sq)=0, flipping is needed for sin_kxdx_over_two_sq>1 but not needed for sin_kxdx_over_two_sq<0.
        ind_flip = find((real(sin_kxdx_over_two_sq)<0 & imag(sin_kxdx_over_two_sq)<0) | (real(sin_kxdx_over_two_sq)>1 & imag(sin_kxdx_over_two_sq)<=0));
        side.kxdx_all(ind_flip) = -side.kxdx_all(ind_flip);
    end
else
    % use the continuous dispersion for homogeneous space: k0dx2_epsilon = kxdx^2 + kydx^2
    kxdx2 = k0dx2_epsilon - kydx_all.^2;
    side.kxdx_all = sqrt(kxdx2);   
    side.ind_prop = find(real(kxdx2) > 0);
    if ~isreal(k0dx2_epsilon)
        ind_flip = find(real(kxdx2)<0 & imag(kxdx2)<0);
        side.kxdx_all(ind_flip) = -side.kxdx_all(ind_flip);
    end
end

% Number of propagating channels
side.N_prop = length(side.ind_prop);
if side.N_prop==0 && length(kydx_all)==1; side.ind_prop = zeros(1,0); end  % a rare possibility, but in this case ind_prop would be zeros(0,0) while it should be zeros(1,0)

% Wave numbers of the propagating channels
side.kxdx_prop = side.kxdx_all(side.ind_prop);
side.kydx_prop = kydx_all(side.ind_prop);
    
% Square root of the normalized longitudinal group velocity, for the propagating channels
%  TM: nu = sin(kxdx)
%  TE: nu = sin(kxdx)/epsilon_bg
% When k0dx2_epsilon is real, sqrt_nu_prop is also real. When k0dx2_epsilon is complex, sqrt_nu_prop is also complex.
if use_TM
    side.sqrt_nu_prop = sqrt(sin(side.kxdx_prop));
else
    side.sqrt_nu_prop = sqrt(sin(side.kxdx_prop)/epsilon_bg);
end

% Permutation that switches one propagating channel with one having a complex-conjugated transverse profile.
if isempty(ka)
    % For Dirichlet and Neumann boundaries, fun_phi is real, so no permutation needed
    side.ind_prop_conj = 1:side.N_prop;
elseif ka == 0
    % For periodic boundary condition, complex conjugation switches ky and -ky
    if ismember(ind_zero_ky, side.ind_prop) || (mod(side.N_prop,2)==0)
        % Simply flip the ordering
        side.ind_prop_conj = side.N_prop:-1:1;
    else
        % The last channel has -ky equal to ky due to aliasing so should not be flipped
        side.ind_prop_conj = [(side.N_prop-1):-1:1, side.N_prop];
    end
    % TODO: implement side.ind_prop_conj when ka == pi
end

end
