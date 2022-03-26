function [A, is_symmetric_A, xPML, yPML] = mesti_build_fdfd_matrix(epsilon, k0dx, xBC, yBC, xPML, yPML, use_UPML)
% TODO: Documentation.

% Check input arguments
% xBC and yBC are checked in build_laplacian_1d(); xPML and yPML are checked in set_PML_params()
if nargin < 4; error('Not enough input arguments.'); end
if ~(ismatrix(epsilon) && isnumeric(epsilon)); error('Input argument epsilon must be a numeric matrix.'); end
if ~(isscalar(k0dx) && isnumeric(k0dx)); error('Input argument k0dx must be a numeric scalar.'); end

% Assign default value for optional inputs
if nargin < 5; xPML = []; end
if nargin < 6; yPML = []; end
if nargin < 7; use_UPML = true; end % use UPML instead of SC-PML by default
if ~(islogical(use_UPML) && isscalar(use_UPML))
    error('Input argument use_UPML must be a logical scalar, if given.');
end

% Number of grid points in y and x
[ny, nx] = size(epsilon);
nxy = nx*ny;

% Nothing to do if the size is zero
if nxy == 0
    A = sparse(0, 0);
    return
end

% Set default values for PML parameters
xPML = set_PML_params(xPML, k0dx, real([mean(epsilon(:,1)), mean(epsilon(:,end))]), 'x');
yPML = set_PML_params(yPML, k0dx, real([mean(epsilon(1,:)), mean(epsilon(end,:))]), 'y');

% Second derivatives in x and y directions
[laplacian_x, s_x] = build_laplacian_1d(nx, xBC, xPML, 'x', use_UPML);
[laplacian_y, s_y] = build_laplacian_1d(ny, yBC, yPML, 'y', use_UPML);

% Use Kronecker outer product to go from 1D to 2D
% A = [-(d/dx)^2 - (d/dy)^2 - k^2*epsilon(x,y)]*(dx^2)
if use_UPML
    % For uniaxial PML (UPML), the whole operator is multiplied by s_d2f
    s_xy = s_y(:).*reshape(s_x, 1, nx);
    A = -kron(laplacian_x, spdiags(s_y,0,ny,ny)) - kron(spdiags(s_x,0,nx,nx), laplacian_y) - spdiags((k0dx^2)*(s_xy(:).*epsilon(:)), 0, nxy, nxy);
else
    % For stretched-coordinate PML (SC-PML), 1/s_df and 1/s_d2f are both multiplied onto the gradient 
    A = -kron(laplacian_x, speye(ny)) - kron(speye(nx), laplacian_y) - spdiags((k0dx^2)*epsilon(:), 0, nxy, nxy);
end

if (isnumeric(xBC) && xBC ~= 0 && xBC ~= pi) || (isnumeric(yBC) && yBC ~= 0 && yBC ~= pi)
    % Bloch periodic boundary condition with ka != 0 or pi breaks the symmetry of A
    is_symmetric_A = false;
elseif (isempty(xPML) && isempty(yPML)) || use_UPML
    is_symmetric_A = true;
else
    % SC-PML also breaks the symmetry of A
    is_symmetric_A = false;
end

end


function [A, s_d2f] = build_laplacian_1d(n, BC, PML, direction, use_UPML)
% Builds dx^2 times the Laplacian on a 1D finite-difference grid, namely [(d/dx)^2]*(dx^2), with the given boundary conditions and with PML via coordinate transformation.
% n: total number of sites (including PML sites)

if ~((ischar(BC) && isrow(BC)) || ((isstring(BC) || isnumeric(BC)) && isscalar(BC)))
    error('Input argument %sBC must be a character vector or string, or numeric scalar (for Bloch periodic BC).', direction);
end

% Handle periodic and Bloch periodic boundary conditions
if strcmpi(BC, 'Bloch')
    error('To use Bloch periodic boundary condition, set BC to k_B*p where k_B is the Bloch wave number and p is the periodicity.');
elseif isnumeric(BC)
    ka = BC;
    BC = 'Bloch';
    if ~isreal(ka)
        warning('k%s_B*a = %g + 1i*%g is a complex number.', direction, real(ka), imag(ka));
    end
elseif strcmpi(BC, 'periodic')
    ka = 0;
    BC = 'Bloch';
end

% f = [f(1), ..., f(n)].'; 

% First derivative of f
if strcmpi(BC, 'Bloch')
    % f(n+1) = f(1)*exp(1i*ka); f(0) = f(n)*exp(-1i*ka)
    % grad_1*f = df = [df(0.5), ..., df(n-0.5)].'
    grad_1 = spdiags([ones(n,1),-ones(n,1),-exp(-1i*ka)*ones(n,1)], [0,-1,n-1], n, n);
elseif strcmpi(BC, 'Dirichlet') || strcmpi(BC, 'PEC') % PEC on both sides
    % f(0) = f(n+1) = 0
    % grad_1*f = df = [df(0.5), ..., df(n+0.5)].'
    grad_1 = spdiags([ones(n,1),-ones(n,1)], [0,-1], n+1, n);
elseif strcmpi(BC, 'Neumann') || strcmpi(BC, 'PMC') % PMC on both sides
    % df(0.5) = df(n+0.5) = 0
    % grad_1*f = df = [df(1.5), ..., df(n-0.5)].'
    grad_1 = spdiags([ones(n-1,1),-ones(n-1,1)], [1,0], n-1, n);
elseif strcmpi(BC, 'DirichletNeumann') || strcmpi(BC, 'PECPMC') % PEC on the low side, PMC on the high side
    % f(0) = 0; df(n+0.5) = 0
    % grad_1*f = df = [df(0.5), ..., df(n-0.5)].'
    grad_1 = spdiags([ones(n,1),-ones(n,1)], [0,-1], n, n);
elseif strcmpi(BC, 'NeumannDirichlet') || strcmpi(BC, 'PMCPEC') % PMC on the low side, PEC on the high side
    % df(0.5) = 0; f(n+1) = 0
    % grad_1*f = df = [df(1.5), ..., df(n+0.5)].'
    grad_1 = spdiags([ones(n,1),-ones(n,1)], [1,0], n, n);
else
    error('Input argument %sBC = %s is not a supported option.', direction, BC);
end

% grad_2*df = d2f = [d2f(1), ..., d2f(n)].'
grad_2 = -ctranspose(grad_1);

% No coordinate transformation if no PML
if isempty(PML)
    A = grad_2*grad_1;
    s_d2f = ones(n,1);
    return
end

npixels = [PML{1}.npixels, PML{2}.npixels];

% Cannot have more PML pixels than number of pixels
sides = {'low', 'high'};
for ii = 1:2
    if npixels(ii) > n
        error('Specified %d PML pixels on the %s side while there are only %d pixels in %s direction.', npixels(ii), sides{ii}, n, direction);
    end
end

% Shouldn't have the two PMLs overlap
% In our implementation, the conductivity goes to zero at one pixel before PML. So there needs to be at least one more pixel in addition to PML
if sum(npixels) >= n
    warning('Total number of PML pixels = %d + %d = %d should be less than the number of pixels = %d in %s direction.', npixels(1), npixels(2), sum(npixels), n, direction);
end

% s-factor s(x) := 1 + 1i*sigma(x)/omega, where sigma(x) is the conductivity profile
n_df = size(grad_1,1); % number of sites for df
s_df = ones(n_df,1); % s factor for df (on half-integer sites)
s_d2f = ones(n,1); % s factor for d2f (on integer sites)

% Below, u(x) is a function that goes linearly from u(x)=0 at one site before PML to u(x)=1 at the end of PML.
% Let index i=0 be one site before PML, i=1 be first site of PML, so u(i=0)=0.
% Note that even thouugh u(i=npixels) is the last site of PML on the lattice we explicitly simulate, we do not set u=1 there.
% For Dirichlet BC, the BC is such that f=0 at i=npixels+1, so we let u(i=npixels+1)=1, with an effective PML thickness of npixels+1 pixels.
% For Neumann BC, the BC is such that df=0 at i=npixels+0.5, so we let u(i=npixels+0.5)=1, with an effective PML thickness of npixels+0.5 pixels.
% For periodic BC, the BC is such that df(n+0.5)=df(0.5)*exp(1i*ka), so we also let u(i=npixels+0.5)=1 so that if PML is used on both sides, we will have u(0.5)=u(n+0.5)=1 so the conductivity profile is continuous across the edge if sigma_max is the same on both sides.

% Construct u(x) and their corresponding indices for df on half-integer sites
if strcmpi(BC, 'Bloch')
    npixels_effective = [npixels(1)+0.5, npixels(2)+0.5];
    if npixels(1)*npixels(2) == 0
        warning('Bloch periodic boundary condition is used with a single-sided PML in %s direction; transmission through PML will only undergo single-pass attenuation.', direction);
    end
    % No s-factor for df(n+0.5) since we only consider df(0.5)
    u_df = {((1:(npixels(1)+1))-0.5)/npixels_effective(1), ((1:npixels(2))-0.5)/npixels_effective(2)};
    ind_df = {fliplr(1:(npixels(1)+1)), (n+1)-fliplr(1:npixels(2))};
elseif strcmpi(BC, 'Dirichlet') || strcmpi(BC, 'PEC') % PEC on both sides
    npixels_effective = [npixels(1)+1, npixels(2)+1];
    u_df = {((1:(npixels(1)+1))-0.5)/npixels_effective(1), ((1:(npixels(2)+1))-0.5)/npixels_effective(2)};
    ind_df = {fliplr(1:(npixels(1)+1)), (n+2)-fliplr(1:(npixels(2)+1))};
elseif strcmpi(BC, 'Neumann') || strcmpi(BC, 'PMC') % PMC on both sides
    npixels_effective = [npixels(1)+0.5, npixels(2)+0.5];
    u_df = {((1:npixels(1))-0.5)/npixels_effective(1), ((1:npixels(2))-0.5)/npixels_effective(2)};
    ind_df = {fliplr(1:npixels(1)), n-fliplr(1:npixels(2))};
elseif strcmpi(BC, 'DirichletNeumann') || strcmpi(BC, 'PECPMC') % PEC on the low side, PMC on the high side
    npixels_effective = [npixels(1)+1, npixels(2)+0.5];
    u_df = {((1:(npixels(1)+1))-0.5)/npixels_effective(1), ((1:npixels(2))-0.5)/npixels_effective(2)};
    ind_df = {fliplr(1:(npixels(1)+1)), (n+1)-fliplr(1:npixels(2))};
elseif strcmpi(BC, 'NeumannDirichlet') || strcmpi(BC, 'PMCPEC') % PMC on the low side, PEC on the high side
    npixels_effective = [npixels(1)+0.5, npixels(2)+1];
    u_df = {((1:npixels(1))-0.5)/npixels_effective(1), ((1:(npixels(2)+1))-0.5)/npixels_effective(2)};
    ind_df = {fliplr(1:npixels(1)), (n+1)-fliplr(1:(npixels(2)+1))};
else
    error('Input argument %sBC = ''%s'' is not a supported option.', direction, BC);
end

% Construct u(x) and their corresponding indices for d2f on integer sites
u_d2f = {(1:npixels(1))/npixels_effective(1), (1:npixels(2))/npixels_effective(2)};
ind_d2f = {fliplr(1:npixels(1)), (n+1)-fliplr(1:npixels(2))};

% Loop over PML on the two sides
for ii = 1:2
    if npixels(ii) > 0
        % s-factor s(u), ie coordinate transformation profile
        % Eq 7.73 of Taflove & Hagness's 2005 FDTD book
        % kappa is real-coordinate stretching, used to accelerate attenuation of evanescent waves.
        % alpha is used for complex frequency shifting (CFS) to suppress reflection of low-frequency components for time-domain simulations.
        % In general, kappa, sigma, alpha can all be arbitrary functions of position.
        % To minimize discretization-induced reflection, kappa should start from 1, and sigma should start from 0.
        % We use a polynomial grading for all of them: Eqs 7.60 and 7.79 of Taflove & Hagness's 2005 FDTD book
        kappa = @(u) 1 + (PML{ii}.kappa_max-1)*(u.^(PML{ii}.power_kappa));
        sigma_over_omega = @(u) PML{ii}.sigma_max_over_omega*(u.^(PML{ii}.power_sigma));
        alpha_over_omega = @(u) PML{ii}.alpha_max_over_omega*((1-u).^(PML{ii}.power_alpha));
        func_s = @(u) kappa(u) + sigma_over_omega(u)./(alpha_over_omega(u) - 1i);

        % Evaluate s(u) for df (on half-integer sites) and for d2f (on integer sites)
        % It is important that s(u) on both integer and half-integer sites are sampled directly from the polynomial functions. If we use the polynomial for one and linear interpolation for the other, the performance will be much worse.
        s_df(ind_df{ii}) = func_s(u_df{ii}).';  % column vector
        s_d2f(ind_d2f{ii}) = func_s(u_d2f{ii}).';  % column vector
    end
end

if use_UPML
    % For uniaxial PML (UPML), s_d2f will be multiplied later
    A = grad_2*spdiags(1./s_df, 0, n_df, n_df)*grad_1;
else
    % For stretched-coordinate PML (SC-PML), 1/s_df and 1/s_d2f are both multiplied onto the gradient 
    A = spdiags(1./s_d2f, 0, n, n)*grad_2*spdiags(1./s_df, 0, n_df, n_df)*grad_1;
end

end


function PML = set_PML_params(PML, k0dx, epsilon_bg, direction)
% Set default values for PML parameters

if isempty(PML)
    return
elseif ~((isstruct(PML) && isscalar(PML)) || (iscell(PML) && numel(PML) == 2))
    error('Input argument %sPML must be a scalar structure or a two-element cell array.', direction)
elseif isstruct(PML)
    % Use the same PML parameters on both sides if only one set of parameters is given
    PML = {PML, PML};
end

% Loop over PML parameters on the two sides
for ii = 1:2
    % Empty means no PML, and we let npixels = 0
    if isempty(PML{ii})
        PML{ii}.npixels = 0;
        continue
    elseif ~(isstruct(PML{ii}) && isscalar(PML{ii}))
        error('%sPML{%d} must be a scalar structure.', direction, ii);
    end

    % Number of PML pixels
    if ~isfield(PML{ii}, 'npixels') || isempty(PML{ii}.npixels)
        error('%sPML{%d} must contain field ''npixels''.', direction, ii);
    else
        temp = PML{ii}.npixels;
        if ~(isreal(temp) && isscalar(temp) && temp >= 0 && round(temp) == temp)
            error('%sPML{%d}.npixels must be a non-negative integer scalar.', direction, ii);
        end
        if temp == 0
            continue
        end
    end

    % Power of polynomial grading for the conductivity sigma
    if ~isfield(PML{ii}, 'power_sigma') || isempty(PML{ii}.power_sigma)
        % From Table 7.1 of Taflove & Hagness's 2005 FDTD book
        PML{ii}.power_sigma = 3;
    else
        temp = PML{ii}.power_sigma;
        if ~(isreal(temp) && isscalar(temp) && temp >= 0)
            error('%sPML{%d}.power_sigma must be a non-negative scalar.', direction, ii);
        end
    end

    % Conductivity at the end of the PML
    if ~isfield(PML{ii}, 'sigma_max_over_omega') || isempty(PML{ii}.sigma_max_over_omega)
        % Eq 7.67 of Taflove & Hagness's 2005 FDTD book
        PML{ii}.sigma_max_over_omega = 0.8*(PML{ii}.power_sigma+1)/(k0dx*sqrt(epsilon_bg(ii)));
    else
        temp = PML{ii}.sigma_max_over_omega;
        if ~(isreal(temp) && isscalar(temp) && temp >= 0)
            error('%sPML{%d}.sigma_max_over_omega must be a non-negative scalar.', direction, ii);
        end
    end

    % Maximal coordinate stretching factor
    if ~isfield(PML{ii}, 'kappa_max') || isempty(PML{ii}.kappa_max)
        % From Table 7.1 of Taflove & Hagness's 2005 FDTD book
        PML{ii}.kappa_max = 15;
    else
        temp = PML{ii}.kappa_max;
        if ~(isreal(temp) && isscalar(temp) && temp >= 1)
            error('%sPML{%d}.kappa_max must be a real scalar that equals or is larger than 1.', direction, ii);
        end
    end

    % Power of polynomial grading for the coordinate stretching factor kappa
    if ~isfield(PML{ii}, 'power_kappa') || isempty(PML{ii}.power_kappa)
        % From Table 7.1 of Taflove & Hagness's 2005 FDTD book
        PML{ii}.power_kappa = 3;
    else
        temp = PML{ii}.power_kappa;
        if ~(isreal(temp) && isscalar(temp) && temp >= 0)
            error('%sPML{%d}.power_kappa must be a non-negative scalar.', direction, ii);
        end
    end

    % Maximal alpha factor for complex frequency shifting (CFS)
    if ~isfield(PML{ii}, 'alpha_max_over_omega') || isempty(PML{ii}.alpha_max_over_omega)
        % CFS is meant to suppress reflection of low-frequency components for time-domain simulations; it is not necessary for frequency-domain simulations
        PML{ii}.alpha_max_over_omega = 0;
    else
        temp = PML{ii}.alpha_max_over_omega;
        if ~(isreal(temp) && isscalar(temp) && temp >= 0)
            error('%sPML{%d}.alpha_max_over_omega must be a non-negative scalar.', direction, ii);
        end
    end

    % Power of polynomial grading for the alpha factor
    if ~isfield(PML{ii}, 'power_alpha') || isempty(PML{ii}.power_alpha)
        % not relevant when alpha_max_over_omega = 0
        PML{ii}.power_alpha = 1;
    else
        temp = PML{ii}.power_alpha;
        if ~(isreal(temp) && isscalar(temp) && temp >= 0)
            error('%sPML{%d}.power_alpha must be a non-negative scalar.', direction, ii);
        end
    end

    % At this point, PML{ii} should contain 7 fields; check that no other fields exist
    if numel(fieldnames(PML{ii})) > 7
        field_names = setdiff(fieldnames(PML{ii}), {'npixels', 'power_sigma', 'sigma_max_over_omega', 'power_kappa', 'kappa_max', 'power_alpha', 'alpha_max_over_omega'});
        for jj = 1:numel(field_names)
            warning('%sPML{%d} contains unrecognized field ''%s''.', direction, ii, field_names{jj});
        end
    end
end

if PML{1}.npixels + PML{2}.npixels == 0
    PML = [];
end

end
