%% Angle dependence of a mm-wide metalens
% In this example, we:
%  - Build a mm-diameter hyperbolic metalens based on the [meta-atom design example](../2d_meta_atom).
%  - Use mesti() to compute its transmission matrix, using compressed input/output matrices (APF-c).
%  - Use angular spectrum propagation to obtain field profile away from the metalens. 
%  - Map out the transmission efficiency and Strehl ratio for all incident angles.
% The details of APF-c are described in Supplementary Section 5 of the APF paper.

clear 

%% System parameters and general setup

% Parameters for the structure
n_air    = 1.0;     % Refractive index of air on the right
n_silica = 1.46;    % Refractive index of silica substrate on the left
n_TiO2   = 2.43;    % Refractive index of TiO2 ridges
wavelength = 0.532; % Vacuum wavelength [micron]
W  = 1000;          % Approximate width of the metalens [micron]
focal_length = 300; % Focal length of the metalens [micron]
dx = wavelength/40; % Discretization grid size [micron]
Lambda = 18*dx;     % Width (i.e., periodicity) of each unit cell [micron]
L  = 0.6;           % Thickness of the metalens [micron]

n_meta_atom = ceil(W/Lambda); % Number of meta-atoms
W = n_meta_atom*Lambda;       % Actual width of the metalens [micron]

% Parameters for the input and output
W_in  = W;                 % Width of the incident wave [micron] 
W_out = W + 40*wavelength; % Width where we sample the transmitted field [micron]
                           % (a larger W_out ensures all transmitted light is captured)

% Whether or not to compress/decompress the input/output matrices (APF-c).
% Set do_compression = false for exact computation without compression.
do_compression = true;

% Parameters for input/output matrix compression
if do_compression
    % The compression procedure is described in supplementary section 5 of the APF paper:
    % 1. Optionally pad extra channels to the input/output.
    % 2. Optionally scale the the channels with weights q.
    % 3. Fourier transform the input/output matrices so they become spatially localized.
    % 4. Truncate the transformed input/output matrices within a window size.
    % These parameters should be chosen based on the desired level of accuracy;
    % see Supplementary Fig 10 of the APF paper.

    % Whether or not to use the Hann window function for the scaling weight q
    % Set use_Hann_window = false to skip scaling.
    % As shown in Supplementary Fig 10, using the Hann window generally reduces
    % the compression error significantly.
    use_Hann_window = true;

    % pad_ratio = (number of extra channels to pad)/(number of channels of interest)
    % Set pad_ratio = 0 to skip padding.
    % When Hann window is used, it is important to also pad extra channels;
    % otherwise the compression error may go up instead of down.
    % But too much padding can increase computing time and memory usage if
    % nnz(B) or numel(S) grows larger than nnz(A). A padding_ratio of 0.5 is
    % a reasonable choice.
    pad_ratio_L = 0.5;
    pad_ratio_R = 0.5;

    % trunc_width_over_lambda = (truncation window width)/wavelength
    % Larger truncation window gives lower compression error but increases
    % nnz(B). The peak width of the Fourier-transformed input/output matrices is
    % inversely proportional to the momentum space range spanned by the
    % input/output plane waves (counting both the target channels and the padded
    % channels), so if fewer input/output channels are considered, a larger
    % truncation width should be used. A window of 10*lambda is typically more
    % accurate than enough.
    trunc_width_over_lambda_L = 10;
    trunc_width_over_lambda_R = 10;
end

syst.length_unit = 'µm';
syst.wavelength = wavelength;
syst.dx = dx;
syst.PML.npixels = 20;  % Number of PML pixels (on all four sides)

%% Build the relative permittivity profile of the system

% Number of pixels
nx_unitcell = ceil(L/dx);              % Number of pixels of one unit cell in the x direction
ny_unitcell = Lambda/dx;               % Number of pixels of one unit cell in the y direction
nx_metalens = nx_unitcell;             % Number of pixels of the metalens in the x direction
ny_metalens = n_meta_atom*ny_unitcell; % Number of pixels of the metalens in the y direction
                       
% Mid-cell position of the meta-atoms in the y direction [micron] 
y_mid_cell = ((0.5:n_meta_atom) - (n_meta_atom/2))*Lambda;

% Target transmission phase shift of each meta-atom, meant to focus normal-incident light
target_phase = (2*pi/wavelength)*(focal_length - sqrt(focal_length^2 + y_mid_cell.^2));

% Load the ridge widths and associated phase shifts from the meta-atom example.
% ridge_width_design_list is a row vector storing the ridge widths of 8 meta-atoms, in nm.
% phi_over_pi_design_list = ((pi/4)*(0:7))/pi are their relative phase shifts.
load('../2d_meta_atom/meta_atom.mat', 'ridge_width_design_list', 'phi_over_pi_design_list');
ridge_width_design_list = ridge_width_design_list/1e3; % convert from nm to micron
ridge_height = L; % Ridge height is the thickness of the metalens
n_phases = numel(ridge_width_design_list); % Number of discrete phase shifts

% epsilon_library(:,:,ii) stores the permittivity profile of the ii-th meta-atom in that list
epsilon_library = zeros(ny_unitcell, nx_unitcell, n_phases);
for ii = 1:n_phases
    % Build the permittivity profile for ii-th meta-atom.
    epsilon_library(:,:,ii)= build_epsilon_meta_atom(dx, n_air, n_TiO2, ...
        ridge_width_design_list(ii), ridge_height, Lambda);
end

% For each meta-atom in the metalens, assign it to be one of those 8 based on its target phase shift
% The 2*pi phase (same as 0 phase) is added for the interpolation
ind_meta_atom = interp1( ...
    [phi_over_pi_design_list, 2], ...
    [1:n_phases, 1], ...
    mod(target_phase, 2*pi)/pi, 'nearest');

% Construct the permittivity profile of the metalens by assigning the permittivity 
% profile of each meta-atom to its position on the metalens
epsilon_metalens = zeros(ny_metalens, nx_metalens);
for ii = 1:n_meta_atom
    ind_y = (ii-1)*ny_unitcell+(1:ny_unitcell); % The range of indices in the y direction for ii-th meta-atom
    epsilon_metalens(ind_y,:) = epsilon_library(:,:,ind_meta_atom(ii));
end

% W_out > W, so we will pad extra pixels
ny_R_extra_half = round((W_out-W)/dx/2);

% Number of pixels in the y direction for the source (on the left) and the projection (on the right)
ny_L = ny_metalens;
ny_R = ny_metalens + 2*ny_R_extra_half;

W_in  = ny_L*dx; % Actual width of the incident wave [micron] 
W_out = ny_R*dx; % Actual width where we sample the transmitted field [micron] 

% Include homogeneous space and PML to the permittivity profile.
epsilon_L = n_silica^2; % Relative permittivity on the left
epsilon_R = n_air^2;    % Relative permittivity on the right and on the top & bottom
nPML = syst.PML.npixels;
% Extra number of pixels for four sides
nx_extra_left = 1 + nPML; % Add one pixel of free space for source and projection
nx_extra_right = nx_extra_left;
ny_extra_low = ny_R_extra_half + nPML;
ny_extra_high = ny_extra_low;

ny_tot = ny_metalens + ny_extra_low + ny_extra_high; % Total number of pixels in the y direction

% Construct the permittivity profile over all simulation domain (including
% metalens, PML, and extra padded pixels)
syst.epsilon = [epsilon_L*ones(ny_tot, nx_extra_left), ...
    [epsilon_R*ones(ny_extra_low,nx_metalens); epsilon_metalens; epsilon_R*ones(ny_extra_high,nx_metalens)], ...
    epsilon_R*ones(ny_tot, nx_extra_right)];

%% Build the compressed input-source matrix B
% Given the very large aspect ratio of the system, the input and output
% matrices B and C would have more nonzero elements than the
% Maxwell-operator matrix A. So, we compress matrices B and C.

time1 = clock;

% Obtain properties of propagating channels on the two sides.
BC = 'periodic'; % Periodic boundary condition means the propagating channels are plane waves
k0dx = 2*pi/wavelength*dx; % Dimensionless frequency k0*dx
use_continuous_dispersion = true; % Use continuous dispersion relation for (kx,ky), i.e. kx^2 + ky^2 = epsilon_bg*k0^2 (epsilon_bg = epsilon_L or epsilon_R) 
channels_L = mesti_build_channels(ny_L, 'TM', BC, k0dx, epsilon_L, [], use_continuous_dispersion);
channels_R = mesti_build_channels(ny_R, 'TM', BC, k0dx, epsilon_R, [], use_continuous_dispersion);

% We use all propagating plane-wave channels on the right, over width W_out.
M_R = channels_R.N_prop; % Number of channels on the right

% For the incident plane waves (over width W_in) in silica, we only take 
% the ones that can propagate in air. 
% Indices of those channels (among the total channels_L.N_prop channels)
ind_in_L = find(abs(channels_L.kydx_prop/dx) < (n_air*2*pi/wavelength));
M_L = numel(ind_in_L); % Number of channels we use on the left

% Build the input matrix on the left and right surfaces.
% Note that even though we only need the transmission matrix with input
% from the left, here we also include input channels from the right (same
% as the output channels of interest). This allows us to make matrix K =
% [A,B;C,0] symmetric. The computing time and memory usage of such partial
% factorization can be reduced when matrix K is symmetric.
if do_compression
    % Half the number of extra channels to pad
    % There is no point in having more channels than the number of spatial pixels
    M_L_pad_half = min([floor(ny_L-M_L)/2, round(M_L*pad_ratio_L/2)]);
    M_R_pad_half = min([floor(ny_R-M_R)/2, round(M_R*pad_ratio_R/2)]);

    % Total number of channels, including extra padded channels
    N_L = M_L + 2*M_L_pad_half;
    N_R = M_R + 2*M_R_pad_half;

    % Truncation window size
    ny_window_L = min([ny_L, 1 + round(trunc_width_over_lambda_L*wavelength/dx)]);
    ny_window_R = min([ny_R, 1 + round(trunc_width_over_lambda_R*wavelength/dx)]);
    
    % Build the compressed input matrix on the left and right surfaces.
    B_L = build_compressed_B(ny_L, N_L, ny_window_L, use_Hann_window);
    B_R = build_compressed_B(ny_R, N_R, ny_window_R, use_Hann_window);
else
    % Total number of channels
    N_L = M_L;
    N_R = M_R;

    % Build the input matrix on the left and right surfaces.
    B_L = channels_L.fun_u(channels_L.kydx_prop(ind_in_L));
    B_R = channels_R.fun_u(channels_R.kydx_prop);
end

% We take out the -2i prefactor of B, so C will equal transpose(B).
% The flux-normalization prefactor sqrt(nu) will be added later.
opts.prefactor = -2i;

% Preallocate the 2-element structure array B_struct
B_struct = struct('pos', {[],[]}, 'data', {[],[]});

% In mesti(), B_struct.pos = [m1, n1, h, w] specifies the position of a
% block source, where (m1, n1) is the index of the smaller-(y,x) corner,
% and (h, w) is the height and width of the block. Here, we put line
% sources (w=1) on the left surface (n1=n_L) and the right surface
% (n1=n_R) with height ny_L and ny_R centered around the metalens.
n_L = nx_extra_left;         % x pixel immediately before the metalens
n_R = n_L + nx_metalens + 1; % x pixel immediately after the metalens
m1_L = ny_extra_low + 1;     % first y pixel of the metalens
m1_R = nPML + 1;             % first y pixel of the output projection window
B_struct(1).pos = [m1_L, n_L, ny_L, 1];
B_struct(2).pos = [m1_R, n_R, ny_R, 1];

% B_struct.data specifies the source profiles inside such block, with
% B_struct.data(:, a) being the a-th source profile.
B_struct(1).data = B_L;
B_struct(2).data = B_R;

clear epsilon_metalens B_L B_R; % These are no longer needed, so clear them to save memory.

% Time spent to build the B
time2 = clock; timing_build_B = etime(time2,time1);

%% Compute the scattering matrix
% Note that with compression, this step will take about 7 GiB of memory (in addition to the memory MATLAB itself uses).

C = 'transpose(B)'; % Specify C = transpose(B)
D = []; % We only need the transmission matrix, for which D=0
opts.clear_syst = true; % syst can be cleared in mesti()
opts.clear_BC = true;   % B can be cleared in mesti()

[S, info] = mesti(syst, B_struct, C, D, opts);

%% Decompression step for APF-c
% Here we undo the compression, as described in supplementary section 5 of
% the APF paper.

time1 = clock;

% The full scattering matrix is S = [r_L,t_R;t_L,r_R]
% Extract the transmission matrix from left to right, t_L
ind_L = 1:N_L;
ind_R = (N_L+1):(N_L+N_R);
t = S(ind_R,ind_L); 
clear S; % It is no longer needed, so clear it to save memory.

% decompress
if do_compression
    % Here we do:
    % (1) centered  fft along dimension 1, equivalent to multiplying F_R on the left, and
    % (2) centered ifft along dimension 2, equivalent to multiplying inv(F_L) on the right.
    t = fftshift((ifft(ifftshift(fftshift((fft(ifftshift(t,1),[],1)),1),2),[],2)),2)/sqrt(N_R/N_L);
    
    % Remove the extra channels we padded earlier.
    ind_L = M_L_pad_half + (1:M_L);
    ind_R = M_R_pad_half + (1:M_R);
    t = t(ind_R,ind_L);
    
    % Undo the diagonal scaling, per Eq (S37) of the APF paper
    if use_Hann_window
        a_L = (-round((M_L-1)/2):round((M_L-1)/2));   % row vector
        a_R = (-round((M_R-1)/2):round((M_R-1)/2)).'; % column vector
        q_inv_L = 2./(1+cos((2*pi/N_L)*a_L)); % row vector
        q_inv_R = 2./(1+cos((2*pi/N_R)*a_R)); % column vector
        t = q_inv_R.*t.*q_inv_L; % use implicit expansion
    end
else
    % Undo the output-channel permutation when we set C = transpose(B)
    % See Supplementary Section 3 of the APF paper
    t = flipud(t);
end

% Multiply the flux-normalized prefactor sqrt(nu)
% Note that M_R includes all propagating channels on the right, but M_L
% only includes propagating channels on the left that can propagate in air
sqrt_nu_L = channels_L.sqrt_nu_prop(ind_in_L); % row vector
sqrt_nu_R = channels_R.sqrt_nu_prop.';         % column vector
t = sqrt_nu_R.*t.*sqrt_nu_L; % use implicit expansion

% Time spent on post-processing
time2 = clock; timing_post_process = etime(time2,time1);

fprintf('          Total elapsed time: %7.3f secs\n', info.timing.total + timing_build_B + timing_post_process);

%% Angular spectrum propagation parameters
% Below, we use angular spectrum propagation (ASP) to obtain field profile
% in the free space after the metalens, as described in Supplementary
% Section 12 of the APF paper.

% System width used for ASP to remove periodic wrapping artifact.
W_ASP_min = 2*W; % Minimal ASP window [micron]

% dx = wavelength/40 is not necessary. Down-sample to a coarser resolution for ASP.
dy_ASP = 5*dx; % ASP grid size [micron]

% If dy_ASP/dx is not a integer, we round it to a integer.
if round(dy_ASP/dx) ~= dy_ASP/dx
    warning(['dy_ASP/dx = %f should be a positive integer to ensure ' ...
        'down-sampling is possible; rounding it to %d.'], ...
        dy_ASP/dx, max([1, round(dy_ASP/dx)]));
    dy_ASP = max([1, round(dy_ASP/dx)])*dx;
end
% If ny_R is an even number and dy_ASP/dx is also even number, 
% we reduce dy_ASP/dx to dy_ASP/dx-1, making it odd number.
if mod(ny_R, 2) == 0 && mod(dy_ASP/dx, 2) == 0
    warning(['When ny_R is an even number, dy_ASP/dx should be an odd ' ...
        'number so down-sampling can be symmetric; reducing it to %d.'], ...
        (dy_ASP/dx)-1);
    dy_ASP = ((dy_ASP/dx)-1)*dx;
end

% fft is more efficient when the length is a power of 2, so we make ny_ASP
% a power of 2.
ny_ASP = 2^nextpow2(round(W_ASP_min/dy_ASP)); % Number of pixels for ASP
W_ASP = ny_ASP*dy_ASP; % Actual ASP window [micron]

% y index of the points we down-sample for ASP.
ind_ASP = 1:(dy_ASP/dx):ny_R;

% Make the sampling points symmetric around the middle.
if ind_ASP(end) ~= ny_R
    % Make sure that ny_R - ind_ASP(end) is an even number.
    if mod(ny_R - ind_ASP(end), 2) ~= 0
        ind_ASP = ind_ASP(1:(end-1));
    end
    ind_ASP = ind_ASP + (ny_R - ind_ASP(end))/2;
end

ny_ASP_pad = ny_ASP - numel(ind_ASP); % Total number of zeros to pad
ny_ASP_pad_low = round(ny_ASP_pad/2); % Number of zeros to pad on the low side
ny_ASP_pad_high = ny_ASP_pad - ny_ASP_pad_low; % Number of zeros to pad on the high side

% y position of the ASP points, including the padded zeros [micron] 
y_ASP = ((0.5:ny_ASP) - 0.5*(ny_ASP + ny_ASP_pad_low - ny_ASP_pad_high))*dy_ASP;

ny_ASP_half = round(ny_ASP/2); % recall that ny_ASP is an even number

% List of (kx,ky) in ASP
% Note that after fft, ky_ASP will be (2*pi/W_ASP)*(0:(ny_ASP-1)), where the
% latter half needs to be unwrapped to negative ky. We can either use fftshift
% and iffshift while centering ky_ASP, or just unwrap ky_ASP itself; we do the
% latter here.
ky_ASP = (2*pi/W_ASP)*[0:(ny_ASP_half-1), -ny_ASP_half:-1].'; % [1/micron]
kx_ASP = sqrt((n_air*2*pi/wavelength)^2 - ky_ASP.^2); % [1/micron]

% We only use the propagating components in ASP
kx_ASP_prop = kx_ASP(abs(ky_ASP) < (n_air*2*pi/wavelength)); % must be a column vector per asp() syntax

% List of incident angles in air [degree]
% n_air*sin(theta_in) = n_silica*sin(theta_silica), 
% where theta_silica = atan(channels_L.kydx/channels_L.kxdx)
theta_in_list = asind(sin(atan(channels_L.kydx_prop(ind_in_L)./channels_L.kxdx_prop(ind_in_L)))*n_silica/n_air); 

% Later we will use ifft to reconstruct field profile immediately after the
% metalens from the transmission matrix, i.e. Eq (S7) of the APF paper,
% and this is the prefactor we need.
prefactor_ifft = sqrt(ny_R)*exp((-2i*pi/ny_R*(M_R-1)/2)*(0:(ny_R-1)).');

%% Intensity profiles of the transmitted light
% Here, we use the transmission matrix and ASP to generate the transmitted
% intensity profile for different incident angles.
% Note that storing these intensity profiles takes 2.7 GiB of memory.

% List of incident angles in air at which to compute the profiles [degree]
%theta_in_list_profiles = -90:5:90;
% We sample angles near normal incidence with a finer spacing, since most
% features are near normal incidence.
theta_in_list_half = 90*linspace(0, 1, 40).^2;
theta_in_list_profiles = [-fliplr(theta_in_list_half(2:end)), theta_in_list_half];

W_plot = W;       % Plotting range in y [micron]
x_plot_max = 600; % Plotting range in x [micron]
dx_plot = 2.0;    % Resolution in x for plotting [micron]

% (x,y) coordinates used for plotting
ind_plot = find(abs(y_ASP) < W_plot/2);
ind_plot = [ind_plot(1)-1, ind_plot, ind_plot(end)+1]; % add one pixel on each side to ensure we fully cover W_plot
y_plot = y_ASP(ind_plot);
x_plot = 0:dx_plot:x_plot_max; % must be a row vector per asp() syntax

% Indices of the incident angles at which to compute the profiles.
% Pad large angles on two sides for interp1.
a_list = interp1([-200, theta_in_list, 200], [0, 1:M_L, inf], ...
    theta_in_list_profiles, 'nearest'); 

% Loop over the incident angles
n_angles_profiles = numel(theta_in_list_profiles);
intensity_profiles = zeros(numel(y_plot), numel(x_plot), n_angles_profiles);
for ii = 1:n_angles_profiles
    % Reconstruct field profile immediately after the metalens, restricted
    % to the propagating components, per Eq. (S7) of the APF paper:
    %    Ez^(a)(x=L,y) = sum_b t_ba*u_b(y)/sqrt(kx_b)
    % where u_b(y) = 1/sqrt(ny_R)*exp(i*ky_b*(y-y_0)), as in Eq. (S21) of the APF paper
    % The summation over plane waves b can be evaluated with ifft as follows:
    Ez0 = circshift(prefactor_ifft.*ifft((1./sqrt_nu_R).*t(:,a_list(ii)), ny_R), -1);

    % Down-sample the field profile prior to ASP.
    % We could have done so directly inside the ifft above to save time,
    % but it's more straightforward this way.
    Ez0_ASP = Ez0(ind_ASP,:); 

    % Obtain Ez(x,y) using ASP.
    Ez_ASP = asp(Ez0_ASP, x_plot, kx_ASP_prop, ny_ASP, ny_ASP_pad_low);

    % Only keep the part within the plotting range
    intensity_profiles(:,:,ii) = abs(Ez_ASP(ind_plot, :)).^2;
end

%% Angle dependence of the transmission efficiency and Strehl ratio
% Here, we map out the transmission efficiency and Strehl ratio for all
% incident angles.

% We define the Strehl ratio by normalizing against the focal intensity of an ideal lens at normal incidence.
% Here, we build Ez(x=L,y) for an ideal lens at normal incidence.
% The overall prefactor does not matter because we will normalize by the transmission.
y_metalens = ((0.5:ny_metalens) - (ny_metalens/2))*dx;
ideal_phase = -(2*pi/wavelength)*(sqrt(focal_length^2 + y_metalens.^2)).';
Ez_ideal = [zeros(ny_R_extra_half,1); exp(1i*ideal_phase); zeros(ny_R_extra_half,1)]; % Outside the metalens, the Ez_ideal is zero.

% Project Ez_ideal onto the propagating channels on the right
temp  = circshift(exp(-1i*2*pi/ny_R*(0:ny_R-1).').*fft(Ez_ideal), floor(channels_R.N_prop/2))/sqrt(ny_R);
t_ideal = (channels_R.sqrt_nu_prop).'.*temp(1:channels_R.N_prop);

% Use ASP to propagate Ez_ideal to the focal plane
Ez0 = circshift(prefactor_ifft.*ifft((1./sqrt_nu_R).*t_ideal,ny_R), -1);
Ez0_ASP = Ez0(ind_ASP,:); 
Ez_focal_ideal = asp(Ez0_ASP, focal_length, kx_ASP_prop, ny_ASP, ny_ASP_pad_low);

% ideal focal intensity normalized by transmission
I0_ideal = max(abs(Ez_focal_ideal))^2 / sum(abs(t_ideal).^2);

% Loop over the incident angles, in batches
SR_list = zeros(1, M_L); % Strehl ratio
T_list  = zeros(1, M_L); % transmission efficiency
batch_size = 100;
a_end = 0;
for ii = 1:ceil(M_L/batch_size)
    a_start = a_end + 1;
    a_end = min([a_start + batch_size - 1, M_L]);
    a = a_start:a_end;

    % transmission efficiency
    T_list(a) = sum(abs(t(:,a)).^2, 1); 
 
    % Use ASP to propagate Ez to the focal plane, same as the intensity profile
    Ez0 = circshift(prefactor_ifft.*ifft((1./sqrt_nu_R).*t(:,a), ny_R), -1);
    Ez0_ASP = Ez0(ind_ASP,:); 
    Ez_focal = asp(Ez0_ASP, focal_length, kx_ASP_prop, ny_ASP, ny_ASP_pad_low);

    % At large angles, the focal intensity may become even lower than those on the wrong side, so we restrict to the correct side.
    Ez_focal_half = zeros(ny_ASP_half, numel(a));
    Ez_focal_half(:, a<round(M_L/2)) = Ez_focal(1:ny_ASP_half,a<round(M_L/2));
    Ez_focal_half(:,a>=round(M_L/2)) = Ez_focal((ny_ASP_half+1):ny_ASP,a>=round(M_L/2));

    % Strehl ratio
    SR_list(a) = (max(abs(Ez_focal_half),[],1).^2./T_list(a)) / I0_ideal;
end

%% Animate the results

% We apply the same normalization to all incident angles. Note that all of
% them already have the same incident flux. We significantly saturate the
% plot near normal incidence, so the intensity at oblique incidence can be
% seen.
norm_factor = 100/max(intensity_profiles, [], 'all');

% Loop through Gaussian beams focused at different locations.
colormap('hot');
deg = char(176);
figure
for ii = 1:n_angles_profiles
    a_ii = a_list(ii);
    theta_in = theta_in_list(a_ii);

    % Plot the intensity profile
    clf
    subplot(1,2,1);
    imagesc(x_plot, y_plot, norm_factor*intensity_profiles(:,:,ii));
    set(gca,'YDir','normal')
    axis image
    xlim([0 x_plot_max])
    ylim(W_plot*[-0.5, 0.5]);
    c = colorbar;
    c.Ticks = 0;
    caxis([0 1])
    xlabel('x (\mum)')
    ylabel('y (\mum)')
    xticks([0 300 600])
    yticks([-500 0 500])
    set(gca,'fontsize',16);
    %title(['\theta_{in} = ', sprintf('%4.1f', theta_in),  deg]);

    % Plot transmission efficiency and Strehl ratio
    subplot(1,2,2);
    yyaxis left
    semilogy(theta_in_list(1:M_L), SR_list, 'linewidth', 1.5)
    hold on
    plot(theta_in, SR_list(a_ii), 'o', 'MarkerSize', 12, 'LineWidth', 2.5)
    xlabel('Incident angle');
    ylabel('Strehl ratio')
    xlim([-90 90])
    ylim([5e-3 1])
    xticks(-90:45:90)
    xticklabels({['-90' deg], ['-45' deg], ['0' deg], ['45' deg], ['90' deg]})
    yyaxis right
    plot(theta_in_list(1:M_L), T_list, 'linewidth', 1.5)
    hold on
    plot(theta_in, T_list(a_ii), 'o', 'MarkerSize', 12, 'LineWidth', 2.5)
    ylabel('Transmission efficiency')
    ylim([0 1])
    yticks(0:0.2:1)
    set(gca,'fontsize',16);
    set(gca,'TickLength',[0.02 0.035])
    set(gca,'linewidth', 1)

    drawnow
end
