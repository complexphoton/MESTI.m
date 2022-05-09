function epsilon = build_epsilon_fp(dx, n_bg, n_slab, thickness)
%BUILD_EPSILON_FP generate the permittivity profile of a Fabry-Pérot etalon with subpixel smoothing.
%
%   === Input Arguments ===
%   dx (numeric scalar, real):
%       Grid size
%   n_bg (numeric scalar, real or complex):
%       Refractive index of background material
%   n_slab (numeric scalar, real or complex):
%       Refractive index of dielectric slab 
%   thickness (numeric scalar, real):
%       Thickness of the dielectric slab
%   === Output Arguments ===
%   epsilon (numeric matrix, real or complex):
%       Discretized relative permittivity profile of the dielectric slab.

% Number of grid points
nx = ceil(thickness/dx);

% Initialize permittivity profile by assuming every pixel is dielectric slab.
epsilon = (n_slab^2)*ones(1,nx); 

% Average permittivity for the last pixel
last_pixel_bg_ratio = ceil((thickness)/dx) - thickness/dx; % Ratio of last pixel in background
last_pixel_slab_ratio = 1-last_pixel_bg_ratio ; % Ratio of last pixel in dielectric slab
epsilon(1,nx) = ((n_slab^2)*last_pixel_slab_ratio+(n_bg^2)*last_pixel_bg_ratio); 
