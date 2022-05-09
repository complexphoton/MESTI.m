function epsilon = build_epsilon_dbr(dx, n_bg, n_1, n_2, d_1, d_2, N_pair)
%BUILD_epsilon generate the permittivity profile of a distributed Bragg reflector (DBR) with subpixel smoothing.
%
%   === Input Arguments ===
%   dx (numeric scalar, real):
%       Grid size
%   n_bg (numeric scalar, real or complex):
%       Refractive index of background material
%   n_1 (numeric scalar, real or complex):
%       Refractive index of material 1
%   n_2 (numeric scalar, real or complex):
%       Refractive index of material 2
%   d_1 (numeric scalar, real):
%       Thickness of the material 1 in a pair
%   d_2 (numeric scalar, real):
%       Thickness of the material 2 in a pair
%   N_pair (numeric scalar, integer):
%       Number of pair in the distributed Bragg reflector 
%   === Output Arguments ===
%   epsilon (numeric matrix, real or complex):
%       Discretized relative permittivity profile of the distributed Bragg reflector.

epsilon_1 = n_1^2;                  % relative permittivity of of material 1
epsilon_2 = n_2^2;                  % relative permittivity of of material 2
epsilon_bg = n_bg^2;                % relative permittivity of of background material 
a = d_2+d_1;                        % thickness of a pair of DBR 
total_thickness = N_pair*(d_1+d_2); % total thickness DBR 

% Number of grid points
nx = ceil(total_thickness/dx); % number of grid points of Ez on integer sites in x

% Initialize permittivity profile by assuming every pixel is material 1.
epsilon = epsilon_1*ones(1,nx); 

% Loop over integer pixel in x to determine the permittivity profile
for ii = 1:nx
    if mod(ii*dx, a) > d_1
        if mod((ii-1)*dx, a) < d_1 % when crossing from material 1 to material 2
            order_of_pair = floor(ii*dx/a)+1; 
            pixel_n_2_ratio = (ii*dx - ((order_of_pair-1)*a+d_1))/dx;
            pixel_n_1_ratio = (1 - pixel_n_2_ratio);
            epsilon(:,ii) = (epsilon_2*pixel_n_2_ratio + epsilon_1*pixel_n_1_ratio); % average permittivity
        else
            epsilon(:,ii) = epsilon_2; % whole pixel in material 2
        end
    elseif mod(ii*dx, a) < d_1 && mod((ii-1)*dx, a) > d_1 % when crossing from material 2 to material 1
            order_of_pair = floor(ii*dx/a)+1; 
            pixel_n_1_ratio = (ii*dx - ((order_of_pair-1)*a))/dx;
            pixel_n_2_ratio = (1 - pixel_n_1_ratio);
            epsilon(:,ii) = (epsilon_2*pixel_n_2_ratio + epsilon_1*pixel_n_1_ratio); % average permittivity
    end
end

% Take care of the last pixel.
if nx*dx > total_thickness % crossing from material 2 to background material
    last_pixel_bg_ratio_integer_grid = nx - total_thickness/dx;
    last_pixel_n_2_ratio_integer_grid = 1-last_pixel_bg_ratio_integer_grid;    
    epsilon(:,nx) = (epsilon_bg*last_pixel_bg_ratio_integer_grid + epsilon_2*last_pixel_n_2_ratio_integer_grid); % average permittivity
end