function epsilon_dbr = build_epsilon_dbr(dx, n_bg, n1, n2, d1, d2, n_pair)

%BUILD_EPSILON_DBR generate relative permittivity profile of 1D distributed Bragg reflector (DBR) with subpixel smoothing.
%
%   === Input Arguments ===
%   dx (numeric scalar, real):
%       Grid size               
%   n_bg (numeric scalar, real or complex):
%       Refractive index of background material
%   n1 (numeric scalar, real or complex):
%       Refractive index of material 1
%   n2 (numeric scalar, real or complex):
%       Refractive index of material 2
%   d1 (numeric scalar, real):
%       Thickness of the material 1 in a pair
%   d2 (numeric scalar, real):
%       Thickness of the material 2 in a pair
%   n_pair (numeric scalar, integer):
%       Number of pair in DBR
%   === Output Arguments ===
%   epsilon_dbr (numeric matrix, real or complex):
%       Discretized relative permittivity profile of DBR

a = d2+d1; % Thickness of a pair of DBR 
total_thickness = n_pair*(d1+d2); % Total thickness DBR 
nx = ceil(total_thickness/dx); % Number of pixels in x direction

% Setting the permittivity profile with subpixel smoothing
epsilon_dbr = n1^2*ones(1, nx); % Assume that every pixel is material 1 as starting permittivity profile.
% Looping over every pixel to determine the permittivity profile
for ii = 1:nx
    if mod(ii*dx, a) > d1
        if mod((ii-1)*dx, a) < d1 % When crossing from material 1 to material 2
            order_of_pair = floor(ii*dx/a)+1; 
            pixel_n2_ratio = (ii*dx - ((order_of_pair-1)*a+d1))/dx;
            pixel_n1_ratio = (1 - pixel_n2_ratio);
            epsilon_dbr(ii) = (n2^2*pixel_n2_ratio + n1^2*pixel_n1_ratio); % Average permittivity           
        else
            epsilon_dbr(ii) = n2^2; % Whole pixel in material 2
        end
    elseif mod(ii*dx, a) < d1 && mod((ii-1)*dx, a) > d1 % When crossing from material 2 to material 1
            order_of_pair = floor(ii*dx/a)+1; 
            pixel_n1_ratio = (ii*dx - ((order_of_pair-1)*a))/dx;
            pixel_n2_ratio = (1 - pixel_n1_ratio);
            epsilon_dbr(ii) = (n2^2*pixel_n2_ratio + n1^2*pixel_n1_ratio); % Average permittivity       
    end
end

if nx*dx > total_thickness % Crossing from material 2 to background material
    pixel_bg_ratio = nx - total_thickness/dx;
    pixel_n2_ratio = 1-pixel_bg_ratio;    
    epsilon_dbr(nx) = (n_bg^2*pixel_bg_ratio + n2^2*pixel_n2_ratio); % Average permittivity
end

end