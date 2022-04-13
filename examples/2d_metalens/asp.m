function f = asp(f0, x, kx, ind_prop, ny_tot, ny_pad_low)
%ASP Use angular spectrum propagation (ASP) to propagate a scalar field
%   from f0 = f(x=0,y) to f(x,y).
%   === Input Arguments ===
%   f0 (numeric column vector or matrix):
%       Initial field profile at x = x0 plane 
%       The column index is index in transverse (y) direction. If f0
%       is a numeric matrix, then different columns of f0 correspond 
%       to different inputs for ASP.
%   x (numeric scalar or row vector):
%       Position field propagating to. If x is a column vector, 
%       different rows of x correspond to different positions
%       for ASP to propagate.
%   kx (ny_tot-by-1 complex column vector):
%       Longitudinal wave number kx for all channels, including both 
%       propagating and evanescent ones   
%   ind_prop (N_prop-by-1 integer column vector):
%       Indices of the N_prop propagating channels among all ny_tot channels
%   ny_tot (numeric scalar, real):
%       Number of grid points in the transverse (y) direction for the ASP
%   ny_pad_low (numeric scalar, real):
%       Number of zeros to pad on one low side.
%   === Output Arguments ===
%   f (numeric column vector or matrix):
%       Field profile at x0+x plane 
%       The column index is index in transverse (y) direction. If f
%       is a numeric matrix, then different columns of f correspond
%       to different outputs for ASP.

if ~isequal(size(kx), [ny_tot,1])
    error('kx must be a column vector with ny_tot = %d elements.', ny_tot);
end

if size(x,2) ~= 1 && size(f0,2) ~= 1
    error(['Either multiple propagation distances or multiple initial fields ' ...
        'can be supported, but cannot input multiple propagation distances ' ...
        'and multiple initial fields']);    
end

% Fourier transform f(x=0,y) to f(x=0,ky), as in Eq. (S43) of the SCSA paper.
% To get a finer spacing in ky, zeros are padded below and above f0.
% The most straightforward implementation is:
% f0_fft = fft([zeros(ny_pad_low,size(f0,2)); f0; zeros(ny_tot-size(f0,1)-ny_pad_low,size(f0,2))]);
% But that comes with some redundant computations on the zeros.
% The following is equivalent but slightly more efficient:
f0_fft = exp((-2i*pi*ny_pad_low/ny_tot)*(0:(ny_tot-1)).').*fft(f0,ny_tot,1);

% Remove the evanescent components of f(x=0,ky), propagate it to f(x,ky),
% and ifft back to f(x,y), as in Eqs. (S41-S42) of the SCSA paper.
% The most straightforward implementation is:
% f = ifft(exp(1i*kx(ind_prop)*x).*f0_fft(ind_prop,:));
% This can be evaluated efficiently as following:
f = exp((-2i*pi/ny_tot*floor(size(ind_prop,1)/2))*(0:(ny_tot-1)).').* ...
        ifft( ...
        circshift(exp(1i*kx(ind_prop)*x).*f0_fft(ind_prop,:), floor(size(ind_prop,1)/2)), ...
        ny_tot,1);
end
