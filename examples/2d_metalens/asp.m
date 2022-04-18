function f = asp(f0, x, kx_prop, ny_tot, ny_pad_low)
%ASP Use angular spectrum propagation (ASP) to propagate a scalar field
%   from f0 = f(x=0,y) to f(x,y).
%   === Input Arguments ===
%   f0 (numeric column vector or matrix):
%       Initial field profile at x = 0 plane, as a column vector.
%       If f0 has more than one column, each column is treated as a
%       distinct initial field profile. size(f0,1) = ny.
%   x (numeric scalar or row vector):
%       Distance to propagate to. When multiple distances are of interest,
%       they can be given as a row vector; the profiles f(x,y) at these
%       distances will be returned as the different columns of f. When f0
%       has more than one column, x must be a scalar (i.e., cannot
%       propagate to multiple distances when multiple initial fields are
%       given).
%   kx_prop (N_prop-by-1 column vector):
%       Longitudinal wave number kx. It can either (1) include all of the
%       ny_tot wave numbers (both propagating and evanescent ones) in which
%       case N_prop = length(kx_prop) must equal ny_tot and all of them
%       will be propagated, or (2) include a subset of the wave numbers up
%       to some cutoff in ky (e.g., the propagating ones) in which case
%       N_prop = length(kx_prop) must be an odd number and only these
%       components will be propagated.
%   ny_tot (integer scalar; optional):
%       Total number of grid points in y direction, including padded
%       zeros. Must be no smaller than ny and N_prop. Defaults to ny.
%   ny_pad_low (integer scalar; optional):
%       Number of zeros to pad on the low side. Defaults to
%       round((ny_tot-ny)/2).
%
%   === Output Arguments ===
%   f (numeric column vector or matrix):
%       Field profile(s) f(x,y). Different columns correspond to
%       different initial profiles (if f0 has more than one column) or
%       different distance x (if x is a row vector).

if nargin < 3
    error('Not enough input arguments.');
end

if ~ismatrix(f0); error('f0 must be a column vector or a matrix.'); end
ny = size(f0,1);

if ~isrow(x); error('x must be a scalar or a row vector.'); end

if size(x,2) ~= 1 && size(f0,2) ~= 1
    error('x must be a scalar when f0 has more than one column.');
end

if ~iscolumn(kx_prop); error('kx_prop must be a column vector.'); end
N_prop = numel(kx_prop);

% no zero-padding unless ny_tot is given
if nargin < 4
    ny_tot = ny;
elseif ny_tot < ny
    error('ny_tot, when given, must be no smaller than size(f0,1) = %d.', ny);
elseif ny_tot < N_prop
    error('ny_tot, when given, must be no smaller than length(kx_prop) = %d.', N_prop);
else

% pad zeros symmetrically by default
if nargin < 5
    ny_pad_low = round((ny_tot-ny)/2);
elseif ny_pad_low + ny > ny_tot
    error('ny_pad_low + ny must be no greater than size(f0,1) = %d.', ny);
end

% Fourier transform f(x=0,y) to f(x=0,ky), as in Eq. (S43) of the APF paper.
% To get a finer spacing in ky, zeros are padded below and above f0.
% The most straightforward implementation is:
%   f0_fft = fft([zeros(ny_pad_low,size(f0,2)); f0; zeros(ny_tot-ny-ny_pad_low,size(f0,2))]);
% But that comes with some redundant computations on the zeros.
% The following is equivalent but slightly more efficient:
f0_fft = exp((-2i*pi*ny_pad_low/ny_tot)*(0:(ny_tot-1)).').*fft(f0,ny_tot,1);

% Remove the evanescent components of f(x=0,ky), propagate it to f(x,ky),
% and ifft back to f(x,y), as in Eqs. (S41-S42) of the APF paper.
% The most straightforward implementation is:
%   f_fft_prop = zeros(size(f0_fft));
%   f_fft_prop(ind_prop,:) = exp(1i*kx_prop.*x).*f0_fft(ind_prop,:);
%   f = ifft(f_fft_prop);
% But that comes with some redundant computations on the zeros.
% The following is equivalent but slightly more efficient:
if N_prop == ny_tot
    f = ifft(exp(1i*kx_prop.*x).*f0_fft);
else
    if mod(N_prop,2) ~= 1
        error('length(kx_prop) = %d must be an odd number when it is not ny_tot.', N_prop);
    end
    a_max = round((N_prop-1)/2);
    ind_prop = [1:(a_max+1), (ny_tot-a_max+1):ny_tot];
    f_fft_prop = exp(1i*kx_prop.*x).*f0_fft(ind_prop,:);
    f = exp((-2i*pi*a_max/ny_tot)*(0:(ny_tot-1)).').* ...
            ifft(circshift(f_fft_prop, a_max), ny_tot, 1);
end

end
