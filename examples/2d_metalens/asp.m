function f = asp(f0, x, kx_prop, ind_prop, ny_tot, ny_pad_low)
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
%       Longitudinal wave number kx. It can include both propagating and
%       evanescent ones (in which also all of them will be propagated), or
%       only propagating ones (in which case evanescent ones will be
%       ignored). length(kx_prop) must be an odd number or ny_tot.
%   ind_prop (integer vector with length N_prop; optional):
%       Indices of the N_prop propagating channels among all ny_tot
%       channels. When N_prop == ny, it doesn't need to be given.
%   ny_tot (integer scalar; optional):
%       Total number of grid points in y direction, including padded
%       zeros. Must be no smaller than ny and N_prop. Defaults to ny.
%   ny_pad_low (integer scalar; optional):
%       Number of zeros to pad on the low side. Defaults to
%       round((ny_tot-ny)/2).
%
%   === Output Arguments ===
%   f (numeric column vector or matrix):
%       Field profile(s) f(x,y). Different columns corresponds to
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

% ind_prop is optional when N_prop == ny
if nargin < 4 && N_prop ~= ny
    error('When length(kx_prop) ~= size(f0,1), ind_prop must be given.');
elseif ~(isvector(ind_prop) && numel(ind_prop) == N_prop)
    error('ind_prop, when given, must be an integer vector with length equal length(kx_prop) = %d.', N_prop);
end

% no zero-padding unless ny_tot is given
if nargin < 5
    ny_tot = ny;
elseif ny_tot < ny
    error('ny_tot, when given, must be no smaller than size(f0,1) = %d.', ny);
elseif ny_tot < N_prop
    error('ny_tot, when given, must be no smaller than length(kx_prop) = %d.', N_prop);
else

% pad zeros symmetrically by default
if nargin < 6
    ny_pad_low = round((ny_tot-ny)/2);
elseif ny_pad_low + ny > ny_tot
    error('ny_pad_low + ny must be no greater than size(f0,1) = %d.', ny);
end

% Fourier transform f(x=0,y) to f(x=0,ky), as in Eq. (S43) of the SCSA paper.
% To get a finer spacing in ky, zeros are padded below and above f0.
% The most straightforward implementation is:
%   f0_fft = fft([zeros(ny_pad_low,size(f0,2)); f0; zeros(ny_tot-ny-ny_pad_low,size(f0,2))]);
% But that comes with some redundant computations on the zeros.
% The following is equivalent but slightly more efficient:
f0_fft = exp((-2i*pi*ny_pad_low/ny_tot)*(0:(ny_tot-1)).').*fft(f0,ny_tot,1);

% Remove the evanescent components of f(x=0,ky), propagate it to f(x,ky),
% and ifft back to f(x,y), as in Eqs. (S41-S42) of the SCSA paper.
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
    f_fft_prop = exp(1i*kx_prop.*x).*f0_fft(ind_prop,:);
    f = exp((-2i*pi*a_max/ny_tot)*(0:(ny_tot-1)).').* ...
            ifft(circshift(f_fft_prop, a_max), ny_tot, 1);
end

end
