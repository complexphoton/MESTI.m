# Distributed Bragg reflector


In this example, we use mesti2s() to compute the reflectance from an alternating sequence of layers of two dielectric material, also called a distributed Bragg reflector (DBR).


For a 1D system at normal incidence, the TM and TE polarizations are equivalent, and we will explicit check this equivalence.

```matlab
% System parameters
n_bg = 1;                     % refractive index of background material (air)
n_1 = 1.5;                    % refractive index of material 1
n_2 = 3;                      % refractive index for material 2
lambda_mid_gap = 550;         % mid-gap wavelength for DBR [nm]
d_1 = lambda_mid_gap/(4*n_1); % thickness of material 1 in a pair [nm]
d_2 = lambda_mid_gap/(4*n_2); % thickness of material 2 in a pair [nm]
N_pair = 5;                   % number of pair in DBR

lambda_min = 300; % minimum vacuum wavelength [nm]
lambda_max = 800; % maximum vacuum wavelength [nm]
delta_lambda = 2.5; % increment of the wavelength [nm]

lambda_list = lambda_min:delta_lambda:lambda_max; % wavelength list
n_lambda = numel(lambda_list); % number of wavelengths
lambda_0 = lambda_list(round((n_lambda+1)/2)); % central wavelength [nm]

syst.epsilon_L = n_bg^2;     % relative permittivity on the left
syst.epsilon_R = n_bg^2;     % relative permittivity on the right
syst.length_unit = 'nm';
syst.dx = lambda_min/n_2/20; % grid size; 20 points per wavelength in highest refractive index material in system
syst.yBC = 'periodic';       % 1D system at normal incidence has periodic boundary in y

opts.verbal = false;         % suppress output information

% Plot refractive index profile of the DBR
dx = lambda_0/n_2/20;
% Build the relative permittivity profile with subpixel smoothing
epsilon_dbr = build_epsilon_dbr(dx, n_bg, n_1, n_2, d_1, d_2, N_pair);

nx = ceil(N_pair*(d_1+d_2)/dx); % number of pixels in in DBR
n_extra_for_plot = 10; % number of pixels on side for plotting
x = (-n_extra_for_plot+0.5:nx+n_extra_for_plot)*dx;

% Plotting
clf
imagesc(x, [], [syst.epsilon_L*ones(1,n_extra_for_plot), epsilon_dbr, 1*syst.epsilon_R*ones(1,n_extra_for_plot)])
colormap(flipud(pink));
xlabel('Position (nm)');
yticks([])
text(-50, 1, 'air', 'FontSize', 15, 'Rotation', 90)
text(320, 1.1, 'material 1', 'FontSize', 15, 'Rotation', 90)
text(390, 1.1, 'material 2', 'FontSize', 15, 'color', 'white', 'Rotation', 90)
text(730, 1, 'air', 'FontSize', 15, 'Rotation', 90)
title('Distributed Bragg reflector', 'FontSize', 15)
set(gca, 'fontsize', 15, 'FontName', 'Arial')
```


![distributed_bragg_reflector_structure.png](distributed_bragg_reflector_structure.png)


```matlab
% For TM polarization, we need to provide epsilon. 
syst_TM = syst; clear syst
syst.polarization = 'TM'; % use TM, which gives the Ez component
syst_TM.epsilon = build_epsilon_dbr(syst_TM.dx, n_bg, n_1, n_2, d_1, d_2, N_pair);

% For TE polarization, we need to provide inv_epsilon instead.
% For a 1D system at normal incidence, syst_TE.inv_epsilon{1} is not necessary because the y derivative vanishes.
% Here, syst_TE.inv_epsilon{2} is simply 1./syst_TM.epsilon because the interfaces' normal vectors are in x direction. In 2D and 3D, subpixel smoothing would require more care, for example see Opt. Lett. 31, 2972 (2006).
syst_TE = rmfield(syst_TM, 'epsilon');
syst_TE.polarization = 'TE';
syst_TE.inv_epsilon{2} = 1./syst_TM.epsilon;

r_list_TM = zeros(1,n_lambda);
r_list_TE = zeros(1,n_lambda);
t_list_TM = zeros(1,n_lambda);
t_list_TE = zeros(1,n_lambda);

% Loop over wavelengths
for ii = 1:n_lambda
    syst_TM.wavelength = lambda_list(ii);
    syst_TE.wavelength = lambda_list(ii);

    % Compute the scattering matrix with input from the left, output to both sides
    % This gives smatrix = [r; t]
    smatrix_TM = mesti2s(syst_TM, {'left'}, {'left', 'right'}, opts);
    smatrix_TE = mesti2s(syst_TE, {'left'}, {'left', 'right'}, opts);

    r_list_TM(ii) = smatrix_TM(1,1);
    r_list_TE(ii) = smatrix_TE(1,1);
    t_list_TM(ii) = smatrix_TM(2,1);
    t_list_TE(ii) = smatrix_TE(2,1);
end

% Check that TM and TE give the same results.
% Note that r_TM = -r_TE because the reflection coefficient is defined based on Ez in TM, Hz in TE.
fprintf(['max(|r_TM + r_TE|) = %6.3g\n', ...
         'max(|t_TM - t_TE|) = %6.3g\n'], ...
        max(abs(r_list_TM + r_list_TE)), ...
        max(abs(t_list_TM - t_list_TE)));
```


```
max(|r_TM + r_TE|) = 5.4e-13
max(|t_TM - t_TE|) = 5.92e-13
```


```matlab
% Check energy conservation
fprintf('max(|1 - T - R|) = %6.3g\n', max(abs(1-abs(r_list_TM).^2-abs(t_list_TM).^2)));
```


```
max(|1 - T - R|) = 3.89e-15
```


```matlab
% Analytic solution
[r_list_analytical, t_list_analytical] = dbr_analytical(n_bg, n_1, n_2, d_1, d_2, N_pair, lambda_list);

% Plotting
clf
plot(lambda_list, abs(r_list_TM).^2, 'o', 'linewidth', 1)
hold on
plot(lambda_list, abs(r_list_TE).^2, 'x', 'linewidth', 1)
plot(lambda_list, abs(r_list_analytical).^2, 'k-', 'linewidth', 1)
xlabel('Wavelength \lambda (nm)')
ylabel('Reflectance{\it R}');
xlim([lambda_min, lambda_max])
ylim([0, 1])
legend('Numeric; TM', 'Numeric; TE', 'Analytic', 'Location', 'south')
set(gca, 'fontsize', 15, 'FontName', 'Arial')
set(gca, 'linewidth', 1)
```


![distributed_bragg_reflector_spectrum.png](distributed_bragg_reflector_spectrum.png)

