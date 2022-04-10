# Meta-atom design for metasurfaces


Example of Meta-atom design for metasurfaces using mesti2s()


Use mesti2s() to 


1. Computing the transmission coefficient of meta-atom with different ridge widths and find meta-atoms satisfying 8 discrete ideal relative phase over [0, 2pi). 


2. Scanning over ridge width and incident angle to get phase and amplitude map of transmission coefficient.


# System parameters


Set up the parameters for the meta-atom system.

```matlab
clear
    
n_air    = 1;    % Refractive index of air
n_silica = 1.46; % Refractive index of silica
n_TiO2   = 2.43; % Refractive index of TiO2
lambda   = 532;  % Free-space wavelength [nm]
dx = lambda/40;  % Discretization grid size [nm]
w  = 18*dx;      % Width of meta-atom cell [nm]
l  = 600;        % Thickness of meta-atom cell [nm]
ridge_height = l; % Ridge height is the thickness of meta-atom cell.
```

# General setup for mesti2s()


Set up input arguments for mesti2s(). 

```matlab
syst.epsilon_L = n_silica^2; % Relative permittivity on the left hand side
syst.epsilon_R = n_air^2;    % Relative permittivity on the right hand side
syst.wavelength = lambda;    % Free-space wavelength [nm]
syst.dx = dx;                % Grid size of system [nm]
syst.length_unit = 'nm';     % Length unit
in = {'left'};               % Specify input channel on the left.
out = {'right'};             % Specify output channel on the right.
opts.verbal = false;         % Suppress output information.
```

# Structure of meta-atom


Plot refractive index profile of a meta-atom with ridge with = 79.8 nm as an illustration.

```matlab
ridge_width = 79.8; % Ridge width of meta-atom [nm]

% Build permittivity for the meta-atom. 
% Please refer to the function build_epsilon_meta_atom.    
epsilon_meta_atom = build_epsilon_meta_atom(dx, n_air, n_TiO2, ridge_width, ridge_height, w);
[ny, nx]= size(epsilon_meta_atom);

n_extra_for_plot = 10; % Extra pixels on side for plotting

% For plotting the space position
y = [0 ny*dx];
x = [-n_extra_for_plot*dx (nx+n_extra_for_plot)*dx]; 

clf
imagesc(x, y, [syst.epsilon_L*ones(ny,n_extra_for_plot), epsilon_meta_atom, 1*syst.epsilon_R*ones(ny,n_extra_for_plot)])
colormap(flipud(pink));
xlabel('{\itx} (nm)');
ylabel('{\ity} (nm)');
set(gca, 'fontsize', 15, 'FontName','Arial')
caxis([1 12])
xlim([x(1), x(2)])
ylim([y(1), y(2)])
text(670,130,'air','FontSize',20,'Rotation',90)
text(240,120,'TiO_2','FontSize',20)
text(-80,140,'silica','FontSize',20,'Rotation',90)
title('Meta-atom','FontSize',20)
```

<img src="meta_atom_design_structure.png" width="672" height="504"> 

# Transmission coefficient of meta-atom at normal incidence


In standard procedure of designing metasurface, people calculate the phase of transmission coefficient of meta-atom with different parameters (for example, ridge width). Then, use this information to design metasurface. Here following the similar procedure, we calculate transmission coefficient of meta-atom looping over different ridge widths.

```matlab
syst.yBC = 'periodic'; % Periodic boundary in y direction
% In periodic boundary, since 2*w/(lambda/n_silica) ~ 1, only one
% propagation channel is on the left. Our ky(a) = a*2*pi/W in periodic boundary 
% and ky = 0 is the only propagation channel whose incident angle is normal. 

ridge_width_list = 40:0.1:200; % List of ridge width: from 40 nm to 200 nm with 0.1 nm increment

t_list = zeros(1,size(ridge_width_list,2)); % Transmission coefficient list

% Loop over different ridge widths
for ii =1:length(ridge_width_list)
    syst.epsilon = build_epsilon_meta_atom(dx, n_air, n_TiO2, ridge_width_list(ii), ridge_height, w);
    % Compute the transmission matrix, which only contains one coefficient (no diffraction) at normal incidence.
    t_list(1,ii) = mesti2s(syst, in, out, opts);
end

rel_phi_over_pi_list = mod(angle(t_list)-angle(t_list(1)), 2*pi)/pi; % Relative phase over different ridge widths

% Plot the relative phase of meta-atom with different ridge widths
clf
plot(ridge_width_list, rel_phi_over_pi_list, '-','linewidth', 2)
xlabel('Ridge width (nm)')
ylabel('Phase (\pi)')
xlim([40 200])
title('$\Phi - \Phi^0$', 'Interpreter','latex')
set(gca, 'fontsize', 20, 'FontName','Arial')
set(gca,'linewidth', 2)
```

<img src="meta_atom_design_relative_phase.png" width="672" height="504"> 

# Finding meta-atoms satisfying 8 discrete ideal relative phases over [0, 2pi)


In the design point of view, practically people would only use meta-atom with discrete parameters to construct the metasurface. Typically, 8 different meta-atoms covering the relative phases equally spacing over [0, 2pi) are used.



```matlab
n_phases = 8;
ideal_rel_phase_over_pi_list = 2*linspace(0, 1-1/n_phases, n_phases); % Have 8 discrete equally spacing ideal relative phases over [0, 2pi).

% Find meta-atoms which are closest to the ideal relative phase through nearest neighbor interpolation.
ind = interp1(rel_phi_over_pi_list,1:length(rel_phi_over_pi_list),ideal_rel_phase_over_pi_list,'nearest');
phi_over_pi_design_list = rel_phi_over_pi_list(ind); 
ridge_width_design_list = ridge_width_list(ind); 

% Print the relative phases and ridge widths.
fprintf(['Relative phase.(pi) %5.2f %5.2f %5.2f %5.2f %5.2f %5.2f %5.2f %5.2f\n' ...
         'Ridge width....(nm) %5.1f %5.1f %5.1f %5.1f %5.1f %5.1f %5.1f %5.1f\n'],...
         phi_over_pi_design_list,ridge_width_design_list);
```


```text:Output
Relative phase.(pi)  0.00  0.25  0.50  0.75  1.00  1.25  1.50  1.75
Ridge width....(nm)  40.0  49.1  60.7  73.1  87.4 107.1 138.2 172.3
```


```matlab
% Save the phase list and the ridge width list.
save('meta_atom.mat','ridge_width_design_list','phi_over_pi_design_list')
```

# Phase and amplitude map of transmission coefficient of meta-atom with different ridge widths and incident angles


In addition to normal incidence, the response of oblique incidence is also important. With Bloch periodic boundary condition, the transmission coefficient with different ridge widths and incident angles is calculated.

```matlab
syst.yBC = 'Bloch'; % Bloch periodic boundary along transverse direction

ridge_width_list = 40:4:200; % List of ridge width: from 40 nm to 200 nm with 4 nm increment
theta_in_list = -89:1:89;    % List of incident angle in air [degree]
n_angles = numel(theta_in_list);

ky_list = (2*pi/lambda)*sind(theta_in_list);  % wave number in y direction

t_list = zeros(n_angles, numel(ridge_width_list));  % Transmission coefficient list 

% Loop over different ridge widths
for ii = 1:length(ridge_width_list)
    syst.epsilon = build_epsilon_meta_atom(dx, n_air, n_TiO2, ridge_width_list(ii), ridge_height, w);

    % Loop over different incident angles
    for jj = 1:n_angles
        syst.ky_B = ky_list(jj); % Bloch wave number in y direction
        [tmatrix, channels] = mesti2s(syst, in, out, opts); % compute the transmission matrix.

        % At large incident angles, there can be more than one channel on the left (i.e., diffraction).
        % We want the incident channel whose ky equals the Bloch wave number we specify (i.e., zeroth-order).
        [~, ind] = min(abs(channels.L.kydx_prop - syst.ky_B*dx));
        t_list(jj,ii) = tmatrix(1, ind);
    end
end

% Plot the phase of the transmission coefficient relative to the first width
clf
imagesc(ridge_width_list,theta_in_list, mod(angle(t_list)-angle(t_list(:,1)), 2*pi))
caxis([0, 2*pi]);
xlabel('Pillar width (nm)')
ylabel('\theta_{in} (degree)')
title('Phase')
colormap(twilight)
colorbar
hcb=colorbar; hcb.Ticks = [0 pi 2*pi]; hcb.TickLabels = {'0','\pi','2\pi'};
```

<img src="meta_atom_design_phase_map.png" width="672" height="504"> 


```matlab
% Plot the amplitude of transmission coefficient.
figure
imagesc(ridge_width_list,theta_in_list, abs(t_list))
caxis([0, 1]);
xlabel('Pillar width (nm)')
ylabel('\theta_{in} (degree)')
title('Amplitude')
colormap('hot')
colorbar
```


<img src="meta_atom_design_amplitude_map.png" width="672" height="504">

