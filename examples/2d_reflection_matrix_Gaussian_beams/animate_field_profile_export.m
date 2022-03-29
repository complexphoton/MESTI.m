fname = 'img/';

theta = linspace(0, 2*pi, 100);
circ_x = cos(theta);
circ_y = sin(theta);

cmap_bluered = colorcet('D09'); % use a blue-white-red colormap from colorcet

% Loop through Gaussian beams focused at different locations 
clf
for ii = 1:M_in
    % Plot the total field profile; exclude PML
    ax1 = subplot(1,2,1);
    imagesc(x, y, real(field_profiles((nPML+1):(ny-nPML), (nPML+1):(nx-nPML), ii)));
    set(gca,'YDir','normal')
    caxis([-1,1]) % center the colorbar axis so zero = white
    colormap(ax1, cmap_bluered);
    hold on
    plot(x_0+r_0*circ_x, y_0+r_0*circ_y, 'k-', 'linewidth', 0.8);
    axis image off
    
    % Plot the reflection matrix
    ax2 = subplot(1,2,2);
    imagesc(y_f, y_f, abs(r).^2)
    set(gca,'YDir','normal')
    colormap(ax2, flipud(gray))
    xlabel('Input position')
    ylabel('Output position')
    hold on
    % Show the current focal position
    plot(y_f(ii)*[1,1], get(gca,'YLim'), 'b-')
    axis image
    xticks([])
    yticks([])
    set(gca, 'FontSize', 18)
    drawnow

    fname_jpg = [fname, sprintf('%02d.jpg', ii)];
    fname_png = [fname, sprintf('%02d.png', ii)];
    %export_fig(fname_jpg);
    export_fig(fname_png);
end
