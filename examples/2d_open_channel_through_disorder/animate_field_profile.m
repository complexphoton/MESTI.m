function animate_field_profile(field_profile, x0_list, y0_list, ...
    r0_list, x, y, nperiod, nframes_per_period)

theta = linspace(0, 2*pi, 100);
circ_x = cos(theta);
circ_y = sin(theta);

colormap(colorcet_D09); % use a blue-white-red colormap from ColorCET
set(gcf,'color','w');
for ii = 1:nperiod
    for jj = 0:(nframes_per_period-1)
        clf
        imagesc(x, y, ...
            real(field_profile*exp(-1i*2*pi*jj/nframes_per_period)));
        caxis([-1 1]);
        set(gca,'YDir','normal')
        hold on
        for kk = 1:numel(r0_list)
            plot(x0_list(kk) + circ_x*r0_list(kk), ...
                 y0_list(kk) + circ_y*r0_list(kk), ...
                 'k-', 'linewidth', 0.8);
        end
        axis off image
        drawnow
    end
end

end
