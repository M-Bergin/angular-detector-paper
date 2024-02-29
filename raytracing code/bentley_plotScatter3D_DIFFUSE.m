function bentley_plotScatter3D_DIFFUSE(all, numRaysScattered, sample_surf, N, deltax, hrms, lambda, init_theta, varargin)
    verbose = false;
    distribution = 'specular';
    for i_=1:2:length(varargin)
        switch varargin{i_}
            case 'verbose'
                verbose = varargin{i_+1};
            case 'distribution'
                distribution = varargin{i_+1};
        end
    end
    
    if ~strcmp(distribution, 'specular') && ~strcmp(distribution, 'cosine')
        error('I dunno what a %s distribution is.', distribution);
    end
    
    figure(2);
    clf;
    
    set(gcf, 'Units', 'pixels');
    pos = get(gcf, 'Position');
    pos(3) = 626;
    pos(4) = 403;
    set(gcf, 'Position', pos);
    
    % All Scattered Rays Numerical PDF
    subplot(8,4,[4*(0:5)+1,4*(0:5)+2]);
    filt = numRaysScattered > 0;
    scatt = all(:,filt);
    xs = sind(scatt(2,:)) .* cosd(scatt(1,:));
    ys = sind(scatt(2,:)) .* sind(scatt(1,:));
    
    x_bins = linspace(-1,1,75);
    y_bins = linspace(-1,1,75);
    c_x = (x_bins(1:74) + x_bins(2:75)) / 2;
    c_y = (y_bins(1:74) + y_bins(2:75)) / 2;
    [CX, CY] = meshgrid(c_x, c_y);
    Ns = hist2d(xs, ys, x_bins, y_bins, 'pdf', 'tile');
    hold off
    plotme = Ns ./ real(cosd(asind(sqrt(CX.^2+CY.^2))));
    imagesc(c_x, c_y, Ns)
    xticks([])
    yticks([])
    caxis([0 max(max(Ns .* (CX.^2+CY.^2 < 0.95)))])
    set(gca, 'YDir', 'normal')
    title('Ray Tracing Total Radiance')
    hold on
    plot(cosd(linspace(0,360)), sind(linspace(0,360)), ':w', 'Linewidth', 2);
    %scatter(sind(20),0,'r', 'Linewidth',2);
    
    % Single Scattered Rays Numerical PDF
    %subplot(15,8,[8*(0:5)+3,8*(0:5)+4]);
    subplot(8,4,[4*(0:5)+3,4*(0:5)+4]);
    filt = numRaysScattered == 1;
    sing = all(:,filt);
    
    xs = sind(sing(2,:)) .* cosd(sing(1,:));
    ys = sind(sing(2,:)) .* sind(sing(1,:));
    x_bins = linspace(-1,1,75);
    y_bins = linspace(-1,1,75);
    
    hist2d(xs, ys, x_bins, y_bins, 'pdf', 'tile');
    xticks([])
    yticks([])
    set(gca, 'YDir', 'normal')
    title('Ray Tracing Single Scattered Radiance')
    hold on
    plot(cosd(linspace(0,360)), sind(linspace(0,360)), ':w', 'Linewidth', 2);
    %scatter(sind(20),0,'r', 'Linewidth',2);
    
    % Numerical Percentage of Rays Scattered
    sing = numRaysScattered == 1;
    mult = numRaysScattered > 1;
    frac = sum(mult(:)) / (sum(sing(:)) + sum(mult(:)));
    xtix = frac*100;
    if frac > 0.03
        xtix = [0, xtix];
    end
    if frac < 0.97
        xtix = [xtix, 100];
    end
    sl1 = subplot(8,4,29:32);
    sl1.Position(4) = 0;
    scatter(frac*100,0, 'filled', 'r')
    xticks(xtix)
    xtickformat('%d%%')
    xlim([0,100])
    title('Fraction of Rays undergoing Multiple Scattering')
    lg = legend({'Ray Tracing'}, 'location', 'eastoutside');
    set(lg, 'Units', 'pixels')
    set(gca, 'Units', 'pixels')
    pos = get(gca, 'Position');
    lg.Position(2) = lg.Position(2) + 10;
    pos(3) = pos(3) - 155;
    set(gca, 'Position', pos);
    
    
    if strcmp(distribution, 'specular')
        sg = 'Specular Scattering from Surface with ';
    elseif strcmp(distribution, 'cosine')
        sg = 'Diffuse Scattering from Surface with ';
    else
        error('We shouldn''t get here.');
    end
    if exist('hrms', 'var')
        sg = sprintf('%sh_{RMS}=%.3f, ', sg, hrms);
    end
    if exist('lambda', 'var')
        sg = sprintf('%s\\lambda=%.3f, ', sg, lambda);
    end
    if exist('init_theta', 'var')
        sg = sprintf('%s\\theta_0=%.1f, ', sg, init_theta);
    end
    
    if endsWith(sg, ', ')
        sg = sg(1:length(sg)-2);
    end
    if endsWith(sg, ' with ')
        sg = sg(1:length(sg)-6);
    end
    sgtitle(sg);

end