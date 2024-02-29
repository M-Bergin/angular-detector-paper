function bentley_plotScatter3D(all, numRaysScattered, sample_surf, N, deltax, hrms, lambda, init_theta, varargin)
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
    pos(3) = 900;
    pos(4) = 640;
    set(gcf, 'Position', pos);
    
    % All Scattered Rays Numerical PDF
    subplot(15,6,[6*(0:5)+1,6*(0:5)+2]);
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
    imagesc(c_x, c_y, Ns ./ real(cosd(asind(sqrt(CX.^2+CY.^2)))))
    xticks([])
    yticks([])
    set(gca, 'YDir', 'normal')
    title('Ray Tracing Total Radiance')
    hold on
    plot(cosd(linspace(0,360)), sind(linspace(0,360)), 'w', 'Linewidth', 2);
    %scatter(sind(20),0,'r', 'Linewidth',2);
    
    % Single Scattered Rays Numerical PDF
    %subplot(15,8,[8*(0:5)+3,8*(0:5)+4]);
    subplot(15,6,[6*(0:5)+1,6*(0:5)+2]+6*7);
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
    plot(cosd(linspace(0,360)), sind(linspace(0,360)), 'w', 'Linewidth', 2);
    %scatter(sind(20),0,'r', 'Linewidth',2);
    
    % Without Shadowing / Masking
    subplot(15,6,[6*(0:5)+5,6*(0:5)+6]);
    
    anal_size = 200;
    xs = linspace(-1,1,anal_size);
    ys = linspace(-1,1,anal_size);
    [X, Y] = meshgrid(xs,ys);
    
    THETADASH = asind(sqrt(X.^2+Y.^2)) .* ((X.^2 + Y.^2)<1) + ((X.^2 + Y.^2)>1) .* pi/2;
    PHIDASH = atan2d(Y,X);
    
    CHI = acosd((cosd(THETADASH) + cosd(init_theta))./sqrt(2*(1 + cosd(init_theta).*cosd(THETADASH) + cosd(PHIDASH).*sind(init_theta).*sind(THETADASH))));
    PSI = acosd( (sind(init_theta) + sind(THETADASH).*cosd(PHIDASH)) ./ sqrt(2 * (1 + cosd(init_theta).*cosd(THETADASH) + cosd(PHIDASH) .* sind(init_theta) .* sind(THETADASH)) - (cosd(init_theta) + cosd(THETADASH)) .^2 ) );
    
   
    p1 = lambda^2/hrms^2 * exp(- lambda^2 * tand(CHI).^2 / 2 / hrms^2) ...               % Surface Gradient Dist
        .* tand(CHI) ./ cosd(CHI).^2 ;% .* ...                                              % Transform to angles
        %(cosd(CHI).*cosd(init_theta) + sind(CHI).*sind(init_theta).*cosd(THETADASH)) ... % Projection
        %.* sind(THETADASH) ./ 2 ./ sqrt(2 * (1 + cosd(init_theta).*cosd(THETADASH) ...   % Outgoing angle transformation Jacobian
        %+ cosd(PHIDASH) .* sind(init_theta) .* sind(THETADASH)));
    
    p1 = p1 .* (sqrt(X.^2+Y.^2)<=1);
    
    %p1 = lambda^2/hrms^2 * exp(- lambda^2 * tand(CHI).^2 / 2 / hrms^2) ...
    %    .* tand(CHI) ./ cosd(CHI).^2;
    
    p1 = p1 .* (cosd(PSI) + 1/tand(init_theta)./tand(CHI) > 0);
    
    % Numerical integration by adding columns
    % TODO is there a better way to do this?
    
    %phidashs = linspace(-180, 180, anal_size);
    %thetadashs = linspace(0, 90, anal_size);
    %[PHIDASH2, THETADASH2] = meshgrid(phidashs, thetadashs);
    
    %norm = sum(p1 * (psis(2) - psis(1)) * (chis(2) - chis(1)), 'all');
    
    %p1 = p1/norm;
    
    biglambda = @(gam) hrms/lambda/sqrt(2*pi)./cotd(abs(gam)) .* exp(-lambda^2 * cotd(gam).^2 /2 /hrms^2) - 0.5.*erfc(lambda*cotd(abs(gam))/sqrt(2)/hrms) ;
    p2 = p1 ./ (1+biglambda(init_theta)) ./ (1+biglambda(THETADASH));
    
    imagesc(xs, ys, flip(p2 .* sind(THETADASH),2));
    xticks([])
    yticks([])
    set(gca, 'YDir', 'normal') 
    title('Our Model with Shad./Mask.')
    
    hold on
    plot(cosd(linspace(0,360)), sind(linspace(0,360)), 'w', 'Linewidth', 2);
    %scatter(sind(20),0,'r', 'Linewidth',2);
    
    subplot(15,6,[6*(0:5)+5,6*(0:5)+6]+6*7);
    imagesc(xs, ys, flip(p1 .* sind(THETADASH),2));
    set(gca, 'YDir', 'normal') 
    title('Our Model w/o Shad./Mask.')
        xticks([])
    yticks([])
    
    hold on
    plot(cosd(linspace(0,360)), sind(linspace(0,360)), 'w', 'Linewidth', 2);
    
    % Paper's results
    %subplot(15,8,[8*(9:14)+6, 8*(9:14)+7, 8*(9:14)+8]);
    %subplot(15,8,[8*(9:14)+7, 8*(9:14)+8]);
    %subplot(15,8,[8*(0:5)+5,8*(0:5)+6]);
    subplot(15,6,[6*(0:5)+3,6*(0:5)+4]);
    biglambda = @(gam) hrms/lambda/sqrt(2*pi)./cotd(abs(gam)) .* exp(-lambda^2 * cotd(gam).^2 /2 /hrms^2) - 0.5.*erfc(lambda*cotd(abs(gam))/sqrt(2)/hrms) ;
    
    anal_size = 200;
    xs = linspace(-1,1,anal_size);
    ys = linspace(-1,1,anal_size);
    [X, Y] = meshgrid(xs,ys);
    THETA = asind(sqrt(X.^2+Y.^2)) .* ((X.^2 + Y.^2)<1) + ((X.^2 + Y.^2)>1) .* pi/2;
    PHI = atan2d(Y,X);
    
    THETA_SPEC = acosd((cosd(init_theta)+cosd(THETA))./...
                sqrt((cosd(PHI).*sind(THETA)+sind(init_theta)).^2+sind(PHI).^2.*sind(THETA).^2+...
                (cosd(init_theta)+cosd(THETA)).^2));
    
    ALPHA = 4.41*abs(PHI) ./ (4.41*abs(PHI) + 1);
    
    P_ILL_VIS = 1./(1+ biglambda(THETA_SPEC).*((THETA_SPEC > init_theta).*ALPHA + (THETA_SPEC <= init_theta)) + biglambda(init_theta).*((THETA_SPEC > init_theta) + (THETA_SPEC <= init_theta).*ALPHA));
    
    LUM = P_ILL_VIS .* exp(-tand(THETA_SPEC).^2*lambda^2/2/hrms^2) ./ cosd(THETA) ./ cosd(THETA_SPEC).^4;
    INTENSE = LUM .* cosd(THETA);
    
    LUM = LUM .* ((X.^2 + Y.^2)<1);
    INTENSE = INTENSE .* ((X.^2 + Y.^2)<1);
    
    LUM_NO_SHAD = exp(-tand(THETA_SPEC).^2*lambda^2/2/hrms^2) ./ cosd(THETA) ./ cosd(THETA_SPEC).^4;
    INTENSE_NO_SHAD = LUM_NO_SHAD .* cosd(THETA);
    
    LUM_NO_SHAD = LUM_NO_SHAD .* ((X.^2 + Y.^2)<1);
    INTENSE_NO_SHAD = INTENSE_NO_SHAD .* ((X.^2 + Y.^2)<1);
    
    imagesc(xs, ys, flip(LUM,2));
    set(gca, 'YDir', 'normal')   
    title('Ginneken1998 Model Radiance')
    xticks([])
    yticks([])
    
    hold on
    plot(cosd(linspace(0,360)), sind(linspace(0,360)), 'w', 'Linewidth', 2);
    
    
    subplot(15,6,[6*(0:5)+3,6*(0:5)+4]+6*7);
    imagesc(xs, ys, flip(LUM,2));
    set(gca, 'YDir', 'normal')   
    title('Ginneken1998 w/o Shad./Mask.')
    xticks([])
    yticks([])
    
    hold on
    plot(cosd(linspace(0,360)), sind(linspace(0,360)), 'w', 'Linewidth', 2);
    
    % Numerical Percentage of Rays Scattered
    sing = numRaysScattered == 1;
    mult = numRaysScattered > 1;
    frac = sum(mult(:)) / (sum(sing(:)) + sum(mult(:)));
    frac2 = (1 - sum(p2,'all') / sum(p1,'all'));
    frac3 = (1 - sum(INTENSE,'all') / sum(INTENSE_NO_SHAD,'all'));
    xtix = frac*100;
    if frac > 0.03 && frac2 > 0.03
        xtix = [0, xtix];
    end
    if frac < 0.97 && frac2 < 0.97
        xtix = [xtix, 100];
    end
    if abs(frac-frac2)>0.03
        xtix = sort([xtix, frac2*100]);
    end
    if abs(frac3-frac2)>0.03 && abs(frac3-frac)>0.03 && frac3 < 0.97 && frac3 > 0.03
        xtix = sort([xtix, frac3*100]);
    end
    sl1 = subplot(15,6,(42:48)+6*7);
    sl1.Position(4) = 0;
    scatter(frac*100,0, 'filled', 'r')
    hold on
    scatter(frac3.*100, 0, 'filled', 'b')
    scatter(frac2.*100, 0, 'filled', 'g')
    xticks(xtix)
    xtickformat('%d%%')
    xlim([0,100])
    title('Fraction of Rays undergoing Multiple Scattering')
    lg = legend({'Ray Tracing', 'Ginneken1998', 'Our Model'}, 'location', 'eastoutside');
    set(lg, 'Units', 'pixels')
    set(gca, 'Units', 'pixels')
    pos = get(gca, 'Position');
    lg.Position(2) = lg.Position(2) + 23;
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