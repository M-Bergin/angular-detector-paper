addpath('functions', 'random_surface_generation', genpath('ParforProgMon'))

hss = [linspace(0.1,1.5,100), linspace(1.5,3.0,20)];
%hss = [1];
%lambda = 1;
%lambdas = 1./linspace(0.0001,3,100);
%lambdas = [1];
lambda = 1;
%hs = 1;
dx = 1/100;

loops = 20;
%loops = 1;

%for i=1:length(hss)
    %hs = hss(i);
for i=1:length(hss)
    
    hs = hss(i);
    
    in_angle = 20;
    
    % Do the ray tracing:
    ths_all = {};
    ns_all = {};
    
    tic
    parfor asdf=1:loops
        [fs, xs] = rsgene1D(20000, 20000/100, hs * sqrt(tanh(dx/lambda)), lambda);
        [ths, ns] = random_scatter1D('hs', fs, 'xs', xs, 'init_angle', in_angle, 'n_rays', 10000, 'scattering', 'specular');
        ths_all{asdf} = ths;
        ns_all{asdf} = ns;
    end
    toc
    
    ths_arr = [];
    ns_arr = [];
    for asdf=1:loops
        ths_arr = [ths_arr, ths_all{asdf}];
        ns_arr = [ns_arr, ns_all{asdf}];
    end
    
    nbins = 20;
    edges = linspace(min(ths_arr),max(ths_arr),nbins+1);
    hc = zeros(nbins,1);
    hc_all = zeros(nbins,1);
    one_scatter = ths_arr(ns_arr==1);
    for j=1:nbins
        thisbin = one_scatter(logical( (one_scatter >= edges(j)) .* (one_scatter < edges(j+1)) ));
        thisbin = 1 ./ cosd(thisbin);
        hc(j) = sum(thisbin);
        
        thisbin_all = ths_arr(logical( (ths_arr >= edges(j)) .* (ths_arr < edges(j+1)) ));
        thisbin_all = 1 ./ cosd(thisbin_all);
        hc_all(j) = sum(thisbin_all);
    end
    hc = hc / (length(one_scatter) .* (edges(2) - edges(1)));
    hc = length(one_scatter) / length(ths_arr) * hc;
    
    hc_all = hc_all / (length(one_scatter) .* (edges(2) - edges(1)));
    
    
    % Step one: the PDF of being seen with flux conservation but without
    % self-shadowing.
    ths = linspace(-89.99+in_angle,in_angle+89.99,100)/2;
    
    p0 = pi / 180 * lambda / sqrt(8*pi) / hs * exp(- lambda^2 /2 / hs^2 * tand(ths).^2) ./ cosd(ths).^2;
    
    v = hs.^2/lambda.^2;
    p1 = exp(-(tand(ths)).^2/2/v) .* cosd(in_angle-ths) ./ cosd(ths).^2;
    p1 = p1 .* (p1 > 0);
    thsnorm = -89.9:0.1:89.9;
    normintegrand = exp(-(tand(thsnorm)).^2/2/v) .* cosd(in_angle-thsnorm) ./ cosd(thsnorm).^2;
    normintegrand = normintegrand .* (normintegrand > 0);
    n = trapz(thsnorm,normintegrand);
    p1 = p1 / n;
    
    % Step two: add shadowing and masking.
    %biglambda_orig = @(gam) tand(gam)/2/pi .* (hs/lambda * exp(-lambda^2 * cotd(gam).^2 /2 /hs^2) - sqrt(pi/2)*erfc(lambda*cotd(gam)/sqrt(2)/hs).*cotd(gam));
    biglambda = @(gam) hs/lambda/sqrt(2*pi)./cotd(gam) .* exp(-lambda^2 * cotd(gam).^2 /2 /hs^2) - 0.5.*erfc(lambda*cotd(gam)/sqrt(2)/hs) ;
    mask_or_shadow = @(th) (1 - 1./( (1+biglambda(abs(in_angle))) .* (1+biglambda(abs(2*th-in_angle))) )) .* (th<0) + (1- 1./(1+biglambda(max(abs(in_angle),abs(2*th-in_angle))))).*(th>=0);
    %mask_or_shadow = @(th) 1 - 1./( (1+biglambda(max(abs(in_angle), abs(2*th-in_angle)))) );
    %mask_or_shadow3 = @(th) 1 - 1./( (1+biglambda_orig(max(abs(in_angle), abs(2*th-in_angle)))) );
    %mask_or_shadow = @(th) 1 - 1./( (1+biglambda(abs(in_angle)) + biglambda(abs(2*th-in_angle))) );
    p2 = p1 .* (1-mask_or_shadow(ths));
    %p4 = p1 .* (1-mask_or_shadow3(ths));


    % Step three: convert to outgoing angle distribution
    alphas = -2*ths + in_angle;
    p3 = p2 /2; % * pi / 360;
    %p5 = p4/2;
    
    figure(1);
    clf;
    hold on
    histogram('BinCounts', hc_all, 'BinEdges', edges, 'FaceColor', [0.9290, 0.6940, 0.1250]*0.6 + [0.4,0.4,0.4], 'FaceAlpha', 1);
    histogram('BinCounts', hc, 'BinEdges', edges, 'FaceColor', [0, 0.4470, 0.7410]*0.6 + [0.4,0.4,0.4], 'FaceAlpha', 1);
    plot(alphas, p1 ./ cosd(alphas) / 2, 'Linewidth', 2, 'Color', [0.4660, 0.6740, 0.1880])
    plot(alphas, p3 ./ cosd(alphas), 'Linewidth', 2, 'color', [0.6350, 0.0780, 0.1840])
    xlim([-90, 90])
    ylim([0, 0.03])
    xlabel('Outgoing Angle')
    ylabel('Radiance / arb.')
    title(sprintf('Single Scattered Radiance with Specular Scattering\nIncoming Angle=%d°, h_{RMS}=%.3f, \\lambda=%g', in_angle, hs, lambda))
    legend({'Multiple Scattered', 'Single Scattered', 'w/o Shad./Mask.', 'with Shad./Mask.',}, 'location', 'northwest')
    xticks([-90, -45, 0, 45, 90])
    xticklabels({'-90°', '-45°', '0°', '45°', '90°'})
    yticks([])
    
    fig = gcf;
    set(fig,'Units','Inches');
    set(fig,'Position',[4,4,4.3,3]);
    pos = get(fig,'Position');
    set(fig,'Paperorientation', 'landscape', 'PaperType', 'a3', 'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
    
    print(strcat('./bentley/specirrad1/spec1dirrad-', num2str(i), '.pdf'), '-dpdf');
%    save(strcat('./bentley/specirrad1/spec1dirrad-', num2str(i), '.mat'), 'ths_arr', 'ns_arr');
    
    disp(i)
end

function out = try_catch(condition, catch_result)
    try
        out = condition();
    catch
        out = catch_result;
    end
end

function out = replace_nan(num, replace)
    if isnan(num)
        out = replace;
    else
        out = num;
    end
end

function varargout = tern(condition, true_action, false_action)
    if condition() % Works for either a value or function handle.
        [varargout{1:nargout}] = true_action();
    else
        [varargout{1:nargout}] = false_action();
    end
end