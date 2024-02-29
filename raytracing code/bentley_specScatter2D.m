addpath('functions', 'random_surface_generation', genpath('ParforProgMon'))

hss = linspace(0.01,0.1,100);

%lambda = 1;
%lambdas = 1./linspace(0.0001,3,100);
lambdas = [1];
hs = 1;
%hs = 1;
dx = 1/100;

loops = 1;
%loops = 1;

%for i=1:length(hss)
    %hs = hss(i);
for i=1:length(lambdas)
    lambda = lambdas(i);
    
    in_angle = 20;
    
    ths_all = [];
    ns_all = [];
    
    for asdf=1:loops
        [fs, xs] = rsgene1D(20000, 20000/100, hs/sqrt(coth(1/100/lambda)), lambda);
        tic;
        [ths, ns] = random_scatter1D('hs', fs, 'xs', xs, 'init_angle', in_angle, 'n_rays', 10000, 'scattering', 'specular');
        toc

        ths_all = [ths_all, ths];
        ns_all = [ns_all, ns];
        
    end
    
    ths = ths_all;
    ns = ns_all;
    
    all_scatter = pdf(fitdist(ths', 'Kernel'), -90:0.1:90);
    try
        one_scatter = pdf(fitdist(ths(ns==1)', 'Kernel'), -90:0.1:90);
    catch 
        one_scatter = zeros(1801,1);
    end
    try
        multi_scatter = pdf(fitdist(ths(ns~=1)', 'Kernel'), -90:0.1:90);
    catch
        multi_scatter = zeros(1801,1);
    end
    
    ratio_one_scatter = length(ths(ns==1)) / length(ths);
    one_scatter = ratio_one_scatter * one_scatter;
    multi_scatter = (1-ratio_one_scatter) * multi_scatter;

    figure(1)
    clf
    plot(-90:0.1:90, all_scatter, 'Linewidth', 2)
    hold on
    plot(-90:0.1:90, one_scatter, 'Linewidth', 2)
    plot(-90:0.1:90, multi_scatter, 'Linewidth', 2)
    
    tic;
    ths = linspace(-90+in_angle,in_angle+90,1000)/2;
    
    % Step one: the PDF of being seen with flux conservation but without
    % self-shadowing.
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
    %mask_or_shadow = @(th) 1 - 1./( (1+biglambda(abs(in_angle))) .* (1+biglambda(abs(2*th-in_angle))) );
    mask_or_shadow = @(th) 1 - 1./( (1+biglambda(max(abs(in_angle), abs(2*th-in_angle)))) );
    %mask_or_shadow3 = @(th) 1 - 1./( (1+biglambda_orig(max(abs(in_angle), abs(2*th-in_angle)))) );
    %mask_or_shadow = @(th) 1 - 1./( (1+biglambda(abs(in_angle)) + biglambda(abs(2*th-in_angle))) );
    p2 = p1 .* (1-mask_or_shadow(ths));
    %p4 = p1 .* (1-mask_or_shadow3(ths));


    % Step three: convert to outgoing angle distribution
    alphas = -2*ths + in_angle;
    p3 = p2 /2; % * pi / 360;
    %p5 = p4/2;
    
    toc
    
    plot(alphas, p3, 'Linewidth', 2)
    %plot(alphas, p5, 'Linewidth', 2)
    xlim([-90, 90])
    ylim([0, 0.03])
    xlabel('Outgoing Angle / degree')
    ylabel('PDF')
    title(sprintf('Specular Scattering: Outgoing Angles PDF\nIncoming Angle=%d°, h=%g, \\lambda=%g', in_angle, hs, lambda))
    legend({'All scattered rays', 'Single Scattered only', 'Multi Scattered only', 'Single Scattered (Analytical)'}, 'location', 'eastoutside')
    
    fig = gcf;
    set(fig,'Units','Inches');
    set(fig,'Position',[7,5,7,4.3]);
    pos = get(fig,'Position');
    set(fig,'Paperorientation', 'landscape', 'PaperType', 'a3', 'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
    print(strcat('./bentley/spec_scatter3/spec-', num2str(i), '.pdf'), '-dpdf');
    
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