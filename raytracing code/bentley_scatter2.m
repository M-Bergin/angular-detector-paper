% Stuff by Bentley
% Needs Sams stuff: https://github.com/slambrick/SHeM-Ray-Tracing-Simulation
% Also needs hist2d: https://au.mathworks.com/matlabcentral/fileexchange/66629-2-d-histogram-plot
%
% Use like dis:
% bentley_scatterRays3D(N, deltax, hrms, lambda, init_theta, n_rays, recompile, verbose, plot);
% Will default to n_rays = 200000.
% Will default to recompile = false.
% Will default to verbose = false.
% Will default to plot = false.
% Returns [az,el] which is a list of all the reflected rays azimuthal angle
% and elevation angle from the horizontal. Note that the length of az and
% el may be less than n_rays.

function [singleScattered,all,numScattersRay,surfImg] = bentley_scatterRays3D(N, deltax, hrms, lambda, init_theta, n_rays, recompile, verbose)
    if (nargin < 5 || nargin > 8)
        error('Can you plz give me right number of arguments thanks x')
    end
    
    if (nargin < 6)
        n_rays = 200000;
    end
    
    if (nargin < 7)
        recompile = false;
    end
    
    if (nargin < 8)
        verbose = false;
    end
    
    if (nargin < 9)
        plot = false;
    end
    
    %% Set Parameters
    % Starting direction of the rays
    init_dir = [sind(init_theta), 0, -cosd(init_theta)];

    % The maximum number of scattering events a ray is allowed to undergo
    maxScatter = 1000;

    % Surface Parameters
    surface_size = N*deltax;
    h = hrms*tanh(deltax/lambda);

    % Should the mex files be recompiled
    loadpath
    if recompile
        mexCompile();
    end

    %% Generate surface
    x = linspace(-surface_size/2,surface_size/2,N); y = linspace(-surface_size/2,surface_size/2,N);
    [X,Y] = meshgrid(x,y); 
    Z = h.*randn(N,N);
    F = exp(-(abs(X)+abs(Y))/lambda);
    f = ifft2(fft2(Z).*fft2(F));
    
    imagesc(x,y,f(1:100:end, 1:100:end)); % TODO remove
    error('stop here') % TODO remove
    
    surfImg  = f;

    vertices = [X(:), Y(:), f(:)];

    % Calculate the triangulation for the lower left triangles.
    nums1 = 1:N^2;
    upper1 = circshift(nums1, -1);
    right1 = circshift(nums1, -N);
    filt = ((mod(nums1,N)~=0) .* (nums1 < N*(N-1))) == 1;
    nums1 = nums1(filt);
    upper1 = upper1(filt);
    right1 = right1(filt);
    fdef1 = [nums1;upper1;right1];

    fnorm1 = [f(right1)-f(nums1); f(upper1)-f(nums1); deltax*ones(size(nums1))]; % Can check by cross product
    fnorm1 = normalize(fnorm1, 1, 'norm', 2);

    % Calculate the triangulation for the upper right triangles.
    nums2 = 1:N^2;
    lower2 = circshift(nums2, 1);
    left2 = circshift(nums2, N);
    filt = ((mod(nums2,N)~=1) .* (nums2 > N)) == 1;
    nums2 = nums2(filt);
    lower2 = lower2(filt);
    left2 = left2(filt);
    fdef2 = [nums2; lower2; left2];

    fnorm2 = [f(nums2)-f(left2); f(nums2)-f(lower2); deltax*ones(size(nums2))];
    fnorm2 = normalize(fnorm2, 1, 'norm', 2);

    fdef = [fdef1, fdef2];
    fnorm = [fnorm1, fnorm2];
    [fmat{1:length(fdef)}] = deal('poop');
    mat = struct('function', 'pure_specular', 'params', [], 'color', [0.8,0.8,1]);
    materials = containers.Map('poop', mat);

    surface = TriagSurface(vertices, fdef', fnorm', fmat', materials);
    surface.moveBy([0, 0, -max(surface.vertices(:,3))]);
%         figure(1);
%         clf;
%         surface.patchPlot(false);

    %% Let's bounce!
    init_dir = repmat(init_dir', 1, n_rays);
    init_pos = [unifrnd(-surface_size/2,surface_size/2, 2, n_rays); zeros(1,n_rays)];
    init_pos = init_pos - 10*init_dir;

    if verbose
        tic
    end
    [~, numScattersRay, ~, final_dir] = ...
        distributionCalc('sample_surface', surface, 'maxScatter', maxScatter, ...
                         'nrays', n_rays, 'start_pos', init_pos, 'start_dir', init_dir);
    if verbose
        fprintf('Ray tracing took %d seconds.\n', uint8(toc));
    end

    [az, el] = cart2sph(final_dir(1,:), final_dir(2,:), final_dir(3,:));

    az = az*180/pi;
    el = el*180/pi;
    
    az_scattered = az(el>0);
    el_scattered = el(el>0);
    
    singleScattered = [az_scattered; 90-el_scattered];
    all = [az; 90-el];
    
    if verbose
        fprintf('%d%% of the rays undergo at least one scattering event.\n', round(100 - length(numScattersRay(numScattersRay == 0)) / n_rays *100));
        fprintf('%d%% of the rays undergo multiple scattering events.\n', round(length(numScattersRay(numScattersRay > 1)) / n_rays *100));
        fprintf('%d%% of the rays have positive elevation.\n', round(length(el(el>0)) / n_rays *100));
    end
end