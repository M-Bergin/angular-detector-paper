% Copyright (c) 2018, Sam Lambrick.
% All rights reserved.
% This file is part of the Sub-beam Ray Tracing Simulation, subject to the  
% GNU/GPL-3.0-or-later.
%
% random_scatter1D.m
%
% Scatters rays uniformally off a randomly rough surface specified by the ratio
% between the rms height and the correlation length. Provide the number of rays
% and surface elements.
%     The surface to be scattered off can either be specified by its statistics
% or by explicity giving the x and y positions of the profile - do not give
% both, the function will error if you do. Specify inputs as name-value pairs:
% ('name', value).
%
% INPUTS:
%  ratio      - The ratio between the 
%  Nelements  - The number of 1D surface elements to model.
%  xs         - Alternativly specify the surface to be scattered off explicitly,
%               specify the x positions of the surfcace
%  hs         - Specuify the height positions of the surface
%  n_rays     - The number of rays to model.
%  init_angle - The incident angle specified in degrees.
%  make_plots - Boolean. Should plots be made? Plots of an example region of
%               surface, the height distribution function of the generated
%               surface, and histograms of the two components of the final
%               directions and the number of scattering events. Optional,
%               defaults to false.
%  scattering - string, the type of scattering to perform
%
% OUTPUTS:
%  thetas       - A vector of the final angles to the normal of the rays after
%                 they have been scattered off the random surface.
%  num_scatters - The number of scattering events that the rays have undergone.
function [thetas, num_scatters] = random_scatter1D(varargin)
    % Get inputs
    for i_=1:2:length(varargin)
        switch varargin{i_}
            case 'ratio'
                ratio = varargin{i_+1};
            case 'Nelements'
                Nelements = varargin{i_+1};
            case 'xs'
                xs = varargin{i_+1};
            case 'hs'
                hs = varargin{i_+1};
            case 'n_rays'
                n_rays = varargin{i_+1};
            case 'init_angle'
                init_angle = varargin{i_+1};
            case 'make_plots'
                make_plots = varargin{i_+1};
            case 'scattering'
                scattering = varargin{i_+1};
            case 'scattering_parameters'
                scattering_parameters = varargin{i_+1};
        end
    end
    
    % Check that inputs exist and default values for non-existant inputs
    if ~exist('n_rays', 'var')
        error('Specify a number of rays.')
    end
    if ~exist('init_angle', 'var')
        error('Specify an incident angle.')
    end
    if ~exist('make_plots', 'var')
        make_plots = false;
    end
    
    % If the type of scattering isn't broad specular then the scattering
    % parameter isn't used
    if ~strcmp(scattering, 'broad specular')
        scattering_parameters = [0]; %#ok<NBRAK>
    elseif ~exist('scattering_parameters', 'var')
        error(['If using the broad specular distribution as the intrinsic ' ... 
            'distribution then you must specify the broadness parameter.'])
    end
    
    % Switch between the cases where we are provided a surface profile and
    % provided the statistics on the surface profile
    if exist('ratio', 'var') && exist('Nelements', 'var')
        if exist('xs', 'var') || exist('hs', 'var')
            error(['Do not specify both a surface profile and statistic for' ...
                'gernerating a surface profile. Do one or the other.'])
        end
        % Generate a random Gaussian surface with a Gaussian height distribution
        % function and an exponential correlation length.
        [hs, xs] = rsgene1D(Nelements+1, (Nelements+1)/100, ratio, 1);
    elseif exist('xs', 'var') && exist('hs', 'var')
        if exist('ratio', 'var') || exist('Nelements', 'var')
            error(['Do not specify both a surface profile and statistic for' ...
                'gernerating a surface profile. Do one or the other.'])
        end
    else
        error(['Provide either a surface profile or the statistics for a' ...
            'profile to be generated.'])
    end
    
    % Make plots if we want
    if make_plots
        % Plot a section of surface with equal axes
        figure
        plot(xs, hs)
        xlim([-ratio*20, ratio*20])
        axis equal
        title('Example surface profile')
    
        % Plot a histogram of the heights overlayed with the Gaussian distribution it
        % was generated from
        figure
        histogram(hs, 'Normalization', 'pdf')
        hold on
        height = linspace(min(hs), max(hs), 1000);
        plot(height, normpdf(height, 0, ratio))
        hold off
        xlabel('Height of surface in units of the correlation length')
        legend('Generated surface', 'Ideal surface')
    end
       
    % The minimum and maxium nominal intersection points
    min_x = min(xs)/2;
    max_x = max(xs)/2;
    
    % The starting positions
    init_y = range(xs);
    init_x = linspace(min_x, max_x, n_rays);
    init_pos = [init_x; repmat(init_y, 1, n_rays)];
    init_pos = init_pos - [tand(init_angle)*range(xs); 0];
    
    % The starting directions
    init_dir = [sind(init_angle); -cosd(init_angle)];
    
    % Calculate the normals to each surface elements
    vertices = [xs; hs];
    normals = [hs(1:end-1) - hs(2:end); abs(xs(2:end) - xs(1:end-1))];
    normals = normals./sqrt(sum(normals.^2, 1));
    
    % Scatter the rays using the ray tracing mex program
    [final_dirs, ~, num_scatters] = scatter_rays({vertices, normals}, ...
        {init_pos, init_dir}, {scattering, scattering_parameters});
    
    % Calculate the angle to the surface normal for the final directions of the rays
    thetas = atand(final_dirs(1,:)./final_dirs(2,:));
    
    if make_plots       
        % Plot a histogram of the final x and y componenets of the directions
        figure
        histogram(final_dirs(1,:), 'Normalization', 'pdf')
        xlabel('x component of direction')
        
        figure
        histogram(final_dirs(2,:), 'Normalization', 'pdf')
        xlabel('y component of direction')
        
        % Plot a histogram of the number of scattering events
        figure
        histogram(num_scatters, 'BinWidth', 1, 'Normalization', 'pdf')
        xlabel('Number of scatters')
    end
end

