function bentley_scatterAndSave3D(N, deltax, hrms, lambda, init_theta, n_rays, varargin)
    DIR = 'bentley_scatterAndSave';
    
    %% Parameters
    if (nargin < 5)
        error('Can you plz give me right number of arguments thanks x')
    end
    
    if (nargin < 6)
        n_rays = 200000;
    end
    
    recompile = false;
    verbose = false;
    plot = false;
    distribution = 'specular';
    for i_=1:2:length(varargin)
        switch varargin{i_}
            case 'recompile'
                recompile = varargin{i_+1};
            case 'plot'
                plot = varargin{i_+1};
            case 'verbose'
                verbose = varargin{i_+1};
            case 'distribution'
                distribution = varargin{i_+1};
        end
    end
    
    if strcmp(distribution, 'specular') && strcmp(distribution, 'cosine')
        error('I dunno what a %s distribution is.', distribution);
    end
    
    loadpath;
    if recompile
        mexCompile();
    end
    
    %% Scatter
    if verbose
        fprintf('Starting scattering for hrms=%.3f, lambda=%.3f.\n', hrms, lambda);
        tic
    end
    [~, all_dirs, numScattersRay, initPos, surfImg] = bentley_scatterRays3D(N, deltax, hrms, lambda, init_theta, n_rays, 'verbose', verbose, 'distribution', distribution);
    if verbose
        fprintf('Finished scattering. Surface generation and scattering took %d seconds.\n', uint16(toc));
        tic
    end
    
    %% Save
    filenameprefix = strrep(sprintf('data-%s-%d-%.5f-%.5f-%.5f-%.2f-', distribution, N, deltax, hrms, lambda, init_theta), '.', 'd');
    
    if ~exist(DIR, 'dir')
        mkdir(DIR);
    end
    
    i=0;
    while isfile(sprintf('%s/%s%d.mat', DIR, filenameprefix, i))
        i = i+1;
    end
    save(sprintf('%s/%s%d.mat', DIR, filenameprefix, i), '-mat');
    
    if verbose
        fprintf('Saving took %d seconds.\n', uint16(toc));
        tic
    end
    
    if verbose && plot
        fprintf('Now plotting.\n');
        tic
    end
    
    %% Plot
    if plot
        if strcmp(distribution, 'specular')
            bentley_plotScatter3D(all_dirs, numScattersRay, surfImg, N, deltax, hrms, lambda, init_theta, 'distribution', distribution);
        end
        if strcmp(distribution, 'cosine')
            bentley_plotScatter3D_DIFFUSE(all_dirs, numScattersRay, surfImg, N, deltax, hrms, lambda, init_theta, 'distribution', distribution);
        end
        drawnow;
        fig = gcf;
        set(fig,'Units','Inches');
        pos = get(fig,'Position');
        set(fig,'Paperorientation', 'landscape', 'PaperType', 'a3', 'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
        print(sprintf('%s/%s%d.pdf', DIR, filenameprefix, i), '-dpdf');
    end
    
    if verbose && plot
        fprintf('Plotting took %d seconds.\n', uint16(toc))
    end
    
end