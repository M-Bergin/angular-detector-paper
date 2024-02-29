clear;

% num_surfaces = 20;
% 
% hrmss = logspace(-3,-0.5,num_surfaces);

num_surfaces = 3;

hrmss = [0,0.05,0.2];%logspace(-3,-0.5,num_surfaces);
lambda = 0.05;
init_theta = 45;

N = 300;
deltax = 1/100;
num_rays = 2000000;

% N = 20;
% deltax = 1/100;
% num_rays = 2000;

tic
for i=1:length(hrmss)
    hrms = hrmss(i);
%     bentley_scatterAndSave3D(N, deltax, hrms, lambda, init_theta, num_rays, 'distribution', 'specular', 'verbose', true, 'plot', true);
    bentley_scatterAndSave3D(N, deltax, hrms, lambda, init_theta, num_rays, 'distribution', 'cosine', 'verbose', true, 'plot', true);
    %bentley_loadAndSave3D(N, deltax, hrms, lambda, init_theta, 'distribution', 'specular', 'verbose', true, 'plot', true);
    fprintf('Done %d of %d.\n', i, num_surfaces);
end
fprintf('Scattering off multiple surfaces took %d seconds.\n', uint16(toc))
