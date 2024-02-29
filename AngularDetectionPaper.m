%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ----------  Angular Detection Paper - All Matlab Figures ---------- %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%clf; clc; clear all; close all;
set(0,'DefaultFigureWindowStyle','normal');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%Colours for Plots
plotblue = [0 0.4470 0.7410];
plotgold = [0.9290 0.6940 0.1250];
plotorange = [0.8500 0.3250 0.0980];
plotgreen = [28/255 120/255 37/255];
plotred = [0.6350 0.0780 0.1840];

%Font Choice and Text Size
plot_font = 'Arial';
tick_font_size = 30;
axes_font_size = 35;

%Colours for Plots (4e)
Scatter1 = [244/255 202/255 0/255];
Scatter2 = [56/255 129/255 161/255];
Scatter3 = [20/255 48/255 141/255];

%% Figure 2a - LiF angular scan

load('data/An004665.mat')

figure('Color','white','Units', 'centimeters','Position',[1 1 35 20],'Resize', 'off');
axis equal
axis tickaligned
box on

% plot(angle_vec/1e6,counts/1000,'.','MarkerSize',24)
plot(angle_vec/1e6,counts/1000,'Marker',"square", 'MarkerFaceColor', 'white','LineStyle','none', 'LineWidth', 2.5, 'MarkerSize', 10)
xlabel('{\theta / ^\circ}','fontname',plot_font,'fontsize',axes_font_size)
ylabel('Intensity x 10^3','fontname',plot_font,'fontsize',axes_font_size)

ax = gca;
ax.XAxis.FontSize = tick_font_size;
ax.YAxis.FontSize = tick_font_size;
ax.LineWidth = 2;
xlim([-10 60])

text(38,7,'$\left(0,0\right)$','FontSize',axes_font_size,'interpreter','latex')
text(16,23,'$\left(\overline{1},\overline{1}\right)$','FontSize',axes_font_size,'interpreter','latex')
text(-2,13,'$\left(\overline{2},\overline{2}\right)$','FontSize',axes_font_size,'interpreter','latex')




% exportgraphics(gcf,'LiF_AngularScan.eps','ContentType','vector')

%% Figure 2b - Plot LiF image and position of the angular scan
scan_position=starting_position;


%Plot an image

% SHeM image
SHeM_data=read_MkII_data('data/23-Aug-2023_002.txt');


% Flip the image around
SHeM_data.image=rot90(flipud(SHeM_data.image),3);

% % Remove spikes
% SHeM_data.image=spike_im_removal(SHeM_data.image,3);
% 
% % Remove the median from each row to flatten?
% SHeM_data.image=median_diff_removal(SHeM_data.image);


SHeM_fig_process_Newcastle(SHeM_data,0,1000, 0)

plot(-scan_position(1)/1000,-scan_position(2)/1000,'r.','MarkerSize',20,'LineWidth',3)
% exportgraphics(gcf,'LiF_Image.pdf','ContentType','vector')

SHeM_img_norm=255*((SHeM_data.image-min(SHeM_data.image(:)))/(max(SHeM_data.image(:)-min(SHeM_data.image(:)))));
% imwrite(uint8(imresize(SHeM_img_norm,10,'nearest')),'LiF_Image_clean.png')





%% Figure 3a - HOPG SHeM Scans

%Load Data from array

storedstructure = load('data/An000081.mat','angle_vec','counts');
angles_1 = storedstructure.angle_vec;
counts_1 = storedstructure.counts;

storedstructure = load('data/An000105.mat','angle_vec','counts');
angles_2 = storedstructure.angle_vec;
counts_2 = storedstructure.counts;


%Plotting Command
figure('Color','white','Units', 'centimeters','Position',[1 1 25 20],'Resize', 'off');
axis equal
axis tickaligned
box on
plot((angles_1/1e+6)+90,counts_2/1E4,'Marker','diamond', 'MarkerEdgeColor', plotred, 'MarkerFaceColor', 'white','LineStyle','none', 'LineWidth', 2.5, 'MarkerSize', 10)
hold on
plot((angles_2/1e+6)+90,counts_1/1E4,'Marker', 'o', 'MarkerEdgeColor', plotblue, 'MarkerFaceColor', 'white','LineStyle','none', 'LineWidth', 2.5, 'MarkerSize', 10)
hold off
legend('235^{\circ}C','30^{\circ}C', 'Location','northwest','fontname',plot_font,'fontsize',tick_font_size)
legend('boxoff')
xlabel('{\theta / ^\circ}','fontname',plot_font,'fontsize',axes_font_size)
ylabel('Intensity (counts x 10^4)','fontname',plot_font,'fontsize',axes_font_size)
xlim([0 70]);
ylim([0 3]);
ax = gca;
ax.XAxis.FontSize = tick_font_size;
ax.YAxis.FontSize = tick_font_size;
ax.LineWidth = 2;

%% Figure 3b,c,d,e - HOPG images
Full_HOPG_plot


%% Figure 4e - Helium Beam Interaction with AFM Grating


%Load Data

load('data\An001886.mat');
angles_grating_parallel = angle_vec;
counts_grating_parallel = counts;

load('data\An001887.mat');
angles_frame_parallel = angle_vec;
counts_frame_parallel = counts;

load('data\An001888.mat');
angles_grating_perpendicular = angle_vec;
counts_grating_perpendicular = counts;


%Plotting Command
figure('Color','white','Units', 'centimeters','Position',[1 1 35 20],'Resize', 'off');
axis equal
axis tickaligned
box on
plot(angles_grating_parallel/1e+6,counts_grating_parallel/1000,'Marker', 'o', 'MarkerEdgeColor', Scatter1, 'MarkerFaceColor', 'white','LineStyle','none', 'LineWidth', 2, 'MarkerSize', 8)
hold on
plot(angles_frame_parallel/1e+6,counts_frame_parallel/1000,'Marker','Square', 'MarkerEdgeColor', Scatter2, 'MarkerFaceColor', 'white','LineStyle','none', 'LineWidth', 2, 'MarkerSize', 8)
hold on
plot(angles_grating_perpendicular/1e+6,counts_grating_perpendicular/1000,'Marker','Diamond', 'MarkerEdgeColor', Scatter3, 'MarkerFaceColor', 'white','LineStyle','none', 'LineWidth', 2, 'MarkerSize', 8)
hold off
legend('Grating (//)','Outer Frame','Grating ({\perp})', 'Location','northeast','fontname',plot_font,'fontsize',tick_font_size)
legend('boxoff')
xlabel('{\theta / ^\circ}','fontname',plot_font,'fontsize',axes_font_size)
ylabel('Intensity (counts x 10^3)','fontname',plot_font,'fontsize',axes_font_size)
xlim([-10 100]);
ylim([0 8]);
ax = gca;
ax.XAxis.FontSize = tick_font_size;
ax.YAxis.FontSize = tick_font_size;
ax.LineWidth = 2;

%% Figure 5 - Scattering Distribution Comparisons

% Load data
files_ind_1 = 2084:2134;
files_ind_2 = 2137:2176;
files_ind = [files_ind_1,files_ind_2];
files_ind_3 = 1986:2030;
files_ind_4 = 2036:2081;
files_ind_5 = [files_ind_3,files_ind_4];
path1 = 'data\TGZ4 scans';

[data, thetas, phis] = load_an_scans(files_ind, path1);
[data_2, thetas_2, ~] = load_an_scans(files_ind_5, path1);

% phis=phis;
phis_3 = 11:4:187;
phis_4 = 189:-4:9;
phis_2 = [phis_3,phis_4];

[theta_sorted,~]=sort(thetas);
[phi_sorted,phi_ind]=sort(phis);

[theta_sorted_2,theta_ind_2]=sort(thetas_2);
[phi_sorted_2,phi_ind_2]=sort(phis_2);


%Remove spikes
data=spike_im_removal(data,4);
data_2=spike_im_removal(data_2,4);


exp_data = data+data_2;

phi_data = phi_sorted_2-99;
theta_data = theta_sorted_2;
theta_ind = theta_ind_2;

load('data\Gratings_output.mat')


tick_font_size = 24;


% Plotting Command

figure('Units', 'centimeters','Position',[1 1 35 25]);
B = tiledlayout(1,2,'TileSpacing','Compact','Padding','Compact');
% B = tiledlayout('flow','TileSpacing','Compact','Padding','Compact');
B.Units = 'centimeters';
% B.InnerPosition = [0.13 0.13 0.775 0.8];

colormap_spectra = parula;

nexttile()
imagesc(theta_data,phi_data,rot90(rescale(exp_data(theta_ind,phi_ind))));
set(gca,'dataAspectRatio',[1 1 1]);
axis equal;
axis tickaligned;
box on;
clim([0 1]);
xlim([-10 90]);
xticks([-10 10 30 50 70 90]);
ylim([-90 90]);
yticks(-90:30:90);
colormap(colormap_spectra);
ax = gca;
ax.XAxis.FontSize = tick_font_size;
ax.YAxis.FontSize = tick_font_size;
ax.LineWidth = 1.25;
ylabel('{\alpha / ^\circ}','fontname',plot_font,'fontsize',tick_font_size);
xlabel('{\theta / ^\circ}','fontname',plot_font,'fontsize',tick_font_size);


nexttile()
imagesc(phi_c,alpha_vec_tot,rescale(I_mat_tot'));
set(gca,'dataAspectRatio',[1 1 1]);
axis equal;
axis tickaligned;
box on;
clim([0 1]);
xlim([-10 90]);
xticks([-10 10 30 50 70 90]);
ylim([-90 90]);
yticks(-90:30:90);
colormap(colormap_spectra);
ax = gca;
ax.XAxis.FontSize = tick_font_size;
ax.YAxis.FontSize = tick_font_size;
ax.LineWidth = 1.25;
ylabel('{\alpha / ^\circ}','fontname',plot_font,'fontsize',tick_font_size);
xlabel('{\theta / ^\circ}','fontname',plot_font,'fontsize',tick_font_size);

%nexttile([1 2])
cb = colorbar(); 
cb.Layout.Tile = 'south'; % Assign colorbar location
cb.Limits = [0 1];
cb.Ticks = [0, 0.25 0.5 0.75, 1];
cb.LineWidth = 1.25;
cb.FontName = 'Arial';
cb.FontSize = 18;
cb.Position = [0.1075 0.10 0.8477 0.0121];


%exportgraphics(B,'Fig2_c+d.pdf','ContentType','vector')



%% Figure 6a - Helium Intensity with Changing Detector Angle

tick_font_size = 30;

%Load Data

figure3a_data = load('data\Fig3a_Data_V2.mat');
figure3b_data = load('data\Fig3b_Data.mat');

Fig3a_Angles = figure3a_data.Fig3a_Data(1,:);
Fig3a_Intensities = figure3a_data.Fig3a_Data(2,:)/1000;
Fig3a_Residuals = figure3a_data.Fig3a_Data(3,:)/1000;
Fig3a_xconf = [Fig3a_Angles, Fig3a_Angles(end:-1:1)];
Fig3a_yconf = [Fig3a_Intensities-Fig3a_Residuals, Fig3a_Intensities(end:-1:1)+Fig3a_Residuals(end:-1:1)];


%Plotting Command
fig = figure('Color','white','Units', 'centimeters','Position',[1 1 35 20],'Resize', 'off');
axis equal
axis tickaligned
box on
xlabel('{\theta / ^\circ}','fontname',plot_font,'fontsize',axes_font_size)
yyaxis right
p1=plot(figure3b_data.FAHY(1,:),figure3b_data.FAHY(2,:)/1000,'-','Color', plotblue, 'LineWidth', 4);
hold on
p2=plot(figure3b_data.FAHY(1,:),figure3b_data.FAHY(3,:)/1000,'-','Color', plotgold, 'LineWidth', 4);
hold on
p3=plot(figure3b_data.FAHY(1,:),figure3b_data.FAHY(4,:)/1000,'-', 'Color', plotred, 'LineWidth', 4);
hold off
ylabel('Simulated rays x 10^3','fontname',plot_font,'fontsize',axes_font_size)
xlim([-100 100]);
ylim([0 15]);

% Setting the plot order
p1.ZData = zeros(size(p1.XData));
p2.ZData = zeros(size(p2.XData));
p3.ZData = zeros(size(p3.XData));

yyaxis left
p4=plot(figure3a_data.Fig3a_Data(1,1:2:end),figure3a_data.Fig3a_Data(2,1:2:end)/1000,'Marker', 'o', 'MarkerEdgeColor', plotgreen, 'MarkerFaceColor', 'white','LineStyle','none', 'LineWidth', 2.5, 'MarkerSize', 10);
ylabel('Intensity (Counts x 10^3)','fontname',plot_font,'fontsize',axes_font_size)
ylim([2.3 13])
hold on
p4.ZData = ones(size(p4.XData));
set(gca, 'SortMethod', 'depth')
legend('Experiment','{ h_{rms} = 0}','{ h_{rms} = \lambda}','{ h_{rms} = 4\lambda}', 'Location','northeast','fontname',plot_font,'fontsize',26)
legend('boxoff')
ax = gca;
ax.XAxis.FontSize = tick_font_size;
ax.YAxis(1).FontSize = tick_font_size;
ax.YAxis(2).FontSize = tick_font_size;
ax.YAxis(1).Color = 'black';
ax.YAxis(2).Color = 'black';
ax.LineWidth = 2;


%% Figure 6b - Isometric Surfaces


% Colours
rgb_blue = [1 1 1; 19/256 163/256 224/256; 0/256 85/256 207/256; 19/256 163/256 224/256; 1 1 1];
rgb_green = [1 1 1; 4/256 201/256 4/256; 0/256 96/256 0/256; 4/256 201/256 4/256; 1 1 1];
rgb_yellow = [1 1 1; 255/256 219/256 74/256; 240/256 156/256 0/256; 255/256 219/256 74/256; 1 1 1];
rgb_red = [1 1 1; 238/256 0/256 4/256; 125/256 6/256 8/256; 238/256 0/256 4/256; 1 1 1];
pts = [1 4.5 5 5.5 9];
custom_blue = colormap(multigradient(rgb_blue, pts));
custom_green = colormap(multigradient(rgb_yellow, pts));
custom_red = colormap(multigradient(rgb_red, pts));


% Load Data - Surface 1
% load('data\bentley_scatterAndSave\data-cosine-300-0d01000-0d00000-0d05000-45d00-0.mat')
load('data\raytrace_surface_ratio_0.mat')


% Surface Parameters
surface_size = N*deltax;
x = linspace(-surface_size/2,surface_size/2,N); y = linspace(-surface_size/2,surface_size/2,N);

% Surface 1 Plotting
fig_1=figure;
rectangle('Position', [-1.5, -1.5, 3, 3]);
view(45, 35.264);
hold on
surf(x,y,surfImg)
grid on
axis off
shading flat
axis tight
axis equal
colormap(custom_blue);
clim([-1, 1]);
hold off
%exportgraphics(fig_1,'surface_0.png','Resolution',750)


% Load Data - Surface 2
% load('data\bentley_scatterAndSave\data-cosine-300-0d01000-0d05000-0d05000-45d00-0.mat')
load('data\raytrace_surface_ratio_1.mat')


% Surface Parameters
surface_size = N*deltax;
x = linspace(-surface_size/2,surface_size/2,N); y = linspace(-surface_size/2,surface_size/2,N);

% Surface 2 Plotting
fig_2=figure;
rectangle('Position', [-1.5, -1.5, 3, 3]);
view(45, 35.264);
hold on
surf(x,y,surfImg)
grid on
axis off
shading flat
axis tight
axis equal
colormap(custom_green);
clim([-1, 1]);
hold off
%exportgraphics(fig_2,'surface_1.png','Resolution',750)


% Load Data - Surface 3
% load('data\bentley_scatterAndSave\data-cosine-300-0d01000-0d20000-0d05000-45d00-0.mat')
load('data\raytrace_surface_ratio_4.mat')
clear plot

% Surface Parameters
surface_size = N*deltax;
x = linspace(-surface_size/2,surface_size/2,N); y = linspace(-surface_size/2,surface_size/2,N);

% Surface 3 Plotting
fig_3=figure;
rectangle('Position', [-1.5, -1.5, 3, 3]);
view(45, 35.264);
hold on
surf(x,y,surfImg)
grid off
axis off
shading flat
axis tight
axis equal
colormap(custom_red);
clim([-1, 1]);
hold off
%exportgraphics(fig_3,'surface_4.png','Resolution',750)


%% Figure 7 - Model scattering distributions

Plot_all_v4

%% Figure B9 - Mask Length Fitting

mask_length_script

%% Figure C10 - Contamination rate

L_plot_loop