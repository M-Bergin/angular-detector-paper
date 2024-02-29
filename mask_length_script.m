%Load the data
load('data\AL000013.mat');

clip_max=100;


%% 
% Code to find steps using inbuilt ischange function

counts_mat_sorted(172,55) = nan;

% Find where the two steps in the data are
[TF,S1] = ischange(counts_mat_sorted,'linear',1,'MaxNumChanges',2);


N=size(counts_mat_sorted);
N_length=min(N(2),clip_max);

upper_bound=NaN*zeros(1,N_length);
lower_bound=NaN*zeros(1,N_length);
for n=1:N_length
    temp_inds=find(TF(:,n));
    if length(temp_inds)==2 && n<clip_max
        upper_bound(n)=temp_inds(1);
        lower_bound(n)=temp_inds(2);
    else
        upper_bound(n)=NaN;
        lower_bound(n)=NaN;
    end
end


%% 
% Determine the positions of the steps as x co-ordinate


upper_x=NaN*zeros(1,N_length);
lower_x=NaN*zeros(1,N_length);

non_nan_inds=find(~isnan(upper_bound));

upper_x(non_nan_inds)=x_sorted(upper_bound(non_nan_inds))/1e6;
lower_x(non_nan_inds)=x_sorted(lower_bound(non_nan_inds))/1e6;

% % Plot the bounds of the mask to check it is sensible
% figure
% imagesc((angle_sorted/1e6),x_sorted/1e6,counts_mat_sorted);
% hold on
% plot(angle_sorted(1:N_length)/1e6,upper_x,'r.','Linewidth',2,'MarkerSize',12)
% plot(angle_sorted(1:N_length)/1e6,lower_x,'r.','Linewidth',2,'MarkerSize',12)
% 
% xlabel('Angle/^\circ')
% ylabel('x/mm')

%Calculate length of the mask and plot
Length_mm_1=abs(upper_x-lower_x);

% figure
% plot(angle_sorted(1:N_length)/1e6,Length_mm_1,'.','Markersize',12);
% xlabel('Angle/^\circ')
% ylabel('Mask length /mm')
% set(gca,'Fontsize',12,'Linewidth',1)

%% 
% Now can fit to the length vector

[xData, yData] = prepareCurveData( angle_sorted(1:N_length)/1e6, Length_mm_1 );

% Set up fittype and options.
ft = fittype( 'h*(abs((alpha-sind(x))/(beta-cosd(x)))+(alpha-sind(x))/(beta-cosd(x)))/2+h2', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.Lower = [-Inf -Inf -Inf -Inf];
opts.StartPoint = [0 0 0.8 0.8];

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );

% % Plot fit with data.
% figure
% theta_plot=linspace(min(xData),max(xData),200);
% fit_plot=fitresult(theta_plot);
% p11= predint(fitresult,theta_plot,0.6827,'observation','off');
% h2= plot(theta_plot,fit_plot,'LineWidth',2,"Color",	[0.8500, 0.3250, 0.0980]);
% hold on
% %h3 = plot(theta_plot,p11,'--','LineWidth',1,"Color",	[0.8500, 0.3250, 0.0980]);
% h = plot(xData, yData,'.','MarkerSize',12,'Color', [0, 0.4470, 0.7410] );
% 
% % Label axes
% xlabel( '\theta/^\circ');
% ylabel( "L/mm");
% set(gca,'Fontsize',12,'Linewidth',1)

theta_1=angle_sorted(1:N_length)/1e6;
z_pos_1=double(starting_position(3))/1e6;






















%% 
% Repeat for next dataset

%Load the data
load('data\AL000014.mat');

clip_max=100;


%% 
% Code to find steps using inbuilt ischange function


counts_mat_sorted(163,3) = nan;
counts_mat_sorted(184,71) = nan;


% Find where the two steps in the data are
[TF,S1] = ischange(counts_mat_sorted,'linear',1,'MaxNumChanges',2);


N=size(counts_mat_sorted);
N_length=min(N(2),clip_max);

upper_bound=NaN*zeros(1,N_length);
lower_bound=NaN*zeros(1,N_length);
for n=1:N_length
    temp_inds=find(TF(:,n));
    if length(temp_inds)==2 && n<clip_max
        upper_bound(n)=temp_inds(1);
        lower_bound(n)=temp_inds(2);
    else
        upper_bound(n)=NaN;
        lower_bound(n)=NaN;
    end
end


%% 
% Determine the positions of the steps as x co-ordinate


upper_x=NaN*zeros(1,N_length);
lower_x=NaN*zeros(1,N_length);

non_nan_inds=find(~isnan(upper_bound));

upper_x(non_nan_inds)=x_sorted(upper_bound(non_nan_inds))/1e6;
lower_x(non_nan_inds)=x_sorted(lower_bound(non_nan_inds))/1e6;

% % Plot the bounds of the mask to check it is sensible
% figure
% imagesc((angle_sorted/1e6),x_sorted/1e6,counts_mat_sorted);
% hold on
% plot(angle_sorted(1:N_length)/1e6,upper_x,'r.','Linewidth',2,'MarkerSize',12)
% plot(angle_sorted(1:N_length)/1e6,lower_x,'r.','Linewidth',2,'MarkerSize',12)
% 
% xlabel('Angle/^\circ')
% ylabel('x/mm')

%Calculate length of the mask and plot
Length_mm_2=abs(upper_x-lower_x);

% figure
% plot(angle_sorted(1:N_length)/1e6,Length_mm_2,'.','Markersize',12);
% xlabel('Angle/^\circ')
% ylabel('Mask length /mm')
% set(gca,'Fontsize',12,'Linewidth',1)






%% 
% Now can fit to the length vector

[xData, yData] = prepareCurveData( angle_sorted(1:N_length)/1e6, Length_mm_2 );

% Set up fittype and options.
ft = fittype( 'h*(abs((alpha-sind(x))/(beta-cosd(x)))+(alpha-sind(x))/(beta-cosd(x)))/2+h2', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.Lower = [-Inf -Inf -Inf -Inf];
opts.StartPoint = [0 0 0.8 0.8];

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );

% % Plot fit with data.
% figure
% theta_plot=linspace(min(xData),max(xData),200);
% fit_plot=fitresult(theta_plot);
% p11= predint(fitresult,theta_plot,0.6827,'observation','off');
% h2= plot(theta_plot,fit_plot,'LineWidth',2,"Color",	[0.8500, 0.3250, 0.0980]);
% hold on
% %h3 = plot(theta_plot,p11,'--','LineWidth',1,"Color",	[0.8500, 0.3250, 0.0980]);
% h = plot(xData, yData,'.','MarkerSize',12,'Color', [0, 0.4470, 0.7410] );
% 
% % Label axes
% xlabel( '\theta/^\circ');
% ylabel( "L/mm");
% set(gca,'Fontsize',12,'Linewidth',1)

theta_2=angle_sorted(1:N_length)/1e6;
z_pos_2=double(starting_position(3))/1e6;



























%% 
% Repeat for next dataset


%Load the data
load('data\AL000015.mat');

clip_max=100;


%% 
% Code to find steps using inbuilt ischange function


counts_mat_sorted(171,10) = nan;
counts_mat_sorted(137,19) = nan;


% Find where the two steps in the data are
[TF,S1] = ischange(counts_mat_sorted,'linear',1,'MaxNumChanges',2);


N=size(counts_mat_sorted);
N_length=min(N(2),clip_max);

upper_bound=NaN*zeros(1,N_length);
lower_bound=NaN*zeros(1,N_length);
for n=1:N_length
    temp_inds=find(TF(:,n));
    if length(temp_inds)==2 && n<clip_max
        upper_bound(n)=temp_inds(1);
        lower_bound(n)=temp_inds(2);
    else
        upper_bound(n)=NaN;
        lower_bound(n)=NaN;
    end
end



%% 
% Determine the positions of the steps as x co-ordinate


upper_x=NaN*zeros(1,N_length);
lower_x=NaN*zeros(1,N_length);

non_nan_inds=find(~isnan(upper_bound));

upper_x(non_nan_inds)=x_sorted(upper_bound(non_nan_inds))/1e6;
lower_x(non_nan_inds)=x_sorted(lower_bound(non_nan_inds))/1e6;

% % Plot the bounds of the mask to check it is sensible
% figure
% imagesc((angle_sorted/1e6),x_sorted/1e6,counts_mat_sorted);
% hold on
% plot(angle_sorted(1:N_length)/1e6,upper_x,'r.','Linewidth',2,'MarkerSize',12)
% plot(angle_sorted(1:N_length)/1e6,lower_x,'r.','Linewidth',2,'MarkerSize',12)
% 
% xlabel('Angle/^\circ')
% ylabel('x/mm')

%Calculate length of the mask and plot
Length_mm_3=abs(upper_x-lower_x);

% figure
% plot(angle_sorted(1:N_length)/1e6,Length_mm_3,'.','Markersize',12);
% xlabel('Angle/^\circ')
% ylabel('Mask length /mm')
% set(gca,'Fontsize',12,'Linewidth',1)







%% 
% Now can fit to the length vector

[xData, yData] = prepareCurveData( angle_sorted(1:N_length)/1e6, Length_mm_3 );

% Set up fittype and options.
ft = fittype( 'h*(abs((alpha-sind(x))/(beta-cosd(x)))+(alpha-sind(x))/(beta-cosd(x)))/2+h2', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.Lower = [-Inf -Inf -Inf -Inf];
opts.StartPoint = [0 0 0.8 0.8];

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );

% % Plot fit with data.
% figure
% theta_plot=linspace(min(xData),max(xData),200);
% fit_plot=fitresult(theta_plot);
% p11= predint(fitresult,theta_plot,0.6827,'observation','off');
% h2= plot(theta_plot,fit_plot,'LineWidth',2,"Color",	[0.8500, 0.3250, 0.0980]);
% hold on
% %h3 = plot(theta_plot,p11,'--','LineWidth',1,"Color",	[0.8500, 0.3250, 0.0980]);
% h = plot(xData, yData,'.','MarkerSize',12,'Color', [0, 0.4470, 0.7410] );

% % Label axes
% xlabel( '\theta/^\circ');
% ylabel( "L/mm");
% set(gca,'Fontsize',12,'Linewidth',1)

theta_3 =angle_sorted(1:N_length)/1e6;
z_pos_3 =double(starting_position(3))/1e6;


























%% 
% Repeat for next dataset


%Load the data
load('data\AL000016.mat');

clip_max=100;


%% 
% Code to find steps using inbuilt ischange function

 
counts_mat_sorted(41,61) = nan;



% Find where the two steps in the data are
[TF,S1] = ischange(counts_mat_sorted,'linear',1,'MaxNumChanges',2);


N=size(counts_mat_sorted);
N_length=min(N(2),clip_max);

upper_bound=NaN*zeros(1,N_length);
lower_bound=NaN*zeros(1,N_length);
for n=1:N_length
    temp_inds=find(TF(:,n));
    if length(temp_inds)==2 && n<clip_max
        upper_bound(n)=temp_inds(1);
        lower_bound(n)=temp_inds(2);
    else
        upper_bound(n)=NaN;
        lower_bound(n)=NaN;
    end
end



%% 
% Determine the positions of the steps as x co-ordinate


upper_x=NaN*zeros(1,N_length);
lower_x=NaN*zeros(1,N_length);

non_nan_inds=find(~isnan(upper_bound));

upper_x(non_nan_inds)=x_sorted(upper_bound(non_nan_inds))/1e6;
lower_x(non_nan_inds)=x_sorted(lower_bound(non_nan_inds))/1e6;

% % Plot the bounds of the mask to check it is sensible
% figure
% imagesc((angle_sorted/1e6),x_sorted/1e6,counts_mat_sorted);
% hold on
% plot(angle_sorted(1:N_length)/1e6,upper_x,'r.','Linewidth',2,'MarkerSize',12)
% plot(angle_sorted(1:N_length)/1e6,lower_x,'r.','Linewidth',2,'MarkerSize',12)
% 
% xlabel('Angle/^\circ')
% ylabel('x/mm')



%Calculate length of the mask and plot
Length_mm_4=abs(upper_x-lower_x);

% figure
% plot(angle_sorted(1:N_length)/1e6,Length_mm_4,'.','Markersize',12);
% xlabel('Angle/^\circ')
% ylabel('Mask length /mm')
% set(gca,'Fontsize',12,'Linewidth',1)






%% 
% Now can fit to the length vector

[xData, yData] = prepareCurveData( angle_sorted(1:N_length)/1e6, Length_mm_4 );

% Set up fittype and options.
ft = fittype( 'h*(abs((alpha-sind(x))/(beta-cosd(x)))+(alpha-sind(x))/(beta-cosd(x)))/2+h2', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.Lower = [-Inf -Inf -Inf -Inf];
opts.StartPoint = [0 0 0.8 0.8];

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );

% % Plot fit with data.
% figure
% theta_plot=linspace(min(xData),max(xData),200);
% fit_plot=fitresult(theta_plot);
% p11= predint(fitresult,theta_plot,0.6827,'observation','off');
% h2= plot(theta_plot,fit_plot,'LineWidth',2,"Color",	[0.8500, 0.3250, 0.0980]);
% hold on
% %h3 = plot(theta_plot,p11,'--','LineWidth',1,"Color",	[0.8500, 0.3250, 0.0980]);
% h = plot(xData, yData,'.','MarkerSize',12,'Color', [0, 0.4470, 0.7410] );
% 
% % Label axes
% xlabel( '\theta/^\circ');
% ylabel( "L/mm");
% set(gca,'Fontsize',12,'Linewidth',1)

theta_4 =angle_sorted(1:N_length)/1e6;
z_pos_4 =double(starting_position(3))/1e6;



























%% 
% Repeat for next dataset


%Load the data
load('data\AL000017.mat');

clip_max=100;


%% 
% Code to find steps using inbuilt ischange function

 
counts_mat_sorted(87,19) = nan;



% Find where the two steps in the data are
[TF,S1] = ischange(counts_mat_sorted,'linear',1,'MaxNumChanges',2);


N=size(counts_mat_sorted);
N_length=min(N(2),clip_max);

upper_bound=NaN*zeros(1,N_length);
lower_bound=NaN*zeros(1,N_length);
for n=1:N_length
    temp_inds=find(TF(:,n));
    if length(temp_inds)==2 && n<clip_max
        upper_bound(n)=temp_inds(1);
        lower_bound(n)=temp_inds(2);
    else
        upper_bound(n)=NaN;
        lower_bound(n)=NaN;
    end
end


%% 
% Determine the positions of the steps as x co-ordinate


upper_x=NaN*zeros(1,N_length);
lower_x=NaN*zeros(1,N_length);

non_nan_inds=find(~isnan(upper_bound));

upper_x(non_nan_inds)=x_sorted(upper_bound(non_nan_inds))/1e6;
lower_x(non_nan_inds)=x_sorted(lower_bound(non_nan_inds))/1e6;

% % Plot the bounds of the mask to check it is sensible
% figure
% imagesc((angle_sorted/1e6),x_sorted/1e6,counts_mat_sorted);
% hold on
% plot(angle_sorted(1:N_length)/1e6,upper_x,'r.','Linewidth',2,'MarkerSize',12)
% plot(angle_sorted(1:N_length)/1e6,lower_x,'r.','Linewidth',2,'MarkerSize',12)
% 
% xlabel('Angle/^\circ')
% ylabel('x/mm')

%Calculate length of the mask and plot
Length_mm_5=abs(upper_x-lower_x);

% figure
% plot(angle_sorted(1:N_length)/1e6,Length_mm_5,'.','Markersize',12);
% xlabel('Angle/^\circ')
% ylabel('Mask length /mm')
% set(gca,'Fontsize',12,'Linewidth',1)







%% 
% Now can fit to the length vector

[xData, yData] = prepareCurveData( angle_sorted(1:N_length)/1e6, Length_mm_5 );

% Set up fittype and options.
ft = fittype( 'h*(abs((alpha-sind(x))/(beta-cosd(x)))+(alpha-sind(x))/(beta-cosd(x)))/2+h2', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.Lower = [-Inf -Inf -Inf -Inf];
opts.StartPoint = [0 0 0.8 0.8];

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );

% % Plot fit with data.
% figure
% theta_plot=linspace(min(xData),max(xData),200);
% fit_plot=fitresult(theta_plot);
% p11= predint(fitresult,theta_plot,0.6827,'observation','off');
% h2= plot(theta_plot,fit_plot,'LineWidth',2,"Color",	[0.8500, 0.3250, 0.0980]);
% hold on
% %h3 = plot(theta_plot,p11,'--','LineWidth',1,"Color",	[0.8500, 0.3250, 0.0980]);
% h = plot(xData, yData,'.','MarkerSize',12,'Color', [0, 0.4470, 0.7410] );
% 
% % Label axes
% xlabel( '\theta/^\circ');
% ylabel( "L/mm");
% set(gca,'Fontsize',12,'Linewidth',1)

theta_5 =angle_sorted(1:N_length)/1e6;
z_pos_5 =double(starting_position(3))/1e6;
















%% 
% Now will repeat using nonlinear fit explicitly

sep=100;
x_shift = 0.8;

theta_comb=[theta_1,theta_2+sep,theta_3+2*sep,theta_4+3*sep,theta_5+4*sep];
% theta_comb=[theta_1,theta_2+sep,theta_3+2*sep,theta_4+3*sep];
%theta_comb=[theta_4+3*sep,theta_5+4*sep];


length_comb=[Length_mm_1,Length_mm_2,Length_mm_3,Length_mm_4,Length_mm_5];
% length_comb=[Length_mm_1,Length_mm_2,Length_mm_3,Length_mm_4];
%length_comb=[Length_mm_4,Length_mm_5];



% Setting up model using 4 seperate data sets with 100 pont steps between each set.
%  b(1)= 1mm h         b(2)=h2        b(3)= predicted x_cor         b(4)= predicted z_cor       
%  b(5)=r         b(6)= 2mm h
modelfun = @(b,x) (x<=max(theta_1)).*(b(1)*(abs(((b(3)+z_pos_1)/b(5)-sind(x))./((b(4)-z_pos_1)/b(5)-cosd(x)))+((b(3)+z_pos_1)/b(5)-sind(x))./((b(4)-z_pos_1)/b(5)-cosd(x)))/2+b(2))...
    + (x>=sep & x<=sep+max(theta_2)).*(b(1)*(abs(((b(3)+z_pos_2)/b(5)-sind(x-sep))./((b(4)-z_pos_2)/b(5)-cosd(x-sep)))+((b(3)+z_pos_2)/b(5)-sind(x-sep))./((b(4)-z_pos_2)/b(5)-cosd(x-sep)))/2+b(2))...
    + (x>=2*sep & x<=2*sep+max(theta_3)).*(b(1)*(abs(((b(3)+z_pos_3)/b(5)-sind(x-2*sep))./((b(4)-z_pos_3)/b(5)-cosd(x-2*sep)))+((b(3)+z_pos_3)/b(5)-sind(x-2*sep))./((b(4)-z_pos_3)/b(5)-cosd(x-2*sep)))/2+b(2))...
    + (x>=3*sep & x<=3*sep+max(theta_4)).*(b(1)*(abs(((b(3)+z_pos_4)/b(5)-sind(x-3*sep))./((b(4)-z_pos_4)/b(5)-cosd(x-3*sep)))+((b(3)+z_pos_4)/b(5)-sind(x-3*sep))./((b(4)-z_pos_4)/b(5)-cosd(x-3*sep)))/2+b(2))...     
    + (x>=4*sep & x<=4*sep+max(theta_5)).*(b(1)*(abs(((b(3)+z_pos_5+x_shift)/b(5)-sind(x-4*sep))./((b(4)-z_pos_5)/b(5)-cosd(x-4*sep)))+((b(3)+z_pos_5+x_shift)/b(5)-sind(x-4*sep))./((b(4)-z_pos_5)/b(5)-cosd(x-4*sep)))/2+b(2));





%Starting point for parameters
%beta0 = [1 1 -z_pos_1 z_pos_1 6, 2];    
beta0 = [1 1 -z_pos_1 z_pos_1 6.5];   


% Fits a nonlinear regression model using the column vector (length_comb) as a response variable and the columns of the matrix (theta_comb) as predictor variables
mdl_value = fitnlm(theta_comb,length_comb,modelfun,beta0);

% 1 standard deviation (68.27%) for confidence interval
level=0.6827;


theta_comb_plot=linspace(min(theta_comb),max(theta_comb),1000);
%predict_weight_value=(gof_value.sse/gof_value.dfe)./std_fit_value_plot.^2; %Predict uses a different defintion of weight


% Height of step and confidence interval
[ypred_value,ci_value] = predict(mdl_value,theta_comb_plot','Alpha',(1-level),'Simultaneous',false,'Prediction','observation');%,'Weights',predict_weight_value);



% figure;plot(theta_comb,length_comb,'.','MarkerSize',12)
% hold on
% plot(theta_comb_plot,ypred_value,'LineWidth',1)
% %plot(conc_plot,y_ci_3,'--','LineWidth',1,'Color', [0.8500, 0.3250, 0.0980])
% plot(theta_comb_plot,ci_value,'--','LineWidth',1,'Color', [0.8500, 0.3250, 0.0980])
% xlabel( '\theta/^\circ');
% ylabel( "L/mm");
% set(gca,'Fontsize',12,'Linewidth',1)



%% 
% Get some parameters and errors

error_params=diff(coefCI(mdl_value,level),1,2)/2;


%% 
% 
%% 
% Generate some nice looking plots


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


% First start by shoving everything back on top of each other

f_h=figure('Color','white','Units', 'centimeters','Position',[1 1 25 20],'Resize', 'off');
axis tickaligned
box on
hold on


theta_3_plot=linspace(min(theta_3),max(theta_3),100)+2*sep;


% Height of step and confidence interval
[ypred_3,ci_3] = predict(mdl_value,theta_3_plot','Alpha',(1-level),'Simultaneous',false,'Prediction','observation');


theta_4_plot=linspace(min(theta_4),max(theta_4),100)+3*sep;


% Height of step and confidence interval
[ypred_4,ci_4] = predict(mdl_value,theta_4_plot','Alpha',(1-level),'Simultaneous',false,'Prediction','observation');


theta_5_plot=linspace(min(theta_5),max(theta_5),100)+4*sep;


% Height of step and confidence interval
[ypred_5,ci_5] = predict(mdl_value,theta_5_plot','Alpha',(1-level),'Simultaneous',false,'Prediction','observation');

xlim([0 60]);
xticks([0 15 30 45 60]);
ylim([1 3]);
yticks([1.0 1.5 2.0 2.5 3.0]);
ytickformat('%.1f');

xlabel('Angle / ^\circ','fontname',plot_font,'fontsize',axes_font_size);
ylabel('Mask length / mm','fontname',plot_font,'fontsize',axes_font_size);
set(gca,'FontSize',30,'LineWidth',2);


%Plot tan for reference
vars=mdl_value.Coefficients.Variables;
h1=vars(1,1);
h2=vars(2,1);

theta_tan_plot=linspace(0,60,100);


plot(theta_3,Length_mm_3,'Marker','diamond', 'MarkerEdgeColor', plotblue,'LineStyle','none', 'MarkerFaceColor', 'white','MarkerSize',8,'LineWidth',2.5)
plot(theta_4,Length_mm_4,'Marker','square', 'MarkerEdgeColor', plotblue,'LineStyle','none', 'MarkerFaceColor', 'white','MarkerSize',8,'LineWidth',2.5)
plot(theta_5,Length_mm_5,'Marker','o', 'MarkerEdgeColor', plotblue,'LineStyle','none', 'MarkerFaceColor', 'white','MarkerSize',8,'LineWidth',2.5)
plot(theta_3_plot-2*sep,ypred_3,'LineWidth', 3,'Color', plotorange)
plot(theta_4_plot-3*sep,ypred_4,'LineWidth', 3,'Color', plotorange)
plot(theta_5_plot-4*sep,ypred_5,'LineWidth', 3,'Color', plotorange)
plot(theta_tan_plot,h1*tand(theta_tan_plot)+h2, 'Color', plotgold, 'LineWidth', 3)

legend('Z=-2.5mm; X=0.0mm','Z=-3.8mm; X=0.0mm','Z= 3.8mm; X=0.8mm',' Model Fits',' Model Fits',' Model Fits',' tan(θ)', 'Location','northwest','fontname',plot_font,'fontsize',25)
legend('boxoff')

% exportgraphics(gcf,'mask_length_data.pdf','ContentType','vector')