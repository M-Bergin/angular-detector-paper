% Script to analyse monitor scans showing peak decay from contamination
% Bergin 23/2/24

%% Setup memory and input files

T_overall=[180,160,140,100,60,80,120];
file_ind=[6,9,7,11,13,15,17];

N_file=length(T_overall);

b_vec=NaN*zeros(N_file,1);
b_vec_error=NaN*zeros(N_file,1);

% Set the amount of data to ignore at the start of the scan
L_ignore=25;

file_prefix='data/Monitor scans/';


torr_mbar = 1.3332236842;

%% Loop over each file
for n_file=1:N_file

    % Import the file
    file_suffix=sprintf('Mn%06i',file_ind(n_file));
    load([file_prefix,file_suffix,'.mat'])

    % Find time and exposure vectors by integration
    t_vec=seconds(time_vec-time_vec(1));
    L=cumtrapz(t_vec,p_sample_vec)/(torr_mbar*1e-6);

    % figure;plot(t_vec,L)
    % figure;plot(L,count_vec,'.')

    % Fit the data with an exponential
    [xData, yData] = prepareCurveData( L, count_vec );
    ft = fittype( 'a*exp(-x/b)+c', 'independent', 'x', 'dependent', 'y' );
    
    % Set the exluded data
    if T_overall(n_file)==100
        % Ignore extra data at the end of the 100C dataset
        excludedByRule = xData < L_ignore | xData > 80;
    else
        excludedByRule = xData < L_ignore;
    end
    
    opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
    opts.Display = 'Off';
    opts.Lower = [0 0 400];
    opts.StartPoint = [20000 100 10000];
    opts.Exclude = excludedByRule;

    % Fit model to data.
    [fitresult, ~] = fit( xData, yData, ft, opts );

    % % Plot fit with data.
    % figure( 'Name', 'untitled fit 1' );
    % plot( fitresult, xData, yData, excludedByRule );
    % % Label axes
    % xlabel( 'L' );
    % ylabel( 'Counts' );

    % Find the decay rate
    errors= confint(fitresult,0.6827);
    b=fitresult.b;
    b_error=diff(errors(:,2))/2;

    % Save data into vector for later use
    b_vec(n_file)=b;
    b_vec_error(n_file)=b_error;

end

%% Overall plots

%Font Choice and Text Size
plot_font = 'Arial';
tick_font_size = 30;
axes_font_size = 35;
width_line=4;

fig = figure('Color','white','Units', 'centimeters','Position',[1 1 35 20],'Resize', 'off');
% errorbar(T_overall,b_vec,b_vec_error,'.','MarkerSize',12,'LineWidth',1);
plot(T_overall,b_vec,'.','MarkerSize',36,'LineWidth',1);
yline(L_ignore/2,'k--','LineWidth',width_line)

ax1=gca;
ax1.YScale='log';


xlabel('T/^\circC')
ylabel('Characteristic exposure/L')

%axis equal
axis tickaligned
box on
ax = gca;
ax.XAxis.FontSize = tick_font_size;
ax.YAxis.FontSize = tick_font_size;
ax.LineWidth = 2;

xlim([50 190])

% set(ax1,'FontSize',12,'LineWidth',1)

% exportgraphics(gcf,'../L_temp.eps','ContentType','vector')


