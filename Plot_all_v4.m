

%Colours for Plots
plotblue = [0 0.4470 0.7410];
plotgold = [0.9290 0.6940 0.1250];
plotorange = [0.8500 0.3250 0.0980];
plotgreen = [28/255 120/255 37/255];
plotred = [0.6350 0.0780 0.1840];
plotpurple = [0.4940, 0.1840, 0.5560];

lightBLUE = [0.356862745098039,0.811764705882353,0.956862745098039];
darkBLUE = [0.0196078431372549,0.0745098039215686,0.670588235294118];

% lightBLUE = [1,0.9294117647058824,0.6274509803921569];
% darkBLUE = [0.9411764705882353,0.23137254901960785,0.12549019607843137];
 
blueGRADIENTflexible = @(i,N) lightBLUE + (darkBLUE-lightBLUE)*((i-1)/(N-1));

%Font Choice and Text Size
plot_font = 'Arial';
tick_font_size = 30;
axes_font_size = 35;
width_line=4;


% %Load the gratings
% load('An001886.mat')
% angle_vec_grating=angle_vec;
% counts_grating=counts;

% load('An004664.mat')


% load('C:\Users\mberg\The University of Newcastle\Centre for Organic Electronics - SHeM\Data\5. Win10 Data\Scan Data\Angular scans\An004204.mat')

fig = figure('Color','white','Units', 'centimeters','Position',[1 1 35 20],'Resize', 'off');
axis equal
axis tickaligned
box on





%% Plot 3D resin

figure3a_data = load('data\Fig3a_Data_V2.mat');
figure3b_data = load('data\Fig3b_Data.mat');


%figure3b_data.FAHY(4,:)/1000;
resin_data_norm=(figure3b_data.FAHY(4,:)-min(figure3b_data.FAHY(4,:)))/(max(figure3b_data.FAHY(4,:))-min(figure3b_data.FAHY(4,:)));


% resin_data_norm=(figure3a_data.Fig3a_Data(2,:)-min(figure3a_data.Fig3a_Data(2,:)))/(max(figure3a_data.Fig3a_Data(2,:))-min(figure3a_data.Fig3a_Data(2,:)));
plot(figure3b_data.FAHY(1,:),movmean(resin_data_norm,2),'Color', plotpurple, 'LineWidth', width_line)
% plot(figure3b_data.FAHY(1,:),figure3b_data.FAHY(4,:)/1000,'-', 'Color', plotred, 'LineWidth', 4)
hold on

%% Plot gratings

theta_plot=linspace(-90,90,500);
plot(theta_plot,cosd(theta_plot),'-','LineWidth',width_line,'Color',plotgold)


%% Plot HOPG

% Simulated HOPG scattering

m_He=4*1.67e-27;
M=72*1.67e-27;
theta_i=45;

k_B=1.380649e-23;
T_s=400;%273+235;
T_source=300;

E_i=(5/2)*k_B*T_source;
u_i_avg=sqrt(2*E_i/(m_He));
u_i_norm=u_i_avg*cosd(theta_i);


mu=m_He/M;

theta_f_vec=linspace(-50,89,1000);
N_theta_f=length(theta_f_vec);

I=NaN*zeros(N_theta_f,1);
B_1_vec=NaN*zeros(N_theta_f,1);
B_2_vec=NaN*zeros(N_theta_f,1);

G= @(x) sqrt(M/(2*pi*k_B*T_s))*exp(-(M/(2*k_B*T_s)).*x.^2);

for n=1:N_theta_f

    theta_f=theta_f_vec(n);
    B_1=((1+mu)/2)*sind(theta_i)*cotd(theta_f)-((1-mu)/2)*cosd(theta_i);
    B_2=((1+mu)/2)*sind(theta_i)*cscd(theta_f).^2;

    B_1_vec(n)=B_1;
    B_2_vec(n)=B_2;

    fun = @(u_i) (cosd(theta_i)+B_1).*B_2.*u_i.^2.*G(B_1*u_i).*normpdf(u_i,u_i_avg,u_i_avg/100);

    I(n)=(cosd(theta_i)+B_1).*B_2.*u_i_avg.^2.*G(B_1*u_i_avg);%integral(fun,0,Inf)./u_i_norm;
    
end

I_tot=sum(I)*(theta_f_vec(2)-theta_f_vec(1));

%%%%% HOPG plotting here %%%%%%%

% figure;
plot(theta_f_vec,I./max(I), 'LineWidth', width_line,'Color',plotorange)




%% Plot LiF scan


% Function to approximate the diffraction pattern from a LiF crystal and
% what signal it would produce in a SHeM

%%%%%%%%%%%% Setup of parameters %%%%%%%%%%%

%Parameters
lambda=6.63e-34/(sqrt(5*4*1.67e-27*1.38e-23*(273+25)));
theta_in=45;

alph=0; %Flux in diffuse component

%Detector aperture dimensions
x_ap_mid=3.5e-3*2;%sqrt(2)*3e-3;
y_ap_mid=0;
r_ap=0.25e-3;

%Offset of aperture from pinhole plate
z_offset=1.5e-3;
%z positions to calculate
z_min=1.5e-3; %z in m
z_max=12e-3;
N_z=200; %Number of points, increase to 1000 for better quality







% Full loading of Boyao data

str=fileread('data\diffrac10001.out');

% Get the start of each row
tkn=regexp(str,'Required number of z grid points');

N_phi=length(tkn);

phi_Boyao_vec=NaN*zeros(N_phi,1);

for n_phi=1:N_phi

    if n_phi==N_phi
        sub_str=str(tkn(n_phi):end);
    else
        sub_str=str(tkn(n_phi):tkn(n_phi+1)-1);
    end

    tkn_n=regexp(sub_str,'n =');
    tkn_n_end=regexp(sub_str(tkn_n:end),'[\n]','once')+tkn_n-2;
    N_rows=str2double(sub_str(tkn_n+3:tkn_n_end));

    tkn_theta=regexp(sub_str,'theta =');
    tkn_theta_end=regexp(sub_str(tkn_theta:end),'[\n]','once')+tkn_theta-2;
    temp=strsplit(sub_str(tkn_theta+6:tkn_theta_end),' ');
    phi_Boyao_vec(n_phi)=str2double(temp{3});


    varnames{n_phi}=matlab.lang.makeValidName(strcat('phi',num2str(phi_Boyao_vec(n_phi))));



    startRow = 7;
    endRow = 8 + N_rows+1;

    % Format for each line of text:
    %   column1: categorical (%C)
    %	column2: double (%f)
    %   column3: double (%f)
    %	column4: double (%f)
    % For more information, see the TEXTSCAN documentation.
    formatSpec = '%1C%7f%6f%f%[^\n\r]';

    % Read columns of data according to the format.
    % This call is based on the structure of the file used to generate this code. If an error occurs for a different file, try regenerating the code from the Import Tool.
    dataArray = textscan(sub_str, formatSpec, endRow-startRow+1, 'Delimiter', '', 'WhiteSpace', '', 'TextType', 'string', 'HeaderLines', startRow-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');

    % Post processing for unimportable data.
    % No unimportable data rules were applied during the import, so no post processing code is included. To generate code which works for unimportable data, select unimportable cells in a file and regenerate the script.

    % Allocate imported array to column variable names
    VarName1 = dataArray{:, 1};
    b1_full.(varnames{n_phi}) = dataArray{:, 2};
    b2_full.(varnames{n_phi}) = dataArray{:, 3};
    I_full.(varnames{n_phi}) = dataArray{:, 4};


    % Clear temporary variables
    clearvars filename startRow endRow formatSpec fileID dataArray ans;

end




% for n_phi=1:N_phi
n_phi=1;
    phi_in=phi_Boyao_vec(n_phi);%45-59;

%%%%%%%%%% Creation of scattering distribution in k space %%%%%%%%%%

%Calculate positions of the diffraction peaks
[k_out,G,theta_out,phi_out,N_eff,N_x,N_y]=diffraction_peak_locations(theta_in,-phi_in,lambda);



%Generate k space grid to plot diffraction pattern on
k_x_max=1.2e11;
k_x_N=1000;

k_y_max=1.2e11;
k_y_N=1000;

k_x_vec=linspace(-k_x_max,k_x_max,k_x_N);
k_y_vec=linspace(-k_y_max,k_y_max,k_y_N);

[k_X,k_Y]=meshgrid(k_x_vec,k_y_vec);

%Set k_Z by energy conservation
k_mag=(2*pi)/lambda;
k_Z=sqrt(k_mag^2-(k_X.^2+k_Y.^2));
imag_inds=imag(k_Z)>0;


% Import Boyao data

% Find the data from library

ind_B=find(phi_Boyao_vec<phi_in+0.5 & phi_Boyao_vec>phi_in-0.5,1);
b1 = b1_full.(varnames{ind_B});
b2 = b2_full.(varnames{ind_B});
I = I_full.(varnames{ind_B});


% Create the diffraction pattern
I_k=zeros(k_x_N,k_y_N);

N_points=size(k_out,1);

%Set the width of the peaks and how they decay
peak_width=3e9;
pattern_width=1.5;

%Main loop to add in the diffraction pattern
for n_point=1:N_points
    %Get intensity
    Boyao_ind= find(b1==N_x(n_point) & b2==N_y(n_point));
    peak_I=I(Boyao_ind);

    if ~isempty(Boyao_ind)
        I_x_temp=normpdf(k_x_vec,k_out(n_point,1),peak_width);
        I_y_temp=normpdf(k_y_vec,k_out(n_point,2),peak_width);

        I_k=I_k+(peak_I*(I_y_temp'*I_x_temp));
    end
end

%Normalise the diffraction pattern contribution
I_k=(I_k/(sum(sum(I_k))))*(1-alph);


%Add in diffuse component
theta_mat=(atand(sqrt(k_X.^2+k_Y.^2)./k_Z));
theta_mat(theta_mat~=real(theta_mat))=NaN;

%Distribution for diffuse scattering
I_diff=~isnan(theta_mat);
%Distribution for the solid angle
%I_diff=1./real(cosd(theta_mat));

%Normalise the diffuse contribution
I_diff=(I_diff./(nansum(nansum(I_diff))))*alph;


%Calculate total diffraction pattern with diffuse component.
I_tot=I_k+I_diff;


%Remove the imaginary parts
I_tot(imag_inds)=0;




%%%%%% LiF plotting here %%%%%%%

plot(asind(k_x_vec/k_mag),cosd(asind(k_x_vec/k_mag)).*I_tot(500,:)./max(I_tot(500,:)),'Color', plotblue, 'LineWidth', width_line)



xlabel('{\theta / ^\circ}','fontname',plot_font,'fontsize',axes_font_size)
ylabel('Normalised counts','fontname',plot_font,'fontsize',axes_font_size)
ax = gca;
ax.XAxis.FontSize = tick_font_size;
ax.YAxis.FontSize = tick_font_size;
ax.LineWidth = 2;
legend('Rough diffuse', 'Cosine','Inelastic','Elastic', 'Location','northwest','fontname',plot_font,'fontsize',tick_font_size)
legend('boxoff')
xlim([-90 90])
ylim([0 1.1])


% exportgraphics(fig,'AngularPaper_comparison_fig_v4.pdf','ContentType','vector')