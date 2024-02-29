% Script to generate HOPG plots

%% Plot the SHeM data

% Import the data
SHeM_data_hot=read_MkII_data('data\19-Aug-2022_001.txt');
SHeM_data_cold=read_MkII_data('data\29-Jul-2022_001.txt');

% Rotate the data to correct orientation
SHeM_data_cold.image=rot90(flipud(SHeM_data_cold.image),3);
SHeM_data_hot.image=rot90(flipud(SHeM_data_hot.image),3);

% Plot the data
c_cold=SHeM_fig_process_Newcastle_colorbar(SHeM_data_cold,0,2000,0,0,1);
c_cold.Ticks=[150,550,950,1350];
% exportgraphics(gcf,'HOPG_cold.pdf','ContentType','vector','Resolution',1000)

c_hot=SHeM_fig_process_Newcastle_colorbar(SHeM_data_hot,0,2000,0,0,1);
c_hot.Ticks=[2000,4000,6000,8000];
% exportgraphics(gcf,'HOPG_hot.pdf','ContentType','vector','Resolution',1000)

%Import final image
SHeM_data_zoom=read_MkII_data('Data\17-Aug-2022_001.txt');
% Rotate the data to correct orientation
SHeM_data_zoom.image=rot90(flipud(SHeM_data_zoom.image),3);
% Plot the data
c_zoom=SHeM_fig_process_Newcastle_colorbar(SHeM_data_zoom,0,100,0,0,0);
% exportgraphics(gcf,'HOPG_zoom.pdf','ContentType','vector','Resolution',1000)
SHeM_img_norm=255*((SHeM_data_zoom.image-min(SHeM_data_zoom.image(:)))/(max(SHeM_data_zoom.image(:)-min(SHeM_data_zoom.image(:)))));
% imwrite(uint8(imresize(SHeM_img_norm,10,'nearest')),'HOPG_zoom_png.png')


%% Plot the raw profilometer data

opts = delimitedTextImportOptions("NumVariables", 2);

% Specify range and delimiter
opts.DataLines = [4, Inf];
opts.Delimiter = " ";

% Specify column names and types
opts.VariableNames = ["y", "z"];
opts.VariableTypes = ["double", "double"];

% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";
opts.ConsecutiveDelimitersRule = "join";
opts.LeadingDelimitersRule = "ignore";




start_ind=585;
end_ind=735;

n_vec=start_ind:end_ind;
N=length(n_vec);

% Import the data
LSData = readtable("data\2nd set\Test585.txt", opts);

N_points=size(LSData,1);
z_mat=NaN*ones(N,N_points);

for n=start_ind:end_ind


    opts = delimitedTextImportOptions("NumVariables", 2);

    % Specify range and delimiter
    opts.DataLines = [4, Inf];
    opts.Delimiter = " ";


    % Specify column names and types
    opts.VariableNames = ["y", "z"];
    opts.VariableTypes = ["double", "double"];

    % Specify file level properties
    opts.ExtraColumnsRule = "ignore";
    opts.EmptyLineRule = "read";
    opts.ConsecutiveDelimitersRule = "join";
    opts.LeadingDelimitersRule = "ignore";

    % Import the data
    LSData = readtable(['data\2nd set\Test',num2str(n),'.txt'], opts);

    clear opts
    n_run=size(LSData,1);

    z_mat(n-start_ind+1,1:n_run)=LSData.z;

end
% 
% figure;
% for n=start_ind:end_ind
%     plot(z_mat(n-start_ind+1,:))
%     pause(0.5)
% end

%figure;imagesc(z_mat); colormap gray

z_norm=255*((z_mat-min(z_mat(:)))/(max(z_mat(:)-min(z_mat(:)))));
% imwrite(uint8(imresize(z_norm,10,'nearest')),'profilometer_test.png')

%% Perform registration of the images to create profilometer image

% Track points manually from the two images

optical_points=[1764, 441;
                1669, 472;
                1643, 590;
                1406, 652;
                1329, 663;
                919, 954;
                1026, 865;
                978, 342;
                1770, 701;
                1698, 829;
                2028, 926;
                2065, 845;
                1321, 1099;
                1981, 1311];

profile_points=[4304, 9;
                3888, 16;
                3794, 34;
                2774, 49;
                2442, 53;
                862, 106;
                1257, 90;
                714, 15;
                4447, 48;
                4219, 68;
                5733, 73;
                5847, 60;
                2791, 117;
                5908, 132];
% 

SHeM_points = [57,8;
                48, 20;
                61, 28;
                57, 41;
                21,31];

SHeM_points2 = [565,72;
                472, 190;
                615, 277;
                567, 400;
                211,302];

profile_points2=[4304, 9;
                3794, 34;
                4447, 48;
                4219, 68;
                2442, 53];

% Import the file
rawData1 = importdata('data\HOPG_004_001_tif.png');

% For some simple files (such as a CSV or JPEG files), IMPORTDATA might
% return a simple array.  If so, generate a structure so that the output
% matches that from the Import Wizard.
[~,name] = fileparts('data\HOPG_004_001_tif.png');
newData1.(matlab.lang.makeValidName(name)) = rawData1;

% Create new variables in the base workspace from those fields.
vars = fieldnames(newData1);
for i = 1:length(vars)
    assignin('base', vars{i}, newData1.(vars{i}));
end

% figure;imagesc(HOPG_004_001_tif)
% hold on
% plot(optical_points(:,1),optical_points(:,2),'x')
% 
% 
% figure;imagesc(z_norm); colormap gray
% hold on
% plot(profile_points(:,1),profile_points(:,2),'x')


% Register the optical and profile images using the control points
tform = fitgeotform2d(profile_points,optical_points,"affine");

% Register the SHeM and profile images using the control points
tform2 = fitgeotform2d(profile_points2,SHeM_points,"affine");
% Register the SHeM and profile images using the control points
tform3 = fitgeotform2d(profile_points2,SHeM_points2,"affine");

Opt_bottom_lvl=5e4;
Opt_g=rgb2gray(HOPG_004_001_tif);
Opt_norm=uint8(255*(double(Opt_g-Opt_bottom_lvl)/double(max(Opt_g(:))-Opt_bottom_lvl)));

Jregistered = imwarp(z_norm,tform,OutputView=imref2d(size(HOPG_004_001_tif)));
J_norm=uint8(255*(double(Jregistered-min(Jregistered(:)))/double(max(Jregistered(:))-min(Jregistered(:)))));

Jregistered2 = imwarp(z_norm,tform2,OutputView=imref2d(size(SHeM_img_norm)));
J_norm2=uint8(255*(double(Jregistered2-min(Jregistered2(:)))/double(max(Jregistered2(:))-min(Jregistered2(:)))));

Jregistered3 = imwarp(z_norm,tform3,OutputView=imref2d(size(imresize(SHeM_img_norm,10))),interp='nearest');
J_norm3=uint8(255*(double(Jregistered3-min(Jregistered3(:)))/double(max(Jregistered3(:))-min(Jregistered3(:)))));


% figure;imshowpair(Opt_norm,J_norm,"Scaling","none")

%figure;imshowpair(SHeM_img_norm,J_norm2,"Scaling","independent")

%figure;imshowpair(imresize(SHeM_img_norm,10,'nearest'),J_norm3,"Scaling","independent")

% Calculate co-ordinates of the corners
dektak_size=size(z_norm);
top_left=tform.A*[1;1;1];
top_right=tform.A*[dektak_size(2);1;1];
bottom_left=tform.A*[1;dektak_size(1);1];
bottom_right=tform.A*[dektak_size(2);dektak_size(1);1];

horizontal_length_px=sqrt((top_left(1)-top_right(1))^2+(top_left(2)-top_right(2))^2);
vertical_length_px=sqrt((top_left(1)-bottom_left(1))^2+(top_left(2)-bottom_left(2))^2);

%Scale bar in image is 229 px long and that corresponds to 200um.
horizontal_length_um=horizontal_length_px *(200/229);
vertical_length_um=vertical_length_px *(200/229);

%Horziontal length should be 1200um and we are getting 1163um which is
%fairly close.

vertical_step_um=vertical_length_um/(dektak_size(1)-1);

%Expected that each pixel would be 7.38um each, but getting 5.9um which is
%further away but at least reasonable.


%Scale the image to correct aspect ratio
vertical_scaled_height=vertical_length_px/horizontal_length_px*dektak_size(2);

scaled_dektak=imresize(z_norm,[vertical_scaled_height,dektak_size(2)],'nearest');

scaled_croped_dektak = scaled_dektak(1:3750,1250:5250);


figure;imagesc(scaled_croped_dektak);axis equal tight off; colormap gray
hold on
Scalebar_length = 6001/12;

quiver(3200,3450,Scalebar_length,0,'ShowArrowHead','off','AutoScale','off','LineWidth',5,'Color','y');

% exportgraphics(gcf,'Figures/Dektak_v2.eps')
% imwrite(uint8(scaled_croped_dektak),'Figures/Dektak_v2_plain.png')


% %% Make a gif varying the transparency of the images
% C = imfuse(double(rgb2gray(HOPG_004_001_tif))/65518,Jregistered/255,'scaling','none');
% imwrite(C,'Figures/Dektak_reg.png')
% 
% im_size=size(C);
% figure_size=700;
% delay=0.15;
% figure_size2=round(figure_size*im_size(2)/im_size(1));
% fraction=linspace(1,0,10);
% fraction=[fraction,flip(fraction)];
% N=length(fraction);
% h = figure;
% h.Position=[100 100 figure_size2 figure_size];
% filename = 'Figures/Dektak_reg.gif';
% % axis equal tight manual % this ensures that getframe() returns a consistent size
% for n=1:N
%     C_temp=C;
%     % C_temp(:,:,1)=C(:,:,1)/2;
%     % C_temp(:,:,3)=C(:,:,3)/2;
%     C_temp(:,:,2)=uint8(double(C(:,:,2))*fraction(n));
%     image(C_temp);axis equal tight off
%     set(gca,'position',[0 0 1 1],'units','normalized')
%       % Capture the plot as an image 
%       frame = getframe(h); 
%       im = frame2im(frame); 
%       [imind,cm] = rgb2ind(im,256); 
%       % Write to the GIF File 
%       if n == 1 
%           imwrite(imind,cm,filename,'gif', 'DelayTime', delay, 'Loopcount',inf); 
%       else 
%           imwrite(imind,cm,filename,'gif', 'DelayTime', delay,'WriteMode','append'); 
%       end 
% end
% 
% for n=1:N
%     C_temp=C;
%     C_temp(:,:,1)=uint8(double(C(:,:,1))*fraction(n));
%     C_temp(:,:,3)=uint8(double(C(:,:,3))*fraction(n));
% %     C_temp(:,:,2)=uint8(double(C(:,:,2))*fraction(n));
%     image(C_temp);axis equal tight off
%     set(gca,'position',[0 0 1 1],'units','normalized')
%       % Capture the plot as an image 
%       frame = getframe(h); 
%       im = frame2im(frame); 
%       [imind,cm] = rgb2ind(im,256); 
%       % Write to the GIF File 
%       imwrite(imind,cm,filename,'gif', 'DelayTime', delay,'WriteMode','append'); 
% end
% 
% 
% 


