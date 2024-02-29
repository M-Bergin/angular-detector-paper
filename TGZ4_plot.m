% Load z data

files_ind_1 = 2084:2134;
files_ind_2 = 2137:2176;
files_ind = [files_ind_1,files_ind_2];

path1 = 'data\TGZ4 scans';
[data, thetas, phis] = load_an_scans(files_ind, path1);




files_ind_3 = 1986:2030;
files_ind_4 = 2036:2081;
files_ind_5 = [files_ind_3,files_ind_4];

[data_2, thetas_2, phis_2] = load_an_scans(files_ind_5, path1);




% phis=phis;
phis_3 = 11:4:187;
phis_4 = 189:-4:9;
phis_2 = [phis_3,phis_4];



[theta_sorted,theta_ind]=sort(thetas);
[phi_sorted,phi_ind]=sort(phis);

[theta_sorted_2,theta_ind_2]=sort(thetas_2);
[phi_sorted_2,phi_ind_2]=sort(phis_2);

phi_sorted_2 = 9:2:189;

%phi_sorted = phi_sorted-99;

%Remove spikes
data=spike_im_removal(data,4);
data_2=spike_im_removal(data_2,4);


% figure;
% imagesc(phi_sorted,theta_sorted,data(theta_ind,phi_ind))
% title('2nd set')
% xlabel('\phi/^\circ')
% ylabel('\theta/^\circ')
% set(gca,'LineWidth',1,'FontSize',12)
% 
% 
% 
% figure;
% imagesc(phi_sorted_2,theta_sorted_2,data_2(theta_ind_2,phi_ind_2))
% title('1st set')
% xlabel('\phi/^\circ')
% ylabel('\theta/^\circ')
% set(gca,'LineWidth',1,'FontSize',12)




data_combined = data+data_2;

figure;
imagesc(phi_sorted_2,theta_sorted_2,data_combined(theta_ind_2,phi_ind))
%title('Combined')
xlabel('\phi/^\circ')
ylabel('\theta/^\circ')
set(gca,'LineWidth',1,'FontSize',12)



%FIG2 = struct('phi',phi_sorted_2,'theta',theta_sorted_2,'data',data_combined,'theta_ind',theta_ind_2,'phi_ind',phi_ind)



% 
% 
% 
% function theta = theta_of_z(z)
%     theta = atand((7 - (z+1.5))./(z+1.5));
% end