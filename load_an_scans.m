function [data, theta_out, alpha_out] = load_an_scans(files_ind,path_name)

% z_zero = 3.41e6;


fname = [path_name,'/An', num2str(files_ind(1), '%06.f'), '.mat'];
load(fname)

theta_out=angle_vec'/1e6;
% alpha_out=alpha_vec;

alpha_out=NaN*zeros(1,length(files_ind));
 

data = zeros(length(files_ind), length(counts));
for ind=1:length(files_ind)
    fname = [path_name,'/An', num2str(files_ind(ind), '%06.f'), '.mat'];
    load(fname)
    data(ind,:) = counts;
    alpha_out(ind)=alpha_vec(ind)/1e6;
%     zs = meas.z_positions;    
end

% zs = (z_zero - zs)*1e-6;
data=data';
end

