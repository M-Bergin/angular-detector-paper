%% Surface generation for gratings



recompile = true;
verbose = true;
distribution = 'cosine';


%% Set Parameters

deltax=1;
N=20;

pitch=2*10;
depth=0.000001;
y_len=N*pitch/2;

alpha_vec=0;%0:1:90;
init_theta=45;
n_rays=1000;

N_alpha=length(alpha_vec);
N_bins=100;
det_R=7;
det_radius=0.25;

% The maximum number of scattering events a ray is allowed to undergo
maxScatter = 1000;

I_mat=NaN*zeros(N_alpha,N_bins);

% Should the mex files be recompiled
loadpath
if recompile
    mexCompile();
end

for n_alpha=1:N_alpha

    clear x_gen x_pt1 x_pt2 x y f vertices fdef fdef_end fnorm initposses

    alpha=alpha_vec(n_alpha);

    % Starting direction of the rays
    init_dir = [sind(init_theta), 0, -cosd(init_theta)];

    % Surface Parameters
    surface_size = y_len;%N*deltax;

    %% Generate surface
    % x = linspace(-surface_size/2,surface_size/2,N); y = linspace(-surface_size/2,surface_size/2,N);
    x_gen=pitch*(-floor(N/2):0.5:ceil(N/2));
    x_pt1=[x_gen(1),x_gen(1)];
    x_pt2=repelem(x_gen(2:end),4);
    x=[x_pt1,x_pt2];
    x(end-1:end)=[];


    y=repmat([-y_len,y_len],1,4*N);

    f=repmat([-depth,-depth,-depth,-depth,0,0,0,0],1,N);

    %Add caps to the end
    N_temp=length(x);
    x_end=repmat([-pitch*floor(N/2),-pitch*floor(N/2),pitch*ceil(N/2),pitch*ceil(N/2)],1,2);
    y_end=[repelem(-y_len,4),repelem(y_len,4)];
    f_end=repmat([0,-depth],1,4);

    x=[x,x_end];
    y=[y,y_end];
    f=[f,f_end];

    %Rotate the co-ordinates by alpha
    A=[cosd(alpha), -sind(alpha);sind(alpha),cosd(alpha)]*[x;y];
    x=A(1,:);
    y=A(2,:);


    % [X,Y] = meshgrid(x,y);
    % f = depth*(mod(X/pitch,2)>1);

    % f = ifft2(fft2(Z).*fft2(F));
    %
    % clear Z F;

    surfImg  = f;

    vertices = [x(:), y(:), f(:)];

    % % Calculate the triangulation for the lower left triangles.
    % nums1 = 1:N^2;
    % upper1 = circshift(nums1, -1);
    % right1 = circshift(nums1, -N);
    % filt = ((mod(nums1,N)~=0) .* (nums1 < N*(N-1))) == 1;
    % nums1 = nums1(filt);
    % upper1 = upper1(filt);
    % right1 = right1(filt);
    % fdef1 = [nums1;upper1;right1];
    %
    % fnorm1 = [f(right1)-f(nums1); f(upper1)-f(nums1); deltax*ones(size(nums1))]; % Can check by cross product
    % fnorm1 = normalize(fnorm1, 1, 'norm', 2);
    %
    % % Calculate the triangulation for the upper right triangles.
    % nums2 = 1:N^2;
    % lower2 = circshift(nums2, 1);
    % left2 = circshift(nums2, N);
    % filt = ((mod(nums2,N)~=1) .* (nums2 > N)) == 1;
    % nums2 = nums2(filt);
    % lower2 = lower2(filt);
    % left2 = left2(filt);
    % fdef2 = [nums2; lower2; left2];
    %
    % fnorm2 = [f(nums2)-f(left2); f(nums2)-f(lower2); deltax*ones(size(nums2))];
    % fnorm2 = normalize(fnorm2, 1, 'norm', 2);
    %
    % fdef = [fdef1, fdef2];
    % fnorm = [fnorm1, fnorm2];


    nums1=1:8*N;
    shift_1=circshift(nums1, -1);
    shift_2=circshift(nums1, -2);

    fdef = [nums1;shift_1;shift_2];

    fdef(:,end-1:end)=[];

    fdef_end=[N_temp+1,N_temp+2,N_temp+3;
        N_temp+2,N_temp+3,N_temp+4;
        N_temp+5,N_temp+6,N_temp+7;
        N_temp+6,N_temp+7,N_temp+8];
    fdef=[fdef,fdef_end'];

    fnorm_unit=[0,0,-1,-1,0,0,1,1;
        0,0,0,0,0,0,0,0;
        1,1,0,0,1,1,0,0];

    fnorm=repmat(fnorm_unit,1,N);

    fnorm(:,end-1:end)=[];

    fnorm_end=[0,0,0,0;
        1,1,-1,-1;
        0,0,0,0];
    fnorm=[fnorm,fnorm_end];

    %Rotate by alpha
    R_rot=rotz(alpha);
    fnorm=R_rot*fnorm;

    lattice=zeros(length(fdef),6);

    clear fmat
    [fmat{1:length(fdef)}] = deal('poop');
    if strcmp(distribution, 'specular')
        mat = struct('function', 'specular', 'params', [], 'color', [0.8,0.8,1]);
    elseif strcmp(distribution, 'cosine')
        mat = struct('function', 'cosine', 'params', [], 'color', [0.8,0.8,1]);
    else
        error('If this error is reached then something bad happened.')
    end
    materials = containers.Map('poop', mat);

    surface = TriagSurface(vertices, fdef', fnorm',lattice, fmat', materials);
    surface.moveBy([0, 0, -max(surface.vertices(:,3))]);
    % figure;surface.patchPlot(false);


    %% Let's bounce!
    % p = gcp();
    % poolsize = p.NumWorkers;
    poolsize=1;
    if verbose
        tic
    end

    nums = {};
    dirs = {};
    initposses = {};
    n = round(n_rays/poolsize);

    for i=1:poolsize
        itdr = repmat(init_dir', 1, n);
        %itpos = [unifrnd(-surface_size/2,surface_size/2, 2, n); zeros(1,n)];
        %itpos = [unifrnd(-surface_size*2/8,surface_size*2/8, 2, n); zeros(1,n)];
        t = 2*pi*rand(1,n);
        r = surface_size*3/8*sqrt(rand(1,n));
        itpos = [r.*cos(t);r.*sin(t); zeros(1,n)];
        initposses{i} = itpos;
        itpos = itpos - 10*itdr;

        [~, nums{i}, ~, dirs{i}] = ...
            distributionCalc('sample_surface', surface, 'maxScatter', maxScatter, ...
            'nrays', n, 'start_pos', itpos, 'start_dir', itdr);

        fprintf('Finish %d.\n', i);
    end

    numScattersRay = zeros(1, n*poolsize);
    final_dir = zeros(3, n*poolsize);
    initPos = zeros(3, n*poolsize);
    for i=1:poolsize
        numScattersRay((i-1)*n+1:i*n) = nums{i};
        final_dir(:,(i-1)*n+1:i*n) = dirs{i};
        initPos(:,(i-1)*n+1:i*n) = initposses{i};
    end

    if verbose
        fprintf('Ray tracing took %d seconds.\n', uint16(toc));
    end

    [az, el] = cart2sph(final_dir(1,:), final_dir(2,:), final_dir(3,:));

    az = az*180/pi;
    el = el*180/pi;

    az_scattered = az(el>0);
    el_scattered = el(el>0);

    singleScattered = [az_scattered; 90-el_scattered];
    all_dirs = [az; 90-el];

    if verbose
        fprintf('%d%% of the rays undergo at least one scattering event.\n', round(length(numScattersRay(numScattersRay > 0)) / length(numScattersRay(:)) *100));
        fprintf('%d%% of the scattered rays undergo multiple scattering events.\n', round(length(numScattersRay(numScattersRay > 1)) / length(numScattersRay(numScattersRay > 0)) *100));
        fprintf('%d%% of the rays have positive elevation.\n', round(length(el(el>0)) / length(numScattersRay(:)) *100));
    end

        %% Plot the result
    
        %plot all the rays
        figure
        filt = numScattersRay > 0 ;
        scatt = all_dirs(:,filt);
        xs = sind(scatt(2,:)) .* cosd(scatt(1,:));
        ys = sind(scatt(2,:)) .* sind(scatt(1,:));
    
        x_bins = linspace(-1,1,75);
        y_bins = linspace(-1,1,75);
        c_x = (x_bins(1:74) + x_bins(2:75)) / 2;
        c_y = (y_bins(1:74) + y_bins(2:75)) / 2;
        [CX, CY] = meshgrid(c_x, c_y);
        Ns = hist2d(xs, ys, x_bins, y_bins, 'tile');
        hold off
        plotme = Ns ./ real(cosd(asind(sqrt(CX.^2+CY.^2))));
        imagesc(c_x, c_y, Ns)
        xticks([])
        yticks([])
        % caxis([0 max(max(Ns .* (CX.^2+CY.^2 < 0.95)))])
        set(gca, 'YDir', 'normal')
        title('Ray Tracing Total Radiance')
        hold on
        plot(cosd(linspace(0,360)), sind(linspace(0,360)), ':w', 'Linewidth', 2);
    
    %Plot single scattered rays
    figure
    filt = numScattersRay == 1;
    sing = all_dirs(:,filt);
    
    xs = sind(sing(2,:)) .* cosd(sing(1,:));
    ys = sind(sing(2,:)) .* sind(sing(1,:));
    x_bins = linspace(-1,1,75);
    y_bins = linspace(-1,1,75);
    
    hist2d(xs, ys, x_bins, y_bins, 'tile');
    xticks([])
    yticks([])
    set(gca, 'YDir', 'normal')
    title('Ray Tracing Single Scattered Radiance')
    hold on
    plot(cosd(linspace(0,360)), sind(linspace(0,360)), ':w', 'Linewidth', 2);
    
    
    figure;surface.patchPlot(false);
    hold on
    plot3(itpos(1,:),itpos(2,:),itpos(3,:),'.')
    
    filt2=final_dir(2,:)<0.1 & final_dir(2,:)>-0.1;
    % figure; histogram(sign(final_dir(2,filt2)).*acosd(final_dir(3,filt2)),100)
    % figure; histogram(acosd(final_dir(3,filt2)))
    figure; histogram((final_dir(1,filt2)),50)
    

    r_cyl=sqrt(final_dir(1,:).^2+final_dir(3,:).^2);
    phi_cyl=mod(atand(final_dir(3,:)./final_dir(1,:))+180,180);
    z_cyl=final_dir(2,:);

    %Apply scale factor to move rays up to the detector
    sf=det_R./r_cyl;
    z_cyl2=sf.*z_cyl;

    filt2=abs(z_cyl2)<det_radius;
    [I_mat(n_alpha,:),phi_c]=histcounts(90-(phi_cyl(filt2)),N_bins);

    %plot all the rays
    figure
    filt = numScattersRay > 0;
    scatt = all_dirs(:,filt);
    
    
    phi_bins = linspace(0,1,75)*180;
    z_bins = linspace(-1,1,75);
    c_x = (phi_bins(1:74) + phi_bins(2:75)) / 2;
    c_y = (z_bins(1:74) + z_bins(2:75)) / 2;
    [CX, CY] = meshgrid(c_x, c_y);
    Ns = hist2d(phi_cyl, z_cyl, phi_bins, z_bins, 'tile');
    % xticks([])
    % yticks([])
    % caxis([0 max(max(Ns .* (CX.^2+CY.^2 < 0.95)))])
    set(gca, 'YDir', 'normal')
    title('Ray Tracing Total')
    % hold on
    % plot(cosd(linspace(0,360)), sind(linspace(0,360)), ':w', 'Linewidth', 2);

        filt2=abs(z_cyl)<0.2;
    %     figure; histogram(sign(final_dir(2,filt2)).*acosd(final_dir(3,filt2)),100)
    %     figure; histogram(acosd(final_dir(3,filt2)))
        figure; histogram(90-(phi_cyl(filt2)),50)
end

% % figure;imagesc(alpha_vec,phi_c,I_mat')
% % 
% % xlabel('\phi/^\circ')
% % ylabel('\theta/^\circ')
% % set(gca,'LineWidth',1,'FontSize',12)
% 
% % ylim([-10 60])
% 
% I_mat_flip=fliplr(I_mat');
% I_mat_tot=[I_mat_flip(:,1:end-1),I_mat'];
% alpha_vec_tot=[-1*fliplr(alpha_vec),alpha_vec(2:end)];
% 
% inds=find(phi_c>-10 & phi_c<65);
% inds_start=min(inds);
% inds_end=max(inds);
% 
% I_mat_crop=I_mat_tot(inds_start:inds_end,:);
% phi_crop=phi_c(inds_start:inds_end);
% 
% 
% figure;imagesc(alpha_vec_tot,phi_crop,I_mat_crop)
% xlabel('\phi/^\circ')
% ylabel('\theta/^\circ')
% set(gca,'LineWidth',1,'FontSize',12)
% 
% 
% % figure;imagesc(alpha_vec_tot,phi_c,I_mat_tot)
% % 
% % xlabel('\phi/^\circ')
% % ylabel('\theta/^\circ')
% % set(gca,'LineWidth',1,'FontSize',12)


