% script: CSA_layer_detector
% By John Paden and Mingze Xu, July 2016

% fn_dir: Directory where files are at
fn_dir = './data/';

% Compile Cpp files
mex -largeArrayDims fuse.cpp
mex -largeArrayDims train_model.cpp
mex -largeArrayDims detect.cpp
mex -largeArrayDims refine.cpp

% img_01: right looking (positive theta)
% img_02: nadir looking (theta == 0)
% img_03: left looking (negative theta)
nadir_img = 2;

%% Automated loading section
% =========================================================================
% Load Data
mdata = {};
for img=1:3
    fn = fullfile(fn_dir,sprintf('Data_img_%02.0f_20140401_03_044.mat',img));
    mdata{img} = load(fn);
end

% Check for twtt DEM data
if ~isfield(mdata{1},'twtt')
    Nx = length(mdata{1}.GPS_time);
    Nsv = length(mdata{1}.param_combine.array_param.theta);
    twtt = NaN*zeros(Nsv,Nx);
end

% Interpolate images onto a common propagation time axis
for img  = [1 3]
    mdata{img}.Data = interp1(mdata{img}.Time,mdata{img}.Data,mdata{nadir_img}.Time);
    mdata{img}.Topography.img = interp1(mdata{img}.Time,mdata{img}.Topography.img,mdata{nadir_img}.Time);
    mdata{img}.Time = mdata{nadir_img}.Time;
end

% Convert Surface and Bottom variables from propagation time into image pixels
Surface_bin = interp1(mdata{1}.Time, 1:length(mdata{1}.Time), mdata{1}.Surface);
Surface_Mult_bin = interp1(mdata{1}.Time, 1:length(mdata{1}.Time), 2 * mdata{1}.Surface);
Bottom_bin = interp1(mdata{1}.Time, 1:length(mdata{1}.Time), mdata{1}.Bottom);
twtt_bin = interp1(mdata{1}.Time, 1:length(mdata{1}.Time), mdata{1}.twtt);
ice_mask = mdata{1}.ice_mask;
Bottom_bin(isnan(Bottom_bin)) = 0;
twtt_bin(isnan(twtt_bin)) = 0;

%% Automated labeling section
% =========================================================================
% Specify which range lines to browse
skip = 10;
rlines = 1:skip:size(mdata{1}.Topography.img,3);

% Training parameters and preparing fusion slices
mu = [];
sigma = [];
dataset = [];
for rline = rlines
    fusion = fuse(double(db(mdata{1}.Topography.img(:,:,rline))), ...
        double(db(mdata{2}.Topography.img(:,:,rline))), ...
        double(db(mdata{3}.Topography.img(:,:,rline))));
    fusion(fusion>27) = 27;
    fusion = reshape(fusion, size(db(mdata{1}.Topography.img(:,:,rline))));
    [m, s] = train_model(fusion, double(twtt_bin(:,rline)));
    dataset(:,:,rline) = fusion;
    mu(:,rline) = m;
    sigma(:,rline) = s;
end

% Plot echogram
figure(2); clf;
imagesc(db(mdata{nadir_img}.Data));
hold on
rline = 1;
ylim([50 500]);
plot(Surface_bin,'k-');
plot(Bottom_bin,'k-');
hplot = plot([rline rline],[1 size(mdata{1}.Data,1)],'k');

% Setup plots
h_image = [];
h_axes = [];
h_title = [];
h_surf_plot = [];
h_bot_plot = [];
h_surf_mult_plot = [];
h_dem_plot= [];

labels = detect(dataset(:,:,rline), double(twtt_bin(:,rline)), double(Bottom_bin(rline)), ...
    [], double(ice_mask(:,rline)), double(mu(:,rline)), double(sigma(:,rline)));

figure(1); clf;
h_image(1) = imagesc(dataset(:,:,rline));
h_axes(1) = gca;
hold on;
h_title(1) = title(sprintf('%d',rline));
h_surf_plot(1) = plot(33,Surface_bin(rline),'kx','MarkerSize',7,'LineWidth',2);
h_bot_plot(1) = plot(33,Bottom_bin(rline),'kx','MarkerSize',7,'LineWidth',2);
h_surf_mult_plot(1) = plot(33,Surface_Mult_bin(rline),'mx','MarkerSize',7,'LineWidth',2);
h_dem_plot(1) = plot(twtt_bin(:,rline),'kx','MarkerSize',7,'LineWidth',2);
h_ice_bed_plot(1) = plot(labels,'b^','MarkerSize',7,'LineWidth',2);
ylim(Surface_bin(rline)+[-25 500]);
caxis([4 27]);
linkaxes(h_axes,'xy');

%% Automated cycling section
% =========================================================================
% 3D surface
bottom_surface = [];

for rline = rlines
    fprintf('rline: %d\n', rline);
    labels = detect(dataset(:,:,rline), double(twtt_bin(:,rline)), double(Bottom_bin(rline)), ...
        [], double(ice_mask(:,rline)), double(mu(:,rline)), double(sigma(:,rline)));
    
    % If empty, assign random values
    if isempty(labels)
        bottom_surface = [bottom_surface; randi([1 700], 1, size(bottom_surface, 2))];
    else
        bottom_surface = [bottom_surface; labels];
    end
    
    set(h_image(1),'CData',dataset(:,:,rline));
    set(h_title(1),'String',sprintf('%d',rline));
    set(h_surf_plot(1),'YData',Surface_bin(rline))
    set(h_bot_plot(1),'YData',Bottom_bin(rline))
    set(h_surf_mult_plot(1),'YData',Surface_Mult_bin(rline))
    set(h_dem_plot(1),'YData',twtt_bin(:,rline))
    set(h_ice_bed_plot(1),'YData',labels)
    ylim(Surface_bin(rline)+[-25 500]);
    set(hplot,'XData',[rline rline]);
    
    keyboard
end
