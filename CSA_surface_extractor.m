% script: CSA_surface_extractor_trws
% By Mingze Xu, July 2016

% fn_dir: Directory where files are at
fn_dir = './data/';

% Compile C++ files
mex -largeArrayDims fuse.cpp
mex -largeArrayDims train_model.cpp
mex -largeArrayDims extract.cpp

% img_01: right looking (positive theta)
% img_02: nadir looking (theta == 0)
% img_03: left looking (negative theta)
nadir_img = 2;

%% Automated loading section
% ======================================================================

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

% Convert surface and bottom variables from propagation time into image pixels
Surface_bin = interp1(mdata{1}.Time, 1:length(mdata{1}.Time), mdata{1}.Surface);
Surface_Mult_bin = interp1(mdata{1}.Time, 1:length(mdata{1}.Time), 2 * mdata{1}.Surface);
Bottom_bin = interp1(mdata{1}.Time, 1:length(mdata{1}.Time), mdata{1}.Bottom);
twtt_bin = interp1(mdata{1}.Time, 1:length(mdata{1}.Time), mdata{1}.twtt);
ice_mask = mdata{1}.ice_mask;
twtt_bin(isnan(twtt_bin)) = 0;
Bottom_bin(isnan(Bottom_bin)) = 0;

%% Automated extracting section
% ======================================================================

% Specify which range lines to browse
skip = 1;
rlines = 1:skip:size(mdata{1}.Topography.img,3);

% Training parameters
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

% TRWS
Extra_bin = [];
bottom_surface = extract(double(dataset), double(twtt_bin), double(Bottom_bin), ...
    double(Extra_bin), double(ice_mask), double(mean(mu,2)), double(mean(sigma,2)));
bottom_surface = reshape(bottom_surface, size(twtt_bin));
filename = 'bottom_surface';
save(filename, 'bottom_surface');

%% Automated correcting section
% ======================================================================

% Plot slices
% rline = 1;
% figure(1); clf;
% h_image(1) = imagesc(dataset(:,:,rline));
% h_axes(1) = gca;
% hold on;
% h_title(1) = title(sprintf('%d',rline));
% h_surf_plot(1) = plot(33,Surface_bin(rline),'kx','MarkerSize',7,'LineWidth',2);
% h_bot_plot(1) = plot(33,Bottom_bin(rline),'kx','MarkerSize',7,'LineWidth',2);
% h_surf_mult_plot(1) = plot(33,Surface_Mult_bin(rline),'mx','MarkerSize',7,'LineWidth',2);
% h_dem_plot(1) = plot(twtt_bin(:,rline),'kx','MarkerSize',7,'LineWidth',2);
% h_ice_bed_plot(1) = plot(bottom_surface(:,rline),'b^','MarkerSize',7,'LineWidth',2);
% ylim(Surface_bin(rline)+[-25 500]);
% caxis([4 27]);
% linkaxes(h_axes,'xy');

% Plot surface
% while(1)
%     figure(2); clf;
%     imagesc(bottom_surface');
%     % Choose which range line
%     [x1, y1] = ginput(1);
%     x1 = floor(x1);
%     y1 = floor(y1);
%
%     if x1 < 0 || x1 > 64 || y1 < 0 || y1 > length(rlines)
%         break;
%     end
%
%     % Plot the chosen slice
%     rline = rlines(y1);
%     figure(1); clf;
%     h_image(1) = imagesc(dataset(:,:,rline));
%     h_axes(1) = gca;
%     hold on;
%     h_title(1) = title(sprintf('%d',rline));
%     h_surf_plot(1) = plot(33,Surface_bin(rline),'kx','MarkerSize',7,'LineWidth',2);
%     h_bot_plot(1) = plot(33,Bottom_bin(rline),'kx','MarkerSize',7,'LineWidth',2);
%     h_surf_mult_plot(1) = plot(33,Surface_Mult_bin(rline),'mx','MarkerSize',7,'LineWidth',2);
%     h_dem_plot(1) = plot(twtt_bin(:,rline),'kx','MarkerSize',7,'LineWidth',2);
%     h_ice_bed_plot(1) = plot(bottom_surface(:,rline),'b^','MarkerSize',7,'LineWidth',2);
%     ylim(Surface_bin(rline)+[-25 500]);
%     caxis([4 27]);
%     linkaxes(h_axes,'xy');
% end
