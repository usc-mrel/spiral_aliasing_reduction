%% add path
addpath ./mfile/

%% add orchestra path
% addpath ./orchestra-sdk-1.9-1.matlab
% addpath ./vdspiral

%% 
ccc

load('./RawData/SI.mat')
if ~isfolder('./figures/')
    mkdir('./figures')
end

%% parameters
iso_center = [384, 330];
im_siz = [1024, 1024];

%% zero pad image
si_siz = size(SI);
im = zeros(im_siz,'single');
im(im_siz(1)/2-iso_center(1)+1:im_siz(1)/2-iso_center(1)+si_siz(1),im_siz(2)/2-iso_center(2)+1:im_siz(2)/2-iso_center(2)+si_siz(2)) = SI;
im = circshift(im, [0, -70]);

%% draw figure 
f = figure;
imagesc(imresize(im(1:876, 74:74+876), [1024, 1024], 'nearest'))
axis image 
axis off
colormap gray
brighten(0.3)
set(gcf, 'Position', [100, 100, 400, 400])
text(512, 51, 'Numerical Phantom','FontSize',20,'VerticalAlignment', 'middle', 'HorizontalAlignment', 'center', 'Color', 'white')
text(51, 51, '(a)','FontSize',20,'VerticalAlignment', 'middle', 'HorizontalAlignment', 'center', 'Color', 'white')
set(gca, 'pos', [0.0, 0.0, 1, 1])
c = colorbar;
c.Position = [0.85 0.12 0.05 0.7];
c.Color = [1,1,1];
c.FontSize = 14;

text(512, 227, '20 cm','FontSize',14,'VerticalAlignment', 'middle', 'HorizontalAlignment', 'center', 'Color', 'white')
hold on

plot([512.5, 512.5] - 233.79/2, [227-20, 227+20], 'Color', 'white', 'LineWidth', 2)
plot([512.5, 512.5] + 233.79/2, [227-20, 227+20], 'Color', 'white', 'LineWidth', 2)

hgexport(f, 'figures/phantom.eps')

%% parameter 
n_spiral = 13;
resolution_scale = 1024/200;

%% Spiral Traj Design
matrix_size = 84;
gmax   = 4;       % maximum gradient [G/cm]
nos    = 13;       % number of spiral interleaves => 8 to 12 interleaves

fov    = 20;      % FOV [cm]
res    = fov*10/matrix_size; % resolution [mm]
ngmax  = 1e5;     % max number of points along gradient
td     = 1;       % Gradient dwell time [usec]
smax   = 15000;   % maximum slew rate [G/cm/sec] (e.g., 150 [T/m/sec])
bw     = 1/(2*td*1e-6)/1000; % receive bandwidth in kHz (BW = [-1/2dt, 1/2dt])
[k_base,g_base,slew,t] = vdsmex(nos, fov, res, gmax, smax, td*1e-6, ngmax); % nt x 2

t      = t*1000; % change unit from [s] to [ms]
k_max  = 1/2/res*10; % k-space max unit in [cm^-1] (1/(2*resolution))

kx = k_base(:,1)/k_max*matrix_size;
ky = k_base(:,2)/k_max*matrix_size;

np = length(t);
kx = k_base(:,1)/k_max*matrix_size/2;
ky = k_base(:,2)/k_max*matrix_size/2;

np = length(t);

%% uniform sampling
theta = (0:12)*2*pi/13;
theta = repmat(theta,[1,ceil(n_spiral/13)]);
theta(:,n_spiral+1:end) = [];
theta = permute(theta,[1,3,2]);
R = [cos(theta),-sin(theta); sin(theta),cos(theta)];

k = [kx,ky];
k = permute(k,[3,2,4,1]);
k = R.*k;
k = sum(k,2);
k = permute(k,[4,1,3,2]);

kx_GA = squeeze(k(:,1,:));
ky_GA = squeeze(k(:,2,:));

w = voronoidens(kx_GA(:)+1i*ky_GA(:)); % nk x ni => nk*ni x 1
w = w / max(w(:));
w = reshape(w,[np,nos]);

%% coil images 
load('./RawData/sens_3xFOV.mat')
cmap = circshift(cmap, [-0, -10]);
nc = size(cmap, 3);

cmap = cmap./max(vec(sos(cmap)));

%% draw coil sensitivity 
f = figure;
coil_show = sos(cmap); 
imagesc(imresize(coil_show(1:876, 74:74+876), [1024, 1024], 'nearest'))
axis image 
axis off
colormap gray
brighten(0.3)
set(gcf, 'Position', [100, 100, 400, 400])
text(512, 51, 'SoS Coil Sensitivity Map','FontSize',20,'VerticalAlignment', 'middle', 'HorizontalAlignment', 'center', 'Color', 'white')
text(51, 51, '(b)','FontSize',20,'VerticalAlignment', 'middle', 'HorizontalAlignment', 'center', 'Color', 'white')
set(gca, 'pos', [0.0, 0.0, 1, 1])
c = colorbar;
c.Position = [0.85 0.12 0.05 0.7];
c.Color = [1,1,1];
c.FontSize = 14;

text(512, 227, '20 cm','FontSize',14,'VerticalAlignment', 'middle', 'HorizontalAlignment', 'center', 'Color', 'white')
hold on

plot([512.5, 512.5] - 233.79/2, [227-20, 227+20], 'Color', 'white', 'LineWidth', 2)
plot([512.5, 512.5] + 233.79/2, [227-20, 227+20], 'Color', 'white', 'LineWidth', 2)

hgexport(f, 'figures/coil.eps')

%% gradient nonlinearity
load('./RawData/GradWrapCoef.mat')
load('./RawData/corners.mat')

% we use a small image to reduce computation
im_small = crop_half_FOV(im, [768, 768]);

corners.UpperLeft = [0 -(768/2) -(768/2)];
corners.UpperRight = [0 768/2 -(768/2)];
corners.LowerLeft = [0 -(768/2) 768/2];

[im_grad_wrap, gw_freq, gw_phase, gw_scale] = GERecon('Gradwarp', im_small, corners,'SphericalHarmonicCoefficients', GradWarpCoef);

[ky, kx] = meshgrid(-512:511); kx = single(kx); ky = single(ky);
[gy, gx] = meshgrid(-384:383); gx = single(gx); gy = single(gy);

mask = im_small ~= 0;
im_coil = im_small .* crop_half_FOV(cmap, [768, 768]);
im_in = reshape(im_coil(repmat(mask, [1, 1, 8])), [114615, 1, 8]);

%% generate distorted k-space data

% This will take a long time (~1.5 hour). Instead, a pre-calculated file is
% loaded
load('./RawData/kSpace_distorted.mat', 'k')
% k = exact_fft2(im_in, ky, kx, gx(mask) + gw_freq(mask) - 1, gy(mask) + gw_phase(mask) - 1);

im_wrap = fftshift2(ifft2(fftshift2(k)));

%% draw destorted image 
f = figure;
im_show = sos(im_wrap);

imagesc(imresize(im_show(1:876, 74:74+876), [1024, 1024], 'nearest'))
axis image 
axis off
colormap gray
brighten(0.3)

set(gcf, 'Position', [100, 100, 400, 400])
text(512, 51, 'SoS Phantom Image with','FontSize',20,'VerticalAlignment', 'middle', 'HorizontalAlignment', 'center', 'Color', 'white')
text(512, 51*2.5, 'Gradient Nonlinearity','FontSize',20,'VerticalAlignment', 'middle', 'HorizontalAlignment', 'center', 'Color', 'white')
text(51, 51, '(c)','FontSize',20,'VerticalAlignment', 'middle', 'HorizontalAlignment', 'center', 'Color', 'white')
set(gca, 'pos', [0, 0, 1, 1])
c = colorbar;
c.Position = [0.85 0.12 0.05 0.7];
c.Color = [1,1,1];
c.FontSize = 14;

text(512, 227, '20 cm','FontSize',14,'VerticalAlignment', 'middle', 'HorizontalAlignment', 'center', 'Color', 'white')
hold on

plot([512.5, 512.5] - 233.79/2, [227-20, 227+20], 'Color', 'white', 'LineWidth', 2)
plot([512.5, 512.5] + 233.79/2, [227-20, 227+20], 'Color', 'white', 'LineWidth', 2)
hgexport(f, 'figures/phantom_wrap.eps')

%% low pass filter
N = NUFFT.init(kx_GA * resolution_scale, ky_GA * resolution_scale, 1, [6,6], im_siz(1));

point = zeros(im_siz);
point(im_siz(1)/2, im_siz(2)/2) = 1;

kSpace = NUFFT.NUFFT(permute(im_wrap, [1, 2, 4, 3]), N);
kSpace_lp = zeros(np, nos, 1, nc);
for i = 1:nos
    for j = 1:nc
        kSpace_lp(:, i, :, j) = lowpass(kSpace(:,i, :, j), 0.25);
    end
end
kSpace_lp = kSpace_lp(1:4:end,:, :, :);

%% low pass filtered image 
N = NUFFT.init(kx_GA(1:4:end,:)*3, ky_GA(1:4:end,:)*3, 1, [6,6], 84*3, 84*3);
N.W = w(1:4:end,:);
test_im_lp = NUFFT.NUFFT_adj(kSpace_lp, N);

im_show = sos(test_im_lp);
scale_im_show = max(im_show(:));
im_show = im_show ./ scale_im_show;

f = figure; imagesc(im_show); colormap gray; axis image; axis off
brighten(0.4)
set(gcf, 'Position', [100, 100, 400, 400])
text(126, 12.6, 'Spiral Trajectory Sampled Phantom','FontSize',20,'VerticalAlignment', 'middle', 'HorizontalAlignment', 'center', 'Color', 'white')
text(126, 12.6*2.5, 'with Low-pass Filter','FontSize',20,'VerticalAlignment', 'middle', 'HorizontalAlignment', 'center', 'Color', 'white')
text(12.6, 12.6, '(e)','FontSize',20,'VerticalAlignment', 'middle', 'HorizontalAlignment', 'center', 'Color', 'white')
set(gca, 'pos', [0, 0, 1, 1])

text(126, 55, '20 cm','FontSize',14,'VerticalAlignment', 'middle', 'HorizontalAlignment', 'center', 'Color', 'white')
hold on
plot([126, 126] - 42, [55-4.9, 55+4.9], 'Color', 'white', 'LineWidth', 2)
plot([126, 126] + 42, [55-4.9, 55+4.9], 'Color', 'white', 'LineWidth', 2)

p = patch([111, 145, 145, 111], [195, 195, 232, 232], [1, 1, 1]);
p.FaceColor = 'none';
p.LineStyle = '--';
p.LineWidth = 2;
p.EdgeColor = [1, 1, 1];
hgexport(f, 'figures/phantom_low_pass.eps')

%% generate kspace 
N = NUFFT.init(kx_GA(1:4:end,:) * resolution_scale, ky_GA(1:4:end,:) * resolution_scale, 1, [6,6], im_siz(1));
kSpace = NUFFT.NUFFT(permute(im_wrap, [1, 2, 4, 3]), N);

%% reconstructed image without low pass filter
N = NUFFT.init(kx_GA(1:4:end,:)*3, ky_GA(1:4:end,:)*3, 1, [6,6], 84*3, 84*3);
N.W = w(1:4:end,:);
test_im = NUFFT.NUFFT_adj(kSpace, N);

f = figure; imagesc(sos(test_im) ./ scale_im_show),colormap gray; axis image, axis off
brighten(0.4)
set(gcf, 'Position', [100, 100, 400, 400])
text(126, 12.6, 'Spiral Trajectory Sampled Phantom','FontSize',20,'VerticalAlignment', 'middle', 'HorizontalAlignment', 'center', 'Color', 'white')
text(126, 12.6*2.5, 'without Low-pass Filter','FontSize',20,'VerticalAlignment', 'middle', 'HorizontalAlignment', 'center', 'Color', 'white')
text(12.6, 12.6, '(d)','FontSize',20,'VerticalAlignment', 'middle', 'HorizontalAlignment', 'center', 'Color', 'white')
set(gca, 'pos', [0, 0, 1, 1])

text(126, 55, '20 cm','FontSize',14,'VerticalAlignment', 'middle', 'HorizontalAlignment', 'center', 'Color', 'white')
hold on
plot([126, 126] - 42, [55-4.9, 55+4.9], 'Color', 'white', 'LineWidth', 2)
plot([126, 126] + 42, [55-4.9, 55+4.9], 'Color', 'white', 'LineWidth', 2)

p = patch([111, 145, 145, 111], [195, 195, 232, 232], [1, 1, 1]);
p.FaceColor = 'none';
p.LineStyle = '--';
p.LineWidth = 2;
p.EdgeColor = [1, 1, 1];

hgexport(f, 'figures/phantom_without_low_pass.eps')