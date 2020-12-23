function spiral_ar_3(data_dir, save_folder)

% save_image = 0;
NFOV = 2.5;
nos_one = 1; % number of spirals per time frame
isize = 84;
%%
imsize = [isize, isize] * NFOV;

load(fullfile(data_dir.folder, data_dir.name));

scale_factor = 1e3 * prod(imsize) / max(abs(kdata(:)));

kSpace = single(permute(kdata, [1, 3, 2])) * scale_factor;

sx = size(kSpace,1);
nc = size(kSpace,3);

kx = real(kloc) * imsize(1);
ky = imag(kloc) * imsize(2);

Nframe = size(kSpace,2)/nos_one;
Nframe = floor(Nframe);
% Nframe = 400;

frame_begin = 1;
% drop the data at the end
kSpace = kSpace(:,(frame_begin-1)*nos_one + 1 : (frame_begin + Nframe - 1)*nos_one,:);
kx = kx(:,(frame_begin-1)*nos_one + 1 : (frame_begin + Nframe - 1)*nos_one);
ky = ky(:,(frame_begin-1)*nos_one + 1 : (frame_begin + Nframe - 1)*nos_one);
nos = Nframe * nos_one;
w = repmat(w,[1,ceil(nos/size(w,2))]);
w(:,Nframe*nos_one+1:end) = [];
kx(:,Nframe*nos_one+1:end) = [];
ky(:,Nframe*nos_one+1:end) = [];
nos = Nframe*nos_one;

%% NUFFT
kSpace = reshape(kSpace,[sx,nos_one,Nframe,nc]);
kx = reshape(kx,[sx,nos_one,Nframe]);
ky = reshape(ky,[sx,nos_one,Nframe]);

t1 = tic;
N = NUFFT.init_new_2(squeeze(kx),squeeze(ky),1,[4,4],imsize(1),imsize(1));
N.W = single(w(:,1));

Image_mean = NUFFT.NUFFT_adj_new_2(permute(kSpace, [1, 3, 2, 4]),N);

%% artifact reduction

% Image_mean = mean(Image, 3);

% display
% Image_disp = abs(rot90(Image_mean));
% Image_disp = Image_disp./max(max(Image_disp));
% Image_disp = reshape(Image_disp, [imsize, 2, 4]);
% Image_disp = permute(Image_disp, [1, 3, 2, 4]);
% Image_disp = reshape(Image_disp, [imsize(1) * 2, imsize(2) * 4]);
% f = figure; imagesc(Image_disp);
% axis image
% axis off
% colormap gray
% set(gcf, 'Position', [100, 100, 800, 400])
% set(gca, 'pos', [0,0,1,1])
% if save_image
%     hgexport(f, './artefact_reduction/temporal_combined_image.eps')
% end
% end

%% manual mask selection
% mask = get_mask(sos(abs(Image_mean)));

%% try to indentify if there is hotspot
Image_sos = sos(Image_mean);

figure
imagesc(Image_sos)
axis image
axis off
colormap gray
brighten(0.3)

max_FOV = max(vec(crop_half_FOV(Image_sos, [isize, isize])));
mask_ = zeros(imsize);
mask_(round(imsize(1)/2) + 1, round(imsize(2)/2) + 1) = 1;
mask_ = bwdist(mask_);
Image_out = Image_sos .* (mask_ > imsize(1)/NFOV/2 * 1.25);
max_out = max(vec(Image_out));

% if there is no hotspot, return
max_ratio = max_out / max_FOV;
if max_ratio < 0.5
%     return
end
%% automatic mask selection
Image_out = imgaussfilt(Image_out, 3);

max_intensity = max(Image_out(:));
mask = Image_out > max_intensity * 0.5;
mask = bwdist(mask);
mask = mask < 2;
CC = bwconncomp(mask);
for ii = 1:CC.NumObjects
    max_int(ii) = max(vec(Image_out(CC.PixelIdxList{ii})));
end
idx = find(max_int == max_intensity);
mask = false(imsize);
mask(CC.PixelIdxList{idx}) = true;


% display
% Image_disp = abs(rot90(Image_mean.*mask));
% Image_disp = Image_disp./max(max(Image_disp));
% Image_disp = reshape(Image_disp, [imsize, 2, 4]);
% Image_disp = permute(Image_disp, [1, 3, 2, 4]);
% Image_disp = reshape(Image_disp, [imsize(1) * 2, imsize(2) * 4]);
% f = figure; imagesc(Image_disp);
% axis image
% axis off
% colormap gray
% set(gcf, 'Position', [100, 100, 800, 400])
% set(gca, 'pos', [0,0,1,1])
% if save_image
%     hgexport(f, './artefact_reduction/masked_image.eps')
% end
% end

k = NUFFT.NUFFT_new_2(Image_mean .* mask, N);

mask_im = NUFFT.NUFFT_adj_new_2(k, N);

A = Image_mean(repmat(mask, [1, 1, 1, nc]));
B = mask_im(repmat(mask, [1, 1, 1, nc]));
scale_k = B \ A;

clear A B mask_im

kSpace = kSpace - permute(k, [1, 3, 2, 4]) * scale_k;
kdata = permute(kSpace, [1, 4, 3, 2]) / scale_factor;
time_ar = toc(t1);

Image_ar = NUFFT.NUFFT_adj_new_2(permute(kSpace, [1, 3, 2, 4]), N);

AR.Image_mean = Image_mean;
AR.Image_ar = Image_ar;
AR.time = time_ar;
AR.mask = mask;
AR.max_ratio = max_ratio;

if nargin == 1
    save(fullfile(data_dir.folder, [data_dir.name(1:end-4), '_ar.mat']), 'kdata', 'kloc', 'spokeindex', 'w', 'AR')
else
    save(fullfile(save_folder, [data_dir.name(1:end-4), '_ar.mat']), 'kdata', 'kloc', 'spokeindex', 'w', 'AR')
end

end