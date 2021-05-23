function [kSpace_ar, AR] = spiral_ar(kSpace, kx, ky, w, para)
t1 = tic;

sx      = size(kSpace, 1);
narm    = size(kSpace, 2);
nof     = size(kSpace, 3);
nc      = size(kSpace, 4);
NFOV    = 2.5;

image_size  = round(para.Recon.matrix_size * NFOV);
kx_ = reshape(kx, [sx, narm * nof]) * image_size(1);
ky_ = reshape(ky, [sx, narm * nof]) * image_size(2);

kSpace_raw = reshape(kSpace, [sx, narm * nof, 1, nc]);

N = NUFFT.init(kx_, ky_, 1, [4, 4], image_size(1), image_size(1));
N.W = w(:, 1);

Image_mean = NUFFT.NUFFT_adj(kSpace_raw, N);

Image_sos = sos(Image_mean);

%% display image
figure
imagesc(Image_sos)
axis image
axis off
colormap gray
brighten(0.3)
title 'Temporally Combined Image'

%% caucluate ratio to determin whether this dataset has artifact
max_FOV = max(vec(crop_half_FOV(Image_sos, para.Recon.matrix_size)));
mask_ = zeros(image_size);
mask_(round(image_size(1)/2) + 1, round(image_size(2)/2) + 1) = 1;
mask_ = bwdist(mask_);
Image_out = Image_sos .* (mask_ > image_size(1)/NFOV/2 * 1.25);
max_out = max(vec(Image_out));

max_ratio = max_out / max_FOV;
if max_ratio < 0.4
    clearvars -except kSpace kx ky para
    kSpace_ar = kSpace;
    AR.m1 = ones(image_size, 'single');
    AR.m2 = ones(image_size, 'single');
    return
end

%% detect the mask automaticlly
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
mask = false(image_size);
mask(CC.PixelIdxList{idx}) = true;

k = NUFFT.NUFFT(Image_mean .* mask, N);

mask_im = NUFFT.NUFFT_adj(k, N);

%% scake_k accounts for NUFFT scaling issue
A = Image_mean(repmat(mask, [1, 1, 1, nc]));
B = mask_im(repmat(mask, [1, 1, 1, nc]));
scale_k = B \ A;

clear A B mask_im

kSpace_ar = kSpace_raw - k * scale_k;

time_ar = toc(t1);

%% display image
Image_ar = NUFFT.NUFFT_adj(kSpace_ar, N);
figure
imagesc(sos(Image_ar))
axis image
axis off
colormap gray
brighten(0.3)
title 'Artifact Reduced Temporally Combined Image'

kSpace_ar = reshape(kSpace_ar, [sx, narm, nof, nc]);

if nargout >1
    AR.Image_mean   = Image_mean;
    AR.Image_ar     = Image_ar;
    AR.time         = time_ar;
    AR.m1           = mask_ < para.Recon.matrix_size(1) * sqrt(2) / 2;
    AR.m2           = mask;
    AR.max_ratio    = max_ratio;
end
end