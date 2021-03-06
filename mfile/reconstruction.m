function Image_recon = reconstruction(para)
%--------------------------------------------------------------------------
%   Image_recon = reconstruction(para)
%--------------------------------------------------------------------------
%   Main reconstruction function used in MRM-20-21688, to reproduce no
%   correction, ES correction, and LF reconstruction for a dataset
%   corresponding to Figure 3. 
%--------------------------------------------------------------------------
%   Please cite this paper if you use this code:
%       [1]     Aliasing Artifact Reduction for Spiral Real-Time MRI. MRM,
%               20.21688
%--------------------------------------------------------------------------
%   Author:
%       Ye Tian
%       E-mail: phye1988@gmail.com
%
%   Copyright:
%       MREL, 2020
%       https://mrel.usc.edu
%--------------------------------------------------------------------------

%% read data
% kdata: kspace data, [nsample, ncoil, narm]
% kloc : real (kx) and imag (ky), [nsample, narm]
% w:     density compensation function [nsample, 13]
load(fullfile(para.dir.raw_file), 'kdata', 'kloc', 'w')

matrix_size = para.Recon.matrix_size;
para.Recon.image_size  = round(matrix_size * para.Recon.FOV);

%% fix dimentions for kspace, kx and ky
scale_factor = 1e3 * prod(para.Recon.image_size) / max(abs(kdata(:)));
kSpace = single(permute(kdata,[1, 3, 2])) * scale_factor;

% correct image orintation (90 degree roation)
kx = -imag(kloc);
ky = real(kloc);

clearvars -except kSpace kx ky w para

[sx, ns, nc] = size(kSpace);

narm = para.Recon.narm;

nof = floor(ns/narm);
kSpace(:, nof*narm+1:end, :) = [];
kx(:, nof*narm+1:end) = [];
ky(:, nof*narm+1:end) = [];

kSpace = reshape(kSpace, [sx, narm, nof, nc]);
kx = reshape(kx, [sx, narm, nof]);
ky = reshape(ky, [sx, narm, nof]);

%% select time frames
if isnumeric(para.Recon.time_frames)
    time_frames = para.Recon.time_frames;
    nof = length(time_frames);
    
    kx = kx(:, :, time_frames);
    ky = ky(:, :, time_frames);
    kSpace = kSpace(:, :, time_frames, :);
end

%% Estimate M2 (aliasing source)
switch para.Recon.method
    case 'ES'
        kSpace = spiral_ar(kSpace, kx, ky, w, para);
    case 'LF'
        [~, AR] = spiral_ar(kSpace, kx, ky, w, para);
        para.Recon.FOV = 2.5;
        mask = AR.m1 | AR.m2;
end

%% normalize kx, ky
matrix_size = para.Recon.matrix_size;
para.Recon.image_size  = round(matrix_size * para.Recon.FOV);

kx = kx * para.Recon.image_size(1);
ky = ky * para.Recon.image_size(2);

%% NUFFT structure
Data.N = NUFFT.init(kx, ky, 1, [4, 4], para.Recon.image_size(1), para.Recon.image_size(1));
Data.N.W = w(:, 1);

%% initial estimate, sensitivity map
Data.kSpace = kSpace;
Data.first_est = NUFFT.NUFFT_adj(Data.kSpace, Data.N);

scale = max(abs(Data.first_est(:)));

Data.sens_map = get_sens_map(Data.first_est, '2D');
switch para.Recon.method
    case 'LF'
        Data.sens_map = Data.sens_map .* mask;
        scale = max(vec(abs(Data.first_est) .* mask));
end
Data.first_est = sum(Data.first_est .* conj(Data.sens_map), 4);

%% set parameters
para.Recon.no_comp = nc;
para.Recon.weight_tTV = scale * para.weight_tTV; % temporal regularization weight
para.Recon.weight_sTV = scale * para.weight_sTV; % spatial regularization weight

clearvars -except Data para

%% conjugate gradient reconstruction
[Image_recon, para] = STCR_conjugate_gradient(Data, para);

%% rotate, crop image
Image_recon = abs(Image_recon);
Image_recon = crop_half_FOV(Image_recon, para.Recon.matrix_size);

%% save reconstruction
save(para.dir.save_recon, 'Image_recon')

end
