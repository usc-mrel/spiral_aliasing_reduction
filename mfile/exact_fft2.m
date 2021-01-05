function kSpace = exact_fft2(image, kx, ky, gx, gy)

[sx, sy, nof, nc] = size(image);

skx = size(kx, 1);
nor = size(kx, 2);

if isa(image,'gpuArray')
    kSpace = gpuArray(single(zeros(skx,nor,nof,nc)));
else
    kSpace = single(zeros(skx,nor,nof,nc));
end


for Nx = 1:skx
    Nx
    tic
    for Nrays = 1:nor
        kx_temp = kx(Nx, Nrays,:);
        ky_temp = ky(Nx, Nrays,:);
        phase = bsxfun(@times, gy, kx_temp)/skx + bsxfun(@times, gx ,ky_temp)/skx;
        phase = exp(-2 * pi * 1i * phase);
        temp = bsxfun(@times, image, phase);
        kSpace(Nx, Nrays, :, :) = sum(sum(temp,1),2);
    end
    toc
end

