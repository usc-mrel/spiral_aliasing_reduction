function show_yt(img, brightness, white_threshold)
%--------------------------------------------------------------------------
%   show_yt(img, brightness)
%--------------------------------------------------------------------------
%   Display dynamic images frame-by-frame
%--------------------------------------------------------------------------
%   Inputs:      
%       - img           [sx, sy, nof, ...]
%       - brightness    [scalar within 0-1]
%
%       - img           dynamic images with at least 3 dimensions
%--------------------------------------------------------------------------
%   Author:
%       Ye Tian
%       E-mail: phye1988@gmail.com
%--------------------------------------------------------------------------

if nargin < 2
    brightness = 0.2;
    white_threshold = 1;
elseif nargin < 3
    white_threshold = 1;
end

if ~isreal(img)
    img = abs(img);
end

figure
img = squeeze(img);
nof = size(img,3);
img = img ./ max(img(:)) / white_threshold;
i = 1;
while i
    j = mod(i,nof);
    if j == 0
        j = nof;
    end 
    imagesc(abs(img(:,:,j)), [0, 1])
    colormap gray
    axis image
    title(num2str(j))
    brighten(brightness)
    drawnow
    i = i+1;
end