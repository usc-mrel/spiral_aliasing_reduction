function show_yt_pause(img, brightness, white_threshold)
if nargin < 2
    brightness = 0.2;
    white_threshold = 1;
elseif nargin < 3
    white_threshold = 1;
end

if isnumeric(img) || islogical(img)
    img = squeeze(img);
    if ~isreal(img)
        img = abs(img);
    end
%     img = img ./ max(img(:)) / white_threshold;
    figure
    for i=1:size(img,3)
        imagesc(img(:,:,i))
        colormap gray
        axis image
        brighten(brightness)
        title(i)
        %drawnow
        pause
    end
elseif iscell(img)
    figure
    for i=1:numel(img)
        if isreal(img{i})
            imagesc(img{i})
        else
            imagesc(abs(img{i}))
        end
        colormap gray
        axis image
        brighten(brightness)
        title(i)
        pause
    end
end
end