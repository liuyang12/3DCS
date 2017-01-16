function [ im ] = imreadallfmt( impath )
%IMREADALLFMT read images of all formats
%   im=IMREADALLFMT(impath) returns the image matrix under such impath.
%   See also IMREAD.
impathfmt = [impath '.jpg'];
if ~exist(impathfmt,'file')
    impathfmt = [impath '.png'];
elseif ~exist(impathfmt,'file')
    impathfmt = [impath '.bmp'];
elseif ~exist(impathfmt,'file')
    impathfmt = [impath '.tif'];
elseif ~exist(impathfmt,'file')
    error('specify image format and use imread instead.\n');
end
im = imread(impathfmt); % read image

end

