function [ im_norm ] = imnorm( im_orig )
%IMNORM min-max normalization of the 2-d image to the range of [0,1]

im_norm = (im_orig - min(min(min(im_orig))))...
          / (max(max(max(im_orig))) - min(min(min(im_orig))));

end

