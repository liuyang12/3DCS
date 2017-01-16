function [  ] = savevideoframes( vpath,fpath,params )
%SAVEVIDEOFRAMES Save all the frames of a video
%   SAVEVIDEOFRAMES(vpath,fpath,params) saves all the video frames from
%   vpath to the directory fpath.

% [0] parameters configuration
width  = params.width;
height = params.height;
format = params.format;

% [1] read the video file from vpath
mov = yuv2mov(vpath,width,height,format); % convert 

% [2] save all the frames
if ~exist(fpath,'dir')
    mkdir(fpath);
end
nframe = length(mov);
for iframe = 1:nframe
    cur = mov(iframe).cdata;
%     if size(cur,3)>1
%         cur = rgb2gray(cur);
%     end
    imwrite(cur,sprintf('%s/frame%04d.png',fpath,iframe));
end

