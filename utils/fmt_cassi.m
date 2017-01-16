function [ fy ] = fmt_cassi( H,y )
%FMT_CASSI Forward model of the coded aperture snapshot spectral imager 
%(CASSI).
%   fy=FMT_CASSI(H,y) returns the result of  transpose of the measurement 
%   matrix multiplies the measurement vector, that is fy=H.*y, 
%   where H is Nx-by-Ny-by-Nt matrix denoting the code or patterns of the 
%   imager, f is the original 3D sample with the same dimensions as H.
%   See also FM_CASSI.
fy  = bsxfun(@times,H,y);
end