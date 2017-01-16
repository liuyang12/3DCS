function [ y ] = fm_cassi( H,f )
%FM_CASSI Forward model of the coded aperture snapshot spectral imager 
%(CASSI).
%   y=FM_CASSI(H,f) returns the coded image of CASSI, that is 
%   y=sum(H.*f,3), where H is Nx-by-Ny-by-Nt matrix denoting the code or 
%   patterns of the imager, f is the original 3D sample with the same
%   dimensions as H.
%   See also FMT_CASSI.
y  = sum(H.*f,3);
end

