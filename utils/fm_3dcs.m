function [ y ] = fm_3dcs( A,x )
%FM_3DCS Forward model of three dimensional (3D) compressive sensing (CS).
%   y=FM_3DCS(A,x) returns the multiple frame vector of forward calculation
%   of CS, that is y_k=A_k*x_k for frame k, where A is the M-by-N-by-F 
%   sensing or measurement matrix, x is the N-by-F object matrix, and y is 
%   the M-by-F measurement matrix.
%   See also FMT_3DCS.
[M,N,F] = size(A);
for k = 1:F
    y(:,k) = A(:,:,k)*x(:,k);
end
% y = sum(A.*repmat(x,[M 1 1]),2);
end

