function [ fy ] = fmt_3dcs( A,y )
%FMT_3DCS Forward model of three dimensional (3D) compressive sensing (CS).
%   fy=FMT_3DCS(A,y) returns the multiple frame vector of forward calculation
%   of CS, that is y_k=A_k*x_k for frame k, where A is the M-by-N-by-F 
%   sensing or measurement matrix, y is the M-by-F measurement matrix and 
%   fy is the N-by-F object matrix.
%   See also FMT_3DCS.
[M,N,F] = size(A);
for k = 1:F
    fy(:,k) = A(:,:,k)'*y(:,k);
end

end

