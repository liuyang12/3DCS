function [ H ] = haarmtx( n )
%HAARMTX Haar wavelet matrix
%   H=HAARMTX(n) returns the Haar wavelet matrix, where n is the k-th power
%   of two (2^k).
k = floor(log2(n));
if 2^k < n
    error('%d is not the power of 2.\n',n);
end

H = [1];
NC = 1/sqrt(2);
LP = [ 1  1];
HP = [ 1 -1];

for i = 1:k
    H = NC*[kron(H,LP);kron(eye(size(H)),HP)];
end
end

