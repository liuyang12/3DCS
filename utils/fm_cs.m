function [ y ] = fm_cs( A,x )
%FM_CS Forward model of compressive sensing (CS).
%   y=FM_CS(A,x) returns the vector of forward calculation of CS, that is 
%   y=A*x, where A is the M-by-N sensing or measurement matrix, x is 
%   the N-by-1 object vector, and y is the M-by-1 measurement vector.
%   See also FMT_CS.
y = A*x;
end

