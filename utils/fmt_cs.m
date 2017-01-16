function [ fy ] = fmt_cs( A,y )
%FMT_CS Transpose of forward model of compressive sensing (CS).
%   fy=FMT_CS(A,y) returns the vector of the transpose of forward 
%   calculation of CS, that is fy=A'*y, where A is the M-by-N sensing or 
%   measurement matrix, y is the M-by-1 measurement vector, and fy is the 
%   N-by-1 result vector.
%   See also FM_CS.
fy = A'*y;
end