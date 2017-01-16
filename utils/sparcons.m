function [ sc ] = sparcons( sig,ref,zerodef )
%SPARCONS  Sparse consistancy served as assessment of the consistancies of 
%the zero part and the nonzero part for sparse 1D signals
%   sc=SPARCONS(sig,ref) returns the sparse consistancy of sig with respect
%   to ref.
if nargin < 3
    zerodef = 1e-10; % defined zero
end
sig = sig(:);
ref = ref(:);
Iz  = find(abs(ref)<zerodef);  % zero-valued index
Inz = find(abs(ref)>=zerodef); % nonzero-valued index
sc  = length(find(abs(sig(Iz))<zerodef))-length(Iz)... % consistancy at zero part
      -norm(sig(Inz)-ref(Inz),1)/std(ref(Inz),1);      % consistancy at nonzero part

end

