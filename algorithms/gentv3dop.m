function [ H, intraH, interH ] = gentv3dop( row,col,nChannel )
%GENTV3DOP generate discrete gradient operator matrix for three dimensional
%total variation (TV) regularization.
%   H -- gradient operator
%   See also TV3D, GENTVOP.

if nargin <3
    nChannel = 1;
end

h1=[-1,1];
h2=[-1;1];
h3=[1,-2,1];
h4=[1;-2;1];
h5=[-1,1;1,-1];
% H1
H1=zeros(col*(row-1),col*row);
for x=1:col
    for y=1:row-1
        H1(x+(y-1)*col,x+(y-1)*col)=h1(1);
        H1(x+(y-1)*col,x+y*col)=h1(2);
    end
end
H1=sparse(H1);
% H2
H2=zeros((col-1)*row,col*row);
for x=1:col-1
    for y=1:row
        H2(x+(y-1)*(col-1),x+(y-1)*col)=h2(1);
        H2(x+(y-1)*(col-1),x+1+(y-1)*col)=h2(2);
    end
end
H2=sparse(H2);
% H3
H3=zeros(col*(row-2),col*row);
for x=1:col
    for y=1:row-2
        H3(x+(y-1)*col,x+(y-1)*col)=h3(1);
        H3(x+(y-1)*col,x+y*col)=h3(2);
        H3(x+(y-1)*col,x+(y+1)*col)=h3(3);
    end
end
H3=sparse(H3);
% H4
H4=zeros((col-2)*row,col*row);
for x=1:col-2
    for y=1:row
        H4(x+(y-1)*(col-2),x+(y-1)*col)=h4(1);
        H4(x+(y-1)*(col-2),x+1+(y-1)*col)=h4(2);
        H4(x+(y-1)*(col-2),x+2+(y-1)*col)=h4(3);
    end
end
H4=sparse(H4);
% H5
H5=zeros((col-1)*(row-1),col*row);
for x=1:col-1
    for y=1:row-1
        H5(x+(y-1)*(col-1),x+(y-1)*col)=h5(1);
        H5(x+(y-1)*(col-1),x+1+(y-1)*col)=h5(2);
        H5(x+(y-1)*(col-1),x+y*col)=h5(3);
        H5(x+(y-1)*(col-1),x+1+y*col)=h5(4);
    end
end
H5=sparse(H5);

% H=[H1;H2;H3;H4;H5];
H = [H1;H2];
% H = kron(eye(nChannel),H);
H = repmat(H, 1 ,nChannel);
intraH = H;

H6 = zeros(col*row*(nChannel-1), col*row*nChannel);
for x = 1:size(H6,1)
    H6(x,x)=-1;
    H6(x,x+col*row)=1;
end
interH = H6;

H = [H;H6];

end