function [ H ] = gentvop( m,n,gradord )
%GENTVOP generate discrete gradient operator matrix for total variation 
%(TV) regularization.
%   H=GENTVOP(m,n,gradord) returns the (m*(n-1)+(m-1)*n)-by-(m*n) gradient 
%   operator matrix for TV regularization, where m,n is the number of rows 
%   and columns of the two-dimensional image, respectively, gradord is the
%   gradient order of the operator.
%   case 1: if gradient_order == 1
%             TV(x)=\sum_{i=0}^{m-1}\sum_{j=0}^{n-1}(|x_{i,j+1}-x_{i,j}|...
%                   +|x_{i+1,j}-x_{i,j}|).
%               H*x=[[x(i,j+1)-x(i,j)](:);[x(i+1,j)-x(i,j)](:);
%   See also TV, GENTV3DOP.
if nargin<3
    gradord=2; % [default] gradient order
end

h=[-1,1]; % horizontal gradient vector
v=[-1;1]; % vertical gradient vector

% first-order horizontal gradient operator matrix 
% ih1 = 1:m*(n-1);
% jh1 = 1:m*(n-1);
% vh1 = h(1)*ones(1,m*(n-1));
% ih2 = 1:m*(n-1);
% jh2 = m+1:m*n;
% vh2 = h(2)*ones(1,m*(n-1));
% Hh = sparse([ih1 ih2],[jh1 jh2],[vh1 vh2]);

Hh=spalloc(m*(n-1),m*n,2*m*(n-1)); % allocate space for sparse matrix
for ih=1:m
    for jh=1:n-1
        Hh(ih+(jh-1)*m,ih+(jh-1)*m)=h(1);
        Hh(ih+(jh-1)*m,ih+jh*m)=h(2);
    end
end
% Hh=sparse(Hh);

% first-order vertical gradient operator matrix 
Hv=spalloc((m-1)*n,m*n,2*(m-1)*n); % allocate space for sparse matrix
for iv=1:m-1
    for jv=1:n
        Hv(iv+(jv-1)*(m-1),iv+(jv-1)*m)=v(1);
        Hv(iv+(jv-1)*(m-1),iv+1+(jv-1)*m)=v(2);
    end
end
% Hv=sparse(Hv);
H=[Hh;Hv];

if gradord==2
    h3=[1,-2,1];
    h4=[1;-2;1];
    h5=[-1,1;1,-1];
    % H3
    H3=zeros(m*(n-2),m*n);
    for i=1:m
        for j=1:n-2
            H3(i+(j-1)*m,i+(j-1)*m)=h3(1);
            H3(i+(j-1)*m,i+j*m)=h3(2);
            H3(i+(j-1)*m,i+(j+1)*m)=h3(3);
        end
    end
    H3=sparse(H3);

    % H4
    H4=zeros((m-2)*n,m*n);
    for i=1:m-2
        for j=1:n
            H4(i+(j-1)*(m-2),i+(j-1)*m)=h4(1);
            H4(i+(j-1)*(m-2),i+1+(j-1)*m)=h4(2);
            H4(i+(j-1)*(m-2),i+2+(j-1)*m)=h4(3);
        end
    end
    H4=sparse(H4);

    % H5
    H5=zeros((m-1)*(n-1),m*n);
    for i=1:m-1
        for j=1:n-1
            H5(i+(j-1)*(m-1),i+(j-1)*m)=h5(1);
            H5(i+(j-1)*(m-1),i+1+(j-1)*m)=h5(2);
            H5(i+(j-1)*(m-1),i+j*m)=h5(3);
            H5(i+(j-1)*(m-1),i+1+j*m)=h5(4);
        end
    end
    H5=sparse(H5);
    H=[H;H3;H4;H5];
end

end