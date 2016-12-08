function [U,D,V,E]=MatDecomp(X,K)
% this function aims at singular matrix decomposition
%
% [U,D,V,E]=MatDecomp(X,K)
%
% X:n*p matrix
% K:#components
%
% U:n*k matrix
% D:k*k matrix whose diagnal elements are singular values
% V:p*k matrix

[n,p]=size(X);
if exist('K','var')
    k=K;
else
k=min(n,p);
end

U=zeros(n,k);
D=zeros(k,k);
V=zeros(p,k);
for i=1:k
   [X,d,u,v]=SinglMD(X); 
   D(i,i)=d;
   U(:,i)=u;
   V(:,i)=v;
end
E=X;
end


function [X,d,u,v]=SinglMD(X)
% this function aims at rank-1 singular value decomposition 
%
% [X,d,u,v]=SinglMD(X)
%
% X:input
%
% d:first singular value
% u,v:first singular vectors

tol=1e-4;
MaxIt=10000;
%initilize
[~,~,v_old]=svds(X,1);

v_old=v_old/norm(v_old,2);
u_old=X*v_old/norm(X*v_old,2);

for i=1:MaxIt
    u=X*v_old/norm(X*v_old,2);
    v=X'*u_old/norm(X'*u_old,2);
    if norm(u_old-u,2)<tol && norm(v_old-v,2)<tol
        break;
    end
    u_old=u;
    v_old=v;
end
d=u'*X*v;
X=X-d*u*v';
end
