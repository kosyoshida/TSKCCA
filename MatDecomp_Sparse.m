function [U,D,V,E]=MatDecomp_Sparse(X,c1,c2,K)
% this function aims singular matrix decomposition with L1 reg.
%
% [U,D,V,E]=MatDecomp(X,K)
%
% X:n*p matrix
% K:#components
% c1,c2:upper bound of L1 norm
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
   [X,d,u,v]=SparseSingleMD(X,c1,c2); 
   D(i,i)=d;
   U(:,i)=u;
   V(:,i)=v;
end
E=X;
[d,id]=sort(diag(D),'descend');
D=diag(d);
U=U(:,id);
V=V(:,id);
end


function  [X,d,u,v]=SparseSingleMD(X,c1,c2)
% this function aims rank-1 singular value decomposition with L1 reg.
%
% [X,d,u,v]=SinglMD(X)
%
% X:input
% c1,c2:upper bound of L1 norm
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
    %update u
    argu=X*v_old;
    lmbd1=BinarySearch(argu,c1);
    su=l1_softth(argu,lmbd1);
    u=su/norm(su,2);
    %update v
    argv=X'*u;
    lmbd2=BinarySearch(argv,c2);
    sv=l1_softth(argv,lmbd2);
    v=sv/norm(sv,2);
    if norm(u_old-u,2)<tol && norm(v_old-v,2)<tol
        break;
    end
    u_old=u;
    v_old=v;
end
d=u'*X*v;
X=X-d*u*v';
end
