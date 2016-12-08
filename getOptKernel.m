function optK=getOptKernel(K,d)
% this function returns weighted summation of base kernels
%
% optK=getOptKernel(K,d,r)
%
% inputs:
% d :weight vector of base kernels

[N1,N2,N3]=size(K);

optK=zeros(N1,N2);
for i=1:N3
    optK=optK+K(:,:,i)*d(i);
end
