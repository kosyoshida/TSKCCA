function M=Fnorm(Ka,Kb)
% this function returns HSIC values
%
% M=Fnorm(Ka,Kb)
% inputs:
% Ka  :a set of p base-kernels
% Kb  :a set of q bese-kernels
% 
% outputs:
% M  :pxq matrix 

N=size(Ka,1);
if size(Kb,1)~=N
    error('the sample size should be same !')
end

% the number of kernels
p=size(Ka,3);
q=size(Kb,3);

% centering
H=eye(N)-ones(N)./N;
if p<=q
    for i=1:p
        Ka(:,:,i)=H*Ka(:,:,i)*H;
    end
else
    for i=1:q
        Kb(:,:,i)=H*Kb(:,:,i)*H;
    end
end
VecKa=reshape(Ka,[N^2,p]);
VecKb=reshape(Kb,[N^2,q]);
M=VecKa'*VecKb./((N-1)^2);

nan=isnan(M);
M(nan)=0;