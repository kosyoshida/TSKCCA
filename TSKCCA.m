function [r,U,V,As,Bs,M]=TSKCCA(Ks,Kc,param)
% this function performs Two-Stage KCCA
%
% [r,U,V,As,Bs,M]=TSKCCA(Ks,Kc,param)
%
% inputs:
% Ks :series of sub-kernel matrices (NxNxMs)
% Kc :series of sub-kernel matrices (NxNxMc)
% param :
%   No_Prj :the number of canonical variables
%   No_Cmb :the number of combinations of kernels
%   c1,c2  :sparse parameter for kernel selection
%   kappa  :regularization parameter for KCCA
% outputs:
% r  :correlation coefficients
% U  :weight of sub-kernels (MsxNo_Cmb)
% V  :weight of sub-kernels (McxNo_Cmb)
% As :coefficient of KCCA (NxNo_Prj
% 
if ~exist('param')
    param.No_Prj=1;
    param.No_Cmb=1;
    param.c1=1.2;
    param.c2=1.2;
    param.kappa=0.02;
end
if ~isfield(param,'No_Prj')
    param.No_Prj=1;
end
if ~isfield(param,'No_Cmb')
    param.No_Cmb=1;
end
if ~isfield(param,'c1')
    param.c1=1.2;
end
if ~isfield(param,'c2')
    param.c2=1.2;
end
if ~isfield(param,'kappa')
    param.kappa=0.02;
end

r=zeros(param.No_Prj,param.No_Cmb);
As=cell(1,param.No_Cmb);
Bs=cell(1,param.No_Cmb);

% HSIC
M=Fnorm(Ks,Kc);
[U,~,V,~]=MatDecomp_Sparse_nonneg(M,param.c1,param.c2,param.No_Cmb);

for Cmb=1:param.No_Cmb
    optKs=getOptKernel(Ks,U(:,Cmb));
    optKc=getOptKernel(Kc,V(:,Cmb));
    [~,A,B]=KCCA(optKs,optKc,param);
    As{Cmb}=A;
    Bs{Cmb}=B;
    r(:,Cmb)=diag(corr(optKs*A,optKc*B));
end