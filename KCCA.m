function [r,A,B,Ua,Vb]=KCCA(Ka,Kb,param)
% this function performs KCCA
%
% [r,A,B,Ua,Vb]=KCCA(Ka,Kb,param)
% inputs:
% Ka,Kb  : Kernel Matrices. Their need to be centered.
% param  :'No_Prj' # of projection vector
%        :'No_Cmb' # of combination vector
%        :'kappa'  regularization parametr for kernel CCA
% outputs:
% A,B    :projection matrices. 
% r      :Kernel canonical correlation
% U,V    :Kernel canonical scores

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

N=size(Ka,1);
if size(Kb,1)~=N
    error('the sample size should be same !')
end

% eigen value decomposition
L=[zeros(N),Ka'*Kb;Kb'*Ka,zeros(N)];
R=[(Ka+N*param.kappa*eye(N)./2)^2,zeros(N);zeros(N),(Kb+N*param.kappa*eye(N)./2)^2];
[V,~]=eigs(L,R,param.No_Prj);

A=V(1:N,:);
B=V(N+1:end,:);
Ua=Ka*A;
Vb=Kb*B;

r=diag(corr(Ua,Vb));
