function [x,z]=Data2(N,dimx,dimz,ep)
% this function returns data. 

% create x and z
x=rand(N,dimx)-0.5*ones(N,dimx);
z=rand(N,dimz)-0.5*ones(N,dimx);

% compute fs
f1=x(:,1);
f2=x(:,2).^2;
f3=abs(x(:,3));
f4=exp(-x(:,4).^2);
f5=sin(pi*x(:,5)/2);
f6=sig(x(:,6),0.2);

% all fs are standardized before computing zs
f1=normal(f1);f2=normal(f2);f3=normal(f3);
f4=normal(f4);f5=normal(f5);f6=normal(f6);

% generate multiple relevant components
z(:,1)=f1+f4+ep*randn(N,1);
z(:,2)=f2+f5+ep*randn(N,1);
z(:,3)=f3+f6+ep*randn(N,1);

function y=sig(x,sigma)
y=1./(1+exp(-x/sigma));
end