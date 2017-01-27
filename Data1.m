function [x,z]=Data1(N,dimx,dimz,ep)
% this function returns data. 

% create x and z
x=rand(N,dimx)-0.5*ones(N,dimx);
z=rand(N,dimz)-0.5*ones(N,dimz);

% single non-linear association
f1=x(:,1).^2;

% all fs are standardized before computing zs
f1=normal(f1);

% generate relevant components
z(:,1)=f1+ep*randn(N,1);


