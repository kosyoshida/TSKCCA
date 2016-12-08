function su=l1_softth(u,lmbd)
%%% this function aims at soft-threshold 
p=u>lmbd;
n=u<-lmbd;
su=zeros(length(u),1);
su(p)=u(p)-lmbd;
su(n)=u(n)+lmbd;
