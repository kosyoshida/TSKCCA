function lmbd=BinarySearch(arg,c,nn)
% this function returns optimal lmbd by binary search
%
% lmbd=BinarySearch(arg,c,nn)
%
% |arg|_1<=c
% |arg|_2=1
%
% Input
% arg :target vector
% c   :upper bound of L1-norm
% nn  :non-negative flag. arg must be nonnegative if nn==1
if nargin==2
    nn=0;
end

if(norm(arg,2)==0 || sum(abs(arg/norm(arg,2)))<=c)
    lmbd=0;
else
    lmbd1=0;
    lmbd2=max(abs(arg))-(1e-5);
    for i=1:150
        su=l1_softth(arg,(lmbd1+lmbd2)/2);
        if nn==1
            su(su<=0)=0;
        end
        if(sum(abs(su/norm(su,2)))<c || norm(su,2)==0)
            lmbd2=(lmbd1+lmbd2)/2;
        else
            lmbd1=(lmbd1+lmbd2)/2;
        end
        if (lmbd2-lmbd1)<1e-6
            break;
        end
    end
    lmbd=(lmbd1+lmbd2)/2;
end
end

        
