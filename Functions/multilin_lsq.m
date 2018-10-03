function [H] = multilin_lsq(X,Y)
%[H] = multilin_lsq(X,Y) 
%multiple linear regression, computes the coefficients to a multiple linear
%regression, of the form Y(t)=H_1*X_1(t)+H_2*X_2(t)+H_3*X_3(t)+error, in 
%the least squares sense; takes as inputs: 
            %X: n timesteps by k inputs matrix
            %Y: n timesteps by 1 output column vector
%output:    %H: k input column vector, estimated coefficients in a least
            %squares sense

if 0, %uncomment to test machinery
    n=1000;            
    X=[randn(n,1) randn(n,1) randn(n,1)]; %build inputs
    Y=5*X(:,1)+2.*X(:,2)-4.*X(:,3)+0.6*randn(n,1);
end

n=size(X,1);
k=size(X,2);
xx=nan(k,k); %input covariances
for aa=1:k %iterates through the number of x inputs
    for bb=aa:k        
        xx(aa,bb)=sum(X(:,aa).*X(:,bb));
        xx(bb,aa)=xx(aa,bb);
    end
end
xx=xx/n; %normalize by number of data points


yxa=repmat(Y,1,k).*X;
yx=sum(yxa)/n; 

H=xx\yx';


end

