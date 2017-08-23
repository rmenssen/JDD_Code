function[out]=intfuncD(x1,N,yi,ri,dr,tau)
out=zeros(length(x1),1);
for i = 1:length(x1)
    D = x1(i);
    z=dr/((pi*D*tau)^(1/2)).*exp(-ri.^2/(4*D*tau));
    denom=prod(sqrt(2.*pi.*N.*z));
    expo=sum((yi-z).^2./(z));
    out(i)=sqrt(2.*pi.*N)/denom.*exp(-N/2.*expo);
    if denom==0
        %in this case, the exponential is zero, but dividing by zero gives
        %a NaN, so you have to right things. 
        out(i)=0;
    end
end
end