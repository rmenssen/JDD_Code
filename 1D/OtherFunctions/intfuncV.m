function[out]=intfuncV(X1,X2,N,yi,ri,dr,tau)
out=zeros(length(X1),length(X2));
for i = 1:length(X1)
    for j = 1:length(X2)
        V=X1(i,j);
        D=X2(i,j);
        z = dr/((4*pi*D*tau)^(1/2)).*exp(-(ri.^2+V^2*tau^2)/(4*D*tau)+ri*V/(2*D));
        denom=prod(sqrt(2*pi*N*z));
        expo=sum((yi-z).^2./(z));
        if exp(-N/2*expo)==0
            out(i,j)=0;
        else
            out(i,j)=sqrt(2*pi*N)/denom*exp(-N/2*expo);
        end
        if denom==0
            %in this case, the exponential is zero, but dividing by zero gives
            %a NaN, so you have to right things.
            out(i,j)=0;
        end
    end
end
end