function[out]=intfuncV3D(X1,X2,N,yi,ri,dr,tau)
out=zeros(length(X1),length(X2));
for i = 1:length(X1)
    for j = 1:length(X2)
        V=X1(i,j);
        D=X2(i,j);
        a = -(ri.^2+V^2*tau^2)/(4*D*tau);
        b = ri*V/(2*D);
        z=dr*ri.^2*4*pi/((4*pi*D*tau)^(3/2)).*exp(a).*sinh(b)./b;
        denom=prod(sqrt(2*pi*N*z));
        expo=sum((yi-z).^2./(z));
        if exp(-N/2*expo)==0
            out(i,j)=0;
        else
            out(i,j)=sqrt(2*pi*N)/denom*exp(-N/2*expo);
            if denom==0
                %in this case, the exponential is zero, but dividing by zero gives
                %a NaN, so you have to right things.
                out(i)=0;
            end
        end
    end
end
end