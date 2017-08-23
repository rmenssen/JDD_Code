function[out]=intfuncA(X1,X2,N,yi,ri,dr,tau)
out=zeros(length(X1),length(X2));
for i = 1:length(X1)
    for j = 1:length(X2)
        %i;
        %         %j;
        Dalpha=X1(i,j);
        alpha=X2(i,j);
        if alpha < 0.5
            min=-500^(.5/alpha); %limits on inverse laplace transform
        else
            min=-500;
        end
        fun=@(p) (exp(1i.*p.*tau)).*(1i.*p)^(alpha/2-1)/(2.*pi).*exp(-ri./(sqrt(Dalpha)).*((1i*p)^(alpha/2)));
        z=dr/(Dalpha)^(1/2).*abs(integral(fun,min,-1*min,'ArrayValued',true,'AbsTol',2e-5));
        denom=prod(sqrt(2*pi*N*z));
        expo=sum((yi-z).^2./(z));
        out(i,j)=sqrt(2*pi*N)/denom*exp(-N/2*expo);
        if denom==0
            %in this case, the exponential is zero, but dividing by zero gives
            %a NaN, so you have to right things.
            out(i,j)=0;
        end
    end
end
end