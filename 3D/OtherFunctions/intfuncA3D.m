function[out]=intfuncA3D(X1,X2,N,yi,ri,dr,tau)
out=zeros(length(X1),length(X2));
for i = 1:length(X1)
    for j = 1:length(X2)
        Dalpha=X1(i,j);
        alpha=X2(i,j);
        if alpha < 0.5 && alpha>0.4
            min=-700; %limits on inverse laplace transform
        elseif alpha<= 0.4 && alpha>0.3
            min=-900;
        elseif alpha<=0.3
            min=-1100;
        else
            min=-700;
        end
        fun=@(p) (exp(1i.*p.*tau)).*(1i.*p)^(alpha-1)/(2.*pi).*(exp(-ri./(sqrt(Dalpha)).*((1i*p)^(alpha/2))));
        intsize=200;
        numint=ceil(abs(min+100)/intsize);
        z=zeros(size(ri));
        for k=1:numint
            z=z+2*dr*ri/(Dalpha).*abs(integral(fun,min+(k-1)*intsize-1i*1e-6,min+k*(intsize)-1i*1e-6,...
                'ArrayValued',true,'AbsTol',1e-6,'RelTol',1e-3));
        end
        z=z+dr*ri/(Dalpha).*abs(integral(fun,-100-1i*1e-6,100-1i*1e-6,...
            'ArrayValued',true,'AbsTol',1e-6,'RelTol',1e-3));
        denom=prod(sqrt(2*pi*N*z));
        expo=sum((yi-z).^2./(z));
        out(i,j)=sqrt(2*pi*N)/denom*exp(-N/2*expo);
        if denom==0
            %in this case, the exponential is zero, but dividing by zero gives
            %a NaN, so you have to right things.
            out(i)=0;
        end
    end
end
end