%Integration function for 2D Anomalous Diffusion
%Rebecca Menssen
%Last Updated: 8/21/18

%This function creates the grid necessary to integrate across a range of
%anomalous diffusion parameters using trapz

%%%%%%%%%%INPUTS%%%%%%%%%%
%x1,x2,x3,x4,x5--equally spaced array of V, DV, Dalpha, alpha,  and fd values to integrate over
%N--the number of trajectories/points in JDD
%yi--vector of the percentage of trajectories in each JDD bin
%ri--vector of the midpoints of the JDD bins
%dr--the width of each bin in the JDD histogram
%tau--the duration of each trajectory (the time lag*dt)

%%%%%%%%%%OUTPUTS%%%%%%%%%%
%out--structure of probability values for a range of Dalpha and alpha
%values

function[out]=intfuncVA2D(X1,X2,X3,X4,X5,N,yi,ri,dr,tau)
out=zeros(length(X1),length(X2),length(X3),length(X4),length(X5));
for i = 1:length(X1)
    V=X1(i);
    for j = 1:length(X2)
        D=X2(j);
        a = -(ri.^2+V^2*tau^2)/(4*D*tau);
        b = ri*V/(2*D);
        part1=(dr*ri/(2*D*tau).*exp(a).*besseli(0, b));
        for k=1:length(X3)
            Dalpha=X3(k);
            for l=1:length(X4)
                %                 [i,j,k,l]
                alpha=X4(l);
                %setting integration limits based on alpha
                if alpha < 0.5
                    min=-300^(.5/alpha); %limits on inverse laplace transform
                else
                    min=-500;
                end
                fun=@(p) (exp(1i.*p.*tau)).*(1i.*p)^(alpha-1)/(2.*pi).*...
                    (besselk(0,ri./(sqrt(Dalpha)).*((1i*p)^(alpha/2))));
                part2=dr*ri/(Dalpha).*abs(integral(fun,min-1i*1e-6,...
                    -1*min-1i*1e-6,'ArrayValued',true,'AbsTol',1e-6,'RelTol',1e-3));
                parfor m=1:length(X5)
                    fd=X5(m);
                    
                    z=fd*part1+...
                        (1-fd)*part2;
                    
                    %calculate P(JDD|model,parameters)
                    denom=prod(sqrt(2*pi*N*z));
                    expo=sum((yi-z).^2./(z));
                    %sqrt(2*pi*N)/denom*exp(-N/2*expo)
                    out(i,j,k,l,m)=sqrt(2*pi*N)/denom*exp(-N/2*expo);
                    
                    %fix a slight numerical bug that can occur
                    if denom==0
                        %in this case, the exponential is zero, but dividing by zero gives
                        %a NaN, so you have make sure it is correctly recorded as a zero.
                        out(i,j,k,l,m)=0;
                    end
                end
            end
        end
    end
end
end