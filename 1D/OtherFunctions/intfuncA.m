%Integration function for 1D Anomalous Diffusion
%Rebecca Menssen
%Last Updated: 4/9/19

%This function creates the grid necessary to integrate across a range of
%anomalous diffusion parameters using trapz

%%%%%%%%%%INPUTS%%%%%%%%%%
%x1,x2--equally spaced array of Dalpha and alpha values to integrate over
%N--the number of trajectories/points in JDD
%yi--vector of the percentage of trajectories in each JDD bin
%ri--vector of the midpoints of the JDD bins
%dr--the width of each bin in the JDD histogram
%tau--the duration of each trajectory (the time lag*dt)

%%%%%%%%%%OUTPUTS%%%%%%%%%%
%out--structure of probability values for a range of Dalpha and alpha
%values

function[out]=intfuncA(X1,X2,N,yi,ri,dr,tau)
out=zeros(length(X1),length(X2));
for i = 1:length(X1)
    i %helpful to know where you are, can comment out
    for j = 1:length(X2) %can make it parfor just be aware of where its being used
        Dalpha=X1(i,j);
        alpha=X2(i,j);
        
        %setting integration limits based on alpha
        if alpha < 0.5
            min=-300^(.5/alpha); %limits on inverse laplace transform
        else
            min=-500;
        end
        
        %calculate predicted probabilities
        fun=@(p) (exp(1i.*p.*tau)).*(1i.*p)^(alpha/2-1)/(2.*pi).*...
            exp(-ri./(sqrt(Dalpha)).*((1i*p)^(alpha/2)));
        z=dr/(Dalpha)^(1/2).*abs(integral(fun,min,-1*min,...
            'ArrayValued',true,'AbsTol',2e-5));
        
        %calculate P(JDD|model,parameters)
        denom=prod(sqrt(2*pi*N*z));
        expo=sum((yi-z).^2./(z));
        out(i,j)=sqrt(2*pi*N)/denom*exp(-N/2*expo);
        
        %fix a slight numerical bug that can occur
        if denom==0
            %in this case, the exponential is zero, but dividing by zero gives
            %a NaN, so you have make sure it is correctly recorded as a zero.
            out(i,j)=0;
        end
         if exp(-N/2.*expo)==0
            %in this case, the exponential is zero, but occasionally things can
            %go a bit haywire.
            out(i,j)=0;
        end
    end
end
end