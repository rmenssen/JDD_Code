%Integration function for 1D Directed Diffusion
%Rebecca Menssen
%Last Updated: 4/9/19

%This function creates the grid necessary to integrate across a range of
%directed diffusion parameters using trapz

%%%%%%%%%%INPUTS%%%%%%%%%%
%x1,x2--equally spaced array of V and Dv values to integrate over
%N--the number of trajectories/points in JDD
%yi--vector of the percentage of trajectories in each JDD bin
%ri--vector of the midpoints of the JDD bins
%dr--the width of each bin in the JDD histogram
%tau--the duration of each trajectory (the time lag*dt)

%%%%%%%%%%OUTPUTS%%%%%%%%%%
%out--structure of probability values for a range of V and Dv values


function[out]=intfuncV(X1,X2,N,yi,ri,dr,tau)
out=zeros(length(X1),length(X2));
for i = 1:length(X1)
    for j = 1:length(X2)
        V=X1(i,j);
        D=X2(i,j);
        
        %calculate predicted probabilities
        z = dr/((4*pi*D*tau)^(1/2)).*exp(-(ri.^2+V^2*tau^2)/(4*D*tau)+ri*V/(2*D))+...
            dr/((4*pi*D*tau)^(1/2)).*exp(-(ri.^2+V^2*tau^2)/(4*D*tau)-ri*V/(2*D));
        
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