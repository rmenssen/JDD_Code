%Integration function for 3D Diffusion
%Rebecca Menssen
%Last Updated: 8/30/17

%This function creates the grid necessary to integrate across a range of
%diffusion parameters using trapz


%%%%%%%%%%INPUTS%%%%%%%%%%
%x1--equally spaced array of D values to integrate over
%N--the number of trajectories/points in JDD
%yi--vector of the percentage of trajectories in each JDD bin
%ri--vector of the midpoints of the JDD bins
%dr--the width of each bin in the JDD histogram
%tau--the duration of each trajectory (the time lag*dt)

%%%%%%%%%%OUTPUTS%%%%%%%%%%
%out--structure of probability values for a range of D values

function[out]=intfuncD3D(x1,N,yi,ri,dr,tau)
out=zeros(length(x1),1);
for i = 1:length(x1)
    D = x1(i);
    
    %calculate predicted probabilities
    z=dr*ri.^2/(2*sqrt(pi)*(D*tau)^(3/2)).*exp(-ri.^2/(4*D*tau));
    
    %calculate P(JDD|model,parameters)
    denom=prod(sqrt(2.*pi.*N.*z));
    expo=sum((yi-z).^2./(z));
    out(i)=sqrt(2.*pi.*N)/denom.*exp(-N/2.*expo);
    
    %fix a slight numerical bug that can occur
    if denom==0
        %in this case, the exponential is zero, but dividing by zero gives
        %a NaN, so you have make sure it is correctly recorded as a zero. 
        out(i)=0;
    end
end
end