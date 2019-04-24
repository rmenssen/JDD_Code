%Integration function for 1D Double Diffusion
%Rebecca Menssen
%Last Updated: 4/9/19

%This function creates the grid necessary to integrate across a range of
%diffusion parameters using trapz

%%%%%%%%%%INPUTS%%%%%%%%%%
%x1,x2,x3--equally spaced array of D1,D2, and fd values to integrate over
%N--the number of trajectories/points in JDD
%yi--vector of the percentage of trajectories in each JDD bin
%ri--vector of the midpoints of the JDD bins
%dr--the width of each bin in the JDD histogram
%tau--the duration of each trajectory (the time lag*dt)

%%%%%%%%%%OUTPUTS%%%%%%%%%%
%out--structure of probability values for a range of D values

function[out]=intfuncDD(x1,x2,x3,N,yi,ri,dr,tau)
out=zeros(length(x1),length(x2),length(x3));
for i = 1:length(x1)
    D1 = x1(i);
    z1=dr/((pi*D1*tau)^(1/2)).*exp(-ri.^2/(4*D1*tau));
    for j=1:length(x2)
        [i,j]
        D2=x2(j);
        z2=dr/((pi*D2*tau)^(1/2)).*exp(-ri.^2/(4*D2*tau));
        parfor k=1:length(x3)
            [i,j]
            fd=x3(k);
            z=fd*z1+(1-fd)*z2;
            %calculate P(JDD|model,parameters)
            denom=prod(sqrt(2.*pi.*N.*z));
            expo=sum((yi-z).^2./(z));
            out(i,j,k)=sqrt(2.*pi.*N)/denom.*exp(-N/2.*expo);
            %fix a slight numerical bug that can occur
            if denom==0
                %in this case, the exponential is zero, but dividing by zero gives
                %a NaN, so you have make sure it is correctly recorded as a zero.
                out(i,j,k)=0;
            end
        end
    end
end
end




