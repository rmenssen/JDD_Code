%Integration function for 2D Anomalous Diffusion
%Rebecca Menssen
%Last Updated: 8/21/18

%This function creates the grid necessary to integrate across a range of
%anomalous diffusion parameters using trapz

%%%%%%%%%%INPUTS%%%%%%%%%%
%x1,x2,x3,x4--equally spaced array of D, V, DV and fd values to integrate over
%N--the number of trajectories/points in JDD
%yi--vector of the percentage of trajectories in each JDD bin
%ri--vector of the midpoints of the JDD bins
%dr--the width of each bin in the JDD histogram
%tau--the duration of each trajectory (the time lag*dt)

%%%%%%%%%%OUTPUTS%%%%%%%%%%
%out--structure of probability values for a range of Dalpha and alpha
%values

function[out]=intfuncDV2D(X1,X2,X3,X4, N,yi,ri,dr,tau)
out=zeros(length(X1),length(X2),length(X3),length(X4));
for i = 1:length(X1)
    D=X1(i);
    part1=(dr*ri/(2*D*tau).*exp(-ri.^2/(4*D*tau)));
    for j = 1:length(X2)
         V=X2(j);
        for k=1:length(X3)
%             [i,j,k]
            DV=X3(k);
            a = -(ri.^2+V^2*tau^2)/(4*DV*tau);
            b = ri*V/(2*DV);
            part2=dr*ri/(2*DV*tau).*exp(a).*besseli(0, b);
            parfor l=1:length(X4)
                
                fd=X4(l);
                
                %calculate predicted probabilities
                
                z=fd*part1+...
                    (1-fd)*part2;
                
                %calculate P(JDD|model,parameters)
                denom=prod(sqrt(2*pi*N*z));
                expo=sum((yi-z).^2./(z));
                %sqrt(2*pi*N)/denom*exp(-N/2*expo)
                out(i,j,k,l)=sqrt(2*pi*N)/denom*exp(-N/2*expo);
                
                %fix a slight numerical bug that can occur
                if denom==0
                    %in this case, the exponential is zero, but dividing by zero gives
                    %a NaN, so you have make sure it is correctly recorded as a zero.
                    out(i,j,k,l)=0;
                end
            end
        end
    end
end
end