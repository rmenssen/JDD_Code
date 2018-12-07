%Integration function for 2D Anomalous Diffusion
%Rebecca Menssen
%Last Updated: 8/21/18

%This function creates the grid necessary to integrate across a range of
%anomalous diffusion parameters using trapz

%%%%%%%%%%INPUTS%%%%%%%%%%
%x1,x2,x3--equally spaced array of D1,D2 and fd values to integrate over
%N--the number of trajectories/points in JDD
%yi--vector of the percentage of trajectories in each JDD bin
%ri--vector of the midpoints of the JDD bins
%dr--the width of each bin in the JDD histogram
%tau--the duration of each trajectory (the time lag*dt)

%%%%%%%%%%OUTPUTS%%%%%%%%%%
%out--structure of probability values for a range of Dalpha and alpha
%values

function[out]=intfuncDD2D(X1,X2,X3,N,yi,ri,dr,tau)
out=zeros(length(X1),length(X2),length(X3));
for i = 1:length(X1)
    D1=X1(i);
    part1=(dr*ri/(2*D1*tau).*exp(-ri.^2/(4*D1*tau)));
    for j = 1:length(X2)
        [i,j]
        D2=X2(j);
        part2=(dr*ri/(2*D2*tau).*exp(-ri.^2/(4*D2*tau)));
        parfor k=1:length(X3)
%             D1=X1(i,j,k);
%             D2=X2(i,j,k);
%             fd=X3(i,j,k);
            fd=X3(k);
            
            z=fd*part1+(1-fd)*part2;
          
            %calculate P(JDD|model,parameters)
            denom=prod(sqrt(2*pi*N*z));
            expo=sum((yi-z).^2./(z));
            %sqrt(2*pi*N)/denom*exp(-N/2*expo)
            out(i,j,k)=sqrt(2*pi*N)/denom*exp(-N/2*expo);
            
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