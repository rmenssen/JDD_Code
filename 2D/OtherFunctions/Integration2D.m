function [prob,method,value]=Integration2D(averages,N,yi,ri,dr,tau,beta,M)
%vector for storing probabilities
prob=zeros(3,1);
%first calculate dbeta for each one
dbeta=10*abs(averages-beta);

%DIFFUSION PROBABILITIES
%find max and min of the region
minD=beta(1)-dbeta(1);
if minD < .006
    minD=.006;
end
maxD=beta(1)+dbeta(1);
if maxD < .0061
    maxD=.0061;
end
%find the length of the integration interval
lengthD = maxD-minD;
x1 = linspace(minD,maxD,100);
out=intfuncD(x1,N,yi,ri,dr,tau);
prob(1)=1/lengthD*trapz(x1,out);
%prob


%DIRECTED MOTION INTEGRATION
minV=beta(2)-dbeta(2);
%unphysical to have V <0
if minV <= 0
    minV=0.01;
end
maxV=beta(2)+dbeta(2);
if maxV <= 0
    maxV=0.02;
end
%avoid numbers too small for matlab
minkv=beta(3)-dbeta(3);
if minkv < 3e-4
    minkv=3e-4;
end
maxkv=beta(3)+dbeta(3);
if maxkv < 3.1e-4
    maxkv=3.1e-4;
end
%integration length
lengthV=maxV-minV;
lengthkv=maxkv-minkv;
x1 = linspace(minV,maxV,100);
x2 = linspace(minkv,maxkv,100);
[X1,X2] = meshgrid(x1,x2);
out=intfuncV(X1,X2,N,yi,ri,dr,tau,M);
prob(2)=1/lengthV*1/lengthkv*trapz(x2,trapz(x1,out,1));
%prob

%DOUBLE DIFFUSION
minD1=beta(5)-dbeta(5);
if minD1 < .006
    minD1=.006;
end
maxD1=beta(5)+dbeta(5);
if maxD1 < .0061
    maxD1=.0061;
end
minD2=beta(6)-dbeta(6);
if minD2 < .006
    minD2=.006;
end
maxD2=beta(6)+dbeta(6);
if maxD2 < .0061
    maxD2=.0061;
end
minfd=beta(4)-dbeta(4);
%cannot have a negative fd
if minfd <= 0
    minfd=0;
end
maxfd=beta(4)+dbeta(4);
if maxfd <= 0
    maxfd=0.01;
end
%cannot have fd >1 unphysical
if maxfd > 1
    maxfd=1;
end
%integrate
lengthDD=(maxD1-minD1)*(maxD2-minD2)*(maxfd-minfd);
x1 = linspace(minD1,maxD1,100);
x2 = linspace(minD2,maxD2,100);
x3 = linspace(minfd,maxfd,100);
[X1,X2,X3] = ndgrid(x1,x2,x3);
out = intfuncDD(X1,X2,X3,N,yi,ri,dr,tau);
prob(3) =1/lengthDD*trapz(x3,trapz(x2,trapz(x1,out,1),2),3);
prob

%NOW CHOSE BEST METHOD
%normalize
total=sum(prob);
%perentages
prob=prob/total;
[value,method]=max(prob);
%as you so decide can analyze the probabilities to choose a cut off
%but for double diffusion, note that you could have a lower probability in
%that it might arrive at two (very similar) diffusion constants,
%essentially the same as a single diffusion constant, so you might need to
%check some of those cases by hand.
end

function[out]=intfuncD(x1,N,yi,ri,dr,tau)
out=zeros(length(x1),1);
for i = 1:length(x1)
    D = x1(i);
    z=dr.*ri./(2.*D.*tau).*exp(-ri.^2/(4.*D.*tau));
    denom=prod(sqrt(2.*pi.*N.*z));
    expo=sum((yi-z).^2./(z));
    out(i)=sqrt(2.*pi.*N)/denom.*exp(-N/2.*expo);
end
end

function[out]=intfuncV(X1,X2,N,yi,ri,dr,tau,M)
out=zeros(length(X1),length(X2));
for i = 1:length(X1)
    for j = 1:length(X2)
        V=X1(i,j);
        kv=X2(i,j);
        z=dr*ri./(kv*M).*exp(-(ri.^2+V.^2*tau^2)/(2*M.*kv)).*besseli(0,ri.*V*tau/(M*kv));
        denom=prod(sqrt(2*pi*N*z));
        expo=sum((yi-z).^2./(z));
        out(i,j)=sqrt(2*pi*N)/denom*exp(-N/2*expo);
    end
end
end

function[out]=intfuncDD(x1,x2,x3,N,yi,ri,dr,tau)
out=zeros(length(x1),length(x2),length(x3));
for i = 1:length(x1)
    for j = 1:length(x2)
        for k = 1:length(x3)
            D1=x1(i,j,k);
            D2=x2(i,j,k);
            fd=x3(i,j,k);
            z=fd*dr*ri/(2*D1*tau).*exp(-ri.^2/(4*D1*tau))+...
                (1-fd)*dr*ri/(2*D2*tau).*exp(-ri.^2/(4*D2*tau));
            denom=prod(sqrt(2*pi*N*z));
            expo=sum((yi-z).^2./(z));
            out(i,j,k)=sqrt(2*pi*N)/denom*exp(-N*expo/2);
        end
    end
end
end
