%1D integration for Bayesian Classifier
%Rebecca Menssen
%This version of the code: 10/24/16

function [prob,value,method]=Integration3D(dbeta,beta,N,yi,ri,dr,tau)
%vector for storing probabilities
prob=zeros(3,1);

%first start with pure diffusion.
%I have found integration to be tricky, so in the future I would like to
%have a method to automatically select the place to start integration, but
%this will have to do for now.

%DIFFUSION INTEGRATION
fprintf('Pure Diffusion\n')
minD=beta(1)-dbeta(1);
maxD=beta(1)+dbeta(1);
%find the length of the integration interval and do the integration.
lengthD = maxD-minD;
x1 = linspace(minD,maxD,1000);
out=intfuncD3D(x1,N,yi,ri,dr,tau);
prob(1)=1/lengthD*trapz(x1,out);

%DIRECTED MOTION INTEGRATION
fprintf('Directed Motion \n')
minV=beta(2)-dbeta(2);
minDv=beta(3)-dbeta(3);
maxV=beta(2)+dbeta(2);
maxDv=beta(3)+dbeta(3);

%integration length
lengthV=maxV-minV;
lengthDv=maxDv-minDv;
x1 = linspace(minV,maxV,100);
x2 = linspace(minDv,maxDv,100);
[X1,X2] = meshgrid(x1,x2);
out=intfuncV3D(X1,X2,N,yi,ri,dr,tau);
prob(2)=1/lengthV*1/lengthDv*trapz(x2,trapz(x1,out,1));

%ANOMALOUS DIFFUSION INTEGRATION
fprintf('Anomalous Diffusion\n')
minDalpha=beta(4)-dbeta(4);
maxDalpha=beta(4)+dbeta(4);
minalpha=beta(5)-dbeta(5);
maxalpha=beta(5)+dbeta(5);
lengthDalpha=maxDalpha-minDalpha;
lengthalpha=maxalpha-minalpha;
x1 = linspace(minDalpha,maxDalpha,100);
x2 = linspace(minalpha,maxalpha,100);
[X1,X2] = meshgrid(x1,x2);
out=intfuncA(X1,X2,N,yi,ri,dr,tau);
prob(3)=1/lengthalpha*1/lengthDalpha*trapz(x2,trapz(x1,out,1));

%Choose the best method
%NOW CHOSE BEST METHOD
%normalize
total=sum(prob);
%perentages
prob=prob/total;
[value,method]=max(prob);
end

