%1D integration for Bayesian Classifier
%Rebecca Menssen
%This version of the code: 10/24/16

function [prob,value,method]=Integration1D(dbeta,beta,N,yi,ri,dr,tau)
%vector for storing probabilities
prob=zeros(3,1);

%first start with pure diffusion.
%I have found integration to be tricky, so in the future I would like to
%have a method to automatically select the place to start integration, but
%this will have to do for now.

%DIFFUSION INTEGRATION
fprintf('Pure Diffusion\n')
minD=beta(1)-dbeta(1);
if minD <=0
    minD=6e-5;
end
maxD=beta(1)+dbeta(1);
if maxD < minD
    maxD=6.1e-5;
end
%find the length of the integration interval and do the integration.
lengthD = maxD-minD;
x1 = linspace(minD,maxD,1000);
out=intfuncD(x1,N,yi,ri,dr,tau);
prob(1)=1/lengthD*trapz(x1,out);

%DIRECTED MOTION INTEGRATION
fprintf('Directed Motion \n')
minV=beta(2)-dbeta(2);
%unphysical to have V <0
if minV <= 0
    minV=0.01;
end
maxV=beta(2)+dbeta(2);
if maxV <= minV
    maxV=0.02;
end
minDv=beta(3)-dbeta(3);
if minDv < 6e-5
    minDv=6e-5;
end
maxDv=beta(3)+dbeta(3);
if maxDv < minDv
    maxDv=6.1e-5;
end
%integration length
lengthV=maxV-minV;
lengthDv=maxDv-minDv;
x1 = linspace(minV,maxV,100);
x2 = linspace(minDv,maxDv,100);
[X1,X2] = meshgrid(x1,x2);
out=intfuncV(X1,X2,N,yi,ri,dr,tau);
prob(2)=1/lengthV*1/lengthDv*trapz(x2,trapz(x1,out,1));

%ANOMALOUS DIFFUSION INTEGRATION
fprintf('Anomalous Diffusion\n')
minDalpha=beta(4)-dbeta(4);
if minDalpha < 6e-5
    minDalpha=6e-5;
end
maxDalpha=beta(4)+dbeta(4);
if maxDalpha < minDalpha
    maxDalpha=6.1e-5;
end
%Now for the alpha parameter
minalpha=beta(5)-dbeta(5);
%unphysical to have alpha <0
if minalpha <= 0
    minalpha=0.01;
end
maxalpha=beta(5)+dbeta(5);
if maxalpha <= minalpha
    maxalpha=0.02;
end
% if minalpha>1
%     minalpha=0.99;
% end
% if maxalpha > 1
%     maxalpha=1;
% end
%integration length
lengthDalpha=maxDalpha-minDalpha;
lengthalpha=maxalpha-minalpha;
x1 = linspace(minDalpha,maxDalpha,100);
x2 = linspace(minalpha,maxalpha,100);
[X1,X2] = meshgrid(x1,x2);
out=intfuncA(X1,X2,N,yi,ri,dr,tau);
%out=intfuncA2(X1,X2,N,yi,ri,dr,tau);
prob(3)=1/lengthalpha*1/lengthDalpha*trapz(x2,trapz(x1,out,1));

%Choose the best method
%NOW CHOSE BEST METHOD
%normalize
total=sum(prob);
%perentages
%prob=prob/total;
[value,method]=max(prob);
end

