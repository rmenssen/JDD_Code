%3D integration for Bayesian Classifier
%Rebecca Menssen
%This version of the code: 4/24/19

%%%%%%%%%%INPUTS%%%%%%%%%%
%dbeta--defines the bounds on integration. Have been using 2 times the
%       standard deviation of the bootstrapped parameters. 
%beta--median or mean parameter values. Center of integration window. 
%N--the number of trajectories/points in JDD
%yi--vector of the percentage of trajectories in each JDD bin
%ri--vector of the midpoints of the JDD bins
%dr--the width of each bin in the JDD histogram
%tau--the duration of each trajectory (the time lag*dt)

%%%%%%%%%%OUTPUTS%%%%%%%%%%
%prob--vector giving normalized probabilities in the order pure diffusion,
%directed diffusion, anomalous diffusion
%value--maximal probability
%method--index giving the method of the maximal probability

function [prob,value,method]=Integration3D(dbeta,beta,N,yi,ri,dr,tau)
%Vector for storing probabilities
prob=zeros(3,1);

%DIFFUSION INTEGRATION
%define ranges on integration. 
minD=beta(1)-dbeta(1);
maxD=beta(1)+dbeta(1);

%Check against negative values. They tend to really mess things up and
%should be ignored.
if minD <=0
    minD=1e-10;
end
if maxD <=0
    maxD=1e-9;
end

%Find length of interval, find probabilities at different values of
%parameters, and then do the integration. 
lengthD = maxD-minD;
x1 = linspace(minD,maxD,1000);
out=intfuncD3D(x1,N,yi,ri,dr,tau);
prob(1)=1/lengthD*trapz(x1,out);

%DIRECTED MOTION INTEGRATION
%define ranges on integration.
minV=beta(2)-dbeta(2);
maxV=beta(2)+dbeta(2);
minDv=beta(3)-dbeta(3);
maxDv=beta(3)+dbeta(3);

%Checks against negative values. 
if minV <= 0
    minV=0.01;
end
if maxV <= 0
    maxV=0.02;
end
if minDv < 0
    minDv=1e-10;
end
if maxDv < 0
    maxDv=1e-9;
end

%Find length of interval, find probabilities at different values of
%parameters, and then do the integration. 
lengthV=maxV-minV;
lengthDv=maxDv-minDv;
x1 = linspace(minV,maxV,100);
x2 = linspace(minDv,maxDv,100);
[X1,X2] = meshgrid(x1,x2);
out=intfuncV3D(X1,X2,N,yi,ri,dr,tau);
prob(2)=1/lengthV*1/lengthDv*trapz(x2,trapz(x1,out,1));

%ANOMALOUS DIFFUSION INTEGRATION
%here we do integrate for all values of alpha. Since alpha should be less
%than 1 for the method, you can place constraints on max or min alpha, or
%only integrate if alpha<1 (see 1D integration for examples of how to do
%this). 

minDalpha=beta(4)-dbeta(4);
maxDalpha=beta(4)+dbeta(4);
minalpha=beta(5)-dbeta(5);
maxalpha=beta(5)+dbeta(5);

%checks against negative values
if minDalpha <=0
    minDalpha=1e-10;
end
if maxDalpha <=0
    maxDalpha=1e-9;
end
if minalpha <= 0
    minalpha=0.01;
end

if maxalpha <= 0
    maxalpha=0.02;
end

%Find length of interval, find probabilities at different values of
%parameters, and then do the integration. 
lengthDalpha=maxDalpha-minDalpha;
lengthalpha=maxalpha-minalpha;
x1 = linspace(minDalpha,maxDalpha,100);
x2 = linspace(minalpha,maxalpha,100);
[X1,X2] = meshgrid(x1,x2);
out=intfuncA3D(X1,X2,N,yi,ri,dr,tau);
prob(3)=1/lengthalpha*1/lengthDalpha*trapz(x2,trapz(x1,out,1));

%Normalize Probabilities
prob=prob/sum(prob);

%Select the best method
[value,method]=max(prob);
end

