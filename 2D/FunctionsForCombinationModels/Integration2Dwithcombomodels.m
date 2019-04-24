%2D integration for Bayesian Classifier
%Rebecca Menssen
%This version of the code: 3/20/2019

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

function [prob,value,method]=Integration2Dwithcombomodels(dbeta,beta,N,yi,ri,dr,tau)
%Vector for storing probabilities
prob=zeros(7,1);

%DIFFUSION INTEGRATION
%define ranges on integration.
minD=beta(1)-dbeta(1); maxD=beta(1)+dbeta(1);

%Check against negative values. They tend to really mess things up and
%should be ignored.
if minD <=0, minD=1e-10; end
if maxD <=0, maxD=1e-9; end

%Find length of interval, find probabilities at different values of
%parameters, and then do the integration.
lengthD = maxD-minD;
x1 = linspace(minD,maxD,1000);
out=intfuncD2D(x1,N,yi,ri,dr,tau);
prob(1)=1/lengthD*trapz(x1,out);

%DIRECTED MOTION INTEGRATION
%define ranges on integration.
minV=beta(2)-dbeta(2); maxV=beta(2)+dbeta(2);
minDv=beta(3)-dbeta(3); maxDv=beta(3)+dbeta(3);

%Checks against negative values.
if minV <= 0, minV=0.01; end
if maxV <= 0, maxV=0.02; end
if minDv < 0, minDv=1e-10; end
if maxDv < 0, maxDv=1e-9; end

%Find length of interval, find probabilities at different values of
%parameters, and then do the integration.
lengthV=maxV-minV; lengthDv=maxDv-minDv;
x1 = linspace(minV,maxV,100); x2 = linspace(minDv,maxDv,100);
[X1,X2] = meshgrid(x1,x2);
out=intfuncV2D(X1,X2,N,yi,ri,dr,tau);
prob(2)=1/lengthV*1/lengthDv*trapz(x2,trapz(x1,out,1));

%ANOMALOUS DIFFUSION INTEGRATION

minDalpha=beta(4)-dbeta(4); maxDalpha=beta(4)+dbeta(4);
minalpha=beta(5)-dbeta(5); maxalpha=beta(5)+dbeta(5);

%checks against negative values
if minDalpha <=0, minDalpha=1e-10;end
if maxDalpha <=0, maxDalpha=1e-9; end
if minalpha <= 0, minalpha=0.01; end
if maxalpha <= 0, maxalpha=0.02; end

%option 1, look at values of alpha and integrate based on that.
if maxalpha>1 && minalpha>1
    prob(3)=0;
elseif maxalpha>1 %just upper bound is above.
    maxalpha=1;
    lengthDalpha=maxDalpha-minDalpha; lengthalpha=maxalpha-minalpha;
    x1 = linspace(minDalpha,maxDalpha,100); x2 = linspace(minalpha,maxalpha,100);
    [X1,X2] = meshgrid(x1,x2);
    out=intfuncA2D(X1,X2,N,yi,ri,dr,tau);
    prob(3)=1/lengthalpha*1/lengthDalpha*trapz(x2,trapz(x1,out,1));
else %both below 1.
    lengthDalpha=maxDalpha-minDalpha; lengthalpha=maxalpha-minalpha;
    x1 = linspace(minDalpha,maxDalpha,100); x2 = linspace(minalpha,maxalpha,100);
    [X1,X2] = meshgrid(x1,x2);
    out=intfuncA2D(X1,X2,N,yi,ri,dr,tau);
    prob(3)=1/lengthalpha*1/lengthDalpha*trapz(x2,trapz(x1,out,1));
end

% %option 2, exclude anything with alpha>1
% %Find length of interval, find probabilities at different values of
% %parameters, and then do the integration.
% lengthDalpha=maxDalpha-minDalpha; lengthalpha=maxalpha-minalpha;
% x1 = linspace(minDalpha,maxDalpha,100); x2 = linspace(minalpha,maxalpha,100);
% [X1,X2] = meshgrid(x1,x2);
% if beta(5)>1
%     prob(3)=0;
% else
% out=intfuncA2D(X1,X2,N,yi,ri,dr,tau);
% prob(3)=1/lengthalpha*1/lengthDalpha*trapz(x2,trapz(x1,out,1));
% end


%Double Diffusion
minD1=beta(6)-dbeta(6);maxD1=beta(6)+dbeta(6);
minD2=beta(7)-dbeta(7);maxD2=beta(7)+dbeta(7);
minFD=beta(8)-dbeta(8);maxFD=beta(8)+dbeta(8);

%checks against negative values
if minD1 <=0 ,minD1=1e-10; end
if maxD1 <=0 , maxD1=1e-9; end
if minD2 <= 0, minD2=1e-10; end
if maxD2 <= 0, maxD2=1e-9; end
if minFD <= 0, minFD=0.01; end
if maxFD <= 0, maxFD=0.02; end

%Find length of interval, find probabilities at different values of
%parameters, and then do the integration.
lengthD1=maxD1-minD1; lengthD2=maxD2-minD2;
lengthfd=maxFD-minFD;
x1 = linspace(minD1,maxD1,20); x2 = linspace(minD2,maxD2,20);
x3 = linspace(minFD,maxFD,10);
[X1,X2,X3] = meshgrid(x1,x2,x3);
out=intfuncDD2D(x1,x2,x3,N,yi,ri,dr,tau);
prob(4)=1/lengthD1*1/lengthD2*1/lengthfd*trapz(x3,(trapz(x2,trapz(x1,out,1),2)));


%diffusion directed motion
minD=beta(9)-dbeta(9);maxD=beta(9)+dbeta(9);
minV=beta(10)-dbeta(10);maxV=beta(10)+dbeta(10);
minDV=beta(11)-dbeta(11);maxDV=beta(11)+dbeta(11);
minFD=beta(12)-dbeta(12);maxFD=beta(12)+dbeta(12);

%checks against negative values
if minD <=0 ,minD=1e-10; end
if maxD <=0 , maxD=1e-9; end
if minV <=0, minV=0.01; end
if maxV <=0, maxV=0.02; end
if minDV <= 0, minDV=1e-10; end
if maxDV <= 0, maxDV=1e-9; end
if minFD <= 0, minFD=0.01; end
if maxFD <= 0, maxFD=0.02; end

%Find length of interval, find probabilities at different values of
%parameters, and then do the integration.
lengthD=maxD-minD; lengthV=maxV-minV;
lengthDV=maxDV-minDV; lengthfd=maxFD-minFD;
x1 = linspace(minD,maxD,20); x2 = linspace(minV,maxV,20);
x3=linspace(minDV,maxDV,20); x4 = linspace(minFD,maxFD,10);
out=intfuncDV2D(x1,x2,x3,x4,N,yi,ri,dr,tau);
prob(5)=1/lengthD*1/lengthV*1/lengthDV*1/lengthfd*trapz(x4,trapz(x3,trapz(x2,trapz(x1,out,1),2),3));

%diffusion anomalous diffusion
%diffusion directed motion
minD=beta(13)-dbeta(13);maxD=beta(13)+dbeta(13);
minDalpha=beta(10)+dbeta(10);maxDalpha=beta(14)-dbeta(14);
minalpha=beta(15)+dbeta(15);maxalpha=beta(15)-dbeta(15);
minFD=beta(16)+dbeta(16);maxFD=beta(16)-dbeta(16);

%checks against negative values
if minD <=0 ,minD=1e-10; end
if maxD <=0 , maxD=1e-9; end
if minDalpha <=0 ,minDalpha=1e-10; end
if maxDalpha <=0 , maxDalpha=1e-9; end
if minalpha <=0, minalpha=0.01; end
if maxalpha <=0, maxalpha=0.02; end
if minFD <= 0, minFD=0.01; end
if maxFD <= 0, maxFD=0.02; end

%Find length of interval, find probabilities at different values of
%parameters, and then do the integration.
lengthDalpha=maxDalpha-minDalpha; lengthalpha=maxalpha-minalpha;
lengthD=maxD-minD;  lengthfd=maxFD-minFD;
x1 = linspace(minD,maxD,20); x2=linspace(minDalpha,maxDalpha,20);
x3=linspace(minalpha,maxalpha,10); x4 = linspace(minFD,maxFD,10);
out=intfuncDA2D(x1,x2,x3,x4,N,yi,ri,dr,tau);
prob(6)=1/lengthDalpha*1/lengthalpha*1/lengthD*1/lengthfd*trapz(x4,trapz(x3,trapz(x2,trapz(x1,out,1),2),3));

%directed anomalous diffusion
%diffusion directed motion
minV=beta(17)-dbeta(17);maxV=beta(17)+dbeta(17);
minDV=beta(18)-dbeta(18);maxDV=beta(18)+dbeta(18);
minDalpha=beta(19)-dbeta(19);maxDalpha=beta(19)+dbeta(19);
minalpha=beta(20)-dbeta(20);maxalpha=beta(20)+dbeta(20);
minFD=beta(21)-dbeta(21);maxFD=beta(21)+dbeta(21);

%checks against negative values
if minDalpha <=0 ,minDalpha=1e-10; end
if maxDalpha <=0 , maxDalpha=1e-9; end
if minalpha <=0, minalpha=0.01; end
if maxalpha <=0, maxalpha=0.02; end
if minV <=0, minV=0.01; end
if maxV <=0, maxV=0.02; end
if minDV <= 0, minDV=1e-10; end
if maxDV <= 0, maxDV=1e-9; end
if minFD <= 0, minFD=0.01; end
if maxFD <= 0, maxFD=0.02; end

%Find length of interval, find probabilities at different values of
%parameters, and then do the integration.
lengthD=maxD-minD; lengthV=maxV-minV;
lengthDV=maxDV-minDV; lengthfd=maxFD-minFD;
x1 = linspace(minV,maxV,20); x2=linspace(minDV,maxDV,20);
x3=linspace(minDalpha,maxDalpha,20); x4=linspace(minalpha,maxalpha,10);
x5 = linspace(minFD,maxFD,20);
out=intfuncVA2D(x1,x2,x3,x4,x5,N,yi,ri,dr,tau);
prob(7)=1/lengthD*1/lengthV*1/lengthDV*1/lengthfd*trapz(x5,trapz(x4,trapz(x3,trapz(x2,trapz(x1,out,1),2),3),4));

%Normalize Probabilities
%prob=prob/sum(prob);

%Select the best method
[value,method]=max(prob);
end

