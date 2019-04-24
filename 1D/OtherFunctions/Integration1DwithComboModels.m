%1D integration for Bayesian Classifier
%Rebecca Menssen
%This version of the code: 3/20/19

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

function [prob,value,method]=Integration1DwithComboModels(dbeta,beta,N,yi,ri,dr,tau)
%Vector for storing probabilities
prob=zeros(6,1);

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
out=intfuncD(x1,N,yi,ri,dr,tau);
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
out=intfuncV(X1,X2,N,yi,ri,dr,tau);
prob(2)=1/lengthV*1/lengthDv*trapz(x2,trapz(x1,out,1));

%ANOMALOUS DIFFUSION INTEGRATION

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

%now the case where alpha>1
if maxalpha>1 && minalpha>1
    prob(3)=0;
elseif maxalpha>1 %just upper bound is above.
    maxalpha=1;
    lengthDalpha=maxDalpha-minDalpha;
    lengthalpha=maxalpha-minalpha;
    x1 = linspace(minDalpha,maxDalpha,100);
    x2 = linspace(minalpha,maxalpha,100);
    [X1,X2] = meshgrid(x1,x2);
    out=intfuncA(X1,X2,N,yi,ri,dr,tau);
    prob(3)=1/lengthalpha*1/lengthDalpha*trapz(x2,trapz(x1,out,1));
else %both below 1.
    lengthDalpha=maxDalpha-minDalpha;
    lengthalpha=maxalpha-minalpha;
    x1 = linspace(minDalpha,maxDalpha,100);
    x2 = linspace(minalpha,maxalpha,100);
    [X1,X2] = meshgrid(x1,x2);
    out=intfuncA(X1,X2,N,yi,ri,dr,tau);
    prob(3)=1/lengthalpha*1/lengthDalpha*trapz(x2,trapz(x1,out,1));
end


%%double diffusion
minD1=beta(6)-dbeta(6); maxD1=beta(6)+dbeta(6); 
minD2=beta(7)-dbeta(7); maxD2=beta(7)+dbeta(7);
minfd=beta(8)-dbeta(8); maxfd=beta(8)+dbeta(8);
if minD1 <=0, minD1=1e-10; end
if maxD1 <=0, maxD1=1e-9; end
if minD2 <=0, minD2=1e-10; end
if maxD2 <=0, maxD2=1e-9; end
if minfd <=0, minfd=1e-10; end
if maxfd <=0, maxfd=1e-9; end
if maxfd > 1, maxfd=1; end
lengthD1=maxD1-minD1; lengthD2=maxD2-minD2; lengthfd=maxfd-minfd;
x1 = linspace(minD1,maxD1,50); x2 = linspace(minD2,maxD2,50); x3 = linspace(minfd,maxfd,24);
[out]=intfuncDD(x1,x2,x3,N,yi,ri,dr,tau);
prob(4)=1/lengthD1*1/lengthD2*1/lengthfd*trapz(x3,(trapz(x2,trapz(x1,out,1),2)));


%%diffusion directed combo motion
minDdv=beta(9)-dbeta(9); maxDdv=beta(9)+dbeta(9); 
minVdv=beta(10)-dbeta(10); maxVdv=beta(10)+dbeta(10);
minDvdv=beta(11)-dbeta(11); maxDvdv=beta(11)+dbeta(11);
minfd=beta(12)-dbeta(12); maxfd=beta(12)+dbeta(12);
if minDdv <=0, minDdv=1e-10; end
if maxDdv <=0, maxDdv=1e-9; end
if minVdv <=0, minVdv=1e-10; end
if maxVdv <=0, maxVdv=1e-9; end
if minDvdv <=0, minDvdv=1e-10; end
if maxDvdv <=0, maxDvdv=1e-9; end
if minfd <=0, minfd=1e-10; end
if maxfd <=0, maxfd=1e-9; end
if maxfd > 1, maxfd=1; end
lengthDdv=maxDdv-minDdv; lengthVdv=maxVdv-minVdv; lengthDvdv=maxDvdv-minDvdv; lengthfd=maxfd-minfd;
x1 = linspace(minDdv,maxDdv,12); x2 = linspace(minVdv,maxVdv,12); 
x3 = linspace(minDvdv,maxDvdv,12); x4 = linspace(minfd,maxfd,12);
out=intfuncDV(x1,x2,x3,x4,N,yi,ri,dr,tau);
prob(5)=1/lengthDdv*1/lengthVdv*1/lengthDvdv*1/lengthfd*trapz(x4,trapz(x3,trapz(x2,trapz(x1,out,1),2),3));


%%diffusion anomalous combo motion
minDda=beta(13)-dbeta(13); maxDda=beta(13)+dbeta(13);
minDalphada=beta(14)-dbeta(14); maxDalphada=beta(14)+dbeta(14);
minalphada=beta(15)-dbeta(15); maxalphada=beta(15)+dbeta(15);
minfd=beta(16)-dbeta(16); maxfd=beta(16)+dbeta(16);
if minDda <=0, minDda=1e-14; end
if maxDda <=0, maxDda=1e-13; end
if minDalphada <=0, minDalphada=1e-14; end
if maxDalphada <=0, maxDalphada=1e-13; end
if minalphada <=0, minalphada=1e-14; end
if maxalphada <=0, maxalphada=1e-13; end
if minfd <=0, minfd=1e-14; end
if maxfd <=0, maxfd=1e-13; end
if maxfd > 1, maxfd=1; end
%now cases where alpha is greater than 1
if maxalphada && minalphada > 1
    prob(6)=0;
elseif maxalphada>1
    maxalphada=1;
    lengthDda=maxDda-minDda; lengthDalphada=maxDalphada-minDalphada; 
    lengthalphada=maxalphada-minalphada; lengthfd=maxfd-minfd;
    x1 = linspace(minDda,maxDda,24); x2 = linspace(minDalphada,maxDalphada,24);
    x3 = linspace(minalphada,maxalphada,24); x4 = linspace(minfd,maxfd,24);
    [out]=intfuncDA(x1,x2,x3,x4,N,yi,ri,dr,tau);
    prob(6)=1/lengthDda*1/lengthDalphada*1/lengthalphada*1/lengthfd*trapz(x4,trapz(x3,trapz(x2,trapz(x1,out,1),2),3));
else
    lengthDda=maxDda-minDda; lengthDalphada=maxDalphada-minDalphada; 
    lengthalphada=maxalphada-minalphada; lengthfd=maxfd-minfd;
    x1 = linspace(minDda,maxDda,12); x2 = linspace(minDalphada,maxDalphada,12);
    x3 = linspace(minalphada,maxalphada,12); x4 = linspace(minfd,maxfd,12);
    [out]=intfuncDA(x1,x2,x3,x4,N,yi,ri,dr,tau);
    prob(6)=1/lengthDda*1/lengthDalphada*1/lengthalphada*1/lengthfd*trapz(x4,trapz(x3,trapz(x2,trapz(x1,out,1),2),3));
end


%Normalize Probabilities (if wanted, not necessary really). 
prob=prob/sum(prob);

%Select the best method
[value,method]=max(prob);
end

