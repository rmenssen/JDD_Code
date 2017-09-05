%Experiment on how things work with 1D and bootstrapping. General ranges of
%parameters and the like

%This experiment is for Pure 1D Diffusion. 
rng('shuffle')
%Parameters and Preliminaries
%tau=.14; %end of simulation
tau=30*10;
%tau=45*7;
%dt=20*10^(-3); %timestep
dt=30;
%dt=45;
points=tau/dt+1; %number of points in each trajectory
N=3000; %number of trajectories total
Nb=30; %number of bins %Nb=round(N/100)
D=0.001; %micro meters^2/s for x direction

%simulate Diffusion
[x]=Diffusion1D(D,points,N,dt);
%do the jdd and binning for x (or all unsure which to do yet)
[jd]=JumpDistance1D(x,N); %jd is a vertical vector
%look at xdirection

[dr, Ni, yi, ri] =  BinningHist(jd, N, Nb,'no');

%parameter fit. 
param = ModelFitting1D(tau, dr, ri, yi, Ni);

%now do boostrapping. 
%store parameters:
beta=[param.D,param.V,param.Dv,param.Dalpha, param.alpha];
%now for some bootstrapping and more parameter things.
rng('shuffle')
numbootstrap=20;
%set up storage
Dboot=zeros(numbootstrap,1);
Vboot=zeros(numbootstrap,1);
Dvboot=zeros(numbootstrap,1);
Daboot=zeros(numbootstrap,1);
Aboot=zeros(numbootstrap,1);

for i=1:numbootstrap
    X = randi(N,N,1);
    jdB=jd(X);
    [drB, NiB, yiB, riB] =  BinningHist(jdB, N, Nb,'no');
    paramB = ModelFitting1D(tau, drB, riB, yiB, NiB);
    Dboot(i)=paramB.D;
    Vboot(i)=paramB.V;
    Dvboot(i)=paramB.Dv;
    Daboot(i)=paramB.Dalpha;
    Aboot(i)=paramB.alpha;
end
avg=[mean(Dboot),mean(Vboot),mean(Dvboot),mean(Daboot),mean(Aboot)];
dbeta=10*(abs(beta-avg));
avg2=[median(Dboot),median(Vboot),median(Dvboot),median(Daboot),median(Aboot)];
stda=[std(Dboot), std(Vboot),std(Dvboot),std(Daboot), std(Aboot)];
dbeta2=10*(abs(beta-avg2));

[prob,value,method]=Integration1D(dbeta,beta,N,yi,ri,dr,tau)
