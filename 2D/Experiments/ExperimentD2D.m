%Experiment to see if data is generated properly
%for diffusion

%Parameters and Preliminaries
tau=.14; %end of simulation
dt=20*10^(-3); %timestep
points=tau/dt+1; %number of points in each trajectory
M=points-1;
N=3000; %number of trajectories total
Nb=30; %number of bins %Nb=round(N/100)
%Nb=N/50; <-experiment with this number
D=0.02; %micro meters^2/s (used in D,DV, and DA)
V=1.2; %micro m/s %For V and DV
kv=0.0008;%micro m^2/s %For V and DV. Related to diffusion constant
Dalpha=0.02; %micro meters^2/s^alpha For A and Da
alpha=0.5; %For A and DA
D1=0.02; %micro meters^2/s %for DD
D2=0.1; %micro meters^2/s for DD
fd=0.5; %fraction of particles that are diffusing.

%simulate Diffusion
[x,y]=Diffusion2D(D,points,N,dt);
%plot(x,y)
%do the jdd and binning
[jd]=JumpDistance2D(x,y,N); %jd is a vertical vector
[dr, Ni, yi, ri] =  BinningHist(jd, N, Nb,'yes');
param = ModelFitting(tau, dr, ri, yi, Ni, N, points);

%plot the predicted on top of it
data=N*dr*ri/(2*D*tau).*exp(-ri.^2/(4*D*tau));
data2=N*dr*ri/(2*param.D*tau).*exp(-ri.^2/(4*param.D*tau));
%plot
hold on 
plot(ri,data,'b')
plot(ri,data2,'r')
averages=0.021;
beta=param.D;
%integration step
%[prob,method]=IntegrationD(averages,N,yi,ri,dr,tau,beta,M)

z1 = -(ri.^2+param.V^2*tau^2)/(2*M*param.kv);
y1 = ri*param.V*tau/(M*param.kv);
data2 = N*dr*ri/(param.kv*M).*exp(z1).*besseli(0, y1);
plot(ri,data2,'g')