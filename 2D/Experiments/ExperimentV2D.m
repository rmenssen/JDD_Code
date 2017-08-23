%Experiment to see if data is generated properly for
%directed motion

%Parameters and Preliminaries
tau=15; %end of simulation
dt=1; %timestep
points=tau/dt+1; %number of points in each trajectory
M=points-1;
N=3000; %number of trajectories total
Nb=30; %number of bins %Nb=round(N/100)
D=0.02; %micro meters^2/s (used in D,DV, and DA)
V=1; %micro m/s %For V and DV
Dv=.5;
Dalpha=0.02; %micro meters^2/s^alpha For A and Da
alpha=0.5; %For A and DA
D1=0.02; %micro meters^2/s %for DD
D2=0.1; %micro meters^2/s for DD
fd=0.5; %fraction of particles that are diffusing.

[x,y]=DirectedMotion2D(V,2*Dv*dt,points,N,dt);
[jd]=JumpDistance2D(x,y,N); %jd is a vertical vector
[dr, Ni, yi, ri] =  BinningHist(jd, N, Nb,'yes');
param = ModelFitting(tau, dr, ri, yi, Ni, N, points);
z = -(ri.^2+V^2*tau^2)/(4*Dv*tau);
y = ri*V/(2*Dv);
data = N*dr*ri/(2*Dv*tau).*exp(z).*besseli(0, y);
%data2 = N*dr*ri/(kv*M).*exp(z).*exp(y)./sqrt(2*pi*y);
%data2 = N*dr*ri/(kv*M).*exp(z).*exp(y).*sqrt(y);
%z1 = -(ri.^2+param.V^2*tau^2)/(2*M*param.kv);
%y1 = ri*param.V*tau/(M*param.kv);
%data2 = N*dr*ri/(param.kv*M).*exp(z1).*besseli(0, y1);
%data3=N*dr*ri/(kv*M).*exp(-(ri-V*tau).^2/(2*M*kv))/(2*pi);
%data3=N*dr*ri/(param.kv*M).*exp(z).*exp(y)/(2*pi);
%data3=N*dr*ri/(param.kv*M).*exp(-(ri.^2+param.V^2*tau^2)/(2*M*param.kv));
hold on 
plot(ri,data,'b')
%plot(ri,data2,'r')
%plot(ri,data3,'g')
%param.V
%param.kv/2/dt