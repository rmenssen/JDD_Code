%Experiment to see if data is generated properly for
%directed motion

%Parameters and Preliminaries
tau=.14; %end of simulation
dt=20*10^(-3); %timestep
points=tau/dt+1; %number of points in each trajectory
M=points-1;
N=3000; %number of trajectories total
Nb=30; %number of bins %Nb=round(N/100)
D=0.02; %micro meters^2/s (used in D,DV, and DA)
V=1.2; %micro m/s %For V and DV


[x,y,z]=DirectedMotion3D(V,2*D*dt,points,N,dt);
[jd]=JumpDistance3D(x,y,z,N); %jd is a vertical vector
[dr, Ni, yi, ri] =  BinningHist(jd, N, Nb,'yes');
param = ModelFitting3D(tau, dr, ri, yi, Ni, N, points,dt,x,y,z)

z = -(ri.^2+V^2*tau^2)/(4*D*tau);
y = ri*V/(2*D);
data= N*dr*ri.^2*4*pi/((4*pi*D*tau)^(3/2)).*exp(z).*sinh(y)./y;
hold on 
plot(ri,data,'b')


%param = ModelFitting3D(tau, dr, ri, yi, Ni, N, points)