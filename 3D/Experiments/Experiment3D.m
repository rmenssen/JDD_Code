%Parameters and Preliminaries

dt=23.5800; %timestep
points=6; %number of points in each trajectory
tau=(points-1)*dt;
N=3600; %number of trajectories total
Nb=36; %number of bins %Nb=round(N/100)
D=1e-3; %micro meters^2/s (used in D,DV, and DA)


%simulate Diffusion
[x,y,z]=Diffusion3D(D,points,N,dt);
%do the jdd and binning
[jd]=JumpDistance3D(x,y,z,N); %jd is a vertical vector
[dr, Ni, yi, ri] =  BinningHist(jd, N, Nb,'yes');
%param = ModelFitting(tau, dr, ri, yi, Ni, N, points);

%plot the predicted on top of it
data=N*dr*ri.^2/(2*sqrt(pi)*(D*tau)^(3/2)).*exp(-ri.^2/(4*D*tau));
param = ModelFitting3D(tau, dr, ri, yi, Ni, N, points,dt,x,y,z)

%plot the actual fit
data2=N*dr*ri.^2/(2*sqrt(pi)*(param.D*tau)^(3/2)).*exp(-ri.^2/(4*param.D*tau));

hold on 
plot(ri,data,'b','LineWidth',1.5)
plot(ri,data2,'r','LineWidth',1.5)
