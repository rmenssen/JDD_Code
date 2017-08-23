%Rebecca Menssen
%Modified: 8/23/17

%This function runs through how to analyze code from beginning to end. 
%you first choose your parameters that guide the simulation. 
%Then it will simulate data, parameter fit, bootstrap and select the
%correct model

%Basic Simulation Parameters
rng(15)
%choose timestep
dt=.2;
%choose number of points
points=9;
tau=(points-1)*dt;
%number of trajectories
N=3000; %number of trajectories total
Nb=round(N/100); %number of bins 

%Anomalous Parameters
Dalpha=.2; %micro meters^2/s^alpha For A and Da
alpha=0.4; %For A and DA

%Simulate Data
[x,y]=AnomalousDiffusion2D(Dalpha,alpha,points,N,dt,tau);
%find the jump distance and bin it into the JDD
[jd]=JumpDistance2D(x,y,N); %jd is a vertical vector
[dr, Ni, yi, ri] =  BinningHist(jd, N, Nb,'no');
%change this last parameter to yes to see the histogram

%Parmeter fit
%parameter fitting, leave out for now until it works....
%param = ModelFitting2D(tau, dr, ri, yi, Ni, N, points);
if alpha < 0.5
    min=-500^(.5/alpha); %limits on inverse laplace transform
else
    min=-500;
end 
fun=@(p) (exp(1i.*p.*tau)).*(1i.*p)^(alpha-1)/(2.*pi).*(besselk(0,ri./(sqrt(Dalpha)).*((1i*p)^(alpha/2))));

data=N*dr*ri/(Dalpha).*abs(integral(fun,min,-1*min,'ArrayValued',true,'AbsTol',1e-6));
hold on 
plot(ri,data,'r','LineWidth',2)
