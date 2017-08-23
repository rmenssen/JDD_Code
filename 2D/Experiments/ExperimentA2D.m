%Parameters and Preliminaries
rng(15)
tau=.14; %end of simulation
dt=20*10^(-3); %timestep
dt=.2;
tau=8*dt;
points=round(tau/dt)+1; %number of points in each trajectory
M=points-1;
N=3000; %number of trajectories total
Nb=round(N/100); %number of bins %Nb=round(N/100)
Dalpha=.2; %micro meters^2/s^alpha For A and Da
alpha=0.4; %For A and DA
%Simulate Data
[x,y]=AnomalousDiffusion2D(Dalpha,alpha,points,N,dt,tau);
%[x,y]=AnomalousDiffusion2(Dalpha,alpha,points,N,dt);
[jd]=JumpDistance2D(x,y,N); %jd is a vertical vector
[dr, Ni, yi, ri] =  BinningHist(jd, N, Nb,'yes');
%parameter fitting, leave out for now until it works....
%param = ModelFitting(tau, dr, ri, yi, Ni, N, points);
if alpha < 0.5
    min=-500^(.5/alpha); %limits on inverse laplace transform
else
    min=-500;
end 
fun=@(p) (exp(1i.*p.*tau)).*(1i.*p)^(alpha-1)/(2.*pi).*(besselk(0,ri./(sqrt(Dalpha)).*((1i*p)^(alpha/2))));
%fun=@(p) (exp(1i.*p.*tau)).*(1i.*p)^(2*alpha/2-1)/(2.*pi).*(exp(-ri./(sqrt(Dalpha)).*((1i*p)^(alpha/2))));
data=N*dr*ri/(Dalpha).*abs(integral(fun,min,-1*min,'ArrayValued',true,'AbsTol',1e-6));
hold on 
plot(ri,data,'r','LineWidth',2)
