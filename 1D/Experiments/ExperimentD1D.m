%Experiment for anisotropic diffusion
%Rebecca Menssen

%Parameters and Preliminaries
%tau=.14; %end of simulation
tau=16*1;
%tau=45*7;
%dt=20*10^(-3); %timestep
dt=1;
%dt=45;
points=tau/dt+1; %number of points in each trajectory
N=3000; %number of trajectories total
Nb=30; %number of bins %Nb=round(N/100)
D=0.5; %micro meters^2/s for x direction


%simulate Diffusion
[x]=Diffusion1D(D,points,N,dt,8);
%do the jdd and binning for x (or all unsure which to do yet)
[jd]=JumpDistance1D(x,N); %jd is a vertical vector
%look at xdirection
figure(1)
[dr, Ni, yi, ri] =  BinningHist(jd, N, Nb,'yes');
%plot predicted equation on top of it:
data=N*dr/((pi*D*tau)^(1/2)).*exp(-ri.^2/(4*D*tau));
hold on
plot(ri,data,'b')

param = ModelFitting1D(tau, dr, ri, yi, Ni,N, points, dt, x)

data2=N*dr/((pi*param.D*tau)^(1/2)).*exp(-ri.^2/(4*param.D*tau));
hold on
plot(ri,data2,'r')

z2 = -(ri.^2+param.V^2*tau^2)/(4*param.Dv*tau);
y2 = ri*param.V/(2*param.Dv);
data2 = N*dr/((4*pi*param.Dv*tau)^(1/2)).*exp(z2+y2);
plot(ri,data2,'y')

alpha=param.alpha;
Dalpha=param.Dalpha;

if alpha < 0.5
    min=-500^(.5/alpha); %limits on inverse laplace transform
else
    min=-500;
end
fun=@(p) (exp(1i.*p.*tau)).*(1i.*p)^(alpha/2-1)/(2.*pi).*exp(-ri./(sqrt(Dalpha)).*((1i*p)^(alpha/2)));
data=2*N*dr/((Dalpha)^(1/2)*2).*abs(integral(fun,min,-1*min,'ArrayValued',true,'AbsTol',1e-6,'RelTol',1e-3));
plot(ri,data,'b')
