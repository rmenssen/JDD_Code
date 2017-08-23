%Experiment to see if data is generated properly for
%directed motion

%Parameters and Preliminaries
%tau=.14; %end of simulation
dt=.1;
tau=21*dt;
%dt=20*10^(-3); %timestep

points=round(tau/dt)+1; %number of points in each trajectory
M=points-1;
N=3000; %number of trajectories total
Nb=30; %number of bins %Nb=round(N/100)
%D=0.02; %micro meters^2/s (used in D,DV, and DA)
V=1; %micro m/s %For V and DV
D=0.5;

[x]=DirectedMotion1D(V,2*D*dt,points,N,dt,10);
[jd]=JumpDistance1D(x,N); %jd is a vertical vector
[dr, Ni, yi, ri] =  BinningHist(jd, N, Nb,'yes');

%drz*2/((4*pi*D3*tau)^(1/2)).
z = -(ri.^2+V^2*tau^2)/(4*D*tau);
y = ri*V/(2*D);
data = N*dr/((4*pi*D*tau)^(1/2)).*exp(z+y);

hold on
plot(ri,data,'b')

% param = ModelFitting1DV(tau, dr, ri, yi, Ni, N, points)
% 
% z2 = -(ri.^2+param.V^2*tau^2)/(4*param.Dv*tau);
% y2 = ri*param.V/(2*param.Dv);
% data2 = N*dr/((4*pi*param.Dv*tau)^(1/2)).*exp(z2+y2);
% plot(ri,data2,'g')

param = ModelFitting1D(tau, dr, ri, yi, Ni,N, points, dt, x)
%param = ModelFitting1DV(tau, dr, ri, yi, Ni, N, points,dt,x)

% data2=N*dr/((pi*param.D*tau)^(1/2)).*exp(-ri.^2/(4*param.D*tau));
% hold on
% plot(ri,data2,'r')
% 
z2 = -(ri.^2+param.V^2*tau^2)/(4*param.Dv*tau);
y2 = ri*param.V/(2*param.Dv);
data2 = N*dr/((4*pi*param.Dv*tau)^(1/2)).*exp(z2+y2);
plot(ri,data2,'y')
% 
 alpha=param.alpha;
 Dalpha=param.Dalpha;
% 
if alpha < 0.5
    min=-300^(.5/alpha); %limits on inverse laplace transform
else
    min=-300;
end
fun=@(p) (exp(1i.*p.*tau)).*(1i.*p)^(alpha/2-1)/(2.*pi).*exp(-ri./(sqrt(Dalpha)).*((1i*p)^(alpha/2)));
data=2*N*dr/((Dalpha)^(1/2)*2).*abs(integral(fun,min,-1*min,'ArrayValued',true,'AbsTol',2e-5));
plot(ri,data,'b')
