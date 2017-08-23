%Parameters and Preliminaries
%rng('shuffle')
rng(5)
%tau=.14; %end of simulation
tau=15*1;
%dt=20*10^(-3); %timestep
dt=1;
points=round(tau/dt+1); %number of points in each trajectory
M=points-1;
N=3000; %number of trajectories total
Nb=round(N/100); %number of bins %Nb=round(N/100)
Dalpha=1; %micro meters^2/s^alpha For A and Da
alpha=0.6; %For A and DA
%Simulate Data
[x]=AnomalousDiffusion1D(Dalpha,alpha,points,N,dt,tau,5);
[jd]=JumpDistance1D(x,N); %jd is a vertical vector
[dr, Ni, yi, ri] =  BinningHist(jd, N, Nb,'yes');
%parameter fitting, leave out for now until it works....
%param = ModelFitting(tau, dr, ri, yi, Ni, N, points);
if alpha < 0.5
    %minb=-500^(.5/alpha); %limits on inverse laplace transform
    minb=-1100;
else
    minb=-500;
end

fun=@(p) (exp(1i.*p.*tau)).*(1i.*p)^(alpha/2-1)/(2.*pi).*exp(-ri./(sqrt(Dalpha)).*((1i*p)^(alpha/2)));
%data=N*dr/((Dalpha)^(1/2)).*abs(integral(fun,minb,-1*minb,'ArrayValued',true,'AbsTol',1e-6,'RelTol',1e-3));
%data=N*dr/((Dalpha)^(1/2)).*abs(integral(fun,minb,10,'ArrayValued',...
%    true,'AbsTol',1e-6,'RelTol',1e-3))+N*dr/((Dalpha)^(1/2)).*...
%    abs(integral(fun,10,-1*minb,'ArrayValued',true,'AbsTol',1e-6,'RelTol',1e-3));

data=N*dr/((Dalpha)^(1/2)).*abs(integral(fun,minb,-100,'ArrayValued',...
    true,'AbsTol',1e-6,'RelTol',1e-3))+N*dr/((Dalpha)^(1/2)).*...
    abs(integral(fun,-100,100,'ArrayValued',true,'AbsTol',1e-6,'RelTol',1e-3))+N*dr/((Dalpha)^(1/2)).*...
    abs(integral(fun,100,-1*minb,'ArrayValued',true,'AbsTol',1e-6,'RelTol',1e-3));


hold on
plot(ri,data,'b','LineWidth',1.5)

%try splitting it up into 200 point intevals
%first find the total number of intevals needed
intsize=200;
numint=ceil(abs(minb+100)/intsize);
testtest=zeros(1,Nb);

for i=1:numint
    testtest=testtest+2*N*dr/((Dalpha)^(1/2)).*abs(integral(fun,...
        minb+(i-1)*intsize,minb+i*intsize,...
        'ArrayValued',true,'AbsTol',1e-6,'RelTol',1e-3));
end
testtest=testtest+N*dr/((Dalpha)^(1/2)).*abs(integral(fun,-100,...
    100,'ArrayValued',true,'AbsTol',1e-6,'RelTol',1e-3));

plot(ri,testtest,'r')

%param1 = ModelFitting1DA(tau, dr, ri, yi, Ni, N, points);
% param2 = ModelFitting1DA2(tau, dr, ri, yi, Ni, N, points);
% 
% alpha=param2.alpha;
% Dalpha=param2.Dalpha;
% 
% fun=@(p) (exp(1i.*p.*tau)).*(1i.*p)^(alpha/2-1)/(2.*pi).*exp(-ri./(sqrt(Dalpha)).*((1i*p)^(alpha/2)));
% data=2*N*dr/((Dalpha)^(1/2)*2).*abs(integral(fun,min,-1*min,'ArrayValued',true,'AbsTol',2e-5));
% plot(ri,data,'r')

%param = ModelFitting1DA(tau, dr, ri, yi, Ni, N, points)
%param = ModelFitting1DA(tau, dr, ri, yi, Ni, N, points, dt, x)
% 
% data2=N*dr/((pi*param.D*tau)^(1/2)).*exp(-ri.^2/(4*param.D*tau));
% hold on
% plot(ri,data2,'r')
% 
% z2 = -(ri.^2+param.V^2*tau^2)/(4*param.Dv*tau);
% y2 = ri*param.V/(2*param.Dv);
% data2 = N*dr/((4*pi*param.Dv*tau)^(1/2)).*exp(z2+y2);
% plot(ri,data2,'y')

% alpha=param.alpha;
% Dalpha=param.Dalpha;
% 
% if alpha < 0.5
%     %minb=-500^(.5/alpha); %limits on inverse laplace transform
%     minb=-700;
% else
%     minb=-500;
% end
% 
% fun=@(p) (exp(1i.*p.*tau)).*(1i.*p)^(alpha/2-1)/(2.*pi).*exp(-ri./(sqrt(Dalpha)).*((1i*p)^(alpha/2)));
% data=2*N*dr/((Dalpha)^(1/2)*2).*abs(integral(fun,minb,-1*minb,'ArrayValued',true,'AbsTol',1e-6));
% plot(ri,data,'r','LineWidth',1.5)


