%Pure Diffusion Demonstration of JDD Method
%Rebecca Menssen
%Last Updated 9/1/17

%This code serves as a way to experiment with the JDD method. It provides a
%demonstration of how the method works from start to finish. Parameters can
%be edited to examine accuracy of the method. The code has five sections:
%Parameters, Simulation, Initial Fitting, Bootstrapping, and Model
%Selection. On a two core computer, Intel i5 with 8 GB of ram takes ~ 2 hrs
%Often it is much faster

%%
%%%%%%%%%%SIMULATION PARAMETERS%%%%%%%%%%

%Diffusion Constant
D=1; %micro meters^2/s 

%Time Step
dt=1;

%Time Lag, points and tau 
timelag=15;
points=timelag+1;
tau=dt*timelag;

%Number of trajectories
N=3000; 

%Number of Bins for fitting
Nb=round(N/100);

%Number of Bootstraps
numboot=50;
%%
%%%%%%%%%%DIFFUSION SIMULATION AND CREATION OF JDD%%%%%%%%%%
%set a seed
seed=randi(1000);

%Simulate Diffusion: 
[x,y]=Diffusion2D(D,points,N,dt,seed);

%Create the Jump Distance
[jd]=JumpDistance2D(x,y,N); 

%Plot the Jump Distance
figure(1)
[dr, Ni, yi, ri]=BinningHist(jd, N, Nb,'yes');

%Plot the predicted JDD on top of it
predictedJDD=N*dr*ri/(2*D*tau).*exp(-ri.^2/(4*D*tau));
hold on
plot(ri,predictedJDD,'k','LineWidth',1.5)

xlabel('Jump Distance')
ylabel('Count')
title('Pure Diffusion Jump Distance Distribution in 2D')

%%
%%%%%%%%%%MODEL FITTING%%%%%%%%%%

param = ModelFitting2D(tau, dr, ri, yi, Ni, N, points,dt,x,y);

%plotting best fit for each model
diffusionbest=N*dr*ri/(2*param.D*tau).*exp(-ri.^2/(4*param.D*tau));
hold on
plot(ri,diffusionbest,'b','LineWidth',1.5)

z1 = -(ri.^2+param.V^2*tau^2)/(4*param.Dv*tau);
y1 = ri*param.V/(2*param.Dv);
directedbest = N*dr*ri/(2*param.Dv*tau).*exp(z1).*besseli(0, y1);
plot(ri,directedbest,'r','LineWidth',1.5)

if param.alpha < 0.5
    min=-300^(.5/param.alpha); 
else
    min=-500;
end
fun=@(p) (exp(1i.*p.*tau)).*(1i.*p)^(param.alpha-1)/(2.*pi).*...
    (besselk(0,ri./(sqrt(param.Dalpha)).*((1i*p)^(param.alpha/2))));
anombest=N*dr*ri/(param.Dalpha).*abs(integral(fun,min-1i*1e-6,-1*min-1i*1e-6,...
    'ArrayValued',true,'AbsTol',1e-6,'RelTol',1e-3));
plot(ri,anombest,'g','LineWidth',1.5)

legend('Jump Distance Distribution',['Predicted Diffusion Fit,D=',num2str(D)],...
    ['Fit Diffusion, D=',num2str(param.D)],...
    ['Fit Directed, V=',num2str(param.V),', D_V=',num2str(param.Dv)],...
    ['Fit Anomalous, \alpha=',num2str(param.alpha),', D_\alpha=',num2str(param.Dalpha)])



%%
%%%%%%%%%%BOOTSTRAPPING%%%%%%%%%%

%Set Up Storage
Dboot=zeros(numboot,1);
Vboot=zeros(numboot,1);
Dvboot=zeros(numboot,1);
Daboot=zeros(numboot,1);
Aboot=zeros(numboot,1);

parfor i=1:numboot
    X = randi(N,N,1);
    jdB=jd(X);
    [drB, NiB, yiB, riB] =  BinningHist(jdB, N, Nb,'no');
    paramB = ModelFitting2D(tau, drB, riB, yiB, NiB,N, points, dt, x,y);
    Dboot(i)=paramB.D;
    Vboot(i)=paramB.V;
    Dvboot(i)=paramB.Dv;
    Daboot(i)=paramB.Dalpha;
    Aboot(i)=paramB.alpha;
end

beta=[param.D,param.V,param.Dv,param.Dalpha,param.alpha];
dbeta=2*[std(Dboot),std(Vboot),std(Dvboot),std(Daboot), std(Aboot)];

%%
%%%%%%%%%%MODEL SELECTION%%%%%%%%%%
[prob,value,method]=Integration2D(dbeta,beta,N,yi,ri,dr,tau);