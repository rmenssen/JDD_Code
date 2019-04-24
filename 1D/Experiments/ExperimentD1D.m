%Pure Diffusion Demonstration of JDD Method
%Rebecca Menssen
%Last Updated: 4/9/19

%This code serves as a way to experiment with the JDD method. It provides a
%demonstration of how the method works from start to finish. Parameters can
%be edited to examine accuracy of the method. The code has five sections:
%Parameters, Simulation, Initial Fitting, Bootstrapping, and Model
%Selection. This code is for non-sliding data, but can be edited for
%a sliding JDD

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

%Number of Bootstraps
numboot=50;

%%
%%%%%%%%%%DIFFUSION SIMULATION AND CREATION OF JDD%%%%%%%%%%
%set a seed
seed=randi(1000);

%Simulate Diffusion:
[x]=Diffusion1D(D,points,N,dt,seed);

%Create the Jump Distance
[jd]=JumpDistance1D(x,N);

%Number of Bins for fitting
%Choose option here. 
Nb=round(1+log2(N)); %Sturges Rule
sigma=sqrt(6*(N-2)/(N+1)/(N+3)); %for Doane's rule
%Nb=round(1+log2(N)+log2(1+abs(skewness(jd))/sigma)); %Doanes Rule
%Nb=round(2*(N^(1/3))); %Rice Rule
%Nb=round(sqrt(N)); %square root guidance
%Nb=round((max(jd)-min(jd))*N^(1/3)/(3.5*std(jd))); %Scott's Normal Reference Rule. 
%Nb=round((max(jd)-min(jd))*N^(1/3)/(2*iqr(jd))); %freedman diaconis Rule 

%Plot the Jump Distance
figure(1)
[dr, Ni, yi, ri] =  BinningHist(jd, N, Nb,'yes');

%Plot the predicted JDD on top of it
predictedJDD=N*dr/((pi*D*tau)^(1/2)).*exp(-ri.^2/(4*D*tau));
hold on
plot(ri,predictedJDD,'k','LineWidth',1.5)

xlabel('Jump Distance')
ylabel('Count')
title('Pure Diffusion Jump Distance Distribution in 1D')

%%
%%%%%%%%%%MODEL FITTING%%%%%%%%%%
param = ModelFitting1D(tau, dr, ri, yi, Ni,N, points, dt, x);

%plotting best fit for each model
diffusionbest=N*dr/((pi*param.D*tau)^(1/2)).*exp(-ri.^2/(4*param.D*tau));
hold on
plot(ri,diffusionbest,'b','LineWidth',1.5)

z2 = -(ri.^2+param.V^2*tau^2)/(4*param.Dv*tau);
y2 = ri*param.V/(2*param.Dv);
directedbest = N*dr/((4*pi*param.Dv*tau)^(1/2)).*exp(z2+y2)+N*dr/((4*pi*param.Dv*tau)^(1/2)).*exp(z2-y2);
plot(ri,directedbest,'r','LineWidth',1.5)

alpha=param.alpha;
Dalpha=param.Dalpha;
if alpha < 0.5
    min=-300^(.5/alpha);
else
    min=-500;
end
fun=@(p) (exp(1i.*p.*tau)).*(1i.*p)^(alpha/2-1)/(2.*pi).*exp(-ri./(sqrt(Dalpha)).*((1i*p)^(alpha/2)));
anombest=2*N*dr/((Dalpha)^(1/2)*2).*abs(integral(fun,min,-1*min,'ArrayValued',true,'AbsTol',1e-6,'RelTol',1e-3));
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
beta=[param.D,param.V,param.Dv,param.Dalpha,param.alpha];
seeds=beta;
parfor i=1:numboot
    X = randi(N,N,1);
    jdB=jd(X);
    %if using Doane's rule can choose to update here.
    %Nb=round(1+log2(N)+log2((1+skewness(jdB))/(sqrt((6*(N-2))/((N+1)*(N+3))))));
    [drB, NiB, yiB, riB] =  BinningHist(jdB, N, Nb,'no');
    paramB = ModelFitting1D(tau, drB, riB, yiB, NiB,N, points, dt, x);
    %if you want to use seeds instead
    %paramB = ModelFitting1DwithSeeds(tau, drB, riB, yiB, NiB,N, points, dt, seeds);
    Dboot(i)=paramB.D; Vboot(i)=paramB.V; Dvboot(i)=paramB.Dv;
    Daboot(i)=paramB.Dalpha; Aboot(i)=paramB.alpha;
end

%beta=[param.D,param.V,param.Dv,param.Dalpha,param.alpha];
dbeta=2*[std(Dboot),std(Vboot),std(Dvboot),std(Daboot), std(Aboot)];

%%
%%%%%%%%%%MODEL SELECTION%%%%%%%%%%
[prob,value,method]=Integration1D(dbeta,beta,N,yi,ri,dr,tau);
