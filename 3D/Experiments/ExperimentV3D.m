%Directed Diffusion Demonstration of JDD Method: 3D
%Rebecca Menssen
%Last Validated 3/14/19

%This code serves as a way to experiment with the JDD method. It provides a
%demonstration of how the method works from start to finish. Parameters can
%be edited to examine accuracy of the method. The code has five sections:
%Parameters, Simulation, Initial Fitting, Bootstrapping, and Model
%Selection. On a two core computer, Intel i5 with 8 GB of ram takes ~ 2 hrs
%Often it is much faster

%%
%%%%%%%%%%SIMULATION PARAMETERS%%%%%%%%%%

%Diffusion Constant
V=5; %micro m/s
Dv=1; %micro meters^2/s 

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

%Simulate Directed Diffusion: 
[x,y,z]=DirectedMotion3D(V,Dv,points,N,dt,seed);

%Create the Jump Distance
[jd]=JumpDistance3D(x,y,z,N); 

%Number of Bins for fitting
%Choose option here. 
Nb=round(1+log2(N)); %Sturges Rule
sigma=sqrt(6*(N-2)/(N+1)/(N+3)); %for Doane's rule
%Nb=round(1+log2(N)+log2(1+abs(skewness(jd))/sigma)); %Doanes Rule
%Nb=round(2*(N^(1/3))); %Rice Rule
%Nb=round(sqrt(N)); %square root guidance
%Nb=round((max(jd)-min(jd))*N^(1/3)/(3.5*std(jd))); %Scott's Normal Reference Rule. 
%Nb=round((max(jd)-min(jd))*N^(1/3)/(2*iqr(jd))); %Freedman Diaconis Rule 

%Plot the Jump Distance
figure(1)
[dr, Ni, yi, ri]=BinningHist(jd, N, Nb,'yes');

%Plot the predicted JDD on top of it
z1 = -(ri.^2+V^2*tau^2)/(4*Dv*tau);
y1 = ri*V/(2*Dv);
predictedJDD= N*dr*ri.^2*4*pi/((4*pi*Dv*tau)^(3/2)).*exp(z1).*sinh(y1)./y1;
hold on
plot(ri,predictedJDD,'k','LineWidth',1.5)

xlabel('Jump Distance')
ylabel('Count')
title('Directed Diffusion Jump Distance Distribution in 3D')

%%
%%%%%%%%%%MODEL FITTING%%%%%%%%%%

param = ModelFitting3D(tau, dr, ri, yi, Ni, N, points,dt,x,y,z);

%plotting best fit for each model
diffusionbest=N*dr*ri.^2/(2*sqrt(pi)*(param.D*tau)^(3/2)).*...
    exp(-ri.^2/(4*param.D*tau));
hold on
plot(ri,diffusionbest,'b','LineWidth',1.5)

z1 = -(ri.^2+param.V^2*tau^2)/(4*param.Dv*tau);
y1 = ri*param.V/(2*param.Dv);
directedbest= N*dr*ri.^2*4*pi/((4*pi*param.Dv*tau)^(3/2)).*exp(z1).*sinh(y1)./y1;
plot(ri,directedbest,'r','LineWidth',1.5)

if param.alpha < 0.5
    min=-300^(.5/param.alpha); 
else
    min=-500;
end
min=-700;
intsize=200;
numint=ceil(abs(min+100)/intsize);
anombest=zeros(1,Nb);
fun2=@(p) (exp(1i.*p.*tau)).*(1i.*p)^(param.alpha-1)/(2.*pi).*...
    (exp(-ri./(sqrt(param.Dalpha)).*((1i*p)^(param.alpha/2))));
for i=1:numint
    anombest=anombest+2*N*dr*ri/(param.Dalpha).*abs(integral(fun2,min+(i-1)*intsize-1i*1e-6,min+i*(intsize)-1i*1e-6,...
        'ArrayValued',true,'AbsTol',1e-6,'RelTol',1e-3));
end
anombest=anombest+N*dr*ri/(param.Dalpha).*abs(integral(fun2,-100-1i*1e-6,100-1i*1e-6,...
    'ArrayValued',true,'AbsTol',1e-6,'RelTol',1e-3));
plot(ri,anombest,'g','LineWidth',1.5)

legend('Jump Distance Distribution',['Predicted Directed Fit, V=',num2str(V),', D_V=',num2str(Dv)],...
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
    %If using Doanes Rule, can consider assessing if a different Nb is
    %needed based on skewness, or can stay with original choice 
    %Nb=round(1+log2(N)+log2((1+skewness(jd))/(sqrt((6*(N-2))/((N+1)*(N+3)))))); %Doanes Rule
    [drB, NiB, yiB, riB] =  BinningHist(jdB, N, Nb,'no');
    paramB = ModelFitting3D(tau, drB, riB, yiB, NiB,N, points, dt, x,y,z);
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
[prob,value,method]=Integration3D(dbeta,beta,N,yi,ri,dr,tau);

