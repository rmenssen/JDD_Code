%Double Diffusion Demonstration of JDD Method
%Rebecca Menssen
%Last Validated 4/9/19

%This code serves as a way to experiment with the JDD method. It provides a
%demonstration of how the method works from start to finish. Parameters can
%be edited to examine accuracy of the method. The code has five sections:
%Parameters, Simulation, Initial Fitting, Bootstrapping, and Model
%Selection.

%%
%%%%%%%%%%SIMULATION PARAMETERS%%%%%%%%%%

%Diffusion Constant
%always make D1 < D2 for consistency (reason being sometimes model fitting
%can flip them, so it's best to be consistent)
D1=1; %micro meters^2/s
D2=5; %micro meters^2/s
fd=0.7; %fraction diffusing at D1

%Time Step
dt=1;

%Time Lag, points and tau
timelag=15;
points=timelag+1;
tau=dt*timelag;

%Number of trajectories
N=3000;
N1=fd*N; %number at D1
N2=N-N1; %number at D2

%Number of Bootstraps
numboot=50;


%%
%%%%%%%%%%DIFFUSION SIMULATION AND CREATION OF JDD%%%%%%%%%%
%set a seed (for reproducibility)
seed1=randi(1000);
seed2=randi(1000);

%Simulate Diffusion:
[x1]=Diffusion1D(D1,points,N1,dt,seed1);
[x2]=Diffusion1D(D2,points,N2,dt,seed2);
x=horzcat(x1,x2);

%Create the Jump Distance
[jd]=JumpDistance1D(x,N);

%Number of Bins for fitting
%Choose option here. 
Nb=round(1+log2(N)); %Sturges Rule
sigma=sqrt(6*(N-2)/(N+1)/(N+3)); %for Doane's rule
%Nb=round(1+log2(N)+log2(1+abs(skewness(jd))/sigma)); %Doanes
%Nb=round(2*(N^(1/3))); %Rice Rule
%Nb=round(sqrt(N)); %square root guidance
%Nb=round((max(jd)-min(jd))*N^(1/3)/(3.5*std(jd))); %Scott's Normal Reference Rule. 
%Nb=round((max(jd)-min(jd))*N^(1/3)/(2*iqr(jd))); %freedman diaconis Rule 

%Plot the Jump Distance
figure(1)
[dr, Ni, yi, ri] =  BinningHist(jd, N, Nb,'yes');

%Plot the predicted JDD on top of it
predictedJDD=fd*N*dr/((pi*D1*tau)^(1/2)).*exp(-ri.^2/(4*D1*tau))+(1-fd)*N*dr/((pi*D2*tau)^(1/2)).*exp(-ri.^2/(4*D2*tau));
hold on
plot(ri,predictedJDD,'k','LineWidth',1.5)

xlabel('Jump Distance')
ylabel('Count')
title('Double Diffusion Jump Distance Distribution in 1D')

%%
%%%%%%%%%%MODEL FITTING%%%%%%%%%%
param = ModelFitting1DwithComboModels(tau, dr, ri, yi, Ni,N, points, dt, x);

%if D1 < D2, change nothing, otherwise flip them....just a precaution. 
if param.D2<param.D1
    test=param.D1;
    param.D1=param.D2;
    param.D2=test;
    param.fdDD=1-param.fdDD;
end

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

%plotting best fit for each model
doublediffusionbest=param.fdDD*N*dr/((pi*param.D1*tau)^(1/2)).*exp(-ri.^2/(4*param.D1*tau))+(1-param.fdDD)*N*dr/((pi*param.D2*tau)^(1/2)).*exp(-ri.^2/(4*param.D2*tau));
hold on
plot(ri,doublediffusionbest,'c','LineWidth',1.5)

%DV
z = -(ri.^2+param.Vdv^2*tau^2)/(4*param.Dvdv*tau);
y = ri*param.Vdv/(2*param.Dvdv);
DVbest=param.fdDV*N*dr/((pi*param.Ddv*tau)^(1/2)).*exp(-ri.^2/(4*param.Ddv*tau))...
    +(1-param.fdDV)*(N*dr/((4*pi*param.Dvdv*tau)^(1/2)).*exp(z+y)+N*dr/((4*pi*param.Dvdv*tau)^(1/2)).*exp(z-y));
hold on
plot(ri,DVbest,'m','LineWidth',1.5)

%Diffusion Anom Combo
fun=@(p) (exp(1i.*p.*tau)).*(1i.*p)^(param.alphada/2-1)/(2.*pi).*exp(-ri./(sqrt(param.Dalphada)).*((1i*p)^(param.alphada/2)));
if param.alphada < 0.5
    min=-300^(.5/param.alphada); 
else
    min=-500;
end
%Plot the predicted JDD on top of it
predictedJDD=fd*N*dr/((pi*param.Dda*tau)^(1/2)).*exp(-ri.^2/(4*param.Dda*tau))+...
    (1-fd)*(2*N*dr/((param.Dalphada)^(1/2)*2).*abs(integral(fun,min,-1*min,'ArrayValued',true,'AbsTol',2e-5)));
hold on
plot(ri,predictedJDD,'y','LineWidth',1.5)

legend('Jump Distance Distribution',['Predicted Doubld Diffusion Fit,D1=',num2str(D1),', D2=',num2str(D2),', fd=',num2str(fd)],...
    ['Fit Diffusion, D=',num2str(param.D)],...
    ['Fit Directed, V=',num2str(param.V),', D_V=',num2str(param.Dv)],...
    ['Fit Anomalous, \alpha=',num2str(param.alpha),', D_\alpha=',num2str(param.Dalpha)],...
    ['Fit Double Diffusion, D1=',num2str(param.D1),', D2=',num2str(param.D2),', fd=',num2str(param.fdDD)],...
    ['Fit DV Fit,D=',num2str(param.Ddv),', V=',num2str(param.Vdv),', Dv=',num2str(param.Dvdv) ', fd=',num2str(param.fdDV)],...
    ['Fit DA Fit,D=',num2str(param.Dda),', D_\alpha=',num2str(param.Dalphada),', \alpha=',num2str(param.alphada) ', fd=',num2str(param.fdDA)])


%%
%%%%%%%%%%BOOTSTRAPPING%%%%%%%%%%

%in this case bootstrapping was used. For this type of system, subsampling
%might be a better alternative. One could even use the median parameter
%value to assess model fitting if it is markedly different. 

%Set Up Storage
Dboot=zeros(numboot,1);
Vboot=zeros(numboot,1); Dvboot=zeros(numboot,1);
Daboot=zeros(numboot,1); Aboot=zeros(numboot,1);
D1boot=zeros(numboot,1); D2boot=zeros(numboot,1); fdDDboot=zeros(numboot,1);
Ddvboot=zeros(numboot,1); Dvdvboot=zeros(numboot,1); Vdvboot=zeros(numboot,1); fdDVboot=zeros(numboot,1);
Ddaboot=zeros(numboot,1); Dalphadaboot=zeros(numboot,1); alphadaboot=zeros(numboot,1); fdDAboot=zeros(numboot,1);

parfor i=1:numboot
    X = randi(N,N,1);
    jdB=jd(X);
    %If using Doanes Rule, can consider assessing if a different Nb is
    %needed based on skewness, or can stay with original choice 
    %Nb=round(1+log2(N)+log2((1+skewness(jdB))/(sqrt((6*(N-2))/((N+1)*(N+3))))));
    [drB, NiB, yiB, riB] =  BinningHist(jdB, N, Nb,'no');
    paramB = ModelFitting1DwithComboModels(tau, drB, riB, yiB, NiB,N, points, dt, x);
    if paramB.D2<paramB.D1
        test=paramB.D1;
        paramB.D1=paramB.D2; paramB.D2=test; paramB.fdDD=1-paramB.fdDD;
    end
    Dboot(i)=paramB.D;
    Vboot(i)=paramB.V; Dvboot(i)=paramB.Dv;
    Daboot(i)=paramB.Dalpha; Aboot(i)=paramB.alpha;
    D1boot(i)=paramB.D1; D2boot(i)=paramB.D2; fdDDboot(i)=paramB.fdDD;
    Ddvboot(i)=paramB.Ddv; Vdvboot(i)=paramB.Vdv; Dvdvboot(i)=paramB.Dvdv; fdDVboot(i)=paramB.fdDV;
    Ddaboot(i)=paramB.Dda; Dalphadaboot(i)=paramB.Dalphada; alphadaboot(i)=paramB.alphada; fdDAboot(i)=paramB.fdDA;
end

beta=[param.D,param.V,param.Dv,param.Dalpha,param.alpha,...
    param.D1,param.D2,param.fdDD,...
    param.Ddv,param.Vdv,param.Dvdv,param.fdDV,...
    param.Dda,param.Dalphada,param.alphada,param.fdDA];
dbeta=2*[std(Dboot),std(Vboot),std(Dvboot),std(Daboot), std(Aboot),...
    std(D1boot),std(D2boot),std(fdDDboot),...
    std(Ddvboot),std(Vdvboot),std(Dvdvboot),std(fdDVboot),...
    std(Ddaboot),std(Dalphadaboot),std(alphadaboot),std(fdDAboot)];

%%
%%%%%%%%%%MODEL SELECTION%%%%%%%%%%
[prob,value,method]=Integration1DwithComboModels(dbeta,beta,N,yi,ri,dr,tau);
