%JDD analysis of data.
%We used this code to do our analysis on JDDs of chemotaxis trajectories.
%The data is not included here, but the general form can be used as a
%template. 

%initially data is inputting, including a time step

x=data1;
y=data2;

timelag=10;
points=timelag+1;
tau=dt*timelag;

%do a sliding JDD
jd=jddsweep2D(x,y,timelag);
jd=jd(:);

N=length(jd);

%Number of Bins for fitting
Nb=50;

%Number of Bootstraps
numboot=100;

figure
[dr, Ni, yi, ri]=BinningHist(jd, N, Nb,'yes');

%%
%%%%%%%%%%MODEL FITTING%%%%%%%%%%
maxtimestep=1000;
xsqavg=zeros(maxtimestep+1,1);
for i = 1:maxtimestep
    i
    jd2=jddsweep2Dupdate(x,y,i);
    jd2=jd2(:);
    jd2 = jd2.^2;
    xsqavg(i+1) = nanmean(jd2);
end

MSD=xsqavg;
N=length(jd);

%if not using combination models
%param = ModelFitting2Dexperimentaldata(tau, dr, ri, yi, Ni, N, maxtimestep,dt,MSD);
%if using combination models but using MSD for seeding
%param = ModelFitting2DexperimentaldataCombo(tau, dr, ri, yi, Ni, N, maxtimestep,dt,MSD);
%if using seeds...
seeds=beta; %we used our seed as the time lag before it. 
param = ModelFitting2DexperimentaldataCombowithSeeds(tau, dr, ri, yi, Ni, N, maxtimestep,dt,MSD,seeds);
    
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

%If putting combination models to the test, plot them here.
%double diffusion
predicted1 = dr*ri/(2*param.D1*tau).*exp(-ri.^2/(4*param.D1*tau));
predicted2 = dr*ri/(2*param.D2*tau).*exp(-ri.^2/(4*param.D2*tau));
DoubleDiffusionBest=N*(param.fdDD*predicted1+(1-param.fdDD)*predicted2);
hold on
plot(ri,DoubleDiffusionBest,'c','LineWidth',1.5)

%diffusion directed motion
predicted1 = dr*ri/(2*param.Dddv*tau).*exp(-ri.^2/(4*param.Dddv*tau));
z = -(ri.^2+param.Vddv^2*tau^2)/(4*param.DVddv*tau);
y = ri*param.Vddv/(2*param.DVddv);
predicted2 = dr*ri/(4*param.DVddv*tau).*exp(z).*besseli(0, y);
DVbest=N*(param.fdDDV*predicted1+(1-param.fdDDV)*predicted2);
plot(ri,DVbest,'m','LineWidth',1.5)

%Diffusion Anomalous Diffusion
if param.Ada < 0.5
    min=-300^(.5/param.Ada); %limits on inverse laplace transform
else
    min=-300;
end
predicted1 = dr*ri/(2*param.Dda*tau).*exp(-ri.^2/(4*param.Dda*tau));
fun=@(p) (exp(1i.*p.*tau)).*(1i.*p)^(param.Ada-1)/(2.*pi).*...
    (besselk(0,ri./(sqrt(param.DAda)).*((1i*p)^(param.Ada/2))));
predicted2=dr*ri/(param.DAda).*abs(integral(fun,min,-1*min,...
    'ArrayValued',true,'AbsTol',1e-6));
DAbest=N*(param.fdDDA*predicted1+(1-param.fdDDA)*predicted2);
plot(ri,DAbest,'y','LineWidth',1.5)

%directed diffusion, anomalous diffusion
if param.Advda < 0.5
    min=-300^(.5/param.Advda); %limits on inverse laplace transform
else
    min=-300;
end
z = -(ri.^2+param.Vdvda^2*tau^2)/(4*param.DVdvda*tau);
y = ri*param.Vdvda/(2*param.DVdvda);
predicted1 = dr*ri/(4*param.DVdvda*tau).*exp(z).*besseli(0, y);
fun=@(p) (exp(1i.*p.*tau)).*(1i.*p)^(param.Advda-1)/(2.*pi).*...
    (besselk(0,ri./(sqrt(param.DAdvda)).*((1i*p)^(param.Advda/2))));
predicted2=dr*ri/(param.DAdvda).*abs(integral(fun,min,-1*min,...
    'ArrayValued',true,'AbsTol',1e-6));
VAbest=N*(param.fdDVDA*predicted1+(1-param.fdDVDA)*predicted2);
plot(ri,VAbest,'k','LineWidth',1.5)

legend('Jump Distance Distribution','Fit Diffusion', 'Fit Directed Motion',...
    'Fit Anomalous', 'Fit Double Diffusion', 'Fit Diffusion-Directed',...
    'Fit Diffusion-Anomalous', 'Fit Directed-Anomalous')


%%
%%%%%%%%%%BOOTSTRAPPING%%%%%%%%%%

%Set Up Storage (so many)
Dboot=zeros(numboot,1);
Vboot=zeros(numboot,1); Dvboot=zeros(numboot,1);
Daboot=zeros(numboot,1); Aboot=zeros(numboot,1);
D1boot=zeros(numboot,1); D2boot=zeros(numboot,1); fdDDboot=zeros(numboot,1);
Dddvboot=zeros(numboot,1); Vddvboot=zeros(numboot,1);
DVddvboot=zeros(numboot,1); fdDDVboot=zeros(numboot,1);

Ddaboot=zeros(numboot,1); DAdaboot=zeros(numboot,1);
Adaboot=zeros(numboot,1); fdDDAboot=zeros(numboot,1);

Vdvdaboot=zeros(numboot,1); DVdvdaboot=zeros(numboot,1); DAdvda=zeros(numboot,1);
Advdaboot=zeros(numboot,1); fdDVDAboot=zeros(numboot,1);

beta=cell2mat(struct2cell(param)); beta=beta';
seeds=beta;
parfor i=1:numboot
    X = randi(N,N,1);
    jdB=jd(X);
    [drB, NiB, yiB, riB] =  BinningHist(jdB, N, Nb,'no');
    %paramB = ModelFitting2DexperimentaldataCombo(tau, drB, riB, yiB, NiB, N, maxtimestep,dt,MSD);
    paramB = ModelFitting2DexperimentaldataCombowithSeeds(tau, drB, riB, yiB, NiB, N, maxtimestep,dt,MSD,seeds);
    Dboot(i)=paramB.D;
    Vboot(i)=paramB.V; Dvboot(i)=paramB.Dv;
    Daboot(i)=paramB.Dalpha; Aboot(i)=paramB.alpha;
   
    D1boot(i)=paramB.D1; D2boot(i)=paramB.D2; fdDDboot(i)=paramB.fdDD;
    
    Dddvboot(i)=paramB.Dddv; Vddvboot(i)=paramB.Vddv; DVddvboot(i)=paramB.DVddv;
    fdDDVboot(i)=paramB.fdDDV;
    
    Ddaboot(i)=paramB.Dda; DAdaboot(i)=paramB.DAda; Adaboot(i)=paramB.Ada; 
    fdDDAboot(i)=paramB.fdDDA;
    
    Vdvdaboot(i)=paramB.Vdvda; DVdvdaboot(i)=paramB.DVdvda; 
    DAdvdaboot(i)=paramB.DAdvda; Advdaboot(i)=paramB.Advda; fdDVDAboot(i)=paramB.fdDVDA;
    
end

dbeta=2*[std(Dboot),std(Vboot),std(Dvboot),std(Daboot), std(Aboot),...
        std(D1boot), std(D2boot), std(fdDDboot),...
        std(Dddvboot), std(Vddvboot), std(DVddvboot), std(fdDDVboot),...
        std(Ddaboot), std(DAdaboot), std(Adaboot), std(fdDDAboot),...
        std(Vdvdaboot), std(DVdvdaboot), std(DAdvdaboot), std(Advdaboot), std(fdDVDAboot)];
combine=vertcat(beta,dbeta); %just for comparison of sizes
combine2=vertcat(beta,dbeta,dbeta./beta);
%%
%%%%%%%%%%MODEL SELECTION%%%%%%%%%%
[prob,value,method]=Integration2Dwithcombomodels(dbeta,beta,N,yi,ri,dr,tau);

%save data
%want to save beta, dbeta, jd and all binning hist stuff. 
save('Results.mat','beta','dbeta','jd','dr','ri','yi','Ni','prob')



