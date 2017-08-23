%Parameters and Preliminaries
rng(15)
dt=23.58;
dt=1;
tau=20*dt;
points=round(tau/dt)+1; %number of points in each trajectory
N=3000; %number of trajectories total
Nb=round(N/100); %number of bins %Nb=round(N/100)
Dalpha=.5; %micro meters^2/s^alpha For A and Da
alpha=0.5; %For A and DA
%Simulate Data
[x,y,z]=AnomalousDiffusion3D(Dalpha,alpha,points,N,dt,tau,randi([1,N],1,1));
%[x,y]=AnomalousDiffusion2(Dalpha,alpha,points,N,dt);
[jd]=JumpDistance3D(x,y,z,N); %jd is a vertical vector
[dr, Ni, yi, ri] =  BinningHist(jd, N, Nb,'yes');
%parameter fitting, leave out for now until it works....
%param = ModelFitting(tau, dr, ri, yi, Ni, N, points);
if alpha < 0.5
    min=-700^(.5/alpha); %limits on inverse laplace transform
else
    min=-500;
end
fun=@(p) (exp(1i.*p.*tau)).*(1i.*p)^(alpha-1)/(2.*pi).*(exp(-ri./(sqrt(Dalpha)).*((1i*p)^(alpha/2))));
min=-700;
hold on

%try splitting it up into 200 point intevals
%first find the total number of intevals needed
intsize=200;
numint=ceil(abs(min+100)/intsize);
testtest=zeros(1,Nb);

for i=1:numint
    testtest=testtest+2*N*dr*ri/(Dalpha).*abs(integral(fun,min+(i-1)*intsize-1i*1e-6,min+i*(intsize)-1i*1e-6,...
        'ArrayValued',true,'AbsTol',1e-6,'RelTol',1e-3));
end
testtest=testtest+N*dr*ri/(Dalpha).*abs(integral(fun,-100-1i*1e-6,100-1i*1e-6,...
    'ArrayValued',true,'AbsTol',1e-6,'RelTol',1e-3));

plot(ri,testtest,'r','LineWidth',1.5)

param = ModelFittingA3D(tau, dr, ri, yi, Ni, N, points,dt,x,y,z)
param2= ModelFitting3D(tau, dr, ri, yi, Ni, N, points,dt,x,y,z)
%param3= ModelFittingA3D2(tau, dr, ri, yi, Ni, N, points,dt,x,y,z)

fun2=@(p) (exp(1i.*p.*tau)).*(1i.*p)^(param.alpha-1)/(2.*pi).*...
    (exp(-ri./(sqrt(param.Dalpha)).*((1i*p)^(param.alpha/2))));
testtest2=zeros(size(testtest));
for i=1:numint
    testtest2=testtest2+2*N*dr*ri/(param.Dalpha).*abs(integral(fun2,min+(i-1)*intsize-1i*1e-6,min+i*(intsize)-1i*1e-6,...
        'ArrayValued',true,'AbsTol',1e-6,'RelTol',1e-3));
end
testtest2=testtest2+N*dr*ri/(param.Dalpha).*abs(integral(fun2,-100-1i*1e-6,100-1i*1e-6,...
    'ArrayValued',true,'AbsTol',1e-6,'RelTol',1e-3));

plot(ri,testtest2,'b','LineWidth',1.5)

legend('Simulated','Fitted')

% %plotting the tails of the integrand.
% %first set the bin that you want to look at
% k=1;
% %set the minimum and maximum for your range and # of points
% mini=300;
% maxi=700;
% points=10000;
% pts=linspace(mini,maxi,points);
% fun3=@(p) (exp(1i.*p.*tau)).*(1i.*p)^(alpha-1)/(2.*pi).*(exp(-ri(k)./(sqrt(Dalpha)).*((1i*p)^(alpha/2))));
% integranddata=zeros(points,1);
% for i=1:points
%     integranddata(i)=fun3(pts(i));
% end
% 
% figure(2)
% hold on
% plot(pts,imag(integranddata),'r')
% axis([mini maxi min(imag(integranddata)) max(imag(integranddata))])


%setpref('Internet','E_mail','beccagobel@gmail.com');
mail = 'beccagobel@gmail.com'; %Your GMail email address
password = 'bigmouth65';  %Your GMail password
setpref('Internet','SMTP_Server','smtp.gmail.com');
setpref('Internet','SMTP_Username',mail);
setpref('Internet','SMTP_Password',password);
props = java.lang.System.getProperties;
props.setProperty('mail.smtp.auth','true');
props.setProperty('mail.smtp.socketFactory.class', 'javax.net.ssl.SSLSocketFactory');
props.setProperty('mail.smtp.socketFactory.port','465');
sendmail('beccagobel@gmail.com','Hello From MATLAB!');



