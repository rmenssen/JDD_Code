%data creation for Brownian Motion
%Rebecca Menssen
%7/6/15

function [x,y,z]=AnomalousDiffusion3D(Dalpha,alpha,points,N,dt,endtime,seed)
rng(seed)
%primed time (for sampling)
%need to sample times for a distribution with a long tail for longer times
%and a shorter tail for earlier times
%once you have the waiting time, steps are sampled from a normal with
%variance 2Dalpha(dtprime)^alpha
dtprime=dt/2000;
%create array to store data
x=zeros(points, N);
y=zeros(points,N);
z=zeros(points,N);
var=2*Dalpha*dtprime^(alpha); %variance of spatial step
%var=2*Dalpha*dt^(alpha);
times=(0:dt:endtime);
rng('shuffle')
xi=dtprime;
%simulate data trajectories
for i=1:N
    timecount=0;
    xcountold=0;%record x position
    ycountold=0;%record y position
    zcountold=0; %record z position
    for j=1:points-1
        while timecount < times(j+1)
            newtime=-xi*log(rand)*(sin(alpha*pi)/tan(alpha*pi*rand)-cos(alpha*pi))^(1/alpha);
            timecount=timecount+newtime;
            %update position
            xcount=xcountold+sqrt(var)*randn;
            ycount=ycountold+sqrt(var)*randn;
            zcount=zcountold+sqrt(var)*randn;
            nextstep=j+1; %check the timecount
            for k=nextstep:points %check the next step until the end
                if timecount >times(k) %could skip over multiple steps so have to consider all points
                    x(k,i)=xcountold; %need old time step since you have skipped over in the next step
                    y(k,i)=ycountold;
                    z(k,i)=zcountold;
                end
            end
            xcountold=xcount;
            ycountold=ycount;
            zcountold=zcount;
        end
    end
end
end