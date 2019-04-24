%Simulation of 2D Anomalous Diffusion
%Rebecca Menssen
%Last edited: 8/30/17-adding final comments and clearing up notation

%This code simulates anomalous diffusion in two dimensions

%%%%%%%%%%INPUTS%%%%%%%%%%
%alpha--anomalous exponent (given in micro m/s)
%Dalpha--diffusion constant (given in micro m^2/s^alpha)
%points--the length of the simulation, including the start point
%N--the total number of trajectories to simulate
%dt--time step (given in seconds)
%endtime--the endtime of the simulation. Should be equal to (points-1)*dt
%seed--seed for the simulation (for repeatability)

%%%%%%%%%%OUTPUTS%%%%%%%%%%
%x,y--arrays of simulated trajectories

function [x,y]=AnomalousDiffusion2D(Dalpha,alpha,points,N,dt,endtime,seed)
%set seed
rng(seed)

%For the continuous time random walk, we need to set a primed time for our
%sampling rate. Essentially, this is the time that is used for setting the
%variance of a spatial step when a step is made. This primed time is also
%used in how we select the time of our spatial step (parameter xi)
dtprime=dt/2000;
xi=dtprime;

%create array to store data
x=zeros(points, N);
y=zeros(points,N);
%once you have the waiting time, steps are sampled from a normal with
%variance 2Dalpha(dtprime)^alpha
var=2*Dalpha*dtprime^(alpha);

%set up an array with the specific times that need to be recorded
times=(0:dt:endtime);

%simulate data trajectories using a continuous time random walk
for i=1:N
    i
    timecount=0; %current time
    xcountold=0; %record x position
    ycountold=0; %record y position
    for j=1:points-1
        while timecount < times(j+1)
            %find the time of the next move (time step)
            newtime=-xi*log(rand)*(sin(alpha*pi)/tan(alpha*pi*rand)-cos(alpha*pi))^(1/alpha);
            %update the time
            timecount=timecount+newtime;
            %update position
            xcount=xcountold+sqrt(var)*randn;
            ycount=ycountold+sqrt(var)*randn;
            %index of next timepoint
            nextstep=j+1; 
            %Now look and see if you have skipped over any time steps and
            %fill in those positions with the old x position. This is how
            %you get positions for each time. 
            for k=nextstep:points 
                if timecount >times(k)
                    x(k,i)=xcountold; 
                    y(k,i)=ycountold;
                end
            end
            %update the old position with the new one
            xcountold=xcount;
            ycountold=ycount;
        end
    end
end
end