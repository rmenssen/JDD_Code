%Simulate 2D Directed Diffusion
%Rebecca Menssen
%Last edited: 8/30/17-cleaning up notation and adding final comments

%This code simulates directed diffusion in a two dimensions.

%%%%%%%%%%INPUTS%%%%%%%%%%
%V--directed motion constant (given in micro m/s)
%Dv--diffusion constant (given in micro m^2/s)
%points--the length of the simulation, including the start point
%N--the total number of trajectories to simulate
%dt--time step (given in seconds)
%seed--seed for the simulation (for repeatability)

%%%%%%%%%%OUTPUTS%%%%%%%%%%
%x,y--arrays of simulated trajectories. N columns, points rows

function [x,y]=DirectedMotion2D(V,Dv,points,N,dt,seed)
%set seed
rng(seed)
%set direction of motion
theta=2*pi*rand;  %normal around zero
%create arrays to store data
x=zeros(points, N);
y=zeros(points,N);
%set spatial variance
var=2*Dv*dt;
%simulate trajectories
for i=1:N
    for j=1:points-1
        x(j+1,i)=x(j,i)+V*dt*cos(theta)+sqrt(var)*randn;
        y(j+1,i)=y(j,i)+V*dt*sin(theta)+sqrt(var)*randn;
    end
end    
end