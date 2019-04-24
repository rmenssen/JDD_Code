%Simulate 1D Directed Diffusion
%Rebecca Menssen
%Last edited: 4/9/19 (comments)

%This code simulates directed diffusion in a single dimension.

%%%%%%%%%%INPUTS%%%%%%%%%%
%V--directed motion constant (given in micro m/s)
%Dv--diffusion constant (given in micro m^2/s)
%points--the length of the simulation, including the start point
%N--the total number of trajectories to simulate
%dt--time step (given in seconds)
%seed--seed for the simulation (for repeatability)

%%%%%%%%%%OUTPUTS%%%%%%%%%%
%x--array of simulated trajectories. N columns, points rows

function [x]=DirectedMotion1D(V,Dv,points,N,dt,seed)
%set seed
rng(seed)
%create array to store data
x=zeros(points, N);
%set spatial variance
var=2*Dv*dt;
%simulate trajectories
for i=1:N
    for j=1:points-1
        x(j+1,i)=x(j,i)+V*dt+sqrt(var)*randn;
    end
end       
end