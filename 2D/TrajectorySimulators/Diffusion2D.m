%Simulate 2D Pure Diffusion
%Rebecca Menssen
%Last edited: 8/30/17-adding final comments and description

%This code simulates pure diffusion in a two dimensions.

%%%%%%%%%%INPUTS%%%%%%%%%%
%D--diffusion constant (given in micro m^2/s)
%points--the length of the simulation, including the start point
%N--the total number of trajectories to simulate
%dt--time step (given in seconds)
%seed--seed for the simulation (for repeatability)

%%%%%%%%%%OUTPUTS%%%%%%%%%%
%x,y--arrays of simulated trajectories. N columns, points rows

function [x,y]=Diffusion2D(D,points,N,dt,seed)
%set seed
rng(seed)
%create array to store data
x=zeros(points, N);
y=zeros(points,N);
%set spatial variance
var=2*D*dt;
%simulate trajectories
for i=1:N
    for j=1:points-1
        x(j+1,i)=x(j,i)+sqrt(var)*randn;
        y(j+1,i)=y(j,i)+sqrt(var)*randn;
    end
end
end