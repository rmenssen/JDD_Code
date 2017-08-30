%Simulate 1D Directed Diffusion
%Rebecca Menssen
%Last edited: 8/30/17-cleaning up notation and adding final comments

%This code simulates directed diffusion in a single dimension.

%%%%%%%%%%INPUTS%%%%%%%%%%
%V--directed motion constant (given in micro m/s)
%Dv--diffusion constant (given in micro m^2/s)
%points--the length of the simulation, including the start point
%N--the total number of trajectories to simulate
%dt--time step (given in seconds)
%seed--seed for the simulation (for repeatability)

%%%%%%%%%%OUTPUTS%%%%%%%%%%
%x,y,z--arrays of simulated trajectories. N columns, points rows

function [x,y,z]=DirectedMotion3D(V,Dv,points,N,dt,seed)
%set seed
rng(seed)
%set direction of motion
theta=2*pi*rand;  %uniform from 0 to 2pi
phi=acos(2*rand-1); %uniform from 0 to pi
%create arrays to store data
x=zeros(points, N);
y=zeros(points,N);
z=zeros(points,N);
%set spatial variance
var=2*Dv*dt;
%simulate trajectories
for i=1:N
    for j=1:points-1
        x(j+1,i)=x(j,i)+V*dt*cos(theta)*sin(phi)+sqrt(var)*randn;
        y(j+1,i)=y(j,i)+V*dt*sin(theta)*sin(phi)+sqrt(var)*randn;
        z(j+1,i)=z(j,i)+V*dt*cos(phi)+sqrt(var)*randn;
    end
end 
end