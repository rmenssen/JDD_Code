%Data Creation for 3D Directed Motion

%data creation for Directed Motion
%Rebecca Menssen
%7/6/15

function [x,y,z]=DirectedMotion3D(V,kv,points,N,dt)
%direction of motion
theta=2*pi*rand;  %uniform from 0 to 2pi
phi=acos(2*rand-1); %uniform from 0 to pi(this is the proper way to distribute points)
%create array to store data
x=zeros(points, N);
y=zeros(points,N);
z=zeros(points,N);
var=kv;
%simulate data trajectories
for i=1:N
    for j=1:points-1
        x(j+1,i)=x(j,i)+V*dt*cos(theta)*sin(phi)+sqrt(var)*randn;
        y(j+1,i)=y(j,i)+V*dt*sin(theta)*sin(phi)+sqrt(var)*randn;
        z(j+1,i)=z(j,i)+V*dt*cos(phi)+sqrt(var)*randn;
    end
end 
end