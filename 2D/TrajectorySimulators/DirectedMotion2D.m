%data creation for Directed Motion
%Rebecca Menssen
%7/6/15

function [x,y]=DirectedMotion2D(V,kv,points,N,dt)
theta=2*pi*rand;  %normal around zero
%create array to store data
x=zeros(points, N);
y=zeros(points,N);
var=kv;
%simulate data trajectories
for i=1:N
    for j=1:points-1
        x(j+1,i)=x(j,i)+V*dt*cos(theta)+sqrt(var)*randn;
        y(j+1,i)=y(j,i)+V*dt*sin(theta)+sqrt(var)*randn;
    end
end    
end