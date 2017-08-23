%data creation for Brownian Motion
%Rebecca Menssen
%7/6/15
function [x,y]=Diffusion2D(D,points,N,dt)
%create array to store data
x=zeros(points, N);
y=zeros(points,N);
var=2*D*dt;
%simulate data trajectories
for i=1:N
    for j=1:points-1
        x(j+1,i)=x(j,i)+sqrt(var)*randn;
        y(j+1,i)=y(j,i)+sqrt(var)*randn;
    end
end
end