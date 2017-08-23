%data creation for Directed Motion
%Rebecca Menssen
%7/6/15

function [x]=DirectedMotion1D(V,kv,points,N,dt,j)
rng(j)
%direction of motion is just going to be set to +x (since it doesn't really
%matter)
%create array to store data
x=zeros(points, N);
var=kv;
%simulate data trajectories
for i=1:N
    for j=1:points-1
        %have V*dt in the positive X direction, and a random component
        %afterwards
        x(j+1,i)=x(j,i)+V*dt+sqrt(var)*randn;
    end
end       
end