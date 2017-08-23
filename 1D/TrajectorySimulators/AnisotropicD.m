%Anisotropic Diffusion
%Rebecca Menssen
function[x,y,z]=AnisotropicD(D1,D2,D3,points,N,dt)
%create array to store data
x=zeros(points, N);
y=zeros(points,N);
z=zeros(points,N);
var1=2*D1*dt;
var2=2*D2*dt;
var3=2*D3*dt;
for i=1:N
    for j=1:points-1
        x(j+1,i)=x(j,i)+sqrt(var1)*randn; %each column is a run
        y(j+1,i)=y(j,i)+sqrt(var2)*randn;
        z(j+1,i)=z(j,i)+sqrt(var3)*randn;
    end
end
end