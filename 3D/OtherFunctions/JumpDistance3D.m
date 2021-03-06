%Calculate Jump Distances
%Rebecca Menssen
%Last Updated: 8/30/17

%This function calculates the jump distance in three dimensions. It does not
%allow for a sliding window, it just takes the first and last points. Thus
%it also assumes all trajectories are of length timelag +1

%%%%%%%%%%INPUTS%%%%%%%%%%
%x,y,z--trajectory for each dimension
%N--number of trajectories

%%%%%%%%%%OUTPUTS%%%%%%%%%%
%jd--list of all jump distances

function [jd]=JumpDistance3D(x,y,z,N)
jd=zeros(N,1);
for i=1:N
    %calculate the jump distance for each trajectory
    jd(i,1)=sqrt((x(end,i)-x(1,i))^2+(y(end,i)-y(1,i))^2+(z(end,i)-z(1,i))^2);
end
end