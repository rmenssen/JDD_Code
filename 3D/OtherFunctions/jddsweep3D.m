%JDD Sweep Function for 3 dimensions
%Rebecca Menssen
%Last Updated: 8/30/17

%This function is to be used when you want to do a sweeping JDD, that is
%not ignore data between points, using a sliding window to calculate the
%jump distance for every set of points a time lag L away. 

%%%%%%%%%%INPUTS%%%%%%%%%%
%x,y,z--the data structured with each trajectory in a single column
%timelag--the timelag that is being used. Split by dimension

%%%%%%%%%%OUTPUTS%%%%%%%%%%
%jd--list of all jump distances

function [jd]=jddsweep3D(x,y,timelag)

numsteps=size(x,1); %number of rows
numtraj=size(x,2); %number of columns

%empty struct for the jump distances
jd=[];

for i=1:numtraj
    for j=1:numsteps-timelag
        %calculate jump distance
        jd=vertcat(jd,sqrt((x(j+timelag,i)-x(j,i))^2+...
            (y(j+timelag,i)-y(j,i))^2+(z(j+timelag,i)-z(j,i))^2));
    end
end
end