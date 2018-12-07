%Creating MSDs and bootstrapped MSDs

%JDD Sweep Function for 2 dimensions
%Rebecca Menssen
%Last Updated: 8/30/17

%This function is to be used when you want to do a sweeping JDD, that is
%not ignore data between points, using a sliding window to calculate the
%jump distance for every set of points a time lag L away.

%%%%%%%%%%INPUTS%%%%%%%%%%
%x,y--the data structured with each trajectory in a single column
%timelag--the timelag that is being used. Split by dimension

%%%%%%%%%%OUTPUTS%%%%%%%%%%
%jd--list of all jump distances

function [jd]=jddsweep2Dintervals(x,y,timelag)
%empty struct for the jump distances
jd=[];
Idx = 1:timelag:size(x,1);
jd = vertcat(jd,...
    sqrt((x(Idx(2:end),:)-x(Idx(1:end-1),:)).^2 ...
    + (y(Idx(2:end),:)-y(Idx(1:end-1),:)).^2));
end