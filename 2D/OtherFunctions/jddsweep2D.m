%Creating MSDs and bootstrapped MSDs

%JDD Sweep Function for 2 dimensions
%Rebecca Menssen
%Last Updated: 8/30/18

%This function is to be used when you want to do a sweeping JDD, that is
%not ignore data between points, using a sliding window to calculate the
%jump distance for every set of points a time lag L away.

%%%%%%%%%%INPUTS%%%%%%%%%%
%x,y--the data structured with each trajectory in a single column
%timelag--the timelag that is being used. Split by dimension

%%%%%%%%%%OUTPUTS%%%%%%%%%%
%jd--list of all jump distances

function [jd]=jddsweep2D(x,y,timelag)
%empty struct for the jump distances
jd=[];
jd = vertcat(jd,sqrt((x(1+timelag:end,:)-x(1:end-timelag,:)).^2 ...
    + (y(1+timelag:end,:)-y(1:end-timelag,:)).^2));
jd=jd(:);
jd(isnan(jd))=[];
end