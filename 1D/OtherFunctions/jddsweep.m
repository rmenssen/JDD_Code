function [jd]=jddsweep(x,timelag)
%function takes in data x, and creates a jdd structure
%first find dimensions of the system

numsteps=size(x,1);
numtraj=size(x,2);

jd=[];
for i=1:numtraj;
    for j=1:numsteps-timelag;
        jd=vertcat(jd,sqrt((x(j+timelag,i)-x(j,i))^2));
    end
end

end