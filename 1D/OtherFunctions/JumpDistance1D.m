%Calculate Jump Distances
%Rebecca Menssen
%7/7/15
function [jd]=JumpDistance1D(x,N)
jd=zeros(N,1);
for i=1:N
    %calculate the jump distance for each trajectory
    jd(i,1)=sqrt((x(end,i)-x(1,i))^2);
end
end