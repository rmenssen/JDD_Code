function [dr, Ni, yi, ri] =  BinningHist(jd, N, Nb, keep)
% [dr, Ni, yi, ri] = BinningHist(jd, N, Nb)

%{
Outputs:
dr = the width of the bins
Ni = vector of the number of trajectories in each bin
yi = vector of the number of trajectories in each bin as a fraction of the
total
ri = vector of the midpoints of the bins

Inputs:
jd = a vector of the jump distance for each trajectory
N = the number of trajectories
Nb = the number of bins to be used
keep = create and save a histogram object of the JDD?

%}

[Ni, edges] = histcounts(jd, Nb);
dr = edges(2)-edges(1);
yi = Ni./N;
ri = edges(1:end-1) + dr/2;

if strcmpi(keep, 'yes')
    JDDhist = histogram(jd);
    JDDhist.NumBins = Nb;
    saveas(JDDhist, 'jdd.fig');
end

end











