%BinningHist
%Way to create the JDD histogram given the jdd data points
%Rebecca Menssen
%Last Updated: 8/30/17

%This function takes the jdd points and turns them into a histogram. 

%%%%%%%%%%INPUTS%%%%%%%%%%
% jd--a vector of the jump distance for each trajectory
% N--the number of trajectories
% Nb--the number of bins to be used
% show--create and save(if wanted) a histogram object of the JDD

%%%%%%%%%%OUTPUTS%%%%%%%%%%
% dr--the width of the bins
% Ni--vector of the number of trajectories in each bin
% yi--vector of the number of trajectories in each bin as a fraction of the
% total
% ri--vector of the midpoints of the bins

function [dr, Ni, yi, ri] =  BinningHist(jd, N, Nb, show)
%find all of the data for the histogram
%counts and edges
[Ni, edges] = histcounts(jd, Nb);
%spacing
dr = edges(2)-edges(1);
%probabilities
yi = Ni./N;
%bin centers
ri = edges(1:end-1) + dr/2;

%show the histogram image
if strcmpi(show, 'yes')
    %plot the histogram
    JDDhist = histogram(jd);
    JDDhist.NumBins = Nb;
    %if you want to save the figure, uncomment this feature
    %saveas(JDDhist, 'jdd.fig');
end
end











