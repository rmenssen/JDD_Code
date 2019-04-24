%Parameter fitting for 1D Jump Distance Distributions
%Rebecca Menssen
%Last Updated: 4/9/19

%This function takes JDD data as the input and outputs the parameters for
%three different models using lsqcurvefit with weighting. This function is
%specifically for the directed diffusion model. 

%%%%%%%%%%INPUTS%%%%%%%%%%
%tau--the duration of each trajectory (the time lag*dt)
%dr--the width of each bin in the JDD histogram
%ri--vector of the midpoints of the JDD bins
%yi--vector of the percentage of trajectories in each JDD bin
%N--the number of trajectories
%points--the number of time points in each trajectory
%dt--the time step
%x1--the 1D trajectory

%%%%%%%%%%OUTPUTS%%%%%%%%%%
%param: strut containing parameters determined by the LSQ fit for each model.

function param = ModelFitting1DV(tau, dr, ri, yi, Ni ,N, points, dt, x1)
%Set up optimization options: can be customized for better fits (espeically
%tolerances)
options = optimoptions('lsqnonlin', 'Algorithm', 'levenberg-marquardt',...
    'MaxFunctionEvaluations',1000,'FunctionTolerance',1e-4,...
    'StepTolerance',1e-4);

%Create struct to hold all the parameters
param = struct('D',NaN,'V', NaN, 'Dv', NaN,'Dalpha',NaN,'alpha',NaN);

%Weighting vector
%We chose to weight by 1/bin count probabilties (Ni/N=bin count
%probabilities). Each bin was given an artificial additional count so no
%bin was empty
wi=1./((Ni+1)/(N+length(Ni)));
weights=wi;
%weights=1;

%Getting initial seeds by using MSD Analysis
%make time vector
t=dt*(0:points-1);

%finding avg(x^2). Have to account for possible missing data. Outliers can
%mess this calculation up, but there isn't too much that can be done in
%that case, since using nanmedian causes other issues. 
xsqavg=nanmean((x1-repmat(x1(1,:),size(x1,1),1)).^2,2);
%standard deviation
%stdev=nanstd(((x1-repmat(x1(1,:),size(x1,1),1)).^2)');

t=t';

%Directed Diffusion Seeding
%[p]=polyfit(t, xsqavg,2);
[p]=polyfitZero(t,xsqavg,2);
if p(1)< 0
    param = fmincon(@(x) norm(x(1)^2*t.^2 + 2*x(2)*t - xsqavg),[0.1 0.1],...
        [],[],[],[],[0 0]);
    msdV = param(1);
    msdDV = param(2);
else
    msdV = sqrt(abs(p(1)));
    msdDV = abs(p(2)/2);
end
%msdV=sqrt(abs(p(1)));
%msdDV=abs(p(2)/2);

%Testing 1D Directed Motion Model
x0 = [msdV, msdDV];
[temp] = lsqnonlin(@NV, x0, [], [], options);
param.V = temp(1);
param.Dv = temp(2);

    function [outvals] = NV(x)
        %Generates the values needed to perform the LSQ non-linear fit for
        %the directed motion model against a given JDD.
        
        %Input Parameters
        V = x(1);
        D = x(2);
        
        %predicted JDD probabilities based on input parameters
        z = -(ri.^2+V^2*tau^2)/(4*D*tau);
        y = ri*V/(2*D);
        predicted = dr/((4*pi*D*tau)^(1/2)).*exp(z+y)+dr/((4*pi*D*tau)^(1/2)).*exp(z-y);
        
        %actual JDD probabilities
        actual = yi;
        
        %Weighted Least Squares Function
        outvals = sqrt(weights) .* (predicted - actual);
    end
end
