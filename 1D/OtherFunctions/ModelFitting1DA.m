%Parameter fitting for 1D Jump Distance Distributions
%Rebecca Menssen
%Last Updated: 4/24/19

%This function takes JDD data as the input and outputs the parameters for
%three different models using lsqcurvefit with weighting. This code is
%specific to anomalous subdiffusion

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

function param = ModelFitting1DA(tau, dr, ri, yi, Ni ,N, points, dt, x1)
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

%Anomalous Diffusion Seeding
[p]=polyfit(log(t(2:end)), log(xsqavg(2:end)),1);
msdA=abs(p(1));
msdDA=abs(exp(p(2)-log(2)+log(gamma(1+msdA))));
% [Param,ResNorm] = lsqnonlin(@(x)Anomalous_MSD1D(x,t(2:end), xsqavg(2:end)),...
%     [msdDA msdA]);
% Res1 = norm(Anomalous_MSD1D([msdDA msdA],t(2:end), xsqavg(2:end)))^2;
% if ResNorm < Res1
%     msdA = Param(2);
%     msdDA = Param(1);
% end


%Testing 1D Anomalous Diffusion Model
x0 = [msdDA, msdA];
%if you want to fully constrain alpha, you do it as below
%[temp]=lsqnonlin(@NA,x0,[0 0],[50000 1],options);
%this adds significant computation time, and for alpha approximately 1 (i.e
%up to 1.05 or something similar, the solution is valid enough, and we can
%deal with it in model selection). 
%otherwise
[temp]=lsqnonlin(@NA,x0,[],[],options);
param.Dalpha=temp(1);
param.alpha=temp(2);

    function [outvals] = NA(x)
        %Generates the values needed to perform the LSQ non-linear fit for
        %the directed motion model against a given JDD.
        
        %Input Parameters
        Dalpha = x(1);
        alpha = x(2);
        
        %Setting integration limits based on alpha
        if alpha < 0.5
            min=-300^(.5/alpha); %limits on inverse laplace transform
        else
            min=-500;
        end
        
        %predicted JDD probabilities based on input parameters
        fun=@(p) (exp(1i.*p.*tau)).*(1i.*p)^(alpha/2-1)/(2.*pi).*...
            exp(-ri./(sqrt(Dalpha)).*((1i*p)^(alpha/2)));
        predicted = dr/((Dalpha)^(1/2)).*abs(integral(fun,min,-1*min,'ArrayValued',true,'AbsTol',1e-5,'RelTol',1e-3));
        
        %actual JDD probabilities
        actual = yi;
        
        %Weighted Least Squares Function
        outvals = sqrt(weights) .* (predicted - actual);
    end
end
