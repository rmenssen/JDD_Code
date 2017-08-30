%Parameter fitting for 2D Jump Distance Distributions
%Rebecca Menssen
%Last Updated: 8/30/17-cleaning up comments, fixing prior fitting issues

%This function takes JDD data as the input and outputs the parameters for
%three different models using lsqcurvefit with weighting.

%%%%%%%%%%INPUTS%%%%%%%%%%
%tau--the duration of each trajectory (the time lag*dt)
%dr--the width of each bin in the JDD histogram
%ri--vector of the midpoints of the JDD bins
%yi--vector of the percentage of trajectories in each JDD bin
%N--the number of trajectories
%points--the number of time points in each trajectory
%dt--the time step
%x1,x2--the 1D trajectory for x and y

%%%%%%%%%%OUTPUTS%%%%%%%%%%
%param: strut containing parameters determined by the LSQ fit for each model.


function param = ModelFitting2D(tau, dr, ri, yi, Ni, N, points,dt,x1,x2)
%Set up optimization options:
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

%Getting initial seeds by using MSD Analysis
%make time vector
t=dt*(0:points-1);

%Finding xsquare average. Have to account for possible missing data and
%crazy outliers, so using nanmedian to do this. Not exactly MSD then, but
%for a simple way to find a seed, it works well.
xsqavg=nanmedian((x1-repmat(x1(1,:),size(x1,1),1)).^2+...
    (x2-repmat(x2(1,:),size(x2,1),1)).^2);
t=t';

%Pure Diffusion Seeding
[p]=polyfit(t, xsqavg,1);
msdD=abs(p(1)/4);

%Directed Diffusion Seeding
[p]=polyfit(t, xsqavg,2);
msdV=sqrt(abs(p(1)));
msdDV=abs(p(2)/4);

%Anomalous Diffusion Seeding
[p]=polyfit(log(t(2:end)), log(xsqavg(2:end)),1);
msdA=abs(p(1));
msdDA=abs(exp(p(2)-log(4)+log(gamma(1+msdA))));

%Testing 1D Diffusion Model
x0=msdD;
[temp]=lsqnonlin(@ND,x0,[],[],options);
param.D=temp(1);

%Testing 1D Directed Motion Model
x0 = [msdV, msdDV];
[temp] = lsqnonlin(@NV, x0, [], [], options);
param.V = temp(1);
param.Dv = temp(2);

%Testing 1D Anomalous Diffusion Model
x0 = [msdDA, msdA];
[temp]=lsqnonlin(@NA,x0,[],[],options);
param.Dalpha=temp(1);
param.alpha=temp(2);

    function outvals = ND(x)
        %Generates the values needed to perform the LSQ non-linear fit for
        %the directed motion model against a given JDD.
        
        %Input Parameter
        D = x(1);
        
        %predicted JDD probabilities based on input parameters
        predicted = dr*ri/(2*D*tau).*exp(-ri.^2/(4*D*tau));
        
        %actual JDD probabilities
        actual = yi;
        
        %Weighted Least Squares Function
        outvals = sqrt(weights) .* (predicted - actual);
    end

    function [outvals] = NV(x)
        %Generates the values needed to perform the LSQ non-linear fit for
        %the directed motion model against a given JDD.
        
        %Input Parameters
        V = x(1);
        D= x(2);
        
        %predicted JDD probabilities based on input parameters
        z = -(ri.^2+V^2*tau^2)/(2*M*kv);
        y = ri*V/(2*D);
        predicted = dr*ri/(4*D*tau).*exp(z).*besseli(0, y);
        
        %actual JDD probabilities
        actual = yi;
        
        %Weighted Least Squares Function
        outvals = sqrt(weights) .* (predicted - actual);
    end

    function outvals=NA(x)
        %Generates the values needed to perform the LSQ non-linear fit for
        %the directed motion model against a given JDD.
        
        %Input Parameters
        Dalpha = x(1);
        alpha = x(2);
        
        %Setting integration limits based on alpha
        if alpha < 0.5
            min=-300^(.5/alpha); %limits on inverse laplace transform
        else
            min=-300;
        end
        
        %predicted JDD probabilities based on input parameters
        fun=@(p) (exp(1i.*p.*tau)).*(1i.*p)^(alpha-1)/(2.*pi).*...
            (besselk(0,ri./(sqrt(Dalpha)).*((1i*p)^(alpha/2))));
        predicted=dr*ri/(Dalpha).*abs(integral(fun,min,-1*min,...
            'ArrayValued',true,'AbsTol',1e-6));
        
        %actual JDD probabilities
        actual = yi;
        
        %Weighted Least Squares Function
        outvals = sqrt(weights) .* (predicted - actual);
    end
end
