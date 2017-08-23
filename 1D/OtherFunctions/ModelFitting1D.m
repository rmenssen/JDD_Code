%Parameter fitting for 1D Jump Distance Distributions
%Rebecca Menssen
%This version of the code: 10/24/16

function param = ModelFitting1D(tau, dr, ri, yi, Ni ,N, points, dt, x)
%Inputs:
%tau = the duration of each trajectory (the time lag)
%dr = the width of each bin in the JDD histogram
%ri = vector of the midpoints of the JDD bins
%yi = vector of the percentage of trajectories in each JDD bin
%N = the number of trajectories
%points = the number of time points in each trajectory

%Outputs:
%param: strut containing parameters determined by the LSQ fit for each model.

%Function:
%Takes a JDD and fits dfferent transport models, using
%lsqcurvefit to optimize for the parameters for each model with least
%squares (non-weighted in this case, weighted, while possible for this data,
%often produces worse results than non-weighted)

%In all the @functions for lsqcurvefit:
%x is a vector that holds the parameters to be optimized
%ri is ri (the data used to generate the predicted counts in each bin

%Optimization options:
options = optimoptions('lsqnonlin', 'Algorithm', 'levenberg-marquardt','MaxFunctionEvaluations',1000,'FunctionTolerance',1e-4,'StepTolerance',1e-4);


% Struct to hold all the parameters
param = struct('D',NaN,'V', NaN, 'Dv', NaN,'Dalpha',NaN,'alpha',NaN);

%weighting. Have some options here. 

%add 1 to those that have nothing
%yi(yi==0)=1/N;

wi = 1./yi;

%zero out weights
wi(wi==Inf) = 0;

%add 1 to all(without augmenting total N)
%wi=1./(yi+1/N);

%add 1 to all(with augmentation)
wi=1./((Ni+1)/(N+length(Ni)));
weights=wi;


%SEEDING USING MEAN SQUARE DISPLACEMENT
%time
t=dt*(0:points-1); %xvalues
xsqavg=sum(x.^2,2)./size(x,2);
t=t';
%diffusion
[p]=polyfit(t, xsqavg,1);
msdD=p(1)/2;
msdD=abs(msdD);

%directed motion
[p]=polyfit(t, xsqavg,2);
msdV=sqrt(abs(p(1)));
msdDV=abs(p(2)/2);

%anomalous diffusion
[p]=polyfit(log(t(2:end)), log(xsqavg(2:end)),1);
msdA=p(1);
msdDA=exp(p(2)-log(2)+log(gamma(1+msdA)));
msdA=abs(msdA);
msdDA=abs(msdDA);

%Testing 1D Diffusion Model
x0 = abs(msdD);
[temp]=lsqnonlin(@ND,x0,[],[],options);
param.D=temp(1);

%Testing 1D Directed Motion Model

x0 = [msdV, abs(msdDV)];
[temp] = lsqnonlin(@NV, x0, [], [], options);
param.V = temp(1);
param.Dv = temp(2);

%Testing 1D Anomalous Diffusion Model
x0 = [abs(msdDA), msdA];
%NA(x0)
[temp]=lsqnonlin(@NA,x0,[],[],options);
param.Dalpha=temp(1);
param.alpha=temp(2);

    function [outvals] = ND(x)
        %Generates the values needed to perform the LSQ non-linear fit for
        %the directed motion model against a given JDD.
        D=x(1);
        predicted=dr/((pi*D*tau)^(1/2)).*exp(-ri.^2/(4*D*tau));
        actual = yi;
        outvals = sqrt(weights) .* (predicted - actual);
        %outvals=predicted-actual;
    end

    function [outvals] = NV(x)
        V = x(1);
        D = x(2);
        z = -(ri.^2+V^2*tau^2)/(4*D*tau);
        y = ri*V/(2*D);
        predicted = dr/((4*pi*D*tau)^(1/2)).*exp(z+y);
        actual = yi;
        outvals = sqrt(weights) .* (predicted - actual);
        %outvals=predicted-actual;
    end

    function [outvals] = NA(x)
        %Generates the values needed to perform the LSQ nonlinear fit for
        %the anomalous diffusion model against a given JDD.
        Dalpha = x(1);
        alpha = x(2);
        %setting integration limits
        if alpha < 0.5
            min=-500^(.5/alpha); %limits on inverse laplace transform
            %min=-2000;
        else
            min=-500;
        end
        fun=@(p) (exp(1i.*p.*tau)).*(1i.*p)^(alpha/2-1)/(2.*pi).*exp(-ri./(sqrt(Dalpha)).*((1i*p)^(alpha/2)));
        predicted = dr/((Dalpha)^(1/2)).*abs(integral(fun,min,-1*min,'ArrayValued',true,'AbsTol',1e-5,'RelTol',1e-3));
        actual = yi;
        outvals = sqrt(weights) .* (predicted - actual);
        %outvals=predicted-actual;
    end
end
