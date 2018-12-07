%Parameter fitting for 1D Jump Distance Distributions
%Rebecca Menssen
%Last Updated: 8/30/17-cleaning up comments

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
%x1--the 1D trajectory

%%%%%%%%%%OUTPUTS%%%%%%%%%%
%param: strut containing parameters determined by the LSQ fit for each model.

function param = ModelFitting1DwithCombinedModels(tau, dr, ri, yi, Ni ,N, points, dt, x1)
%Set up optimization options:
options = optimoptions('lsqnonlin', 'Algorithm', 'levenberg-marquardt',...
    'MaxFunctionEvaluations',1000,'FunctionTolerance',1e-4,...
    'StepTolerance',1e-4);

%Create struct to hold all the parameters
param = struct('D',NaN,'V', NaN, 'Dv', NaN,'Dalpha',NaN,'alpha',NaN,'D1',NaN,'D2',NaN,'fd',NaN);

%Weighting vector
%We chose to weight by 1/bin count probabilties (Ni/N=bin count
%probabilities). Each bin was given an artificial additional count so no
%bin was empty
wi=1./((Ni+1)/(N+length(Ni)));
weights=wi;

%Getting initial seeds by using MSD Analysis
%make time vector
t=dt*(0:points-1);

%finding avg(x^2). Have to account for possible missing data. Outliers can
%mess this calculation up, but there isn't too much that can be done in
%that case, since using nanmedian causes other issues. 
xsqavg=nanmean((x1-repmat(x1(1,:),size(x1,1),1)).^2,2);

t=t';

%Pure Diffusion Seeds
[p]=polyfit(t, xsqavg,1);
msdD=abs(p(1)/2);

%Directed Diffusion Seeding
[p]=polyfit(t, xsqavg,2);
msdV=sqrt(abs(p(1)));
msdDV=abs(p(2)/2);

%Anomalous Diffusion Seeding
[p]=polyfit(log(t(2:end)), log(xsqavg(2:end)),1);
msdA=abs(p(1));
msdDA=abs(exp(p(2)-log(2)+log(gamma(1+msdA))));

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

%Testing Double Diffusion Model
x0 = [msdD, msdD,.5];
[temp]=lsqnonlin(@NDD,x0,[],[],options);
param.D1=temp(1);
param.D2=temp(2);
param.fd=temp(3);

    function [outvals] = ND(x)
        %Generates the values needed to perform the LSQ non-linear fit for
        %the directed motion model against a given JDD.
        
        %Input Parameter
        D=x(1);
        
        %predicted JDD probabilities based on input parameters
        predicted=dr/((pi*D*tau)^(1/2)).*exp(-ri.^2/(4*D*tau));
        
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


    function [outvals] = NDD(x)
        %Generates the values needed to perform the LSQ non-linear fit for
        %the directed motion model against a given JDD.
        
        %Input Parameter
        D1=x(1);
        D2=x(2);
        fd=x(3);
        
        %predicted JDD probabilities based on input parameters
        predicted=fd*dr/((pi*D1*tau)^(1/2)).*exp(-ri.^2/(4*D1*tau))+...
            (1-fd).*dr/((pi*D2*tau)^(1/2)).*exp(-ri.^2/(4*D2*tau));
        
        %actual JDD probabilities
        actual = yi;
        
        %Weighted Least Squares Function
        outvals = sqrt(weights) .* (predicted - actual);
    end

%     function [outvals] = NV(x)
%         %Generates the values needed to perform the LSQ non-linear fit for
%         %the directed motion model against a given JDD.
%         
%         %Input Parameters
%         V = x(1);
%         D = x(2);
%         
%         %predicted JDD probabilities based on input parameters
%         z = -(ri.^2+V^2*tau^2)/(4*D*tau);
%         y = ri*V/(2*D);
%         predicted = dr/((4*pi*D*tau)^(1/2)).*exp(z+y);
%         
%         %actual JDD probabilities
%         actual = yi;
%         
%         %Weighted Least Squares Function
%         outvals = sqrt(weights) .* (predicted - actual);
%     end
% 
% 
%     function [outvals] = NA(x)
%         %Generates the values needed to perform the LSQ non-linear fit for
%         %the directed motion model against a given JDD.
%         
%         %Input Parameters
%         Dalpha = x(1);
%         alpha = x(2);
%         
%         %Setting integration limits based on alpha
%         if alpha < 0.5
%             min=-300^(.5/alpha); %limits on inverse laplace transform
%         else
%             min=-500;
%         end
%         
%         %predicted JDD probabilities based on input parameters
%         fun=@(p) (exp(1i.*p.*tau)).*(1i.*p)^(alpha/2-1)/(2.*pi).*...
%             exp(-ri./(sqrt(Dalpha)).*((1i*p)^(alpha/2)));
%         predicted = dr/((Dalpha)^(1/2)).*abs(integral(fun,min,-1*min,'ArrayValued',true,'AbsTol',1e-5,'RelTol',1e-3));
%         
%         %actual JDD probabilities
%         actual = yi;
%         
%         %Weighted Least Squares Function
%         outvals = sqrt(weights) .* (predicted - actual);
%     end

end
