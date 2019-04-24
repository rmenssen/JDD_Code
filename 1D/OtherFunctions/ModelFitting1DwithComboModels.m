%Parameter fitting for 1D Jump Distance Distributions
%Rebecca Menssen
%Last Updated: 4/9/19

%This function takes JDD data as the input and outputs the parameters for
%three different models using lsqcurvefit with weighting.
%our model fitting method is only valid for alpha <1. In the case that
%alpha >1, you can ignore fits, or deal with fitting in the model fitting
%step. As an alternative, you can impose bound constraints on any fits with
%anomalous diffusion.
%This function also allows for combination models, such as double
%diffusion, increasing computation time, but giving more options for fits. 

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

function param = ModelFitting1DwithComboModels(tau, dr, ri, yi, Ni ,N, points, dt, x1)
%Set up optimization options: can be customized for better fits (espeically
%tolerances)
options = optimoptions('lsqnonlin', 'Algorithm', 'levenberg-marquardt',...
    'MaxFunctionEvaluations',1000,'FunctionTolerance',1e-4,...
    'StepTolerance',1e-4,'Display','off');

%Create struct to hold all the parameters
param = struct('D',NaN,'V', NaN, 'Dv', NaN,'Dalpha',NaN,'alpha',NaN,...
    'D1', NaN,'D2',NaN,'fdDD',NaN,...
    'Ddv', NaN,'Vdv', NaN,'Dvdv', NaN,'fdDV', NaN,...
    'Dda', NaN,'Dalphada', NaN,'alphada', NaN,'fdDA', NaN);

%Weighting vector
%We chose to weight by 1/bin count probabilties (Ni/N=bin count
%probabilities). Each bin was given an artificial additional count so no
%bin was empty
wi=1./((Ni+1)/(N+length(Ni)));
weights=wi;
%if you are using a sliding construction for the jdd then you do NOT want
%weights. 
%weights=1; %for a correlated or sliding JDD. 

%Getting initial seeds by using MSD Analysis
%make time vector
t=dt*(0:points-1);

%finding avg(x^2). Have to account for possible missing data. Outliers can
%mess this calculation up, but there isn't too much that can be done in
%that case, since using nanmedian causes other issues.
xsqavg=nanmean((x1-repmat(x1(1,:),size(x1,1),1)).^2,2);

t=t';

%Pure Diffusion Seeds
[p]=polyfitZero(t, xsqavg,1);
msdD=abs(p(1)/2);

%Directed Diffusion Seeding
options2=optimset('Display','off');
[p]=polyfitZero(t,xsqavg,2);
if p(1)< 0
    paramF = fmincon(@(x) norm(x(1)^2*t.^2 + 2*x(2)*t - xsqavg),[0.1 0.1],...
        [],[],[],[],[0 0],[],[],options2);
    msdV = paramF(1);
    msdDV = paramF(2);
else
    msdV = sqrt(abs(p(1)));
    msdDV = abs(p(2)/2);
end

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
%we deal with any anomalous problems (i.e. alpha>1 in model fitting)
x0 = [msdDA, msdA];
%if you want to fully constrain alpha, you do it as below
%[temp]=lsqnonlin(@NA,x0,[0 0],[500 1],options);
%this will change the method used for curve fitting to the trust-region
%reflective method, and can add computation time. We prefer to just exclude
%models that fit over alpha=1. 
[temp]=lsqnonlin(@NA,x0,[],[],options);
param.Dalpha=temp(1);
param.alpha=temp(2);

%Testing 1D Double Diffusion Model
x0=[0.5*msdD,1.5*msdD,.5]; %so they start out different
[temp]=lsqnonlin(@NDD,x0,[],[],options);
param.D1=temp(1);
param.D2=temp(2);
param.fdDD=temp(3);

%Testing 1D Diffusion Directed Combo Model
x0=[msdD,msdV, msdDV,.5]; %so they start out different
[temp]=lsqnonlin(@NDV,x0,[],[],options);
param.Ddv=temp(1);
param.Vdv=temp(2);
param.Dvdv=temp(3);
param.fdDV=temp(4);

%Testing 1D Diffusion-Anomalous Combo Model
x0=[msdD,msdDA,msdA,.5]; %so they start out different
[temp]=lsqnonlin(@NDA,x0,[],[],options);
param.Dda=temp(1);
param.Dalphada=temp(2);
param.alphada=temp(3);
param.fdDA=temp(4);

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
        predicted1=dr/((pi*D1*tau)^(1/2)).*exp(-ri.^2/(4*D1*tau));
        predicted2=dr/((pi*D2*tau)^(1/2)).*exp(-ri.^2/(4*D2*tau));
        predicted=fd*predicted1+(1-fd)*predicted2;
        %actual JDD probabilities
        actual = yi;
        
        %Weighted Least Squares Function
        outvals = sqrt(weights) .* (predicted - actual);
    end

    function [outvals] = NDV(x)
        %Generates the values needed to perform the LSQ non-linear fit for
        %the directed motion model against a given JDD.
        
        %Input Parameters
        D=x(1);
        V=x(2);
        Dv=x(3);
        fd=x(4);
        
        %predicted JDD probabilities based on input parameters
        z = -(ri.^2+V^2*tau^2)/(4*Dv*tau);
        y = ri*V/(2*Dv);
        predicted=fd*dr/((pi*D*tau)^(1/2)).*exp(-ri.^2/(4*D*tau))+...
            (1-fd)*(dr/((4*pi*Dv*tau)^(1/2)).*exp(z+y)+dr/((4*pi*Dv*tau)^(1/2)).*exp(z-y));
        
        %actual JDD probabilities
        actual = yi;
        
        %Weighted Least Squares Function
        outvals = sqrt(weights) .* (predicted - actual);
    end

    function [outvals] = NDA(x)
        %Generates the values needed to perform the LSQ non-linear fit for
        %the directed motion model against a given JDD.
        
        %Input Parameters
        D=x(1);
        Dalpha = x(2);
        alpha = x(3);
        fd=x(4);
        
        fun=@(p) (exp(1i.*p.*tau)).*(1i.*p)^(alpha/2-1)/(2.*pi).*exp(-ri./(sqrt(Dalpha)).*((1i*p)^(alpha/2)));
        if alpha < 0.5
            min=-300^(.5/alpha);
        else
            min=-500;
        end
        %Plot the predicted JDD on top of it
        predicted=fd*dr/((pi*D*tau)^(1/2)).*exp(-ri.^2/(4*D*tau))+...
            (1-fd)*(dr/((Dalpha)^(1/2)).*abs(integral(fun,min,-1*min,'ArrayValued',true,'AbsTol',2e-5)));
        
        %actual JDD probabilities
        actual = yi;
        
        %Weighted Least Squares Function
        outvals = sqrt(weights) .* (predicted - actual);
    end


end
