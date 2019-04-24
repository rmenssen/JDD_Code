%Parameter fitting for 1D Jump Distance Distributions
%Rebecca Menssen
%Last Updated: 4/9/19

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
%x1,x2,x3--the 1D trajectory for x, y, and z

%%%%%%%%%%OUTPUTS%%%%%%%%%%
%param: strut containing parameters determined by the LSQ fit for each model.

function param = ModelFitting3D(tau, dr, ri, yi, Ni, N, points,dt,x1,x2,x3)
%Set up optimization options: can be customized for better fits (espeically
%tolerances)
options = optimoptions('lsqnonlin', 'Algorithm', 'levenberg-marquardt',...
    'MaxFunctionEvaluations',1000,'FunctionTolerance',1e-4,...
    'StepTolerance',1e-4,'Display','off');

%Create struct to hold all the parameters
param = struct('D',NaN,'V', NaN, 'Dv', NaN,'Dalpha',NaN,'alpha',NaN);

%Weighting vector
%We chose to weight by 1/bin count probabilties (Ni/N=bin count
%probabilities). Each bin was given an artificial additional count so no
%bin was empty
wi=1./((Ni+1)/(N+length(Ni)));
weights=wi;
%if using the sliding method to construct the JDD, we do not recommend
%using weights, so set weights to 1. 
%weights=1;

t=dt*(0:points-1); %xvalues
%need to change this once we have non-zero starting values. Check this
%later. %changed to nanmedian to get rid of the effect of crazy outliers
xsqavg=nanmean((x1-repmat(x1(1,:),size(x1,1),1)).^2+...
    (x2-repmat(x2(1,:),size(x2,1),1)).^2+(x3-repmat(x3(1,:),size(x3,1),1)).^2,2);
t=t';

%Pure Difusion Seeding
[p]=polyfitZero(t, xsqavg,1);
msdD=abs(p(1)/6);

%Directed Diffusion Seeding
options2=optimset('Display','off');
[p]=polyfitZero(t, xsqavg,2);
if p(1)< 0
    paramF = fmincon(@(x) norm(x(1)^2*t.^2 + 6*x(2)*t - xsqavg),[0.1 0.1],...
        [],[],[],[],[0 0],[],[],options2);
    msdV = paramF(1);
    msdDV = paramF(2);
else
    msdV = sqrt(abs(p(1)));
    msdDV = abs(p(2)/6);
end

%Anomalous Diffusion Seeding
[p]=polyfit(log(t(2:end)), log(xsqavg(2:end)),1);
msdA=abs(p(1));
msdDA=abs(exp(p(2)-log(6)+log(gamma(1+msdA))));

%Testing 3D Pure Diffusion Model
x0 = msdD;
[temp]=lsqnonlin(@ND,x0,[],[],options);
param.D=temp(1);

%Testing 1D Directed Motion Model
x0 = [msdV, msdDV];
[temp] = lsqnonlin(@NV, x0, [0,0], [], options);
param.V = temp(1);
param.Dv = temp(2);

%Testing 1D Anomalous Diffusion Model
%we do not constrain alpha here even though alpha should be less than one.
%Instead, we treat it in the model selection step or exclude the model
%entirely if that is the case. To add in constraints, you can switch the
%method to the trust-region-reflective method by editing the bound
%constraints below. 
x0 = [msdDA, msdA];
[temp]=lsqnonlin(@NA,x0,[],[],options);
param.Dalpha=temp(1);
param.alpha=temp(2);

    function [outvals] = ND(x)
        %Generates the values needed to perform the LSQ non-linear fit for
        %the directed motion model against a given JDD.
        
        %Input Parameter
        D=x(1);
        
        %predicted JDD probabilities based on input parameters
        predicted=dr*ri.^2/(2*sqrt(pi)*(D*tau)^(3/2)).*exp(-ri.^2/(4*D*tau));
        
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
        predicted=dr*ri.^2*4*pi/((4*pi*D*tau)^(3/2)).*exp(z).*sinh(y)./y;
        
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
        if alpha < 0.5 && alpha>0.4
            min=-700; %limits on inverse laplace transform
        elseif alpha<= 0.4 && alpha>0.3
            min=-900;
        elseif alpha<=0.3
            min=-1100;
        else
            min=-700;
        end
        
        %predicted JDD probabilities based on input parameters
        %For 3D anomalous diffusion, we found that splitting up the
        %integration into chunks led to the most accurate result. 
        fun=@(p) (exp(1i.*p.*tau)).*(1i.*p)^(alpha-1)/(2.*pi).*...
            (exp(-ri./(sqrt(Dalpha)).*((1i*p)^(alpha/2))));
        intsize=200;
        numint=ceil(abs(min+100)/intsize);
        predicted=zeros(size(ri));
        for i=1:numint
            predicted=predicted+2*dr*ri/(Dalpha).*abs(integral(fun,min+(i-1)...
                *intsize-1i*1e-6,min+i*(intsize)-1i*1e-6,...
                'ArrayValued',true,'AbsTol',1e-6,'RelTol',1e-3));
        end
        predicted=predicted+dr*ri/(Dalpha).*abs(integral(fun,-100-1i*1e-6,...
            100-1i*1e-6,'ArrayValued',true,'AbsTol',1e-6,'RelTol',1e-3));
        
        %actual JDD probabilities
        actual = yi;
        
        %Weighted Least Squares Function
        outvals = sqrt(weights) .* (predicted - actual);
    end
end
