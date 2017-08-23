function param = ModelFitting3D(tau, dr, ri, yi, Ni, N, points,dt,x1,x2,x3)
%Inputs:
%tau = the duration of each trajectory (the time lag)
%dr = the width of each bin in the JDD histogram
%ri = vector of the midpoints of the JDD bins
%yi = vector of the percentage of trajectories in each JDD bin
%N = the number of trajectories
%points = the number of time points in each trajectory

%Optimization options:
options = optimoptions('lsqnonlin', 'Algorithm', 'levenberg-marquardt','MaxFunctionEvaluations',1000,'FunctionTolerance',1e-4,'StepTolerance',1e-4);

%wi = 1./yi;
%wi = wi.^2;
%wi(wi==Inf) = 0;

wi=1./((Ni+1)/(N+length(Ni)));
weights=wi;

t=dt*(0:points-1); %xvalues
%need to change this once we have non-zero starting values. Check this
%later. %changed to nanmedian to get rid of the effect of crazy outliers
xsqavg=nanmedian((x1-repmat(x1(1,:),size(x1,1),1)).^2+...
    (x2-repmat(x2(1,:),size(x2,1),1)).^2+(x3-repmat(x3(1,:),size(x3,1),1)).^2,2);
t=t';

%diffusion
[p]=polyfit(t, xsqavg,1);
msdD=p(1)/6;

%directed motion
[p]=polyfit(t, xsqavg,2);
msdV=sqrt(abs(p(1)));
msdDV=p(2)/6;

%anomalous diffusion
[p]=polyfit(log(t(2:end)), log(xsqavg(2:end)),1);
msdA=p(1);
msdDA=exp(p(2)-log(6)+log(gamma(1+msdA)));

% Struct to hold all the parameters
param = struct('D',NaN,'V', NaN, 'Dv', NaN,'Dalpha',NaN,'alpha',NaN);


%Testing 3D Diffusion Model
x0 = msdD;
[temp]=lsqnonlin(@ND,x0,[],[],options);
param.D=temp(1);

%Testing 3D Directed Motion Model

x0 = [msdV, msdDV];
[temp] = lsqnonlin(@NV, x0, [0 0], [], options);
param.V = temp(1);
param.Dv = temp(2);

%Testing 3D Anomalous Diffusion Model
x0 = [msdDA, msdA];
[temp]=lsqnonlin(@NA,x0,[],[],options);
param.Dalpha=temp(1);
param.alpha=temp(2);
%param.Dalpha=0.0087;
%param.alpha=0.5;

    function [outvals] = ND(x)
        %Generates the values needed to perform the LSQ non-linear fit for
        %the directed motion model against a given JDD.
        D=x(1);
        predicted=dr*ri.^2/(2*sqrt(pi)*(D*tau)^(3/2)).*exp(-ri.^2/(4*D*tau));
        actual = yi;
        outvals = sqrt(weights) .* (predicted - actual);
        %outvals=predicted-actual;
    end

    function [outvals] = NV(x)
        %{
        Generates the values needed to perform the LSQ nonlinear fit for
        the directed motion model against a given JDD.
        x = [V, Dv]
        %}
        V = x(1);
        D = x(2);
        z = -(ri.^2+V^2*tau^2)/(4*D*tau);
        y = ri*V/(2*D);
        predicted=dr*ri.^2*4*pi/((4*pi*D*tau)^(3/2)).*exp(z).*sinh(y)./y;
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
        if alpha < 0.5 && alpha>0.4
            min=-700; %limits on inverse laplace transform
        elseif alpha<= 0.4 && alpha>0.3
            min=-900;
        elseif alpha<=0.3
            min=-1100;
        else
            min=-700;
        end
        fun=@(p) (exp(1i.*p.*tau)).*(1i.*p)^(alpha-1)/(2.*pi).*(exp(-ri./(sqrt(Dalpha)).*((1i*p)^(alpha/2))));
        intsize=200;
        numint=ceil(abs(min+100)/intsize);
        predicted=zeros(size(ri));
        
        for i=1:numint
            predicted=predicted+2*dr*ri/(Dalpha).*abs(integral(fun,min+(i-1)*intsize-1i*1e-6,min+i*(intsize)-1i*1e-6,...
                'ArrayValued',true,'AbsTol',1e-6,'RelTol',1e-3));
        end
        predicted=predicted+dr*ri/(Dalpha).*abs(integral(fun,-100-1i*1e-6,100-1i*1e-6,...
            'ArrayValued',true,'AbsTol',1e-6,'RelTol',1e-3));
        
        %predicted = dr*ri/(Dalpha).*abs(integral(fun,min,-1*min,'ArrayValued',true,'AbsTol',1e-6,'RelTol',1e-3));
        actual = yi;
        outvals = sqrt(weights) .* (predicted - actual);
        %outvals=predicted-actual;
    end
end
