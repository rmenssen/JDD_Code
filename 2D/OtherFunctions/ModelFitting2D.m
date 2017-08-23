function param = ModelFitting(tau, dr, ri, yi, Ni, N, points)
%Inputs:
%tau = the duration of each trajectory (the time lag)
%dr = the width of each bin in the JDD histogram
%ri = vector of the midpoints of the JDD bins
%yi = vector of the percentage of trajectories in each JDD bin
%N = the number of trajectories
%points = the number of time points in each trajectory

%Outputs:
%param1,2,3,4,5,6 = vectors containing the set of optimized parameters
%determined by the LSQ fit for each of the 6 models.

%Function:
%Takes a JDD and fits each of the six models (D,V,A,DD,DV,DA), using
%lsqcurvefit to optimize for the parameters for each model with a weighted
%least-squares algorithm.
%The weight of each bin is equal to 1/yi; since
%lsqcurvefit does not natively support weighted least-squares, the data are
%weighted before being passed.

%In all the @functions for lsqcurvefit:
%x is a vector that holds the parameters to be optimized
%ri is ri (the data used to generate the predicted counts in each bin


M=points-1;
wi = 1./yi;
%wi = wi.^2;
wi(wi==Inf) = 0; %CURRENTLY all bins with no entries are given weight of 0 (NOT A CORRECT THING TO DO BECAUSE THEY SHOULD HAVE HIGHESET WEIGHT)

%Optimization options:
options = optimoptions('lsqnonlin', 'Algorithm', 'levenberg-marquardt');

%Find the peak of the JDD - rough, because max may not be peak, but close
%enough for a seed value for optimization. Probably.
maxmid = maxmidp(Ni, ri);

% Struct to hold all the parameters
param = struct('D', NaN, 'V', NaN, 'kv', NaN,'Dalpha',NaN,...
    'alpha',NaN,'fdDD', NaN, 'D1', NaN, 'D2', NaN,...
    'fdDV', NaN, 'Ddv', NaN, 'Vdv', NaN, 'kvdv', NaN);

% Testing Diffusion Model, 1 Parameter: D
% x = [D]
x0 = (maxmid^2)/tau;
[temp] = lsqnonlin(@ND, x0, [], [], options);
%[temp, resnorm, residual, exitflag, output]
param.D = temp(1);

% Testing Directed Motion Model, 2 Parameter: V,kv
% x = [V,kv]
x0 = [maxmid/tau, 2*param.D*tau/M];
[temp] = lsqnonlin(@NV, x0, [], [], options);
param.V = temp(1);
param.kv = temp(2);

% Testing Anomalous Diffusion Model
% x = [Dalpha,alpha]
% x0 = [(maxmid^2)/tau^.5, .5];
% [temp] = lsqnonlin(@NA, x0, [], [], options);
% param.Dalpha = temp(1);
% param.alpha = temp(2);

% Testing DD, 3 Parameters: D1,D2,fd
% x = [D1 D2 fd]
x0 = [(maxmid^2)/tau, 2*(maxmid^2)/tau, 0.5];
[temp] = lsqnonlin(@NDD, x0, [], [], options);
param.D1 = temp(1);
param.D2 = temp(2);
param.fdDD = temp(3);

% Testing DV, 4 Parameters: D, V, kv, fd
% x = [D V kv fd]
x0 = [(maxmid^2)/tau, maxmid/tau, 2*param.D*tau/M,0.5];
[temp] = lsqnonlin(@NDV, x0, [], [], options);
param.Ddv=temp(1);
param.Vdv=temp(2);
param.kvdv=temp(3);
param.fdDV=temp(4);

%}
%%%FUNCTION DEFINTION
    function outvals = ND(x)
        %{
        Generates the values needed to perform the LSQ nonlinear fit for
        the pure diffusion model against a given JDD.
        x = [D]
        %}
        D = x(1);
        predicted = dr*ri/(2*D*tau).*exp(-ri.^2/(4*D*tau));
        weights = wi;
        %weights=1./predicted;
        %weights(weights==Inf) = 0;
        actual = yi;
        outvals = sqrt(weights) .* (predicted - actual);
        %outvals=predicted-actual;
    end

    function [outvals] = NV(x)
        %{
        Generates the values needed to perform the LSQ nonlinear fit for
        the directed motion model against a given JDD.
        x = [V, kv]
        %}
        V = x(1);
        kv = x(2);
        z = -(ri.^2+V^2*tau^2)/(2*M*kv);
        y = ri*V*tau/(M*kv);
        predicted = dr*ri/(kv*M).*exp(z).*besseli(0, y);
        %predicted = dr*ri/(kv*M).*exp(z).*exp(y)./sqrt(2*pi*y);
        %weights = wi;
        %weights=1./predicted;
        %weights(weights==Inf) = 0;
        actual = yi;
        %outvals = sqrt(weights) .* (predicted - actual);
        outvals=predicted-actual;
    end

    function outvals=NA(x)
        Dalpha = x(1);
        alpha = x(2);
        if alpha < 0.5
            min=-300^(.5/alpha); %limits on inverse laplace transform
        else
            min=-300;
        end
        fun=@(p) (exp(1i.*p.*tau)).*(1i.*p)^(alpha-1)/(2.*pi).*(besselk(0,ri./(sqrt(Dalpha)).*((1i*p)^(alpha/2))));
        predicted=dr*ri/(Dalpha).*abs(integral(fun,min,-1*min,'ArrayValued',true,'AbsTol',1e-6));
        weights = wi;
        %weights=1./predicted;
        %weights(weights==Inf) = 0;
        actual = yi;
        outvals = sqrt(weights) .* (predicted - actual);
    end

    function outvals = NDD(x)
        %{
        Generates the values needed to perform the LSQ nonlinear fit for
        the pure double diffusion model against a given JDD.
        x = [D1, D2, fd]
        %}
        D1 = x(1);
        D2 = x(2);
        fd = x(3);
        predicted1 = dr*ri/(2*D1*tau).*exp(-ri.^2/(4*D1*tau));
        predicted2 = dr*ri/(2*D2*tau).*exp(-ri.^2/(4*D2*tau));
        composite  = fd*predicted1 + (1-fd)*predicted2;
        weights = wi;
        %weights=1./composite;
        %weights(weights==Inf) = 0;
        actual = yi;
        outvals = sqrt(weights) .* (composite - actual);
        %outvals=composite-actual;
    end

    function outvals=NDV(x)
        %{
        Generates the values needed to perform the LSQ nonlinear fit for
        the Diffusion/Directed Motion model against a given JDD.
        x = [D,V,kv, fd]
        %}
        D=x(1);
        V=x(2);
        kv=x(3);
        fd=x(4);
        predicted1 = dr*ri/(2*D*tau).*exp(-ri.^2/(4*D*tau));
        z = -(ri.^2+V^2*tau^2)/(2*M*kv);
        y = ri*V*tau/(M*kv);
        predicted2 = dr*ri/(kv*M).*exp(z).*besseli(0, y);
        composite  = fd*predicted1 + (1-fd)*predicted2;
        %weights=1./composite;
        %weights(weights==Inf) = 0;
        %weights = wi;
        actual = yi;
        %outvals = sqrt(weights) .* (composite - actual);
        outvals=composite-actual;
    end
end
