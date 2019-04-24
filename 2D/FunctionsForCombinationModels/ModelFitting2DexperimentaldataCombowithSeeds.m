%Parameter fitting for 2D Jump Distance Distributions
%Rebecca Menssen
%Last Updated: 3/20/19

%This function takes JDD data as the input and outputs the parameters for
%three different models using lsqcurvefit with weighting.

%%%%%%%%%%INPUTS%%%%%%%%%%
%tau--the duration of each trajectory (the time lag*dt)
%dr--the width of each bin in the JDD histogram
%ri--vector of the midpoints of the JDD bins
%yi--vector of the percentage of trajectories in each JDD bin
%N--the number of trajectories
%seeds--Seeding values for each model. These could be based on a prior fit,
    %msd fits, whatever necessary. Order that they should be inputted in is as
    %follows [D, V, DV, Dalpha, alpha, D1, D2, fdDD, D(from DV model), V(from DV model), 
    %DV (from DV model), fdDV, D (from DA model), Dalpha (from DA model), 
    %alpha (from DA model), fdDA, V (from VA model), DV (from VA model), 
    %Dalpha (from VA model), alpha (from VA model), fdDA]

%%%%%%%%%%OUTPUTS%%%%%%%%%%
%param: strut containing parameters determined by the LSQ fit for each model.


function param = ModelFitting2DexperimentaldataCombowithSeeds(tau, dr, ri, yi, Ni, N,seeds)
%Set up optimization options: can be customized for better fits (espeically
%tolerances)
options = optimoptions('lsqnonlin', 'Algorithm', 'levenberg-marquardt',...
    'MaxFunctionEvaluations',1000,'FunctionTolerance',1e-6,...
    'StepTolerance',1e-6);

%Create struct to hold all the parameters
param = struct('D',NaN,'V', NaN, 'Dv', NaN,'Dalpha',NaN,'alpha',NaN, ...
    'D1',NaN,'D2',NaN,'fdDD',NaN,'Dddv',NaN,'Vddv',NaN,'DVddv',NaN, ...
    'fdDDV',NaN,'Dda',NaN,'DAda',NaN,'Ada',NaN,'fdDDA',NaN,...
    'Vdvda',NaN,'DVdvda',NaN,'DAdvda',NaN,'Advda',NaN,'fdDVDA',NaN);

%Weighting vector
%We chose to weight by 1/bin count probabilties (Ni/N=bin count
%probabilities). Each bin was given an artificial additional count so no
%bin was empty
wi=1./((Ni+1)/(N+length(Ni)));
weights=wi;

%Testing 1D Diffusion Model
x0=seeds(1);
[temp]=lsqnonlin(@ND,x0,[],[],options);
param.D=temp(1);

%Testing 1D Directed Motion Model
x0 = [seeds(2), seeds(3)];
%suggested limits are given here. To use them, you need to use the trust
%region reflective algorithm and thus change the options. otherwise it
%defaults to levenberg marquart
[temp] = lsqnonlin(@NV, x0, [0,0], [], options);
param.V = temp(1);
param.Dv = temp(2);

%Testing 1D Anomalous Diffusion Model
x0 = [seeds(4), seeds(5)];
[temp]=lsqnonlin(@NA,x0,[],[],options);
param.Dalpha=temp(1);
param.alpha=temp(2);

%Double Diffusion
x0 = [seeds(6),seeds(7),seeds(8)];
%suggested limits are given here. To use them, you need to use the trust
%region reflective algorithm and thus change the options. otherwise it
%defaults to levenberg marquart
[temp]=lsqnonlin(@NDD,x0,[0,0,0],[1000,1000,1],options);
param.D1=temp(1);
param.D2=temp(2);
param.fdDD=temp(3);

%Diffusion Directed
x0 = [seeds(9),seeds(10),seeds(11),seeds(12)];
%suggested limits are given here. To use them, you need to use the trust
%region reflective algorithm and thus change the options. otherwise it
%defaults to levenberg marquart
[temp]=lsqnonlin(@NDDV,x0,[0,0,0,0],[1000,1000,1000,1],options);
param.Dddv=temp(1);
param.Vddv=temp(2);
param.DVddv=temp(3);
param.fdDDV=temp(4);

%diffusion anomalous
x0 = [seeds(13),seeds(14),seeds(15),seeds(16)];
%suggested limits are given here. To use them, you need to use the trust
%region reflective algorithm and thus change the options. otherwise it
%defaults to levenberg marquart
[temp]=lsqnonlin(@NDDA,x0,[0,0,0,0],[1000,1000,1000,1],options);
param.Dda=temp(1);
param.DAda=temp(2);
param.Ada=temp(3);
param.fdDDA=temp(4);

%dircted anomalous
x0 = [seeds(17),seeds(18),seeds(19),seeds(20),seeds(21)];
%suggested limits are given here. To use them, you need to use the trust
%region reflective algorithm and thus change the options. otherwise it
%defaults to levenberg marquart
[temp]=lsqnonlin(@NDVDA,x0,[0,0,0,0,0],[1000,1000,1000,1000,1],options);
param.Vdvda=temp(1);
param.DVdvda=temp(2);
param.DAdvda=temp(3);
param.Advda=temp(4);
param.fdDVDA=temp(5);

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
        z = -(ri.^2+V^2*tau^2)/(4*D*tau);
        y = ri*V/(2*D);
        predicted = dr*ri/(2*D*tau).*exp(z).*besseli(0, y);
        
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

    function outvals = NDD(x)
        %Generates the values needed to perform the LSQ non-linear fit for
        %the directed motion model against a given JDD.
        
        %Input Parameter
        D1 = x(1);
        D2 = x(2);
        fd = x(3);
        
        %predicted JDD probabilities based on input parameters
        predicted1 = dr*ri/(2*D1*tau).*exp(-ri.^2/(4*D1*tau));
        predicted2 = dr*ri/(2*D2*tau).*exp(-ri.^2/(4*D2*tau));
        
        predicted=fd*predicted1+(1-fd)*predicted2;
        
        %actual JDD probabilities
        actual = yi;
        
        %Weighted Least Squares Function
        outvals = sqrt(weights) .* (predicted - actual);
    end

    function outvals = NDDV(x)
        %Generates the values needed to perform the LSQ non-linear fit for
        %the directed motion model against a given JDD.
        
        %Input Parameter
        Dddv=x(1);
        Vddv=x(2);
        DVddv=x(3);
        fdDDV=x(4);
        
        %predicted JDD probabilities based on input parameters
        predicted1 = dr*ri/(2*Dddv*tau).*exp(-ri.^2/(4*Dddv*tau));
        z = -(ri.^2+Vddv^2*tau^2)/(4*DVddv*tau);
        y = ri*Vddv/(2*DVddv);
        predicted2 = dr*ri/(2*DVddv*tau).*exp(z).*besseli(0, y);
        predicted2(isnan(predicted2))=0;
        predicted=fdDDV*predicted1+(1-fdDDV)*predicted2;
        
        %actual JDD probabilities
        actual = yi;
        
        %Weighted Least Squares Function
        outvals = sqrt(weights) .* (predicted - actual);
    end

    function outvals = NDDA(x)
        %Generates the values needed to perform the LSQ non-linear fit for
        %the directed motion model against a given JDD.
        
        %Input Parameter
        Dda=x(1);
        DAda=x(2);
        Ada=x(3);
        fdDDA=x(4);
        
        if Ada < 0.5
            min=-300^(.5/Ada); %limits on inverse laplace transform
        else
            min=-300;
        end
        
        %predicted JDD probabilities based on input parameters
        predicted1 = dr*ri/(2*Dda*tau).*exp(-ri.^2/(4*Dda*tau));
        fun=@(p) (exp(1i.*p.*tau)).*(1i.*p)^(Ada-1)/(2.*pi).*...
            (besselk(0,ri./(sqrt(DAda)).*((1i*p)^(Ada/2))));
        predicted2=dr*ri/(DAda).*abs(integral(fun,min,-1*min,...
            'ArrayValued',true,'AbsTol',1e-6));
        
        predicted=fdDDA*predicted1+(1-fdDDA)*predicted2;
        
        %actual JDD probabilities
        actual = yi;
        
        %Weighted Least Squares Function
        outvals = sqrt(weights) .* (predicted - actual);
    end

    function outvals = NDVDA(x)
        %Generates the values needed to perform the LSQ non-linear fit for
        %the directed motion model against a given JDD.
        
        %Input Parameter
        Vdvda=x(1);
        DVdvda=x(2);
        DAdvda=x(3);
        Advda=x(4);
        fdDVDA=x(5);
        
        if Advda < 0.5
            min=-300^(.5/Advda); %limits on inverse laplace transform
        else
            min=-300;
        end
        
        %predicted JDD probabilities based on input parameters
        
        z = -(ri.^2+Vdvda^2*tau^2)/(4*DVdvda*tau);
        y = ri*Vdvda/(2*DVdvda);
        predicted1 = dr*ri/(2*DVdvda*tau).*exp(z).*besseli(0, y);
        
        fun=@(p) (exp(1i.*p.*tau)).*(1i.*p)^(Advda-1)/(2.*pi).*...
            (besselk(0,ri./(sqrt(DAdvda)).*((1i*p)^(Advda/2))));
        predicted2=dr*ri/(DAdvda).*abs(integral(fun,min,-1*min,...
            'ArrayValued',true,'AbsTol',1e-6));
        
        predicted=fdDVDA*predicted1+(1-fdDVDA)*predicted2;
        
        %actual JDD probabilities
        actual = yi;
        
        %Weighted Least Squares Function
        outvals = sqrt(weights) .* (predicted - actual);
    end
end
