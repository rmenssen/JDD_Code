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
%maxtimestep--the length of the overall trajectory (or the length of the
%MSD)
%dt--the time step
%xsqavg--the basic MSD (can be as simple as needed)

%%%%%%%%%%OUTPUTS%%%%%%%%%%
%param: strut containing parameters determined by the LSQ fit for each model.

function param = ModelFitting2DexperimentaldataCombo(tau, dr, ri, yi, Ni, N, maxtimestep,dt,xsqavg)
%Set up optimization options: can be customized for better fits (espeically
%tolerances)
options = optimoptions('lsqnonlin', 'Algorithm', 'levenberg-marquardt',...
    'MaxFunctionEvaluations',1000,'FunctionTolerance',1e-4,...
    'StepTolerance',1e-4);

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

%Getting initial seeds by using MSD Analysis
%make time vector
t=dt*(0:maxtimestep);
t=t';

%raw MSD fits to get basic seeds (you could greatly improve these as
%needed, in general the basic fits do well enough for a decent seed)
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

%testing models
%Testing 1D Diffusion Model
x0=msdD;
[temp]=lsqnonlin(@ND,x0,[],[],options);
param.D=temp(1);

%Testing 1D Directed Motion Model
x0 = [msdV, msdDV];
%suggested limits are given here. To use them, you need to use the trust
%region reflective algorithm and thus change the options. otherwise it
%defaults to levenberg marquart
[temp] = lsqnonlin(@NV, x0, [0,0], [], options);
param.V = temp(1);
param.Dv = temp(2);

%Testing 1D Anomalous Diffusion Model
x0 = [msdDA, msdA];
[temp]=lsqnonlin(@NA,x0,[],[],options);
param.Dalpha=temp(1);
param.alpha=temp(2);

%Double Diffusion
x0 = [msdD, msdD/2,.6];
%suggested limits are given here. To use them, you need to use the trust
%region reflective algorithm and thus change the options. otherwise it
%defaults to levenberg marquart
[temp]=lsqnonlin(@NDD,x0,[0,0,0],[1000,1000,1],options);
param.D1=temp(1);
param.D2=temp(2);
param.fdDD=temp(3);

%Diffusion Directed
x0 = [msdD, msdV, msdDV,.5];
%suggested limits are given here. To use them, you need to use the trust
%region reflective algorithm and thus change the options. otherwise it
%defaults to levenberg marquart
[temp]=lsqnonlin(@NDDV,x0,[0,0,0,0],[1000,1000,1000,1],options);
param.Dddv=temp(1);
param.Vddv=temp(2);
param.DVddv=temp(3);
param.fdDDV=temp(4);

%diffusion anomalous
x0 = [msdD, msdDA, msdA,.5];
%suggested limits are given here. To use them, you need to use the trust
%region reflective algorithm and thus change the options. otherwise it
%defaults to levenberg marquart
[temp]=lsqnonlin(@NDDA,x0,[0,0,0,0],[1000,1000,1000,1],options);
param.Dda=temp(1);
param.DAda=temp(2);
param.Ada=temp(3);
param.fdDDA=temp(4);

%dircted anomalous
x0 = [msdV, msdDV, msdDA, msdA, 0.5];
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
