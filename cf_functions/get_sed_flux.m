function [dvol,delV] = get_sed_flux(inp,slr,isplt)
%
%-------function help------------------------------------------------------
% NAME
%   get_sed_flux.m
% PURPOSE
%   single element ASMITA model to compute amount of sedimentation
% USAGE
%   [dvol,delV] = get_sed_flux(inp,slr,isplt)
% INPUTS
%   inp - struct containing the following fields
%         inp.Volume;            %element volume at start of run (m^3)
%         inp.SurfaceArea;       %element surface area (m^2)
%         inp.Prism;             %tidal prism of channel (m^3)
%         inp.VerticalExchange;  %vertical exchange (m/s)
%         inp.HorizontalExchange;%horizontal exchange (m/s)
%         inp.RiverDischarge;    %river discharge (m^3/s) +ve downstream
%         inp.TransportCoeff;    %transport coefficient n (3-5)
%         inp.EqConc;            %equilibrium concentration (-)
%         inp.RiverConc;         %river load imported by advection (-)
%         inp.BedConc;           %concentration of bed (-)
%         inp.TimeInt;           %time increment in sem model (years)
%         inp.y2s;               %factor to convert from years to seconds
%   slr - rate of sea level rise (m/yr)
%   isplt - flag, if true generate a plot
% OUTPUT
%   dvol - change in morphological volume in a year (m3)
%   delV - change in water volume (S*slr) in a year (m3)
% SEE ALSO
%   used in ChannelForm model (eg CF_SediData)
%
% Author: Ian Townend
% CoastalSEA (c) Jan 2022
%--------------------------------------------------------------------------
%
    %
    %input is class object, modelUI object and rate of sea level rise(m/yr)
    %returns 
    if nargin<3
        isplt = false;
    end
    V = inp.Volume;            %element volume at start of run (m^3)
    S = inp.SurfaceArea;       %element surface area (m^2)
    P = inp.Prism;             %tidal prism of channel (m^3)
    w = inp.VerticalExchange;  %vertical exchange (m/s)
    d = inp.HorizontalExchange;%horizontal exchange (m/s)
    q = inp.RiverDischarge;    %river discharge (m^3/s) +ve downstream
    n = inp.TransportCoeff;    %transport coefficient n (3-5)
    %
    cE = inp.EqConc;           %equilibrium concentration (-)
    cr = inp.RiverConc;        %river load imported by advection (-)
    cb = inp.BedConc;          %concentration of bed (-)
    k = cr/cE;
    %
    dt = inp.TimeInt;    %time increment in sem model (years)
    if dt==0
        dt = 1;          %single year (years) was RunDuration
    end
    y2s = inp.y2s;       %factor to convert from years to seconds
    

    C = 1/cb*(w*cE*S)/(d+q+w*S);%constant term in dV/dt
    alpha = ((k*q+d)/(q+d))^(1/n);

    Vf(1) = V; Vm(1) = V;  aa = V/P;
    %             if ~ists
    delV = S*slr;                      %water volume change (m3/yr)
    Ve = aa*(P);
    Gam = (alpha*Ve/(Vm+delV))^n;
    dvol = C*((q+d)*Gam-(k*q+d))*y2s;  %morphological change (m3/yr)
    %             else
%     tlen = inp.RunDuration;         %length of run in years
%     if tlen<100
%         tlen = 100;                 %use a minimum of 100 years
%     end
%     tlen = 1000;
%     dtt = 1.0;                      %sem time step
%     ti = 0:dtt:tlen;                %time (yrs)
%     r = SLR/dt;                     %rate of slr (m/yr)
%     delV = S*r*dtt;
%     Ve(1) = aa*(P);
%     for i = 2:length(ti)
%         Gam = (alpha*Ve(i-1)/Vm(i-1))^n;
%         dV(i) = C*((q+d)*Gam-(k*q+d))*dtt*y2s;
%         Vf(i) = Vf(i-1)+dV(i);
%         Vm(i) = Vm(i-1)+dV(i)+delV;
%         Ve(i) = aa*(P);
%     end
%     %dvol = mean(diff(Vf))*dt/dtt;
%     %dvol = mean(dV)*dt/dtt;
%     dvol = mean(dV(end-100:end))*dt/dtt;
%     dvol = dV(end);
    if isplt
        figure;
        plot(ti,Vf,ti,Ve,ti,Vm);
        legend('Fixed','Equilibirum','Moving')
        %plot(ti,dV);
    end