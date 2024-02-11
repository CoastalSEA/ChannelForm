function [dw,yw,Aw] = ckfa_wave_profile(inp,Uw,Fch,rhoc,df,ds)
%
%-------header-------------------------------------------------------------
% NAME
%   ckfa_wave_profile.m [NB modified from version in CKFA model}
% PURPOSE
%   Depth and width of wave formed profile d=dw*(1-y/yw)^2/3
% USAGE
%   [dw yw] = ckfa_wave_profile(inp,Uw,Fch,rhoc,df,ds)
% INPUTS
%   inp is a struct with fields 
%       rhow = density of water (kg/m^3)
%       zw  = elevation of wind speed (m) - default is 10m
%       taucr= critical threshold bed shear stress (Pa)
%       d50  = median sediment grain size diameter (m)
%       ws   = sediment fall velocity (m/s)
%       me   = erosion rate coeficient (kg/N/s)
%       g    = acceleration due to gravity (m/s2)
%   Uw  - wind speed (m/s)
%   Fch - fetch length (m)
%   rhoc - suspended sediment concentration (kg/m^3)
%   df - average depth oveer fetch (m)
%   ds - depth at site (m)
% OUTPUTS
%   dw = depth at outer edge of wave profile
%   yw = width of wave profile
%   Aw - cross-sectional area of profile
% NOTES
%   Equates rate of erosion due to waves with rate of deposition to find
%   depth at which there is a balance
%   Calls internal function fun_dw, which calculates the sedimentation
%   balance for a given depth and wave conditions
%   v2 - now same as wave_form_solver.m
% SEE ALSO
%   ckfa_form_model.m and ckfa_form_properties.m
%
% Author: Ian Townend
% CoastalSEA (c) Jan 2022
%--------------------------------------------------------------------------
%

    %wind speed is zero
    dw = 0; yw = 0; Aw = 0;
    if Uw==0 || ds<0.01, return; end
    % wave parameters
    [Hs, Tp, ~] = tma_spectrum(Uw,inp.zw,Fch,df,ds);
    if Hs==0 || ds<0.01, return; end
    Hrms = Hs/sqrt(2);

    fd = @(dwe) fun_dw(inp,Hrms,Tp,rhoc,dwe);
    options = optimset('TolX',1e-3);
    dw = fminsearch(fd,eps,options);         %depth of profile
    % width of wave formed profile
    [Cd,~,~] = frictionCoefficient(inp,Hrms,Tp,dw);
    yw = 3*pi*dw^2/4/Cd/Hrms;                %width of profile 
    Aw = 3/5*dw*yw;                          %area of profile
end
%%
function fd = fun_dw(inp,Hrms,Tp,rhoc,dwe)
    % Function used to find the erosion depth of wave profile based on a 
    % balance of eorsion and deposition with Uw determined from the hydraulics 
    [Cd,~,~] = frictionCoefficient(inp,Hrms,Tp,dwe);
    %orbital wave amplitude
    Uw = Hrms/2*sqrt(inp.g/dwe);
    %bed shear stress
    tau = inp.rhow*Cd*Uw^2;
    if  Uw<=0 || dwe<0
        ero = NaN;
    elseif tau<inp.taucr
        ero=0;
    else
        ero = (tau-2*inp.taucr)*(pi/2-asin(sqrt(inp.taucr/tau)));
        ero = inp.me/pi*(ero+sqrt(inp.taucr*(tau-inp.taucr)));
    end
    %
    sed = rhoc*inp.ws;  %NB -mass concentration to be consistent with erosion
    %
    fd = abs(ero-sed);
end

%%
function [Cd,fs,fr] = frictionCoefficient(inp,Hrms,Tp,dw)
    %max of smooth and rough wave friction factor
    g = inp.g;
    ks = 2.5*inp.d50;
    fs  = 0.095*(Hrms^2*g*Tp/dw/inp.visc)^-0.187;
    fr  = 0.237*(Hrms*sqrt(g*dw)*Tp/4/pi/dw/ks)^-0.52;
    Cd  = 0.5*max(fs,fr);
end