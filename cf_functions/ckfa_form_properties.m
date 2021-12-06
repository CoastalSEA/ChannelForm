function formdims = ckfa_form_properties(hme,params)
%
%-------function help------------------------------------------------------
% NAME
%   ckfa_form_properties.m 
% PURPOSE
%   find the hydraulic depth, e-folding length for width and the peak tidal
%   amplitude for defined tidal range and sediment properties
% USAGE
%   formdims = ckfa_form_properties(hme,params)
% INPUTS
%   hme = intial guess of hydraulic depth at mouth (m)
%   params is a struct with fields:
%       a    = tidal amplitude (m)
%       Tp   = tidal period (s)
%       rhow = density of water (kg/m^3)
%     	rhoc = suspended sediment concentration (kg/m^3)
%       taucr= critical threshold bed shear stress (Pa)
%       d50  = median sediment grain size diameter (m)
%       ws   = sediment fall velocity (m/s)
%       me   = erosion rate coeficient (kg/N/s)
% OUTPUTS
%   formdims is a struct containing:
%       La = e-folding length for area (m)
%       Lw = e-folding length for width (m)
%       Am = MTL Area at mouth (m2)
%       Wm = MTL Width at mouth (m)
%       Ucr = peak tidal amplitude (m/s)
% NOTES
%   replaces earlier version that used fzero to find root and called
%   channel_form.m
% SEE ALSO
%   ckfa_wave_profile.m and ckfa_form_model.m
%
% Author: Ian Townend
% CoastalSEA (c)June 2018 updated 2020. converted for ChannelForm App 2022
%--------------------------------------------------------------------------
%
    formdims = [];
    params.omega= 2*pi()/params.Tp;  %angular frequency (1/s)
    fhL = @(hme) fun_hm(hme,params);
    %unconstrained optimisation
    options = optimset('MaxIter',500,'TolFun',1e-9,'TolX',1e-3);
    [hm,~,ok] = fminsearch(fhL,hme,options);    %find minimum of fhL
    if ok<1, return; end       
    Lw = convergenceLength(hm,params.a,params.omega);
    Wm = params.Wr*exp(params.Le/Lw);
    Am = Wm*hm;
    La = -params.Le/log(params.Ar/Am);
    Ucr= params.a*params.omega*La*Wm/Am;
    formdims = struct('La',La,'Lw',Lw,'Am',Am,'Wm',Wm,'Ucr',Ucr);
end
%%
function La = convergenceLength(hm,a,omega)
    %approximate convergence length for given depth  
    g  = 9.81;              %acceleration due to gravity (m/s^2)
    k  = omega/sqrt(g*hm);  %wave number
    eH = pi()*a/4/hm;
    La = 2/k*atan(eH);      %width and area covergence length
    La(La<0) = 0;
end
%%
function fy = fun_hm(hme,params)
% Function used to find the hydraulic depth based on a balance of eorsion
% and deposition with Uc determined from the hydraulics nd Lw from equating
% hydraulic and geometric prism estimates
    p = params;
    [~,Cd] = ucrit(hme,p.d50,p.rhow,p.taucr);    
    %
    La = convergenceLength(hme,p.a,p.omega);
    %
    Uc  = p.a*p.omega*La./hme;
    idu = Uc<0;
    Uc(idu) = 0;
    %
    tau = p.rhow*Cd.*Uc.^2;
    ero = zeros(size(hme));
    ide = tau>p.taucr;
    ero(ide) = (tau(ide)-2*p.taucr).*(pi/2-asin(sqrt(p.taucr./tau(ide))));
    ero(ide) = p.me/pi.*(ero(ide)+ sqrt(p.taucr*(tau(ide)-p.taucr)));               
    %
    sed = p.rhoc*p.ws;
    %
    fy = abs(ero-sed);
end