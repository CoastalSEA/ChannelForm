function formdims = ckfa_form_solver(hme,params)
%
%-------function help------------------------------------------------------
% NAME
%   ckfa_form_solver.m 
% PURPOSE
%   find the hydraulic depth, e-folding length for width and the peak tidal
%   amplitude for defined tidal range and sediment properties
% USAGE
%   formdims = ckfa_form_solver(hme,params)
% INPUTS
%   hme = intial guess of hydraulic depth at mouth (m)
%   params is a struct with fields (*)=not used in this function:
%       am    = tidal amplitude (m)
%       tp   = tidal period (s)
%       Le   = estuary length (m)
%       Qr   = river discharge (m3/s)...(*)
%       hrv  = hydraulic depth of regime river channel (m)...(*)
%       Wrv  = width of regime river channel (m)
%       Arv  = CSA of regime river channel (m2)...(*)
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
%       hm = MTL hydraulic depth at mouth (m)
%       Wm = MTL Width at mouth (m)
%       Am = MTL CSA at mouth (m2)
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
    params.omega= 2*pi()/params.tp;  %angular frequency (1/s)
    fhL = @(hme) fun_hm(hme,params);
    %unconstrained optimisation
    options = optimset('MaxIter',500,'TolFun',1e-9,'TolX',1e-3);
    [hm,~,ok] = fminsearch(fhL,hme,options);    %find minimum of fhL
    if ok<1, return; end       
    Lw = convergenceLength(hm,params.am,params.omega);
    Wm = params.Wrv*exp(params.Le/Lw);
    Am = Wm*hm;
    La = -params.Le/log(params.Arv/Am);
    Ucr= params.am*params.omega*La*Wm/Am;
    formdims = table(La,Lw,hm,Wm,Am,Ucr);
end
%%
function La = convergenceLength(hm,am,omega)
    %approximate convergence length for given depth  
    g  = 9.81;              %acceleration due to gravity (m/s^2)
    k  = omega/sqrt(g*hm);  %wave number
    eH = pi()*am/4/hm;
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
    La = convergenceLength(hme,p.am,p.omega);
    %
    Uc  = p.am*p.omega*La./hme;
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