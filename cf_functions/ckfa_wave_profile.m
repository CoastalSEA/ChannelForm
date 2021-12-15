function [dw,yw] = ckfa_wave_profile(Uw,zw,Fch,cnc,df,ds,d50,tau,me,ws)
%
%-------header-------------------------------------------------------------
% NAME
%   ckfa_wave_profile.m [NB modified from version in CKFA model}
% PURPOSE
%   Depth and width of wave formed profile d=dw*(1-y/yw)^2/3
% USAGE
%   [dw yw] = ckfa_wave_profile(Uw,zw,Fch,cnc,df,ds,d50,tau,me,ws)
% INPUTS
%   Uw  = wind speed (m/s)
%   zw  = elevation of wind speed (m) - default is 10m
%   Fch = fetch length (m)
%   cnc = suspended sediment concentration
%   df  = average depth over fetch (m)
%   ds  = depth at site (m) - taken as depth at edge of profile
%   d50 = median sediment grain size diameter (m)
%   tau = critical threshold bed shear stress (Pa)
%   me  = erosion rate coeficient (kg/N/s)
%   ws  = sediment fall velocity (m/s)
% OUTPUTS
%   dw = depth at outer edge of wave profile
%   yw = width of wave profile
% NOTES
%   Equates rate of erosion due to waves with rate of deposition to find
%   depth at which there is a balance
%   Calls internal function fun_dw, which calculates the sedimentation
%   balance for a given depth and wave conditions
% SEE ALSO
%   ckfa_form_model.m and ckfa_form_properties.m
%
% Author: Ian Townend
% CoastalSEA (c) Jan 2022
%--------------------------------------------------------------------------
%
g  = 9.81;
ks = 2.5*d50;
visc = 1.34*10^-6; %viscosity of water (m2/s)
rhos = 2650;       %density of sediment (kg/m3)
rhow = 1025;       %density of water (kg/m3)
rhoc = cnc*rhos;   %density of suspended sediment (kg/m3)
%
if Uw==0
    dw = 0;
    yw = 0;
    return
end
% Adjust wind speed to 10m using power law profile
% U  = Uw*(10/zw)^(1/7);
% depth at which waves can erode bed (initial <high> estimate)
% kcon = 0.5*10^-3;  %concentration coefficient
%dw1 = 2.89*10^-6*sqrt(kcon/cnc)*g*Fch^0.44*U^2.4;

% wave parameters
[Hs, Tp, ~] = tma_spectrum(Uw,zw,Fch,df,ds);
Hrms = Hs/sqrt(2);

%the convergence algorithm is different to the one used in the CKFA model
dw1 = eps;
fdw1= fun_dw(dw1,Hrms,Tp,rhoc,rhow,tau,me,ws,ks,visc);
if fdw1<0
    dw=0; yw=0;
    sprintf('Warning: no lower bound in wave_profile.m, ds= %g', ds)
    return   %if no lower bound then no wave profile
end
dw2 = ds;
fdw2= fun_dw(dw2,Hrms,Tp,rhoc,rhow,tau,me,ws,ks,visc);
while fdw2>0                         %force upper bound to have f(x)<0
    dw2=dw2*2;
    fdw2= fun_dw(dw2,Hrms,Tp,rhoc,rhow,tau,me,ws,ks,visc);
end

% function f(hm)=(ero-depo)=>equilib
fd = @(dw) fun_dw(dw,Hrms,Tp,rhoc,rhow,tau,me,ws,ks,visc);
% function to find root of f(dw)
options = optimset('TolX',1e-3);
dw = fzero(fd,[dw1 dw2],options); %built-in function for root of f(x)
%
% width of wave formed profile
fs  = 0.095*(Hrms^2*g*Tp/dw/visc)^-0.187;
fr  = 0.237*(Hrms*sqrt(g*dw)*Tp/4/pi/dw/ks)^-0.52;
cd  = 0.5*max(fs,fr);
yw = 3*pi*dw^2/4/cd/Hrms;
%
%--------------------------------------------------------------------------
%
function fd = fun_dw(dwe,Hrms,Tp,rhoc,rhow,tau,me,ws,ks,visc)
%
% Function used to find the erosion depth of wave profile based on a 
% balance of eorsion and deposition with Uw determined from the hydraulics 
%
    g = 9.81;
    %
    %friction factor
    fs  = 0.095*(Hrms^2*g*Tp/dwe/visc)^-0.187;
    fr  = 0.237*(Hrms*sqrt(g*dwe)*Tp/4/pi/dwe/ks)^-0.52;
    Cd  = 0.5*max(fs,fr);
    %
    Uw = Hrms/2*sqrt(g/dwe);
    %   
    ust = rhow*Cd*Uw^2;
    if  Uw<=0 || dwe<0
        ero = NaN;
    elseif ust-tau<0
        ero=0;
    else
        ero1= 2*ust/Uw*sqrt((ust-tau)/ust)*sqrt(tau/rhow/Cd);
        ero2= (ust-2*tau)*(pi-2*asin(sqrt(tau/rhow/Cd)/Uw));
        ero = me*(ero1+ero2)/2/pi;
    end
    %
    sed = rhoc*ws;
    %
    fd = ero-sed;
%