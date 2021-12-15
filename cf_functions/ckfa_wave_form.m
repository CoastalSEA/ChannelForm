function [dhw,yhw,dlw,ylw] = ckfa_wave_form(am,Uw,d50,tau,me,ws,cnc,hm,Wlw,Whw)
%
%-------header-------------------------------------------------------------
% NAME
%   ckfa_wave_form.m [NB modified from version in CKFA model}*********
% PURPOSE
%   Obtain depth and width of wave formed profile at high and low water
% USAGE
%   [dw yw] = ckfa_wave_form(Uw,zw,Fch,cnc,df,ds,d50,tau,me,ws)
% INPUTS
%   am   = tidal amplitude at the mouth (m)
%   Uw  = wind speed (m/s)
%   d50 = median sediment grain size diameter (m)
%   tau = critical threshold bed shear stress (Pa)
%   me  = erosion rate coeficient (kg/N/s)
%   ws  = sediment fall velocity (m/s)
%   cnc = suspended sediment concentration
%   hm  = hydraulic depth to mtl (m)
%   Wlw = width at low water (m)
%   Whw = width at high watere (m)
% OUTPUTS
%   dhw = high water depth at outer edge of wave profile
%   yhw = high water width of wave profile
%   dlw = low water depth at outer edge of wave profile
%   ylw = low water width of wave profile
% SEE ALSO
%   ckfa_form_model.m and ckfa_form_properties.m
%
% Author: Ian Townend
% CoastalSEA (c) Jan 2022
%--------------------------------------------------------------------------
%
    zw   = 10;                    %elevation of wind speed (m)
    
    %width and depth at high water
    Fhw = sqrt(2)*Whw;            %fetch at high water
    hhw = hm+am;                  %hydraulic depth at high water
    %conc = 5*0.0001*rhow/rhos*Ucr;  %SSC based on eqn of Prandle etal, GRL, 2005
    cnc = cnc*hm/hhw;             %modified concentration at high water
    %depth and width of wave profile
    [dhw,yhw] = ckfa_wave_profile(Uw,zw,Fhw,cnc,hhw,2*am,d50,tau,me,ws);
    
    %width and depth at low water
    Flw = sqrt(2)*Wlw;            %fetch at low water
    if hm>am
    hlw = hm-am;                  %hydraulic depth at low water
    else
    hlw = 0.5;                    %minimum value for drainage channel
    end
    cnc = cnc*hm/hlw;             %modified concentration at low water
    %depth and width of wave profile
    [dlw,ylw] = ckfa_wave_profile(Uw,zw,Flw,cnc,hlw,hlw,d50,tau,me,ws);
end