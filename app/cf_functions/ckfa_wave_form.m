function [dhw,yhw,dlw,ylw] = ckfa_wave_form(inp,hm,Wlw,Whw)
%
%-------header-------------------------------------------------------------
% NAME
%   ckfa_wave_form.m [NB modified from version in CKFA model}*********
% PURPOSE
%   Obtain depth and width of wave formed profile at high and low water
% USAGE
%   [dw yw] = ckfa_wave_form(Uw,zw,Fch,cnc,df,ds,d50,tau,me,ws)
% INPUTS
%   inp is a struct with fields
%       amp    = tidal amplitude (m)
%       omega =  angular frequency, 2pi/Tp (1/s)
%       Uw  = wind speed (m/s)
%       zw  = elevation of wind speed (m) - default is 10m
%       rhow = density of water (kg/m^3)
%     	rhoc = suspended sediment mass concentration (kg/m^3)
%       taucr= critical threshold bed shear stress (Pa)
%       d50  = median sediment grain size diameter (m)
%       ws   = sediment fall velocity (m/s)
%       me   = erosion rate coeficient (kg/N/s)
%       g    = acceleration due to gravity (m/s2)
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
    amp = inp.amp;
    %width and depth at high water
    Fhw = sqrt(2)*Whw;             %fetch at high water
    hhw = hm+amp;                  %hydraulic depth at high water
    cnc = inp.rhoc*hm/hhw;         %modified concentration at high water
    %depth and width of wave profile
    [dhw,yhw] = ckfa_wave_profile(inp,inp.Uw,Fhw,cnc,hhw,2*amp);
    
    %width and depth at low water
    Flw = sqrt(2)*Wlw;             %fetch at low water
    if hm>amp
    hlw = hm-amp;                  %hydraulic depth at low water
    else
    hlw = 0.5;                     %minimum value for drainage channel
    end
    cnc = inp.rhoc*hm/hlw;         %modified concentration at low water
    %depth and width of wave profile
    [dlw,ylw] = ckfa_wave_profile(inp,inp.Uw,Flw,cnc,hlw,hlw);
end