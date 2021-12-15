function [hav,Wrv,Arv] = get_river_regime(obj,tr,Qr)
%
%-------function help------------------------------------------------------
% NAME
%   get_river_regime.m
% PURPOSE
%   get river profile using regime theoretical profile (Cao & Knight,1996)
% USAGE
%   [hrv,bh,mur] = get_river_regime(mobj,tr,Qr)
% INPUTS
%   obj - CF_FormModel class instance
%   tr - tidal range (m)
%   Qr - river discharge (m3/s) - optional
% OUTPUT
%   hav - hydraulic depth of river regime channel (m)
%   Wrv - width of regime channel (m)
%   Arv - CSA of regime channel (m2)
% SEE ALSO
%   used in channel_form_models and pr_form_model
%
% Author: Ian Townend
% CoastalSEA (c) Jan 2022
%--------------------------------------------------------------------------
%
    hydobj = obj.RunParam.CF_HydroData;
    sedobj = obj.RunParam.CF_SediData;
    cnsobj = muiConstants.Evoke;
    %assign properties
    Le = hydobj.xTidalLimit;       %length of channel (m)
    
    d50river = sedobj.d50river;    %sediment grain size, D50 (m)
    tauriver = sedobj.tauriver;    %critical bed shear stress (Pa)
    rhos = cnsobj.SedimentDensity;
    rhow = cnsobj.WaterDensity;  
    
    if nargin<3
        Qr = hydobj.Qr;  %river discharge used in form model (m3/s)
    end                  %otherwise use input value

    Sr  = tr/Le;                  %mean water surface slope in estuary
    if Qr>0
        [hav,Wrv,~] = river_regime(Qr,Sr,d50river,tauriver,rhos,rhow);
    else
        hav = 0; Wrv = 0;         %no channel if no flow
    end
    Arv = hav*Wrv;
end