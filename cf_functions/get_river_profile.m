function [hrv,bh,mur] = get_river_profile(yi,mobj)
%
%-------function help------------------------------------------------------
% NAME
%   get_river_profile.m
% PURPOSE
%   get river profile using regime theoretical profile (Cao & Knight,1996)
% USAGE
%   [hrv,bh,mur] = get_river_profile(yi,mobj)
% INPUTS
%   yi - y co-ordinates
%   mobj - ChannelForm model UI instance
% OUTPUT
%   hrv - depths of river regime channelfrom 0-bh (ie half-width)
%   bh  - half width of regime channel 
%   mur - effective submerged static coefficient of Coulomb friction 
%         (estimated from geometry)
% SEE ALSO
%   used in channel_form_models and pr_form_model
%
% Author: Ian Townend
% CoastalSEA (c) Jan 2022
%--------------------------------------------------------------------------
%
    expobj = getClassObj(mobj,'Inputs','CF_ExpData');
    hydobj = getClassObj(mobj,'Inputs','CF_HydroData');
    sedobj = getClassObj(mobj,'Inputs','CF_SediData');
    cnsobj = mobj.Constants;    
    %assign properties
    Le = expobj.TotalLength;         %total length of channel (m)
    nc = expobj.ChannelShapeParam;   %channel shape parameter (-)
    d50river = sedobj.d50river;      %sediment grain size, D50 (m)
    tauriver = sedobj.tauriver;      %critical bed shear stress (Pa)
    rhos = cnsobj.SedimentDensity;
    rhow = cnsobj.WaterDensity;   
    Qr = hydobj.RiverDischarge;    

    tr = (hydobj.zhw(1)-hydobj.zlw(1)); %tidal range at mouth
    Sr  = tr/Le;                %mean water surface slope in estuary
    if Qr>0
        [hav,Wrv,~] = river_regime(Qr,Sr,d50river,tauriver,rhos,rhow);
    else
        hav = 1; Wrv = 8;  %dummy head
    end

    bh = Wrv/2;            %half-width at head when mean water level = HW
    hc = hav*(nc+1)/nc;
    mur = nc*hc/bh(1);%submerged static coefficient of Coulomb friction (estimated from geometry)
    hrv = hc*(1-(yi/bh(1)).^nc).*(yi<bh(1));
end