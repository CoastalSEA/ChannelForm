function [hrv,bh,mur] = get_river_profile(obj,tr,yi)
%
%-------function help------------------------------------------------------
% NAME
%   get_river_profile.m
% PURPOSE
%   get river profile using regime theoretical profile (Cao & Knight,1996)
% USAGE
%   [hrv,bh,mur] = get_river_profile(obj,tr,yi)
% INPUTS
%   obj - CF_FormModel instance with RunParam to use
%   tr - tidal range (m)
%   yi - y co-ordinates
% OUTPUT
%   hrv - depths of river regime channel from 0-bh (ie half-width)
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
    frmobj = obj.RunParam.CF_FormData;
    nc = frmobj.ChannelShapeParam;           %channel shape parameter (-)
    [hav,Wrv,~] = get_river_regime(obj,tr);  %Cao & Knight, 1996
    if hav>0
        bh = Wrv/2;         %half-width at head when mean water level = HW
        hc = hav*(nc+1)/nc; %depth at centre-line of channel section
        mur = nc*hc/bh(1);  %submerged static coefficient of Coulomb friction (estimated from geometry)
        hrv = hc*(1-(yi/bh(1)).^nc).*(yi<bh(1));
    else
        hrv = 0; bh = 0; mur = 0;
    end
end