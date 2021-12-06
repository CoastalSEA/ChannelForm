function [xi,yi,zgrd,yz] = ckfa_form_model(obj,mobj,isfull)
%
%-------function help------------------------------------------------------
% NAME
%   ckfa_form_model.m
% PURPOSE
%   construct idealised channel form using 3D CKFA model
% USAGE
%   [xi,yi,zgrd,yz] = ckfa_form_model(obj,inputs,isfull)
% INPUTS
%   obj - CF_FormModel class instance
%   mobj - ChannelForm model UI instance
%   isfull - true returns full grid, false half-grid
% OUTPUTS
%   xi - x co-ordinate (m)
%   yi - y co-ordinate (m)
%   zgrd - bed elevation grid (m)
%   yz - width at hw,mt,lw (m)
% NOTES
%   CKFA cross-section comprises a channel using Cao&Knight section and an
%   intertidal using the profile proposed by Friedrichs & Aubrey (tide only)
% SEE ALSO
%   used in CF_FormModel as part of ChannelForm model
%
% Author: Ian Townend
% CoastalSEA (c) Jan 2022
%--------------------------------------------------------------------------
%
    if nargin<3
        isfull = true;
    end
    %get the required input parameter classes
    expobj = getClassObj(mobj,'Inputs','CF_ExpData');
    grdobj = getClassObj(mobj,'Inputs','CF_GridData');
    hydobj = getClassObj(mobj,'Inputs','CF_HydroData');
    sedobj = getClassObj(mobj,'Inputs','CF_SediData');
    
    %model run parameters
    Lt = diff(grdobj.XaxisLimits);   %length of model domain (m)
    nintx = grdobj.Xint;             %no of intervals in the x direction
    bt = diff(grdobj.YaxisLimits)/2; %half width of model domain (m)   
    ninty = grdobj.Yint/2;           %no of intervals in the y direction
    
    %set-up co-ordinate system
    delx = Lt/nintx;
    dely = bt/ninty;
    xi = 0:delx:Lt;
    yi = 0:dely:bt;   yi(1) = 0.01;      %the offset ensures no duplicates
                                         %values when matrix mirrored     
    zi = zeros(length(xi),length(yi));
    yz = zeros(length(xi),3);
    
    %water level properties based on amplitude+mtl or CST model (mAD)
    am0 = (hydobj.zhw(1)-hydobj.zlw(1))/2;     %tidal amplitude at mouth
    if isscalar(hydobj.zhw)
        am = am0*((xi+(Lt-Ll))/Lt);            %increasing tidal amplitude 
        zHWxi = ones(size(xi))*hydobj.zhw;     %assume constant HW surface
        zLWxi = zHWxi-2*am;                    %reduce tidal range upstream
    else
        zHWxi = hydobj.zhw; 
        zLWxi = hydobj.zlw; 
    end
    
    for ix=1:length(xi)
        zhw = zHWxi(ix); zlw = zLWxi(ix);
        ax = (zhw-zlw)/2;               %tidal amplitude(m) at x(j)
        mtl = zhw-ax+zo;                %mean tide level(mOD)
    end
    
    
    %generate complete 3D channel form by mirroring half section
    if isfull                        %return full grid
        zgrd = cat(2,fliplr(zi),zi);
        yi  = [-fliplr(yi), yi];
    else                             %return half grid
        zgrd = zi;
    end
end