function [obj,ok] = cf_set_hydroprops(obj,iscst)
%
%-------function help------------------------------------------------------
% NAME
%   cf_set_hydroprops.m
% PURPOSE
%   set water levels for form model using either the surface defined by the
%   cst_model or high water at the mouth and a reducing tidal amplitude
% USAGE
%   obj = cf_set_hydroprops(obj,iscst)
% INPUTS
%   obj - CF_FormModel class instance 
%   iscst - logical flag, true if CSTmodel to be used to define water levels
% OUTPUTS
%   obj - CF_FormModel class instance updated with water levels
%   ok - flag to indicate whether cst_model found a solution
% SEE ALSO
%   used in channel_form_model and ckfa_form_model as part of ChannelForm model
%
% Author: Ian Townend
% CoastalSEA (c) Jan 2022
%--------------------------------------------------------------------------
%
    hydobj = obj.RunParam.CF_HydroData; 
    if iscst
        %use cst_model to set-up water levels for model           
        [resX,~,~,xyz] = runHydroModel(hydobj,obj);
        if isempty(resX), ok = 0; return; end
        %interpolate CSTmodel results onto model grid + reverse for xmin @ head
        obj.RunParam.CF_HydroData = cf_interp2grid(obj,resX,xyz{:},true);
    else
        %use constant high water level and linear reducing low water
        grdobj = obj.RunParam.GD_GridProps;
        Lt = diff(grdobj.XaxisLimits);   %length of model domain (m)
%         nintx = grdobj.Xint;             %no of intervals in the x direction
%         xi = 0:Lt/nintx:Lt;
        xi = getGridDimensions(grdobj);
        zhw = hydobj.zhw(end); zlw = hydobj.zlw(end);
        amp0 = (zhw-zlw)/2;                    %tidal amplitude at mouth
        ampx = amp0*(xi/Lt);
        zHWxi = ones(size(xi))*zhw;            %assume constant HW surface
        obj.RunParam.CF_HydroData.zhw = zHWxi;        %high water
        obj.RunParam.CF_HydroData.zmt = zHWxi-ampx;   %mean tide level
        obj.RunParam.CF_HydroData.zlw = zHWxi-2*ampx; %low water
        obj.RunParam.CF_HydroData.cstres = [];
    end
    ok = 1;
end
