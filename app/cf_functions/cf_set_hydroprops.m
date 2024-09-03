function [obj,ok] = cf_set_hydroprops(obj)
%
%-------function help------------------------------------------------------
% NAME
%   cf_set_hydroprops.m
% PURPOSE
%   set water levels for form model using either the surface defined by the
%   cst_model or high water at the mouth and a reducing tidal amplitude
% USAGE
%   obj = cf_set_hydroprops(obj)
% INPUTS
%   obj - CF_FormModel class instance with grid and water level model
%         inputs defined in RunParam versions of GD_GridProps and CF_HydroData
%         obj.Selection.wlflag - flag to indicate type of water surface to use
%            0=CSTmodel used to define water levels
%            1=constant HW tapering LW 
%            2=constant HW & LW
% OUTPUTS
%   obj - CF_FormModel class instance updated with water levels
%   ok - flag to indicate whether cst_model found a solution
% NOTE
%   x-axis landward and zero at mouth (ie xM=0)
% SEE ALSO
%   used in channel_form_model and ckfa_form_model as part of ChannelForm model
%
% Author: Ian Townend
% CoastalSEA (c) Jan 2022
%--------------------------------------------------------------------------
%
    ok = 0; 
    hydobj = obj.RunParam.CF_HydroData; 
    grdobj = obj.RunParam.GD_GridProps;
    xi = getGridDimensions(grdobj);
    nx = length(xi);
    
    if obj.Selection.wlflag==0
        %use cst_model to set-up water levels for model           
        [res,xyz,~] = runHydroModel(hydobj,obj);
        if isempty(res), return; end
        %interpolate CSTmodel results onto model grid + reverse for xmin @ head
        obj.RunParam.CF_HydroData = cf_cst2grid(obj,res,xyz{:},false);
    else
        %use constant high water level and linear reducing low water
        Lt = obj.RunParam.CF_HydroData.xTidalLimit; %distance from mouth to tidal limit 
        zhw = hydobj.zhw(1); 
        zlw = hydobj.zlw(1);
        zHWxi = ones(1,nx)*zhw;      %assume constant HW surface
        amp0 = (zhw-zlw)/2;          %tidal amplitude at mouth
        if obj.Selection.wlflag==1            
            ampx = amp0*(1-xi'/Lt);  %linear reducing low water
            ampx(xi>Lt) = 0;
        else
            ampx = amp0;             %use constant high and low water level 
        end
        obj.RunParam.CF_HydroData.zhw = zHWxi;        %high water
        obj.RunParam.CF_HydroData.zmt = zHWxi-ampx;   %mean tide level
        obj.RunParam.CF_HydroData.zlw = zHWxi-2*ampx; %low water
        obj.RunParam.CF_HydroData.cstres = [];
    end
    ok = 1;
end
