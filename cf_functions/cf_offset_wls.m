function obj = cf_offset_wls(obj)
%
%-------function help------------------------------------------------------
% NAME
%   cf_offset_wls.m
% PURPOSE
%   translate water levels landwards when mouth is some distance from x=0 
%   and pad the water level vectors with values at the mouth
% USAGE
%   obj = cf_offset_wls(obj)
% INPUTS
%   obj - CF_FormModel or CF_TransModel class instance with grid and 
%         water level model inputs defined in RunParam versions of 
%         GD_GridProps and CF_HydroData
% OUTPUTS
%   obj - CF_FormModel or CF_TransModel class instance updated with water levels
% NOTE
%   x-axis increases landward
%   NB 'cstres' is NOT adjusted in code below
% SEE ALSO
%   used in CF_TransModel as part of ChannelForm model
%
% Author: Ian Townend
% CoastalSEA (c) Jan 2022
%--------------------------------------------------------------------------
%
    if isprop(obj,'Grid') && ~isempty(obj.Grid) && obj.Grid.xM>0
        %coast is not at x=0        
        [~,ixM] = gd_basin_indices(obj.Grid);      %x-index of coast
        %delx = abs(obj.Grid.x(2)-obj.Grid.x(1));
        %ixM = floor(obj.Grid.xM/delx)+1;           %x-index of coast
        if ~isempty(ixM) && ixM>1
            zhw = obj.RunParam.CF_HydroData.zhw;   %high water
            zmt = obj.RunParam.CF_HydroData.zmt;   %mean tide level
            zlw = obj.RunParam.CF_HydroData.zlw;   %low water

            mask = ones(1,ixM);
            zHWxM = [mask*zhw(1),zhw(1:end-ixM)];
            zMTxM = [mask*zmt(1),zmt(1:end-ixM)];     
            zLWxM = [mask*zlw(1),zlw(1:end-ixM)];
            obj.RunParam.CF_HydroData.zhw = zHWxM; %high water
            obj.RunParam.CF_HydroData.zmt = zMTxM; %mean tide level
            obj.RunParam.CF_HydroData.zlw = zLWxM; %low water
        end
    end
end