function hydobj = cf_cst2grid(obj,resX,cstx,isflip)
%
%-------function help------------------------------------------------------
% NAME
%   cf_cst2grid.m
% PURPOSE
%   interpolate the CST water levels onto the model grid
% USAGE
%   cstvar = cf_cst2grid(obj,resX,cstx,isflip)
% INPUTS
%   obj - CF_FormModel class instance
%   resX - along-channel results from cst_model (part of CSTmodel App)
%   cstx - distances used in cst_model
%   isflip - logical true if results are to be reversed (optional, default is false)
% OUTPUTS
%   hydobj - updated version of CF_HydroData thaat is assigned to obj.RunParam
% SEE ALSO
%   used in channel_form_model, pr_form_model and ckfa_form_model as part 
%   of ChannelForm model
%
% Author: Ian Townend
% CoastalSEA (c) Jan 2022
%--------------------------------------------------------------------------
%
    if nargin<4
        isflip = false;
    end
    hydobj = obj.RunParam.CF_HydroData;
    grdobj = obj.RunParam.GD_GridProps;
     
    %set-up co-ordinate system
    xi = getGridDimensions(grdobj);   %x co-ordinates         
    
    if isflip
        resX = cellfun(@fliplr,resX,'UniformOutput',false);
    end
    
    z = interp1(cstx,resX{1},xi,'linear','extrap'); %water level
    a = interp1(cstx,resX{2},xi,'linear','extrap'); %elevation amplitude
    U = interp1(cstx,resX{3},xi,'linear','extrap'); %tidal velocity amplitude
    v = interp1(cstx,resX{4},xi,'linear','extrap'); %river velocity 
    d = interp1(cstx,resX{5},xi,'linear','extrap'); %hydraulic depth
    hydobj.cstres = struct('z',z,'a',a,'U',U,'v',v,'d',d);
    
    hydobj.zhw = z+a;  %mean tide level+tidal amplitude
    hydobj.zmt = z;    %mean tide level
    hydobj.zlw = z-a;  %mean tide level-tidal amplitude 
end