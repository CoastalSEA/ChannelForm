function hydobj = cf_interp2grid(obj,resX,cstx,isflip)
%
%-------function help------------------------------------------------------
% NAME
%   cf_interp2grid.m
% PURPOSE
%   interpolate the CST water levels onto the model grid
% USAGE
%   cstvar = cf_interp2grid(obj,resX,cstx,isflip)
% INPUTS
%   obj - CF_FormModel class instance
%   resX - along-channel results from cst_model (part of CSTmodel App)
%   cstx - distances used in cst_model
%   isflip - logical true if results are to be reversed (optional, default is false)
% OUTPUTS
%   cstvar - struct containing the  
% SEE ALSO
%   used in channel_form_model and ckfa_form_model as part of ChannelForm model
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
%     Lt = diff(grdobj.XaxisLimits);   %length of model domain (m)
%     nintx = grdobj.Xint;             %no of intervals in the x direction
%     delx = Lt/nintx;                 %x interval
%     xi = 0:delx:Lt;                  
    
    if isflip
        resX = cellfun(@fliplr,resX,'UniformOutput',false);
    end

    cst_hw = resX{1}+resX{2};        %mean tide level+tidal amplitude
    cst_mt = resX{1};                %mean tide level
    cst_lw = resX{1}-resX{2};        %mean tide level-tidal amplitude

    hydobj.zhw = interp1(cstx,cst_hw,xi,'linear','extrap'); %high water
    hydobj.zmt = interp1(cstx,cst_mt,xi,'linear','extrap'); %mean water
    hydobj.zlw = interp1(cstx,cst_lw,xi,'linear','extrap'); %low water
    U = interp1(cstx,resX{3},xi,'linear','extrap'); %tidal velocity amplitude
    v = interp1(cstx,resX{4},xi,'linear','extrap'); %river velocity 
    d = interp1(cstx,resX{5},xi,'linear','extrap'); %hydraulic depth
    hydobj.cstres = struct('U',U,'v',v,'d',d);
end