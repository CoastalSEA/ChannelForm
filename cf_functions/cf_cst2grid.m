function hydobj = cf_cst2grid(obj,res,cstx,isflip)
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
%   res -  struct of results from cst_model (part of CSTmodel App)
%   cstx - distances used in cst_model
%   isflip - logical true if results are to be reversed (optional, default is false)
% OUTPUTS
%   hydobj - updated version of CF_HydroData thaat is assigned to obj.RunParam
% SEE ALSO
%   used in cf_exp_model, cf_pow_model and ckfa_form_model as part 
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
    xi = getGridDimensions(grdobj)';   %x co-ordinates         
    
    if isflip
        res.X = cellfun(@fliplr,res.X,'UniformOutput',false);
        res.F = cellfun(@fliplr,res.F,'UniformOutput',false);        
    end
    
    z = interp1(cstx,res.X{1},xi,'linear','extrap'); %water level
    a = interp1(cstx,res.X{2},xi,'linear','extrap'); %elevation amplitude
    U = interp1(cstx,res.X{4},xi,'linear','extrap'); %tidal velocity amplitude
    v = interp1(cstx,res.X{5},xi,'linear','extrap'); %river velocity 
    d = interp1(cstx,res.F{2},xi,'linear','extrap'); %hydraulic depth
    hydobj.cstres = struct('z',z,'a',a,'U',U,'v',v,'d',d);
    
    hydobj.zhw = z+a;  %mean tide level+tidal amplitude
    hydobj.zmt = z;    %mean tide level
    hydobj.zlw = z-a;  %mean tide level-tidal amplitude 
end