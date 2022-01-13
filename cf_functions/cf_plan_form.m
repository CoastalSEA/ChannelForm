function yz = cf_plan_form(grid,zwl)
%
%-------function help------------------------------------------------------
% NAME
%   cf_pan_form.m
% PURPOSE
%   compute planform variation along the x-axis at specified planar levels 
% USAGE
%   yz = cf_plan_form(grid,zwl)
% INPUTS
%   grid - struct of x, y, z (eg as used in getGrid in the GDinterface)
%   zwl - struct or cell array of levels, constant value for each level
% OUTPUTS
%   yz - cell array of widths for each level 
% SEE ALSO
%   called by cf_valley_model as part of ChannelForm model and
%   addFormProperties in GDinterface
%
% Author: Ian Townend
% CoastalSEA (c) Jan 2022
%--------------------------------------------------------------------------
%
    if isstruct(zwl)
        zwl = struct2cell(zwl);
    end
    nlevels = length(zwl);
    xi = grid.x;
    nx = length(xi);
    
    dely = grid.y(2)-grid.y(1);  %grid interval
    yz = zeros(nx,nlevels);
    
    zi = grid.z;
    if grid.ishead  %orientation of x-axis, x=0 is nearest the head
        zi = flipud(zi);
    end
    
    zmax = max(zi,[],'all');
    if zmax<zwl{1}(1)
        %grid does not extend to zhw
        warndlg('Grid does not does not extend to high water. Plan form not set');
        yz = num2cell(yz',2)';      %formatted to load into dstable
        return;
    end

    for j=1:nlevels
        zij = zi;
        zij(zij>zwl{j}(1)) = NaN;             %set values above zwl to NaN
        for i=1:nx
            yz(i,j) = sum(~isnan(zij(i,:))*dely); %width at zwl
        end
    end
    yz = num2cell(yz',2)';      %formatted to load into dstable
end