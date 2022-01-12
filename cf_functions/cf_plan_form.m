function yz = cf_plan_form(zwl,grid)
%
%-------function help------------------------------------------------------
% NAME
%   cf_pan_form.m
% PURPOSE
%   compute planform variation along the x-axis at specified levels 
% USAGE
%   yz = cf_plan_form(zwl,grid)
% INPUTS
%   zwl - struct or cell array of levels, constant value for each level
%   grid - struct of x, y, z (eg as used in getGrid in the GDinterface)
%          NB: at the moment this can only handle a half-grid ****
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
    xi = grid.x;
    nx = length(xi);
    if isstruct(zwl)
        zwl = struct2cell(zwl);
    end
    nlevels = length(zwl);
%     if isvector(zwl)              
%         %single value for each level
%         nlevels = length(zwl);
%     else
%         yz=[]; return;
%         % %levels specified along the x-axis of the grid
%         % [m,n] = size(zwl);
%         % yz = [];
%         % if m==n
%         %     warndlg('Unable to determine number of levels in cf_plan_form')
%         %     return;
%         % elseif m==nx
%         %     nlevels = n;
%         % elseif n==nx
%         %     nlevels = m;
%         % else
%         %     warndlg('Levels variable, zwl, does not match x-dimension in cf_plan_form')
%         %     return;
%         % end
%     end
    
    zmax = max(grid.z,[],'all');
    if zmax<zwl{1}
        %grid does not rxtend to zhw
        warndlg('Grid does not extend to high water level. Plan form not set');
        yz = [];
        return;
    end
    hf = figure('Tag','PlotFig','Visible','off');
    ax = axes(hf);
    yz = zeros(nx,nlevels);
    for j=1:nlevels
        M = contour(ax,grid.y,xi,grid.z,[zwl{j},zwl{j}]);
        xj = M(2,2:end);
        yj = M(1,2:end); 
        [xj,idx] = unique(xj,'stable');
        yj = yj(idx);               
        yz(:,j) = interp1(xj,yj,xi,'linear');
    end
    delete(hf);
    yz = num2cell(yz',2)';      %formatted to load into dstable
end