function obj = cf_offset_wls(obj,isextend)
%
%-------function help------------------------------------------------------
% NAME
%   cf_offset_wls.m
% PURPOSE
%   translate water levels landwards when mouth is some distance from x=0 
%   and pad the water level vectors with values at the mouth
% USAGE
%   obj = cf_offset_wls(obj,isextend)
% INPUTS
%   obj - CF_FormModel or CF_TransModel class instance with grid and 
%         water level model inputs defined in RunParam versions of 
%         GD_GridProps and CF_HydroData
%         NB: grid is in obj.Data.Grid for saved models, whereas the grid
%         is in obj.Grid in given time step (i.e. transient).
%   isextend - true extends the water level vectors by adding values,
%              whereas false maintains the size of the vectors by removing 
%              an equivalent number of records from the upstream end            
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
    if isprop(obj,'Grid') && ~isempty(obj.Grid)
        %called from CF_TransModel and using RunParam.CF_HydroProps
        grid = obj.Grid;
        hydobj = obj.RunParam.CF_HydroData;
        istransient = true;
    elseif isprop(obj,'Data') && isfield(obj.Data,'Grid') && ...
                                                   ~isempty(obj.Data.Grid)
       %
       grid = getGrid(obj,1);
       hydobj = obj.Data.WaterLevels;
       istransient = false;
    else
        grid = [];
    end
    
    if  ~isempty(grid) && grid.xM>0
        %coast is not at x=0        
        [~,ixM] = gd_basin_indices(grid);      %x-index of coast
        %delx = abs(obj.Grid.x(2)-obj.Grid.x(1));
        %ixM = floor(obj.Grid.xM/delx)+1;           %x-index of coast
        if ~isempty(ixM) && ixM>1
            zhw = hydobj.zhw;   %high water
            zmt = hydobj.zmt;   %mean tide level
            zlw = hydobj.zlw;   %low water

            mask = ones(1,ixM-1);
            xrec = length(zhw);
            ix0 = 1;                %index of first wls at zwl(1)
            if isextend             %extend wl vector seaward
                xrec = min([xrec,length(grid.x)-ixM+1]);
                range = 1:xrec;              
            else                    %wl vector correct length
                if isnan(zhw(end))        %add seaward and remove landward values
                    range = 1:xrec-ixM+1;
                elseif isnan(zhw(1))
                    range = ixM:xrec;     %overwrite seaward values
                    ix0 = ixM;            %first wls are at zwl(ixM)
%                 elseif xrec~=length(grid.x)
%                     range = 1:xrec;
                else
                    range = 1:xrec-ixM+1; %add seaward and remove landward values
                end
            end

            zHWxM = [mask*zhw(ix0),zhw(range)];
            zMTxM = [mask*zmt(ix0),zmt(range)];     
            zLWxM = [mask*zlw(ix0),zlw(range)];            
            
            if istransient
                obj.RunParam.CF_HydroData.zhw = zHWxM; %high water
                obj.RunParam.CF_HydroData.zmt = zMTxM; %mean tide level
                obj.RunParam.CF_HydroData.zlw = zLWxM; %low water
            else
                obj.Data.WaterLevels.zhw = zHWxM; %high water
                obj.Data.WaterLevels.zmt = zMTxM; %mean tide level
                obj.Data.WaterLevels.zlw = zLWxM; %low water                
            end
        end
    end
end