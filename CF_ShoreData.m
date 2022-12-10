classdef CF_ShoreData < muiPropertyUI             
%
%-------class help---------------------------------------------------------
% NAME
%   CF_ShoreData.m
% PURPOSE
%   Class for shoreline form parameters used in ChannelForm model
% USAGE
%   obj = CF_ShoreData.setInput(mobj); %mobj is a handle to Main UI
% SEE ALSO
%   inherits muiPropertyUI
%
% Author: Ian Townend
% CoastalSEA (c) Jan 2021
%--------------------------------------------------------------------------
%      
    properties (Hidden)
        %abstract properties in muiPropertyUI to define input parameters
        PropertyLabels = {'Width of shore to depth of closure (m)',...
                          'Depth of shore closure from mtl (m)',...
                          'Offshore bed slope (1:obs)',...
                          'Upper beach slope (1:ubs)'...
                          'Beach crest height above mean sea level (m)'};
        %abstract properties in muiPropertyUI for tab display
        TabDisplay   %structure defines how the property table is displayed 
    end
    
    properties        
        ShoreWidth = 5000    %shore normal width to depth of closure (m) 
        ShoreDepth = 8       %closure depth of profile from mtl (m) 
        OffshoreBS = 5000    %offshore bed slope - defines bed alope of channel
        UpperBeachBS = 20    %upper beach slope - defines slope above msl
        UpperBeachHeight = 5 %beach crest height above mean sea level (m)
    end    

%%   
    methods (Access=protected)
        function obj = CF_ShoreData(mobj)      
            %constructor code:            
            %TabDisplay values defined in UI function setTabProperties used to assign
            %the tabname and position on tab for the data to be displayed
            obj = setTabProps(obj,mobj);  %muiPropertyUI function
        end 
    end
%%  
    methods (Static)  
        function obj = setInput(mobj,editflag)
            %gui for user to set Parameter Input values
            classname = 'CF_ShoreData';             
            obj = getClassObj(mobj,'Inputs',classname);
            if isempty(obj)
                obj = CF_ShoreData(mobj);         
            end
            
            %use muiPropertyUI function to generate UI
            if nargin<2 || editflag
                %add nrec to limit length of props UI (default=12)
                obj = editProperties(obj);  
                %add any additional manipulation of the input here
            end
            setClassObj(mobj,'Inputs',classname,obj);
        end     
    end       
%%  
    methods
        function grid = setShoreline(obj,grid,z0,isfull)
            %create shoreline strip based on equilibrium profile 
            % z0 - elevation of mean sea level to model datum
            % isfull - flag return coastal strip added to input grid if
            % true, otherwise just returns coastal strip
            %function intended for model grids where xM=0 in the source grid
            %and this is modified by the width of the shore (xM=Shorewidth)
            [~,ixM] = gd_basin_indices(grid); %nearest grid point
            delx = abs(grid.x(2)-grid.x(1));
            
            %parameters for beach profile
            obs = 1/obj.OffshoreBS;     %offshore bed slope
            ubs = 1/obj.UpperBeachBS;   %beach slope above msl
            zBC = obj.UpperBeachHeight; %beach crest height above msl
            nint = round(obj.ShoreWidth/delx);
            xS = nint*delx;             %width of shore along x-axis
            newx = 0:delx:(xS-delx);    %shore x-axis
            zdc = z0-obj.ShoreDepth;    %elevation of closure depth
            z1km = [xS,zdc];   %distance and depth of offshore limit
            %beach shore profile
            [~,zp,~] = deanbeachprofile(newx,zBC,z1km,ubs,false);
            zp = flipud(zp(1:nint));  %function adds upper beach so clip to required size
            shorez = repmat(zp,1,length(grid.y));
            
            %slope seaward from estuary mouth 
            zM = grid.z(ixM,:); %elevation at mouth
            estz = zM-(xS-newx(1:nint)')*obs; %bed sloping at 1:obs from grid levels at shoreline
            shorez = min(shorez,estz);            
            if isfull
                %return original grid with coastal strip added
                if grid.x(1)~=0   %grid co-ordinates are being used 
                    dx = grid.x(2)-grid.x(1);
                    newx = grid.x(ixM)-sign(dx)*xS:dx:(grid.x(ixM)-dx);
                    grid.x = [newx';grid.x(ixM:end)];
                    if sign(dx)>0
                        grid.xM = grid.x(1)+xS;    %distance to shore from x(1)
                    else
                        %grid xM does not change if axis is descending
                    end
                else       %model grid so move origin to include shoreface
                    grid.x = [newx';grid.x(ixM:end)+xS]; 
                    grid.xM = xS;                  %distance to shore from x=0
                end
                grid.z = [shorez;grid.z(ixM:end,:)];                
            else
                %return just the coastal strip
                grid.x = newx';
                grid.z = shorez;
                grid.xM = xS;                     %distance to shore from x=0
            end
        end
    end  
end