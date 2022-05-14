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
                          'Upper beach slope (1:ubs)'};
        %abstract properties in muiPropertyUI for tab display
        TabDisplay   %structure defines how the property table is displayed 
    end
    
    properties        
        ShoreWidth = 5000    %shore normal width to depth of closure (m) 
        ShoreDepth = 8       %closure depth of profile from mtl (m) 
        OffshoreBS = 5000    %offshore bed slope - defines bed alope of channel
        UpperBeachBS = 20    %upper beach slope - defines slope above 0mAD
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
        function grid = setShoreline(obj,grid,isfull)
            %create shoreline strip based on equilibrium profile 
            % isfull - flag return coastal strip added to input grid if
            % true, otherwise just returns coastal strip
            obs = 1/obj.OffshoreBS;
            ubs = 1/obj.UpperBeachBS;
            
            delx = grid.x(2)-grid.x(1);
            nint = floor(obj.ShoreWidth/delx);
            xM = nint*delx;
            newx = 0:delx:xM-delx;
            zBC = grid.z(1,1);
            z1km = [xM,-obj.ShoreDepth];
            [~,zp,~] = deanbeachprofile(newx,zBC,z1km,ubs,false);
            zp = flipud(zp(1:nint));  %function adds upper beach so clip to required size
            shorez = repmat(zp,1,length(grid.y));
            
            %slope seaward from estuary mouth
            zM = grid.z(1,:); %elevation at mouth
            estz = zM-(xM-newx(1:nint)')*obs;
            shorez = min(shorez,estz);
            if isfull
                %return original grid with coastal strip added
                if grid.x(1)~=0   %grid co-ordinates are being used 
                    newx = grid.x(1)-xM:delx:grid.x(1)-delx;
                    grid.x = [newx';grid.x];
                else       %model grid so move origin to include shoreface
                    grid.x = [newx';grid.x+xM];
                end
                grid.z = [shorez;grid.z];                
            else
                %return just the coastal strip
                grid.x = newx';
                grid.z = shorez;
            end
            grid.xM = xM;
        end
    end  
end