classdef CF_GridData < muiPropertyUI    
%
%-------class help---------------------------------------------------------
% NAME
%   CF_GridData.m
% PURPOSE
%   Class for grid and time step parameters used in ChannelForm model
% USAGE
%   obj = CF_GridData.setInput(mobj); %mobj is a handle to Main UI
% SEE ALSO
%   inherits muiPropertyUI
%
% Author: Ian Townend
% CoastalSEA (c) Jan 2021
%--------------------------------------------------------------------------
%      
    properties (Hidden)
        %abstract properties in muiPropertyUI to define input parameters
        PropertyLabels = {'X-axis definition [x0 xN]',...
                          'No. of intervals in the x direction',...
                          'X-coordinate of mouth',...
                          'Y-axis definition [y0 yN]',...
                          'No. of intervals in the y direction',...
                          'Y-coordinate of centre-line',...
                          'Vertical resolution for hypsometry'};
        %abstract properties in muiPropertyUI for tab display
        TabDisplay   %structure defines how the property table is displayed 
    end
    
    properties
        XaxisLimits = [0 100000]   %definition of the x-axis [x0 xN]         
        Xint = 500                 %no of intervals in the x direction
        Xmouth = 0                 %x co-ordinate of mouth
        YaxisLimits = [-5000 5000] %definition of the y-axis [y0 yN]
        Yint = 500                 %no of intervals in the y direction
        Ycentre = 0                %y co-ordinate of centre-line
        histint = 0.1;             %vertical resolution for hypsometry histogram
    end    

%%   
    methods (Access=protected)
        function obj = CF_GridData(mobj)          
            %constructor code:            
            %TabDisplay values defined in UI function setTabProperties used to assign
            %the tabname and position on tab for the data to be displayed
            obj = setTabProps(obj,mobj);  %muiPropertyUI function
            
            %to use non-numeric entries then one can either pre-assign 
            %the values in the class properties defintion, above, or 
            %specify the PropertyType as a cell array here in the class 
            %constructor, e.g.:
            % obj.PropertyType = [{'datetime','string','logical'},...
            %                                       repmat({'double'},1,8)];
        end 
    end
%%  
    methods (Static)  
        function obj = setInput(mobj,editflag)
            %gui for user to set Parameter Input values
            classname = 'CF_GridData';        
            obj = getClassObj(mobj,'Inputs',classname);
            if isempty(obj)
                obj = CF_GridData(mobj);         
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
        %add other functions to operate on properties as required   
end