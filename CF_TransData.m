classdef CF_TransData < muiPropertyUI                
%
%-------class help---------------------------------------------------------
% NAME
%   CF_TransData.m
% PURPOSE
%   Class for transgression parameters used in ChannelForm model
% USAGE
%   obj = CF_TransData.setInput(mobj); %mobj is a handle to Main UI
% SEE ALSO
%   inherits muiPropertyUI
%
% Author: Ian Townend
% CoastalSEA (c) Jan 2021
%--------------------------------------------------------------------------
%      
    properties (Hidden)
        %abstract properties in muiPropertyUI to define input parameters
        PropertyLabels = {'Rate of sea level rise (m/yr)',...
                          'Linear or exponential slr (1 or 0)',...
                          'Change in tidal amplitude (m)',...
                          'Open coast transgression ratio (-)',...
                          'Integration distance from mouth (m)',...
                          'High water constraint (1 or 0)',...
                          'Start of channel bed constraints (m)',...
                          'End of channel bed constraints (m)',...
                          'Constant sediment flux (m^3/yr) if reqd'};
        %abstract properties in muiPropertyUI for tab display
        TabDisplay   %structure defines how the property table is displayed 
    end
    
    properties
        dSLR                     %rate of sea level rise (m/yr)
        isLinSLR = 1             %linear or exponential slr (1 or 0)
        delTidalAmp = 0          %change in tidal amplitude (m) - linear increase
        BruunRatio = 0           %open coast erosion ratio (Bruun scaling of dx = slr.L/h)
        IntDist                  %distance to use for volume integration (m)
                                 %uses estuary length if empty
        isConstrained = 0        %flag to include a constraint at high water
        StConstraints            %start and end distance from mouth for 
        NdConstraints            %a set of geological constraints   
        SedFlux                  %user specified value of sediment exchange
    end                          %calculated using get_sed_flux if empty   

%%   
    methods (Access=protected)
        function obj = CF_TransData(mobj)         
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
            classname = 'CF_TransData';           
            obj = getClassObj(mobj,'Inputs',classname);
            if isempty(obj)
                obj = CF_TransData(mobj);              
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