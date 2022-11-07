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
        PropertyLabels = {'Include flood plain in mass balance (1 or 0)',...                                                   
                          'High water constraint (1 or 0)',...
                          'Geological constraint (1 or 0)',...
                          'Start of geological constraints (m)',...
                          'End of geological constraints (m)',...
                          'Integration distance from mouth (m)',...
                          'Flood plain offset above high water (m)',...
                          'Open coast transgression ratio (L/h)',...
                          'Include migrating meander (NaN, 1, or 0)',...
                          'Constant sediment flux (m^3/yr): 0,NaN,or value',...
                          'Current (0) or Source-model (1) water levels'};
        %abstract properties in muiPropertyUI for tab display
        TabDisplay   %structure defines how the property table is displayed 
    end
    
    properties
        inclFloodPlain = false        
        inclHWConstraint = false %flag to include a constraint at high water
        inclGeoConstraint = false%flag to include geological constaints
        StConstraints = 0        %start and end distance from mouth for 
        NdConstraints = 0        %a set of geological constraints         
        IntDist = 0              %distance to use for volume integration (m)
        FPoffset = 0.2           %flood plain offset - default consistent with form models (m)
        BruunRatio = 0           %open coast erosion ratio, L/h (Bruun scaling of dx = slr.L/h)  
        isMeander = NaN          %include meander - NaN=exclude,1=migrate,0=fixed
        SedFlux = 0              %user specified value of sediment exchange
        isModelWL = true         %flag to select water levels to use - true uses model values
    end                            

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
            if length(obj.StConstraints)~=length(obj.NdConstraints)
                warndlg('Number of start and end constraints must be the same')
            end
            setClassObj(mobj,'Inputs',classname,obj);
        end     
    end
%%        
        %add other functions to operate on properties as required   
end