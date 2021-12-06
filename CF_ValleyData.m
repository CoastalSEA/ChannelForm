classdef CF_ValleyData < muiPropertyUI             
%
%-------class help---------------------------------------------------------
% NAME
%   CF_ValleyData.m
% PURPOSE
%   Class for valley form parameters used in ChannelForm model
% USAGE
%   obj = CF_ValleyData.setInput(mobj); %mobj is a handle to Main UI
% SEE ALSO
%   inherits muiPropertyUI
%
% Author: Ian Townend
% CoastalSEA (c) Jan 2021
%--------------------------------------------------------------------------
%      
    properties (Hidden)
        %abstract properties in muiPropertyUI to define input parameters
        PropertyLabels = {'Width & depth exponents at high water [m1 n1]',...
                          'Width & depth exponents at low water [m2 n2]',...
                          'Cut-off valley elevation (mAD)',...
                          'Valley width at mouth (m)',...
                          'Valley depth at mouth (mAD)',...
                          'Distance to tidal limit (m)',...
                          'Elevation of tidal limit (mAD)',...
                          'Distance to head of valley (m)',...
                          'Elevation of head of valley (mAD)'};
        %abstract properties in muiPropertyUI for tab display
        TabDisplay   %structure defines how the property table is displayed 
    end
    
    properties
        numu = 1        %width exponent at high water (-) can be two values nu and mu*
        nlml = 1        %width exponent at low water (-) can be two values nl and ml*
        ValleyEle       %cut-off valley elevation (mAD)
        ValleyWidth     %valley mouth at mouth (m) 
        ValleyDepth     %valley depth at mouth (mAD)
        xTidalLimit     %distance to tidal limit(m)
        zTidalLimit     %elevation of tidal limit (mAD)
        xValleyHead     %distance to head of valley (m)
        zValleyHead     %elevation of head of valley  (mAD)
    end    

%%   
    methods (Access=protected)
        function obj = CF_ValleyData(mobj)      
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
            classname = 'CF_ValleyData';             
            obj = getClassObj(mobj,'Inputs',classname);
            if isempty(obj)
                obj = CF_ValleyData(mobj);         
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