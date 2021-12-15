classdef CF_PowerData < muiPropertyUI            
%
%-------class help---------------------------------------------------------
% NAME
%   CF_PowerData.m
% PURPOSE
%   Class for power form parameters for use in ChannelForm model
% USAGE
%   obj = CF_PowerData.setInput(mobj); %mobj is a handle to Main UI
% SEE ALSO
%   inherits muiPropertyUI
% NOTES
%   Notation used in model, with exobj an instance of this class:
%     Le = pwrobj.ChannelLength;     %total length of channel (m)
%     bu = pwrobj.HWmouthWidth/2;    %half-width of mouth at high water(m)
%     bl = pwrobj.LWmouthWidth/2;    %half-width of mouth at low water level(m)
%     nu = pwrobj.HWwidthExponent;   %width exponent at high water (-)
%     nl = pwrobj.LWwidthExponent;   %width exponent at low water (-)
%     mu = pwrobj.HWdepthExponent;   %depth exponent at high water (-)
%     ml = pwrobj.LWdepthExponent;   %depth exponent at low water (-)
%     zm = pwrobj.zMouthInvert;      %thalweg bed level at mouth to zero datum (m)
%
% Author: Ian Townend
% CoastalSEA (c) Jan 2022
%--------------------------------------------------------------------------
%      
    properties (Hidden)
        %abstract properties in muiPropertyUI to define input parameters
        PropertyLabels = {'Width of mouth at high water (m)',...
                          'Width of mouth at low water level (m)',...
                          'Width exponent at high water (-)',...
                          'Width exponent at low water (-)',...  
                          'Depth exponent at high water (-)',...
                          'Depth exponent at low water (-)',...
                          'Thalweg bed level at mouth (mOD)'};
        %abstract properties in muiPropertyUI for tab display
        TabDisplay   %structure defines how the property table is displayed 
        ChannelShapeParam = 2 %used to define river profile
    end
    
    properties           
        HWmouthWidth        %bu - width of mouth at high water(m)
        LWmouthWidth        %bl - width of mouth at low water level(m)
        HWwidthExponent     %nu - width exponent at high water (-)
        LWwidthExponent     %nl - width exponent at low water (-)
        HWdepthExponent     %mu - depth exponent at high water (-)
        LWdepthExponent     %ml - depth exponent at low water (-)
        zMouthInvert        %zm - thalweg bed level at mouth to zero datum (m)]
    end    

%%   
    methods (Access=protected)
        function obj = CF_PowerData(mobj)            
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
            classname = 'CF_PowerData';        
            obj = getClassObj(mobj,'Inputs',classname);
            if isempty(obj)
                obj = CF_PowerData(mobj);  
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