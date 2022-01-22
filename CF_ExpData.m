classdef CF_ExpData < muiPropertyUI             
%
%-------class help---------------------------------------------------------
% NAME
%   CF_ExpData.m
% PURPOSE
%   Class for exponential form parameters for use in ChannelForm model
% USAGE
%   obj = CF_ExpData.setInput(mobj); %mobj is a handle to Main UI
% SEE ALSO
%   inherits muiPropertyUI
% NOTES
%   Notation used in model, with expobj an instance of this class:
%     Le = expobj.ChannelLength;     %total length of channel (m)
%     bu = expobj.HWmouthWidth/2;    %half-width of mouth at high water(m)
%     bl = expobj.LWmouthWidth/2;    %half-width of mouth at low water level(m)
%     nc = expobj.ChannelShapeParam; %channel shape parameter (-)
%     Lwu = expobj.HWwidthELength;   %width convergence length at high water (m)
%     Lwl = expobj.LWwidthELength;   %width convergence length at low water (m)
%     nu = expobj.HWwidthPower;      %width exponent at high water (-)
%     nl = expobj.LWwidhPoser;       %width exponent at low water (-)
%     zm = expobj.zMouthInvert;      %thalweg bed level at mouth to zero datum (m)
%     ki = expobj.FlatShapeParam;    %intertidal shape parameter[ki*100; range:0.01-0.5]
%
% Author: Ian Townend
% CoastalSEA (c) Jan 2022
%--------------------------------------------------------------------------
%      
    properties (Hidden)
        %abstract properties in muiPropertyUI to define input parameters
        PropertyLabels = {'Width of mouth at high water (m)',...
                          'Width of mouth at low water level (m)',...
                          'Channel shape parameter',...
                          'Width convergence length at high water (m)',...
                          'Width convergence length at low water (m)',... 
                          'Width exponent at high water (-)',...
                          'Width exponent at low water (-)',...
                          'Thalweg bed level at mouth (mOD)',...
                          'Intertidal shape parameter (range:0.01-0.5)'};
        %abstract properties in muiPropertyUI for tab display
        TabDisplay   %structure defines how the property table is displayed 
    end
    
    properties
        HWmouthWidth        %bu - width of mouth at high water(m)
        LWmouthWidth        %bl - width of mouth at low water level(m)
        ChannelShapeParam   %nc - channel shape parameter (-)
        HWwidthELength      %Lwu - width convergence length at high water (m)
        LWwidthELength      %Lwl- width convergence length at low water (m)
        HWwidthPower        %nu - width exponent at high water (-)
        LWwidthPower        %nl - width exponent at low water (-)
        zMouthInvert        %zm - thalweg bed level at mouth to zero datum (m)
        FlatShapeParam      %ki - intertidal shape parameter[ki*100; range:0.01-0.5]
    end    

%%   
    methods (Access=protected)
        function obj = CF_ExpData(mobj)       
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
            classname = 'CF_ExpData';               
            obj = getClassObj(mobj,'Inputs',classname);
            if isempty(obj)
                obj = CF_ExpData(mobj);         
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