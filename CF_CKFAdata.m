classdef CF_CKFAdata < GDinterface        
%
%-------class help---------------------------------------------------------
% NAME
%   CF_CKFAmodel.m
% PURPOSE
%   Class for exponential form parameters for use in ChannelForm model
% USAGE
%   obj = CF_CKFAmodel.setInput(mobj); %mobj is a handle to Main UI
% SEE ALSO
%   inherits muiPropertyUI
% NOTES
%   Notation used in model, with exobj an instance of this class:
%     Le = expobj.ChannelLength;       %total length of channel (m)
%     bu = expobj.HWmouthWidth;      %width of mouth at high water(m)
%     bl = expobj.LWmouthWidth;      %width of mouth at low water level(m)
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
    properties
        %inherits Data, RunParam, MetaData and CaseIndex from muiDataSet
        %Additional properties:   
    end

%      
%     properties (Hidden)
%         %abstract properties in muiPropertyUI to define input parameters
%         PropertyLabels = {'Total length of channel (m)',...
%                           'Width of mouth at mean tide level (m)',...
%                           'Width convergence length at mean tide level (m)',...
%                           'Area of mouth at mean tide level (m^2)',...
%                           'Area convergence length at mean tide level (m^2)'};
                          
%                           'Channel shape parameter',...
%                           
%                           'Width convergence length at low water (m)',... 
%                           'Width exponent at high water (-)',...
%                           'Width exponent at low water (-)',...
%                           'Thalweg bed level at mouth (mOD)',...
%                           'Intertidal shape parameter (range:0.01-0.5)'
        %abstract properties in muiPropertyUI for tab display
%         TabDisplay   %structure defines how the property table is displayed 
%     end
    
    properties
        ChannelLength       %Le - total length of channel (m)
        MouthWidth          %bu - width of mouth at high water(m)
        WidthELength        %Lwu - width convergence length at high water (m)
        MouthCSA            %bu - CSA of mouth at high water(m2)
        AreaELength         %Lwl- area convergence length mean tide level (m2)
%         ChannelShapeParam   %nc - channel shape parameter (-)
        
%         
%         HWwidthPower        %nu - width exponent at high water (-)
%         LWwidhPoser         %nl - width exponent at low water (-)
%         zMouthInvert        %zm - thalweg bed level at mouth to zero datum (m)
%         FlatShapeParam      %ki - intertidal shape parameter[ki*100; range:0.01-0.5]
    end    
    
    methods (Access = private)
        function obj = CF_CKFAdata()                   
            %class constructor
        end
    end   
%%   
%     methods (Access=protected)
%         function obj = CF_CKFAmodel(mobj)       
%             %constructor code:            
%             %TabDisplay values defined in UI function setTabProperties used to assign
%             %the tabname and position on tab for the data to be displayed
%             obj = setTabProps(obj,mobj);  %muiPropertyUI function
%             
%             %to use non-numeric entries then one can either pre-assign 
%             %the values in the class properties defintion, above, or 
%             %specify the PropertyType as a cell array here in the class 
%             %constructor, e.g.:
%             % obj.PropertyType = [{'datetime','string','logical'},...
%             %                                       repmat({'double'},1,8)];
%         end 
%     end
%%  
%     methods (Static)  
%         function obj = setInput(mobj,editflag)
%             %gui for user to set Parameter Input values
%             classname = 'CF_CKFAmodel';               
%             obj = getClassObj(mobj,'Inputs',classname);
%             if isempty(obj)
%                 obj = CF_CKFAmodel(mobj);         
%             end
%             
%             %use muiPropertyUI function to generate UI
%             if nargin<2 || editflag
%                 %add nrec to limit length of props UI (default=12)
%                 obj = editProperties(obj);  
%                 %add any additional manipulation of the input here
%             end
%             setClassObj(mobj,'Inputs',classname,obj);
%         end     
%     end
%%        
        %add other functions to operate on properties as required   
end