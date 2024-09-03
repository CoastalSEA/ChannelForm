classdef CF_InletData < muiPropertyUI             
%
%-------class help---------------------------------------------------------
% NAME
%   CF_InletData.m
% PURPOSE
%   Class for exponential form parameters for use in ChannelForm model
% USAGE
%   obj = CF_InletData.setInput(mobj); %mobj is a handle to Main UI
% SEE ALSO
%   inherits muiPropertyUI
% NOTES
%   Notation used in model, with expobj an instance of this class:
%     bu = expobj.HWmouthWidth/2;    %half-width of mouth at high water(m)
%     bl = expobj.LWmouthWidth/2;    %half-width of mouth at low water level(m)
%     bt = expobj.HWbasinWidth/2;    %half-width of basin (m)
%     Lwl = expobj.LWwidthELength;   %width convergence length at low water (m)
%     nl = expobj.LWwidthPower;       %width exponent at low water (-)
%     zm = expobj.zMouthInvert;      %thalweg bed level at mouth to zero datum (m)
%     Li = xInletLength;             %length of mouth inlet (m)
%     nc = expobj.ChannelShapeParam; %channel shape parameter (-)
%     ki = expobj.FlatShapeParam;    %intertidal shape parameter[ki*100; range:0.01-0.5]
%
% Author: Ian Townend
% CoastalSEA (c) June 2024
%--------------------------------------------------------------------------
%      
    properties (Hidden)
        %abstract properties in muiPropertyUI to define input parameters
        PropertyLabels = {'Width of mouth at high water (m)',...
                          'Width of mouth at low water level (m)',...
                          'Width of basin at high water (m)',...
                          'Width convergence length at low water (m)',... 
                          'Width exponent at low water (-)',...
                          'Depth at mouth to MTL (m)',...
                          'Length of mouth inlet (m)',...
                          'Length of LW channel (m)',...
                          'Length of basin at high water (m)',...
                          'Channel shape parameter',...                                   
                          'Intertidal shape parameter (range:0.01-0.5)'};
        %abstract properties in muiPropertyUI for tab display
        TabDisplay   %structure defines how the property table is displayed 
    end
    
    properties
        HWmouthWidth        % width of mouth at high water(m)
        LWmouthWidth        %width of mouth at low water level(m)
        HWbasinWidth        %maximum width of basin at high water (m)        
        LWwidthELength      %width convergence length at low water (m)
        LWwidthPower        %width exponent at low water (-)
        MTmouthDepth        %depth at mouth to MTL (m)
        xInletLength        %length of mouth inlet (m)
        xLWchannel;         %distance from mouth to end of LW channel (m)
        xHWbasinLength;     %distance from mouth to tidal limit (m)
        ChannelShapeParam   %channel shape parameter (-)        
        FlatShapeParam      %intertidal shape parameter[ki*100; range:0.01-0.5] for L&M profile
    end    
    
    %Note zMouthInvert - thalweg bed level at mouth to zero datum (m) is
    %now a Dependent property in CF_FormModel
%%   
    methods (Access=protected)
        function obj = CF_InletData(mobj)       
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
            classname = 'CF_InletData';               
            obj = getClassObj(mobj,'Inputs',classname);
            if isempty(obj)
                obj = CF_InletData(mobj);         
            end
            
            %use muiPropertyUI function to generate UI
            if nargin<2 || editflag
                %add nrec to limit length of props UI (default=12)
                obj = editProperties(obj);  
                %add any additional manipulation of the input here
            end
            setClassObj(mobj,'Inputs',classname,obj);

            %set-up co-ordinate system
            grdobj = getClassObj(mobj,'Inputs','GD_GridProps');
            if ~isempty(grdobj)
                [xi,yi,delx,dely] = getGridDimensions(grdobj);
                Lt = obj.xHWbasinLength;
                Bt = obj.HWbasinWidth/2;
                Xmax = xi(end)-2*delx;
                Ymax = yi(end)-2*dely;
                if Lt>Xmax
                    warndlg('Length of basin at high water is longer than grid domain')
                elseif Bt>Ymax
                    warndlg('Width of basin at high water is wider than grid domain')
                end
            end
        end     
    end      
%%
        %add other functions to operate on properties as required   
end