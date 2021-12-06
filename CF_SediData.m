classdef CF_SediData < muiPropertyUI             
%
%-------class help---------------------------------------------------------
% NAME
%   CF_SediData.m
% PURPOSE
%   Class for sediment parameters for use in ChannelForm model
% USAGE
%   obj = CF_SediData.setInput(mobj); %mobj is a handle to Main UI
% SEE ALSO
%   inherits muiPropertyUI
%
% Author: Ian Townend
% CoastalSEA (c) Jan 2022
%--------------------------------------------------------------------------
%      
    properties (Hidden)
        %abstract properties in muiPropertyUI to define input parameters
        PropertyLabels = {'Sediment grain size, d50 (m)',...
                          'Critical bed shear stress (Pa)',...
                          'Equilibrium sediment density (kg/m^3)', ...
                          'Sediment load in river (kg/m^3)',...
                          'Bed density (kg/m^3)',...
                          'Transport coefficient (+/-n)',...
                          'Average depth over marsh (m)',...
                          'Maximum marsh depth (m), 0=no marsh',...
                          'Sediment grain size in river, d50 (m)',...
                          'Critical bed shear stress in river (Pa)'};
        %abstract properties in muiPropertyUI for tab display
        TabDisplay   %structure defines how the property table is displayed 
    end
    
    properties
        SedimentSize         %sediment grain size, D50 (m)
        CritBedShear         %critical bed shear stress (Pa)
        EqDensity = 0        %equilibrium concentration density(kg/m^3)
        RiverDensity = 0     %river load imported by advection (kg/m^3)
        BedDensity = 0       %bed density (kg/m^3)
        TransportCoeff = 0   %transport coefficient n (3-5) 
        AvMarshDepth = 0     %average depth of marsh surface
        MaxMarshDepth = 0    %maximum depth of saltmarsh   
        d50river             %sediment grain size in river, D50 (m)
        tauriver             %critical bed shear stress in river (Pa)
    end    
    properties (Dependent)
        EqConcentration      %equilibrium concentration (-)
        BedConcentration     %concentration of bed (-)
    end
%%   
    methods (Access=protected)
        function obj = CF_SediData(mobj)          
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
            classname = 'CF_SediData';             
            obj = getClassObj(mobj,'Inputs',classname);
            if isempty(obj)
                obj = CF_SediData(mobj);      
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
        function eqConc = get.EqConcentration(obj)
            %dependent property derived from EqRhoCoarse
            cn = muiConstants.Evoke;
            eqConc = obj.EqDensity/cn.SedimentDensity;
        end
%%
        function eqConc = get.BedConcentration(obj)
            %dependent property derived from EqRhoCoarse
            cn = muiConstants.Evoke;
            eqConc = obj.BedDensity/cn.SedimentDensity;
        end
    end
end