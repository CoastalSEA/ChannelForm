classdef ChannelForm < muiModelUI                    
%
%-------class help---------------------------------------------------------
% NAME
%   ChannelForm.m
% PURPOSE
%   Main GUI for a generic model interface, which implements the 
%   muiModelUI abstract class to define main menus.
% SEE ALSO
%   Abstract class muiModelUI.m and tools provided in muitoolbox
%
% Author: Ian Townend
% CoastalSEA (c) Jan 2021
%--------------------------------------------------------------------------
% 
    properties  (Access = protected)
        %implement properties defined as Abstract in muiModelUI
        vNumber = '2.0'
        vDate   = 'Jan 2022'
        modelName = 'ChannelForm'                    
        %Properties defined in muiModelUI that need to be defined in setGui
        % ModelInputs  %classes required by model: used in isValidModel check 
        % DataUItabs   %struct to define type of muiDataUI tabs for each use                         
    end
    
    methods (Static)
        function obj = ChannelForm                 
            %constructor function initialises GUI
            obj = setMUI(obj);             
        end
    end
%% ------------------------------------------------------------------------
% Definition of GUI Settings
%--------------------------------------------------------------------------  
    methods (Access = protected)
        function obj = setMUI(obj)
            %initialise standard figure and menus    
            %classes required to run model, format:
            %obj.ModelInputs.<model classname> = {'Param_class1',Param_class2',etc}
            defaultprops = {'CF_HydroData','CF_SediData','GD_GridProps'};
            obj.ModelInputs.CF_FormModel = defaultprops;
            obj.ModelInputs.CF_ValleyModel = [defaultprops,{'CF_ValleyData'}];
            obj.ModelInputs.CF_HydroData = {'CF_HydroData','CF_SediData'};
            obj.ModelInputs.CF_TransModel = defaultprops;
            %tabs to include in DataUIs for plotting and statistical analysis
            %select which of the options are needed and delete the rest
            %Plot options: '2D','3D','4D','2DT','3DT','4DT'
            obj.DataUItabs.Plot = {'2D','3D','4D','2DT','3DT','4DT'};  
            %Statistics options: 'General','Timeseries','Taylor','Intervals'
            obj.DataUItabs.Stats = {'General','Timeseries','Taylor','Intervals'};  
            
            modelLogo = 'ChannelForm_logo.jpg';  %default splash figure - edit to alternative
            initialiseUI(obj,modelLogo); %initialise menus and tabs   
            
            
        end    
        
%% ------------------------------------------------------------------------
% Definition of Menu Settings
%--------------------------------------------------------------------------
        function menu = setMenus(obj)
            %define top level menu items and any submenus
            %MenuLabels can any text but should avoid these case-sensitive 
            %reserved words: "default", "remove", and "factory". If label 
            %is not a valid Matlab field name this the struct entry
            %is modified to a valid name (eg removes space if two words).
            %The 'gcbo:' Callback text triggers an additional level in the 
            %menu. Main menu labels are defined in sequential order and 
            %submenus in order following each brach to the lowest level 
            %before defining the next branch.         
                                    
            MenuLabels = {'File','Tools','Project','Setup','Utilities',...
                                                'Run','Analysis','Help'};
            menu = menuStruct(obj,MenuLabels);  %create empty menu struct
            %
            %% File menu --------------------------------------------------
             %list as per muiModelUI.fileMenuOptions
            menu.File.List = {'New','Open','Save','Save as','Exit'};
            menu.File.Callback = repmat({@obj.fileMenuOptions},[1,5]);
            
            %% Tools menu -------------------------------------------------
            %list as per muiModelUI.toolsMenuOptions
            menu.Tools(1).List = {'Refresh','Clear all'};
            menu.Tools(1).Callback = {@obj.refresh, 'gcbo;'};  
            
            % submenu for 'Clear all'
            menu.Tools(2).List = {'Model','Figures','Cases'};
            menu.Tools(2).Callback = repmat({@obj.toolsMenuOptions},[1,3]);

            %% Project menu -----------------------------------------------
            menu.Project(1).List = {'Project Info','Cases','Export/Import'};
            menu.Project(1).Callback = {@obj.editProjectInfo,'gcbo;','gcbo;'};
            
            %list as per muiModelUI.projectMenuOptions
            % submenu for Scenarios
            menu.Project(2).List = {'Edit Description','Edit Data Set',...
                                    'Save Data Set','Delete Case','Reload Case',...
                                    'View Case Settings'};                                               
            menu.Project(2).Callback = repmat({@obj.projectMenuOptions},[1,6]);
            
            % submenu for 'Export/Import'                                          
            menu.Project(3).List = {'Export Case','Import Case'};
            menu.Project(3).Callback = repmat({@obj.projectMenuOptions},[1,2]);
            
            %% Setup menu -------------------------------------------------
            menu.Setup(1).List = {'Form Parameters','System Parameters',...
                                  'Run Parameters','Import Grid',...
                                  'Grid Tools','Add Properties',...
                                  'Delete Properties','Model Constants'}; 
            menu.Setup(1).Callback = [repmat({'gcbo;'},[1,5]),...
                                      repmat({@obj.gridMenuOptions},[1,2]),...
                                      {@obj.setupMenuOptions}];                  
            %add separators to menu list (optional - default is off)
            menu.Setup(1).Separator = [repmat({'off'},[1,3]),{'on'},...
                                       repmat({'off'},[1,3]),{'on'}]; %separator preceeds item
            menu.Setup(2).List = {'Exp Form Parameters','Power Form Parameters',...
                                  'CKFA Form Parameters','Valley Parameters'};
            menu.Setup(2).Callback = repmat({@obj.setupMenuOptions},[1,4]); 
            
            menu.Setup(3).List = {'Hydraulic Parameters','Sediment Parameters',...
                                  'Transgression Parameters','Morphological Modifications'};
            menu.Setup(3).Callback = [{'gcbo;'},repmat({@obj.setupMenuOptions},[1,3])];  
            
            menu.Setup(4).List = {'Tidal Forcing','Hydraulic Model'};
            menu.Setup(4).Callback = repmat({@obj.setupMenuOptions},[1,2]); 
            
            menu.Setup(5).List = {'Grid Parameters','Run Time Parameters'};           
            menu.Setup(5).Callback = repmat({@obj.setupMenuOptions},[1,2]);
            
            menu.Setup(6).List = {'Load','Add','Delete'};
            menu.Setup(6).Callback = repmat({@obj.loadMenuOptions},[1,3]);
            
            menu.Setup(7).List = {'Translate Grid','Rotate Grid',...
                                  'Re-Grid','Sub-Grid',...
                                  'Combine Grids','Add Surface','Export xyz Grid'};                                                                        
            menu.Setup(7).Callback = repmat({@obj.gridMenuOptions},[1,7]);
            menu.Setup(7).Separator = [repmat({'off'},[1,6]),{'on'}]; %separator preceeds item
            %% Utilities menu ---------------------------------------------------
            menu.Utilities(1).List = {'Hydraulic Model',...
                                      'Add Form to Valley',...
                                      'Add Modifications',...
                                      'River Dimensions',...
                                      'Valley Thalweg',...
                                      'Area of Flood Plain',...
                                      'CKFA Channel Dimensions',...
                                      'Morphological Timescale',...
                                      'Channel-Valley Sub-Plot'};
            menu.Utilities(1).Callback = repmat({@obj.utilsMenuOptions},[1,9]);
             menu.Utilities(1).Separator = [repmat({'off'},[1,3]),{'on'},...
                                            repmat({'off'},[1,5])]; %separator preceeds item
            
%             %% Process menu ---------------------------------------------------
%             menu.Process(1).List = {'Add Form to Valley','Add Modifications',...
%                                     'Restore Form'};
%             menu.Process(1).Callback = repmat({@obj.procMenuOptions},[1,3]);
%             
            %% Run menu ---------------------------------------------------
            menu.Run(1).List = {'Exponential form model','Power form model',...
                                'CKFA form model','Valley form model',...
                                'Hydro-Form model','Transgression model',...
                                'Derive Output'};
            menu.Run(1).Callback = repmat({@obj.runMenuOptions},[1,7]);
            menu.Run(1).Separator = [repmat({'off'},[1,6]),{'on'}]; %separator preceeds item
            
            %% Plot menu --------------------------------------------------  
            menu.Analysis(1).List = {'Plots','Statistics'};
            menu.Analysis(1).Callback = repmat({@obj.analysisMenuOptions},[1,2]);
            
            %% Help menu --------------------------------------------------
            menu.Help(1).Callback = {@obj.Help}; %make model specific?
            
        end
        
%% ------------------------------------------------------------------------
% Definition of Tab Settings
%--------------------------------------------------------------------------
        function [tabs,subtabs] = setTabs(obj)
            %define main tabs and any subtabs required. struct field is 
            %used to set the uitab Tag (prefixed with sub for subtabs). 
            %Order of assignment to struct determines order of tabs in figure.
            %format for tabs: 
            %    tabs.<tagname> = {<tab label>,<callback function>};
            %format for subtabs: 
            %    subtabs.<tagname>(i,:) = {<subtab label>,<callback function>};
            %where <tagname> is the struct fieldname for the top level tab.
            tabs.Cases  = {'   Cases  ',@obj.refresh};   
            
            tabs.Form = {'  Form  ',''};
            subtabs.Form(1,:) = {' Exponential ',@obj.InputTabSummary};
            subtabs.Form(2,:) = {'  Power  ',@obj.InputTabSummary};
            subtabs.Form(3,:) = {'  Valley  ',@obj.InputTabSummary};
            tabs.Settings = {'  Settings  ',''};
            subtabs.Settings(1,:) = {' Forcing ',@obj.InputTabSummary}; 
            subtabs.Settings(2,:) = {' Sediments ',@obj.InputTabSummary};
            subtabs.Settings(3,:) = {' Transgression ',@obj.InputTabSummary};
            subtabs.Settings(4,:) = {' Modifications ',@obj.InputTabSummary};
            subtabs.Settings(5,:) = {' Run Parameters ',@obj.InputTabSummary};
            tabs.HydroProps = {' Hydro-Props ',''};
            subtabs.HydroProps(1,:) = {' Water Levels ',@obj.setCFMtabs};
            subtabs.HydroProps(2,:)   = {' Hydraulics ',@obj.setCFMtabs};
            tabs.FormProps = {' Form-Props ',@obj.getTabData};
            tabs.Plot   = {'  Q-Plot  ',@obj.getTabData};
            tabs.Stats = {'   Stats   ',@obj.setTabAction};
        end
       
%%
        function props = setTabProperties(~)
            %define the tab and position to display class data tables
            %props format: {class name, tab tag name, position, ...
            %               column width, table title}
            % position and column widths vary with number of parameters
            % (rows) and width of input text and values. Inidcative
            % positions:  top left [0.95,0.48];    top right [0.95,0.97]
            %         bottom left [0.45, 0.48]; bottom right [0.45,0.97]                                              
            props = {...                                    
                'CF_ExpData','Exponential',[0.90,0.60],{220,70}, 'Exponential model parameters:';...  
                'CF_PowerData','Power',[0.90,0.60],{220,70}, 'Power model parameters:';...
                'CF_ValleyData','Valley',[0.90,0.60],{220,70}, 'Valley model parameters:';... 
                'CF_ModsData','Modifications',[0.90,0.96],{240,260},'Morphological modifications:';...
                'WaterLevels','Forcing',[0.90,0.48],{160,80},'Hydraulic forcing parameters:';...
                'CF_HydroData','Forcing',[0.90,0.96],{160,80},'Hydraulic model parameters:';...
                'CF_SediData','Sediments',[0.90,0.50],{180,70}, 'Sediment parameters:';... 
                'RunProperties','Run Parameters',[0.40,0.50],{170,80},'Run time parameters:';...
                'GD_GridProps','Run Parameters',[0.90,0.50],{160,90}, 'Grid parameters:';...
                'CF_TransData','Transgression',[0.90,0.50],{180,70}, 'Transgression parameters:'}; 
        end    
 %%
        function setTabAction(obj,src,cobj)
            %function required by muiModelUI and sets action for selected
            %tab (src)
            switch src.Tag                                   
                case 'Plot' 
                     tabPlot(cobj,src);
                case 'Stats'
                    lobj = getClassObj(obj,'mUI','Stats');
                    if isempty(lobj), return; end
                    tabStats(lobj,src);           
                case 'FormProps'   
                    cf_model_tabs(cobj,src); %the function cf_form_tabs selects Properties
%                 case 'Transgression'
%                     channel_transgression(gobj,inp);
                case 'Export Grid'
                    readwrite_grid(gobj,inp);      
            end
        end     
%%
        function setCFMtabs(obj,src,~)
            %update the Settings tabs 
            switch src.Tag
%                 case 'Saltmarsh'
%                     InputTabSummary(obj,src,evt)
%                     msgtxt = 'Saltmarsh parameters have not been defined';
%                     cobj = getClassObj(obj,'Inputs','Saltmarsh',msgtxt);                    
                case 'Water Levels'
                    msgtxt = 'Water level parameters have not been defined';
                    cobj = getClassObj(obj,'Inputs','WaterLevels',msgtxt);
                case 'Hydraulics'
                    msgtxt = 'CSTmodel parameters have not been defined';
                    cobj = getClassObj(obj,'Inputs','CF_HydroData',msgtxt);
            end
            %
            if isempty(cobj), return; end
            tabPlot(cobj,src,obj);
        end  
%% ------------------------------------------------------------------------
% Callback functions used by menus and tabs
%-------------------------------------------------------------------------- 
        %% File menu ------------------------------------------------------
        %use default menu functions defined in muiModelUI
            
        %% Tools menu -----------------------------------------------------
        %use default menu functions defined in muiModelUI
                
        %% Project menu ---------------------------------------------------
        %use default menu functions defined in muiModelUI           

        %% Setup menu -----------------------------------------------------
        function setupMenuOptions(obj,src,~)
            %callback functions for data input            
            switch src.Text
                case 'Exp Form Parameters'
                    CF_ExpData.setInput(obj);  
                    %update tab display with input data
                    tabsrc = findobj(obj.mUI.Tabs,'Tag','Exponential');
                    InputTabSummary(obj,tabsrc);
                case 'Power Form Parameters'
                    CF_PowerData.setInput(obj);  
                    %update tab display with input data
                    tabsrc = findobj(obj.mUI.Tabs,'Tag','Power');
                    InputTabSummary(obj,tabsrc);
                case 'CKFA Form Parameters'
                    CF_CKFAdata.setInput(obj);  
                    %update tab display with input data
                    tabsrc = findobj(obj.mUI.Tabs,'Tag','CKFA');
                    InputTabSummary(obj,tabsrc);
                case 'Valley Parameters'    
                    CF_ValleyData.setInput(obj);  
                    %update tab display with input data
                    tabsrc = findobj(obj.mUI.Tabs,'Tag','Valley');
                    InputTabSummary(obj,tabsrc);
                case 'Tidal Forcing'
                    WaterLevels.setInput(obj);                    
                    %update tab display with input data
                    tabsrc = findobj(obj.mUI.Tabs,'Tag','Forcing');
                    InputTabSummary(obj,tabsrc);
                case 'Hydraulic Model'
                    CF_HydroData.setInput(obj);                      
                    %update tab display with input data
                    tabsrc = findobj(obj.mUI.Tabs,'Tag','Forcing');
                    InputTabSummary(obj,tabsrc);   
                case 'Sediment Parameters'   
                    CF_SediData.setInput(obj);  
                    %update tab display with input data
                    tabsrc = findobj(obj.mUI.Tabs,'Tag','Sediments');
                    InputTabSummary(obj,tabsrc);
                case 'Transgression Parameters'
                    CF_TransData.setInput(obj);  
                    %update tab display with input data
                    tabsrc = findobj(obj.mUI.Tabs,'Tag','Run Parameters');
                    InputTabSummary(obj,tabsrc);
                case 'Morphological Modifications'
                    CF_ModsData.setInput(obj);  
                    %update tab display with input data
                    tabsrc = findobj(obj.mUI.Tabs,'Tag','Run Parameters');
                    InputTabSummary(obj,tabsrc);
                case 'Grid Parameters'
                    GD_GridProps.setInput(obj);  
                    %update tab display with input data
                    tabsrc = findobj(obj.mUI.Tabs,'Tag','Run Parameters');
                    InputTabSummary(obj,tabsrc);
                case 'Run Time Parameters'
                    RunProperties.setInput(obj);  
                    %update tab display with input data
                    tabsrc = findobj(obj.mUI.Tabs,'Tag','Run Parameters');
                    InputTabSummary(obj,tabsrc);
                case 'Model Constants'
                    obj.Constants = setInput(obj.Constants);
            end
        end  
%%
        function loadMenuOptions(obj,src,~)
            %callback functions to import data
            classname = 'GD_ImportData';
            switch src.Text
                case 'Load'
                    fname = sprintf('%s.loadData',classname);
                    callStaticFunction(obj,classname,fname); 
                case 'Add'
                    useCase(obj.Cases,'single',{classname},'addData');
                case 'Delete'
                    useCase(obj.Cases,'single',{classname},'deleteGrid');
            end
            DrawMap(obj);
        end   

        %% Utilities menu -------------------------------------------------------
        function utilsMenuOptions(obj,src,~)
            %callback functions to run model
            switch src.Text      
                case 'Hydraulic Model'
                    msgtxt = ('Hydraulic properties have not been defined');
                    cobj = getClassObj(obj,'Inputs','CF_HydroData',msgtxt);
                    if isempty(cobj), return; end
                    runModel(cobj,obj);
                case 'Add Form to Valley'
                    CF_ValleyModel.addForm2Valley(obj);
                case 'Restore Form'
                    CF_ValleyModel.restoreForm(obj);     
                case 'River Dimensions'
                    CF_HydroData.displayRiverDims(obj);
                case 'Valley Thalweg'
                    CF_ValleyModel.checkValleyThalweg(obj);
                case 'Area of Flood Plain'
                    CF_ValleyModel.checkFloodPlainArea(obj);
                case 'CKFA Channel Dimensions'
                    ckfa_dimensions(obj);
                case 'Morphological Timescale'
                    CF_SediData.displayMorphTime(obj);
                case 'Channel-Valley Sub-Plot'
                    CF_ValleyModel.componentsPlot(obj);
            end            
        end  


%         %% Process menu -------------------------------------------------------
%         function procMenuOptions(obj,src,~)
%             %callback functions to run model
%             switch src.Text                   
%                 case 'Add Form to Valley'
%                     CF_ValleyModel.addForm2Valley(obj);
%                 case 'Restore Form'
%                     CF_ValleyModel.restoreForm(obj); 
%                 case 'Add Modifications'
%                     CF_FormModel.addMorphMods(obj);
%                 
%             end            
%         end  
        %%
        function gridMenuOptions(obj,src,~)
            %callback functions for grid tools options
            gridclasses = {'CF_FormModel','CF_ValleyModel','GD_ImportData'};
            CF_FormModel.gridMenuOptions(obj,src,gridclasses);
            DrawMap(obj);
        end

        %% Run menu -------------------------------------------------------
        function runMenuOptions(obj,src,~)
            %callback functions to run model
            switch src.Text                   
                case 'Exponential form model'
                    CF_FormModel.runModel(obj,'Exponential');
                case 'Power form model'
                    CF_FormModel.runModel(obj,'Power');
                case 'CKFA form model'
                    CF_FormModel.runModel(obj,'CKFA');
                case 'Valley form model'
                    CF_ValleyModel.runModel(obj);
                case 'Hydro-Form model'
                case 'Transgression model'
                case 'Derive Output'
                    obj.mUI.Manip = muiManipUI.getManipUI(obj);
            end            
        end               
            
        %% Analysis menu ------------------------------------------------------
        function analysisMenuOptions(obj,src,~)
            switch src.Text
                case 'Plots'
                    obj.mUI.PlotsUI = muiPlotsUI.getPlotsUI(obj);
                case 'Statistics'
                    obj.mUI.StatsUI = muiStatsUI.getStatsUI(obj);
            end            
        end

        %% Help menu ------------------------------------------------------
        function Help(~,~,~)
            docsearch ChannelForm 
        end
%% ------------------------------------------------------------------------
% Overload muiModelUI.MapTable to customise Tab display of records (if required)
%--------------------------------------------------------------------------     
%         function MapTable(obj,ht)
%             %create tables for Record display tabs - called by DrawMap
%             % ht - tab handle
%         end
    end
end    
    
    
    
    
    
    
    
    
    
    