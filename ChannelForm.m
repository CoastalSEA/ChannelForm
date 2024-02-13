classdef ChannelForm < muiModelUI                    
%
%-------class help---------------------------------------------------------
% NAME
%   ChannelForm.m
% PURPOSE
%   Main UI for the ChannelForm interface, which implements the 
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
        vNumber = '3.3'
        vDate   = 'Jan 2024'
        modelName = 'ChannelForm'                    
        %Properties defined in muiModelUI that need to be defined in setGui
        % ModelInputs  %classes required by model: used in isValidModel check 
        % DataUItabs   %struct to define type of muiDataUI tabs for each use                         
    end
    
    methods (Static)
        function obj = ChannelForm                 
            %constructor function initialises GUI
            isok = check4muitoolbox(obj);
            if ~isok, return; end
            %
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
            obj.ModelInputs.CF_FormModel = [defaultprops,{'WaterLevels'}];
            obj.ModelInputs.CF_ValleyModel = [defaultprops,{'CF_ValleyData'}];
            obj.ModelInputs.CF_HydroData = {'CF_HydroData','CF_SediData'};
            obj.ModelInputs.CF_TransModel = [defaultprops,...
                           {'CF_TransData','RunProperties','WaterLevels'}];
            %tabs to include in DataUIs for plotting and statistical analysis
            %select which of the options are needed and delete the rest
            %Plot options: '2D','3D','4D','2DT','3DT','4DT'
            obj.DataUItabs.Plot = {'2D','3D','2DT','3DT'};  
            %Statistics options: 'General','Timeseries','Taylor','Intervals'
            obj.DataUItabs.Stats = {'General','Timeseries','Taylor'};  
            
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
                                  'Delete Properties','Edit Inlet Definition',...
                                  'Model Constants'}; 
            menu.Setup(1).Callback = [repmat({'gcbo;'},[1,5]),...
                                      repmat({@obj.formMenuOptions},[1,3]),...
                                      {@obj.loadMenuOptions}];                  
            %add separators to menu list (optional - default is off)
            menu.Setup(1).Separator = [repmat({'off'},[1,3]),{'on'},...
                                       repmat({'off'},[1,4]),{'on'}]; %separator preceeds item
            menu.Setup(2).List = {'Exp Form Parameters','Power Form Parameters',...
                                        'Valley Parameters','Shore Parameters'};
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
                                  'Combine Grids','Add Surface',...
                                  'To curvilinear','From curvilinear',... 
                                  'Display Dimensions','Difference Plot',...
                                  'Plot Sections','Digitise Line',...
                                  'Export xyz Grid','User Function'};                                                                          
            menu.Setup(7).Callback = repmat({@obj.gridMenuOptions},[1,14]);
            menu.Setup(7).Separator = [repmat({'off'},[1,6]),...
                             {'on','off','on','off','off','on','on','on'}];%separator preceeds item
            %% Utilities menu ---------------------------------------------------
            menu.Utilities(1).List = {'Hydraulic Model',...
                                      'Add Form to Valley',...
                                      'Add Meander',...
                                      'Add Shoreline',...
                                      'Add Modifications',...
                                      'Add Thalweg to Valley',...
                                      'River Dimensions',...
                                      'Valley Thalweg',...
                                      'Area of Flood Plain',...
                                      'CKFA Channel Dimensions',...
                                      'Morphological Timescale',...
                                      'Channel-Valley Sub-Plot',...
                                      'Centre-line Plot'};
            menu.Utilities(1).Callback = [repmat({@obj.utilsMenuOptions},[1,3]),...
                                {'gcbo;'},repmat({@obj.utilsMenuOptions},[1,9])];
            menu.Utilities(1).Separator = [repmat({'off'},[1,6]),{'on'},...
                                            repmat({'off'},[1,6])]; %separator preceeds item
            
            menu.Utilities(2).List = {'Model Shore','Extrapolate Shore'};
            menu.Utilities(2).Callback = repmat({@obj.utilsMenuOptions},[1,2]);

            %% Run menu ---------------------------------------------------
            menu.Run(1).List = {'Exponential form model','Power form model',...
                                'CKFA form model','Valley form model',...
                                'Transgression model','Derive Output'};                                
            menu.Run(1).Callback = repmat({@obj.runMenuOptions},[1,6]);
            menu.Run(1).Separator = [repmat({'off'},[1,5]),{'on'}]; %separator preceeds item
            
            %% Plot menu --------------------------------------------------  
            menu.Analysis(1).List = {'Plots','Statistics','Transgression Plots'};
            menu.Analysis(1).Callback = [repmat({@obj.analysisMenuOptions},[1,2]),{'gcbo;'}];
            menu.Analysis(1).Separator = {'off','off','on'};
            
            menu.Analysis(2).List = {'Grid Plot','Change Plot','Section Plot',...
                                    'Transgression Plot','Animation','Export Tables'};
            menu.Analysis(2).Callback = repmat({@obj.analysisPlotOptions},[1,6]);
            menu.Analysis(2).Separator = [repmat({'off'},[1,5]),{'on'}]; 
            
            %% Help menu --------------------------------------------------
            menu.Help.List = {'Documentation','Manual'};
            menu.Help.Callback = repmat({@obj.Help},[1,2]);
            
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
            subtabs.Form(4,:) = {'  Shore  ',@obj.InputTabSummary};
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
                'CF_ShoreData','Shore',[0.90,0.60],{220,70}, 'Shore model parameters:';... 
                'CF_ModsData','Modifications',[0.90,0.96],{240,260},'Morphological modifications:';...
                'WaterLevels','Forcing',[0.90,0.48],{160,80},'Hydraulic forcing parameters:';...
                'CF_HydroData','Forcing',[0.90,0.96],{160,80},'Hydraulic model parameters:';...
                'CF_SediData','Sediments',[0.90,0.60],{220,70}, 'Sediment parameters:';... 
                'RunProperties','Run Parameters',[0.40,0.50],{170,80},'Run time parameters:';...
                'GD_GridProps','Run Parameters',[0.90,0.50],{160,90}, 'Grid parameters:';...
                'CF_TransData','Transgression',[0.90,0.70],{220,120}, 'Transgression parameters:'}; 
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
                    tabtxt = 'Exponential';
                case 'Power Form Parameters'
                    CF_PowerData.setInput(obj);  
                    tabtxt = 'Power';
                case 'Valley Parameters'    
                    CF_ValleyData.setInput(obj); 
                    tabtxt = 'Valley';
                case 'Shore Parameters'
                    CF_ShoreData.setInput(obj);  
                    tabtxt = 'Shore';
                case 'Tidal Forcing'
                    WaterLevels.setInput(obj);    
                    tabtxt = 'Forcing';
                case 'Hydraulic Model'
                    CF_HydroData.setInput(obj);      
                    tabtxt = 'Forcing'; 
                case 'Sediment Parameters'   
                    CF_SediData.setInput(obj);  
                    tabtxt = 'Sediments';
                case 'Transgression Parameters'
                    CF_TransData.setInput(obj);  
                    tabtxt = 'Transgression';
                case 'Morphological Modifications'
                    CF_ModsData.setInput(obj);  
                    tabtxt = 'Modifications';
                case 'Grid Parameters'
                    GD_GridProps.setInput(obj);  
                    tabtxt = 'Run Parameters';
                case 'Run Time Parameters'
                    RunProperties.setInput(obj);  
                    tabtxt = 'Run Parameters';
            end
            %update tab display with input data
            tabsrc = findobj(obj.mUI.Tabs,'Tag',tabtxt);            
            InputTabSummary(obj,tabsrc);            
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
                case 'Model Constants'
                    obj.Constants = setInput(obj.Constants);
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
                    DrawMap(obj);
                case 'Add Meander'
                    gridclasses = {'CF_FormModel','CF_TransModel'};
                    gd.Text = 'To curvilinear';
                    CF_FormModel.gridMenuOptions(obj,gd,gridclasses);
                    DrawMap(obj);
                case 'Model Shore'
                    CF_FormModel.addShoreline(obj);
                    DrawMap(obj);
                case 'Extrapolate Shore'
                    formMenuOptions(obj,src,[])
                case 'Add Modifications'
                    CF_FormModel.addMorphMods(obj);   
                case 'Add Thalweg to Valley'
                    formMenuOptions(obj,src,[]);
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
                case 'Centre-line Plot'
                    cf_plot_centreline();
            end            
        end   
        %%
        function gridMenuOptions(obj,src,~)
            %callback functions for grid tools options
            gridclasses = {'CF_FormModel','CF_ValleyModel',...
                                          'CF_TransModel','GD_ImportData'};
            %CF_FromModel inherits GDinterface, which includes the grid tools
            GDinterface.gridMenuOptions(obj,src,gridclasses);
            DrawMap(obj);
        end
        %%
        function formMenuOptions(obj,src,~)
            %callback functions for grid tools options
            gridclasses = {'CF_FormModel','CF_ValleyModel',...
                                          'CF_TransModel','GD_ImportData'};
            %CF_FromModel inherits GDinterface, which includes the grid tools
            FGDinterface.formMenuOptions(obj,src,gridclasses);
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
                case 'Transgression model'
                    CF_TransModel.runModel(obj);
                case 'Derive Output'
                    obj.mUI.ManipUI = muiManipUI.getManipUI(obj);
            end            
        end               
            
        %% Analysis menu ------------------------------------------------------
        function analysisMenuOptions(obj,src,~)
            %callback functions for plot and stats menus
            switch src.Text
                case 'Plots'
                    obj.mUI.PlotsUI = muiPlotsUI.getPlotsUI(obj);
                case 'Statistics'
                    obj.mUI.StatsUI = muiStatsUI.getStatsUI(obj);
            end            
        end
%%
        function analysisPlotOptions(obj,src,~)
            %callback functions for summary transgression plots
            if strcmp(src.Text,'Export Tables')
                promptxt = 'Select a Case to export';                 
            else
                promptxt = 'Select a Case to plot'; 
            end
            
            cobj = selectCaseObj(obj.Cases,[],{'CF_TransModel'},promptxt);
            if isempty(cobj)
                getdialog('No transgression model results available');
                return; 
            end
            cobj.cns.y2s = obj.Constants.y2s;
            switch src.Text
                case 'Grid Plot'
                    cf_summarygridplot(cobj);
                case 'Change Plot'
                    cf_changeplot(cobj);
                case 'Section Plot'
                    cf_sectionplot(cobj);
                case'Transgression Plot'
                    cf_summarytransplot(cobj); 
                case 'Animation'
                    cf_animation(cobj);
                case 'Export Tables'
                    cf_writetable(cobj);
            end
        end
        
        %% Help menu ------------------------------------------------------
        function Help(~,src,~)
            %menu to access online documentation and manual pdf file
            switch src.Text
                case 'Documentation'
                    doc channelform   %must be name of html help file  
                case 'Manual'
                    cfm_open_manual;
            end
        end 
        %% Check that toolboxes are installed------------------------------
        function isok = check4muitoolbox(~)
            %check that dstoolbox and muitoolbox have been installed
            fname = 'dstable.m';
            dstbx = which(fname);
        
            fname = 'muiModelUI.m';
            muitbx = which(fname);
        
            if isempty(dstbx) && ~isempty(muitbx)
                warndlg('dstoolbox has not been installed')
                isok = false;
            elseif ~isempty(dstbx) && isempty(muitbx)
                warndlg('muitoolbox has not been installed')
                isok = false;
            elseif isempty(dstbx) && isempty(muitbx)
                warndlg('dstoolbox and muitoolbox have not been installed')
                isok = false;
            else
                isok = true;
            end
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
    
    
    
    
    
    
    
    
    
    