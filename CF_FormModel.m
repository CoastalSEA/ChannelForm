classdef CF_FormModel < FGDinterface  
%
%-------class help---------------------------------------------------------
% NAME
%   CF_FormModel.m
% PURPOSE
%   Class for exponential, power plan and ckfa form model used in the 
%   ChannelForm App
% NOTES
%   Co-ordinate convention is that the x-axis is landward from the mouth to
%   the tidal limit. The y-axis is the cross-channel axis, mirrored about the
%   centre-line.
%   If a hydraulic surface can be based on the CSTmodel, a tidal amplitude
%   that decays linearly to the tidal limit, or constant high and low water 
%   levels. The surrounding surface is at high water level + a small offset
%
%   Inherits FGDInterface which is a subclass of GDinterface and muiDataSet 
%   and is part of the grid tools used by mui Apps.
% SEE ALSO
%   muiDataSet and FGDinterface
%
% Author: Ian Townend
% CoastalSEA (c) Jan 2022
%--------------------------------------------------------------------------
%     
    properties
        %inherits Data, RunParam, MetaData and CaseIndex from muiDataSet
        %Additional properties:   
        Selection %struct for plan, channel and intertidal form selection
                  %and type of along-channel water level model        
    end
 
    properties (Transient)
        CSTparams %struct for model summary parameters used when calling 
                  %cst model in CF_HydroData, or get_sed_flux in
                  %CF_TransModel
    end    
    
    properties (Dependent, SetAccess=private)
        zMouthInvert        %zm - thalweg bed level at mouth to zero datum (m)
    end    
    
    methods
        function obj = CF_FormModel()                   
            %class constructor
            %model being run, selections for cf_exp_model, wlflag for wl 
            %selection, flag to indicate if modifiaction have been added
            obj.Selection = struct('modeltype','','planform',0,'intertidalform',0,...
                                   'channelform',0,'wlflag',0,'incmods',false);
        end
    end      
%%
    methods (Static)        
%--------------------------------------------------------------------------
% Model implementation
%--------------------------------------------------------------------------         
        function obj = runModel(mobj,option)
            %function to run a simple 2D diffusion model
            obj = CF_FormModel;             
            %now check that the input data has been entered
            %isValidModel checks the InputHandles defined in ChannelForm
            if ~isValidModel(mobj, metaclass(obj).Name)  
                warndlg('Use Setup to define model input parameters');
                return;
            end
            
%--------------------------------------------------------------------------
% Model code
%--------------------------------------------------------------------------
            hydobj = getClassObj(mobj,'Inputs','CF_HydroData');
            setTransHydroProps(hydobj,mobj); %initialise transient properties
            
            %assign the run parameters to the model instance           
            setRunParam(obj,mobj); %assigns a copy of Input classes to obj
            %add water level definition to the run parameters
            [obj.Selection.wlflag,wlstxt] = setWaterLevels(obj);
            
            obj.Selection.modeltype = option;
            switch option
                case 'Exponential'                    
                    %prompt user to select plan form and x-sect form 
                    %model selection assigned to obj.Selection struct
                    obj.RunParam.CF_FormData = getClassObj(mobj,'Inputs','CF_ExpData');
                    hf = setFormSelection(obj); 
                    waitfor(hf);
                    [xi,yi,zi,yz] = cf_exp_models(obj);
                    sel = obj.Selection;
                    meta.data = sprintf('%s plan form, %s intertidal, %s channel, %s',...
                                        sel.planform,sel.intertidalform,...
                                        sel.channelform,wlstxt);
                case 'Power'
                    obj.RunParam.CF_FormData = getClassObj(mobj,'Inputs','CF_PowerData');
                    [xi,yi,zi,yz] = cf_pow_model(obj);
                    meta.data = sprintf('PR power form, %s',wlstxt);
                case 'CKFA'
                    [xi,yi,zi,yz] = ckfa_form_model(obj);
                    meta.data = sprintf('CKFA exogenous form, %s',wlstxt);
            end
            if isempty(xi), return; end
            
            griddata = reshape(zi,1,length(xi),length(yi));  
            %check that x and y are 1xn and yz is 1xnx3
            if size(yz,2)~=3, yz = yz'; end
            %now assign results to object properties  
            gridyear = years(0);  %durataion data for rows 
            
            %could use function to find tidal limit???? Lt ******
            
            dims = struct('x',xi,'y',yi,'t',gridyear,'xM',0,...
                                          'Lt',max(xi),'ishead',false);
%--------------------------------------------------------------------------
% Assign model output to dstable using GDinterface.setGrid
%--------------------------------------------------------------------------                   
            %assign metadata about model and save grid
            meta.source = sprintf('%s(%s)',metaclass(obj).Name,option);
            obj = setGrid(obj,{griddata},dims,meta);
%--------------------------------------------------------------------------
% Add property dstables in function FGDinterface.setProperties
%--------------------------------------------------------------------------              
            %zwl and Wz are empty and resolved in setProperties
            %limits=0 to use grid to determine hypsometry limits        
            histint = obj.RunParam.GD_GridProps.histint;
            obj = setProperties(obj,[],[],0,histint);  
%--------------------------------------------------------------------------
% Save results
%--------------------------------------------------------------------------             
            setCase(mobj.Cases,obj,'form_model');
            getdialog('Run complete');
            DrawMap(mobj);
        end
%%
        function addShoreline(mobj)
            %add shoreline strip based on the idealised beach profile
            muicat = mobj.Cases;
            ftxt = 'Select Form Model to use:';
            obj = selectCaseObj(muicat,[],{'CF_FormModel','CF_ValleyModel','GD_ImportData'},ftxt);
            if isempty(obj), return; end
 
            grid = getGrid(obj,1);
            %add shoreline strip based on equilibrium profile 
            obj.RunParam.CF_ShoreData = getClassObj(mobj,'Inputs','CF_ShoreData');
            if isempty(obj.RunParam.CF_ShoreData)
                warndlg('No shoreline data available')
                return
            end
            grid = setShoreline(obj.RunParam.CF_ShoreData,grid,true);

            %create new grid dstable and update  
            formdst = copy(obj.Data.Grid);
            formdst.Dimensions.X = grid.x;
            sz = num2cell(size(grid.z));
            formdst.DataTable.Z = reshape(grid.z,1,sz{:}); 
            formdst.UserData.xM = grid.xM;

            %save as a new case
            setGridObj(obj,muicat,formdst); 
        end
%%
        function addMorphMods(mobj)
            %add a prismatic dredge channel to a channel form
            muicat = mobj.Cases;
            ftxt = 'Select Form Model to use:';
            obj = selectCaseObj(muicat,[],{'CF_FormModel'},ftxt);
            if isempty(obj), return; end
 
            grid = getGrid(obj,1);
            %apply the modifications defined in CF_ModsData to define new grid
            obj.RunParam.CF_ModsData = getClassObj(mobj,'Inputs','CF_ModsData');
            new_z = setMorphMods(obj.RunParam.CF_ModsData,grid);
            obj.Selection.incmods = true;

            answer = questdlg('Add or update existing?','Add Mods','Add','Update','Add');
            
            if strcmp(answer,'Add')
                %create new record
                fdst = copy(obj.Data.Grid);
                fdst.Z(1,:,:) = new_z;
                setGridObj(obj,muicat,fdst); %copies form property table to new instance
            else            
                %overwrite exisitng form data set with new form             
                obj.Data.Grid.Z(1,:,:) = new_z;  
                classrec = classRec(muicat,caseRec(muicat,obj.CaseIndex));
                updateCase(muicat,obj,classrec,true);
            end 
        end 
    end        
%%
    methods
        function tabPlot(obj,src) %abstract class method for muiDataSet
            %generate plot for display on Q-Plot tab
            cf_model_tabs(obj,src);
        end
%%
        function zmouth = get.zMouthInvert(obj)
            %zm - thalweg bed level at mouth to zero datum (m)
            %dependent property needs mean sea level and depth at mouth
            frmobj = obj.RunParam.CF_FormData;
            hydobj = obj.RunParam.CF_HydroData;
            if isempty(obj.Data)   %new grid not yet defined
                ixM = 1;
            else                   %grid may include an offset
                grid = getGrid(obj,1);            %irow=1 assumed
                [~,ixM] = gd_basin_indices(grid); %account for offset to mouth
            end
            zmouth = hydobj.zmt(ixM)-frmobj.MTmouthDepth;
        end
    end 
%%    
    methods (Access = private) 
        function [wlflag,mtxt] = setWaterLevels(~)
            %set water levels for form model using either the surface
            %defined by the cst_model, or high water at the mouth and a 
            %reducing tidal amplitude    
            % wlflag - flag to indicate type of water surface to use
            % mtxt - text to define type of water surface
            answer = questdlg('Use which hydro surface?','Select hydro',...
                              'CST surface','Constant HW',...
                              'Constant HW&LW','CST surface');

            if strcmp(answer,'CST surface')
                wlflag = 0;
                mtxt = 'hydraulic surfaces from CSTmodel';
            elseif strcmp(answer,'Constant HW')
                wlflag = 1; 
                mtxt = 'constant HW surface, tapered LW surface';
            else
                wlflag = 2;
                mtxt = 'constant HW and LW surfaces';
            end                              
        end
%%
        function hf = setFormSelection(obj)   
            %initialise ui for selection of form options
            hf = figure('Name','Channel Form', ...
                'NumberTitle','off', ...
                'MenuBar','none', ...
                'Units','normalized', ...
                'Position',[0.2 0.5 0.2 0.2], ...
                'Resize','on','HandleVisibility','on', ...
                'Tag','ChannelForm');
            axes('Parent',hf, ...
                'Color',[0.94,0.94,0.94], ...
                'Position',[0 0.002 1 1], ...
                'XColor','none', ...
                'YColor','none', ...
                'ZColor','none', ...
                'Tag','Gui'); 
            
            vartitle = {'Plan form:','Channel form:','Intertidal form:'};
            varorder = {'PlanUI','ChannelUI','IntertidalUI'};
            uioptions{1} = {'Exponential','Power'};
            uioptions{2} = {'Parabolic','Rectangular'};
            uioptions{3} = {'Linear','Rectangular','Stepped','Uniform Shear','L&M muddy shore'};

            for i=1:3
                height = 1-i/5;
                uicontrol('Parent',hf, 'Style','text',...
                    'String',vartitle{i},...
                    'HorizontalAlignment', 'right',...
                    'Units','normalized', ...
                    'Position',[0.01 height-0.05 0.3 0.09]);
                uicontrol('Parent',hf, ...
                    'Style','popupmenu', ...
                    'Units','normalized', ...
                    'Position',[0.32 height 0.58 0.05], ...
                    'String',uioptions{i},...
                    'Tag',varorder{i}, ...
                    'Value',int16(1)); %max list length is 32767
            end            
            uicontrol('Parent',hf,...  %callback button
                'Style','pushbutton',...
                'String', 'Accept selection',...
                'Units','normalized', ...
                'Position', [0.3,0.1,0.4,0.1],...
                'Callback',@obj.pushbutton_callback);  
        end
%%
        function pushbutton_callback(obj,src,~)
            %callback function for the accept button on the selection UI
            hf = src.Parent;
            temp = findobj(hf.Children,'Tag','PlanUI');
            obj.Selection.planform = temp.String{temp.Value};
            temp = findobj(hf.Children,'Tag','ChannelUI');
            obj.Selection.channelform = temp.String{temp.Value};
            temp = findobj(hf.Children,'Tag','IntertidalUI');
            obj.Selection.intertidalform = temp.String{temp.Value}; 
            close(hf);
        end
    end    
end