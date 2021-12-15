classdef CF_FormModel < GDinterface  
%
%-------class help---------------------------------------------------------
% NAME
%   CF_FormModel.m
% PURPOSE
%   Class for exponential and power plan form model used in ChannelForm App
% NOTES
%   Co-ordinate convention is that the x-axis is seaward from the head or
%   tidal limit. The y-axis is the cross-channel axis mirrored about the
%   centre-line.
%   If a hydraulic surface (based on CSTmodel) is not included the
%   surrounding surface is at high water level +a small offset and the 
%   tidal amplitude is assumed to decay linearly to the tidal limit.
% SEE ALSO
%   muiDataSet and GDinterface
%
% Author: Ian Townend
% CoastalSEA (c) Jan 2022
%--------------------------------------------------------------------------
%     
    properties
        %inherits Data, RunParam, MetaData and CaseIndex from muiDataSet
        %Additional properties:   
        ModelType %type of form model being used (Exp, Power, CKFA)
        Selection %struct for plan, channel and intertidal form selection
        CKFAform  %struct for CKFA model property tables for:
                  %form, flow and wave
        Channel   %struct for Channel model summary parameters
    end
    
    methods (Access = private)
        function obj = CF_FormModel()                   
            %class constructor
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
            %isValidModel checks the InputHandles defined in ModelUI
            if ~isValidModel(mobj, metaclass(obj).Name)  
                warndlg('Use Setup to define model input parameters');
                return;
            end
            
%--------------------------------------------------------------------------
% Model code
%--------------------------------------------------------------------------
%             wlvobj = getClassObj(mobj,'Inputs','WaterLevels');
            hydobj = getClassObj(mobj,'Inputs','CF_HydroData');
            setTransHydroProps(hydobj,mobj); %initialise transient properties
            
            %assign the run parameters to the model instance           
            setRunParam(obj,mobj); 
            %add water level definition to the run parameters
            [iscst,wlstxt] = setWaterLevels(obj);
            
            obj.ModelType = option;
            switch option
                case 'Exponential'                    
                    %prompt user to select plan form and x-sect form 
                    %model selection assigned to obj.Selection struct
                    obj.RunParam.CF_FormData = getClassObj(mobj,'Inputs','CF_ExpData');
                    hf = setFormSelection(obj); 
                    waitfor(hf);
                    [xi,yi,zi,yz] = channel_form_models(obj,iscst);
                    sel = obj.Selection;
                    meta.data = sprintf('%s plan form, %s intertidal, %s channel, %s',...
                                        sel.planform,sel.intertidalform,...
                                        sel.channelform,wlstxt);
                case 'Power'
                    obj.RunParam.CF_FormData = getClassObj(mobj,'Inputs','CF_PowerData');
                    [xi,yi,zi,yz] = pr_form_model(obj,iscst);
                    meta.data = sprintf('PR power form, %s',wlstxt);
                case 'CKFA'
                    [xi,yi,zi,yz] = ckfa_form_model(obj,iscst);
                    meta.data = sprintf('CKFA exogenous form, %s',wlstxt);
%                 case 'Valley'
%                     [xi,yi,zi,yz] = cf_valley_model(obj,iscst);
%                     meta.data = sprintf('Valley form, %s',wlstxt);
            end
            if isempty(xi), return; end
            
            griddata = reshape(zi,1,length(xi),length(yi));  
            %check that x and y are 1xn and yz is 1xnx3
            if size(yz,2)~=3, yz = yz'; end
            results = [{griddata},yz];
%             xydata = {xi,yi};     %assign xyz coordinates as row vectors
            %now assign results to object properties  
            gridyear = years(0);  %durataion data for rows 
            dims = struct('x',xi,'y',yi,'t',gridyear);
%--------------------------------------------------------------------------
% Assign model output to dstable using the defined dsproperties meta-data
%--------------------------------------------------------------------------                   
            %assign metadata about model and save grid
            meta.source = metaclass(obj).Name;
            dst = setGridOuput(obj,results,dims,meta);
%--------------------------------------------------------------------------
% Add property dstables in function GDinterface.setFormProperties
%--------------------------------------------------------------------------  
            dst = setFormProps(obj,dst,1,meta);            
%--------------------------------------------------------------------------
% Save results
%--------------------------------------------------------------------------             
            setDataSetRecord(obj,mobj.Cases,dst,'form_model');
            getdialog('Run complete');
            DrawMap(mobj);
        end
    end
%%
    methods
        function tabPlot(obj,src) %abstract class method for muiDataSet
            %generate plot for display on Q-Plot tab
            cf_model_tabs(obj,src);
        end
    end 
%%    
    methods (Access = private) 
        function [iscst,mtxt] = setWaterLevels(obj)
            %set water levels for form model using either the surface
            %defined by the cst_model, or high water at the mouth and a 
            %reducing tidal amplitude            
            answer = questdlg('Include hydro surface?','Select hydro',...
                              'CST surface','Constant HW','CST surface');
            
            if strcmp(answer,'Constant HW')
                iscst = false; 
                mtxt = 'hydraulic surfaces use HW at mouth';
            else
                iscst = true;
                mtxt = 'hydraulic surfaces from CSTmodel';
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