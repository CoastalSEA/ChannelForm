classdef CF_FormModel < GDinterface  
%
%-------class help---------------------------------------------------------
% NAME
%   CF_FormModel.m
% PURPOSE
%   Class for exponential and power plan form model used in ChannelForm App
%
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
        Selection %structure for plan, channel and intertidal form selection
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
            %assign the run parameters to the model instance
            %may need to be after input data selection to capture caserecs
            %prompt user to select plan form and x-sect form  
            setRunParam(obj,mobj); 
%--------------------------------------------------------------------------
% Model code
%--------------------------------------------------------------------------
            wlvobj = getClassObj(mobj,'Inputs','WaterLevels');
            hydobj = getClassObj(mobj,'Inputs','CF_HydroData');
            [hydobj.zhw,hydobj.zmt,hydobj.zlw] = newWaterLevels(wlvobj,0,0);
            switch option
                case 'Exponential'
                    %model selection assigned to obj.Selection struct
                    hf = setFormSelection(obj); 
                    waitfor(hf);
                    [xi,yi,zi,yz] = channel_form_models(obj,mobj);
                    sel = obj.Selection;
                    meta.data = sprintf('%s plan form, %s intertidal, %s channel',...
                          sel.planform,sel.intertidalform,sel.channelform);
                case 'Power'
                    [xi,yi,zi,yz] = pr_form_model(obj,mobj);
                    meta.data = 'PR power form';
                case 'CKFA'
                    [xi,yi,zi,yz] = ckfa_properties(obj,mobj);
                    meta.data = 'CKFA ideal form';
            end

            if isempty(xi), return; end
            griddata = reshape(zi,1,length(xi),length(yi));  
            results = [{griddata},yz'];
            xydata = {xi',yi'};      %assign xyz coordinates
            %now assign results to object properties  
            gridyear = years(0);  %durataion data for rows 
%--------------------------------------------------------------------------
% Assign model output to dstable using the defined dsproperties meta-data
%--------------------------------------------------------------------------                   
            %assign metadata about model and save grid
            meta.source = metaclass(obj).Name;
            dst = setGridOuput(obj,results,xydata,gridyear,meta);
%--------------------------------------------------------------------------
% Add property dstables in function GDinterface.setFormProperties
%--------------------------------------------------------------------------  
            grdobj = getClassObj(mobj,'Inputs','CF_GridData');
            dst = setFormProps(obj,dst,1,grdobj,hydobj);            
%--------------------------------------------------------------------------
% Save results
%--------------------------------------------------------------------------             
            setDataSetRecord(obj,mobj.Cases,dst,'model');
            getdialog('Run complete');
        end
    end
%%
    methods
        function tabPlot(obj,mobj,src) %abstract class method for muiDataSet
            %generate plot for display on Q-Plot tab
            cf_model_tabs(obj,mobj,src);
        end
    end 
%%    
    methods (Access = private)       
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