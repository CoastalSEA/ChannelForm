classdef Model_template < muiDataSet                         % << Edit to classname
%
%-------class help---------------------------------------------------------
% NAME
%   Model_template.m
% PURPOSE
%   Class description - Class for Model XXXXX to be run as a muitoolbox app
%
% SEE ALSO
%   muiDataSet
%
% Author: Ian Townend
% CoastalSEA (c) Jan 2021
%--------------------------------------------------------------------------
%     
    properties
        %inherits Data, RunParam, MetaData and CaseIndex from muiDataSet
        %Additional properties:     
    end
    
    methods (Access = private)
        function obj = Model_template()                      % << Edit to classname
            %class constructor
        end
    end      
%%
    methods (Static)        
%--------------------------------------------------------------------------
% Model implementation
%--------------------------------------------------------------------------         
        function obj = runModel(mobj)
            %function to run a simple 2D diffusion model
            obj = Model_template;                            % << Edit to classname
            dsp = modelDSproperties(obj);
            
            %now check that the input data has been entered
            %isValidModel checks the InputHandles defined in ModelUI
            if ~isValidModel(mobj, metaclass(obj).Name)  
                warndlg('Use Setup to define model input parameters');
                return;
            end
            muicat = mobj.Cases;
            %assign the run parameters to the model instance
            %may need to be after input data selection to capture caserecs
            setRunParam(obj,mobj); 
%--------------------------------------------------------------------------
% Model code  <<INSERT MODEL CODE or CALL MODEL>>
%--------------------------------------------------------------------------
            inp = mobj.Inputs.PropsInput_template_1;         % <<Edit to suit model
            run = mobj.Inputs.PropsInput_template_2;
            [results,xy,modeltime] = your_model(inp,run);
            %now assign results to object properties  
            modeltime = seconds(modeltime);  %durataion data for rows 
            %code for datetime format
            % modeldate = datetime(rnpobj.StartYear,1,1,0,0,0);
            % modeltime = modeldate + seconds(modeltime);
            % modeltime.Format = dsp.Row.Format;
            %code for duration format
            % modeltime = years(StartYear)+seconds(modeltime);
            % modeltime.Format = dsp.Row.Format;
            %code for calendatDuration format 
            % modeldate = datetime(StartYear,1,1,0,0,0);
            % modeltime = modeldate + seconds(modeltime);
            % modeltime.Format = dsp.Row.Format;
            % modeltime = calendarDuration(datevec(modeltime));
%--------------------------------------------------------------------------
% Assign model output to a dstable using the defined dsproperties meta-data
%--------------------------------------------------------------------------                   
            %each variable should be an array in the 'results' cell array
            %if model returns single variable as array of doubles, use {results}
            dst = dstable(results,'RowNames',modeltime,'DSproperties',dsp);
            dst.Dimensions.X = xy{:,1};     %grid x-coordinate
            dst.Dimensions.Y = xy{:,2};     %grid y-coordinate
%--------------------------------------------------------------------------
% Save results
%--------------------------------------------------------------------------                        
            %assign metadata about model
            dst.Source = metaclass(obj).Name;
            dst.MetaData = 'Any additional information to be saved';
            %save results
            setDataSetRecord(obj,muicat,dst,'model');
            getdialog('Run complete');
        end
    end
%%
    methods
        function tabPlot(obj,src) %abstract class for muiDataSet
            %generate plot for display on Q-Plot tab
            
            %add code to define plot format or call default tabplot using:
            tabDefaultPlot(obj,src);
        end
    end 
%%    
    methods (Access = private)
        function dsp = modelDSproperties(~) 
            %define a dsproperties struct and add the model metadata
            dsp = struct('Variables',[],'Row',[],'Dimensions',[]); 
            %define each variable to be included in the data table and any
            %information about the dimensions. dstable Row and Dimensions can
            %accept most data types but the values in each vector must be unique
            
            %struct entries are cell arrays and can be column or row vectors
            dsp.Variables = struct(...                       % <<Edit metadata to suit model
                'Name',{'Var1','Var2'},...
                'Description',{'Variable 1','Variable 2'},...
                'Unit',{'m/s','m/s'},...
                'Label',{'Plot label 1','Plot label 2'},...
                'QCflag',{'model','model'}); 
            dsp.Row = struct(...
                'Name',{'Time'},...
                'Description',{'Time'},...
                'Unit',{'s'},...
                'Label',{'Time (s)'},...
                'Format',{'s'});        
            dsp.Dimensions = struct(...    
                'Name',{'X','Y'},...
                'Description',{'X co-ordinate','Y co-ordinate'},...
                'Unit',{'m','m'},...
                'Label',{'X co-ordinate (m)','Y co-ordinate (m)'},...
                'Format',{'-','-'});  
        end
    end    
end