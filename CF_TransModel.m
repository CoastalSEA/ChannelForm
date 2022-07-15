classdef CF_TransModel < GDinterface  
%
%-------class help---------------------------------------------------------
% NAME
%   CF_TransModel.m
% PURPOSE
%   Class for transgression model to compute the movement of a
%   channel within a valley in response to sea level rise using a simple
%   kinetic model.
% NOTES
%   Co-ordinate convention is that the x-axis is landward from the mouth to
%   the tidal limit. The y-axis is the cross-channel axis, mirrored about the
%   centre-line.
%   If a hydraulic surface can be based on the CSTmodel, a tidal amplitude
%   that decays linearly to the tidal limit, or constant high and low water 
%   levels. The surrounding surface is at high water level + a small offset
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
        Selection %struct for plan, channel and intertidal form selectio
    end
%     
    properties (Transient)
        CSTparams      %struct for model summary parameters used when calling 
                       %cst model in CF_HydroData
        Channel        %handle to instance of channel model being used
        Valley         %handle to instance of valley model being used
        Time = 0       %time elapsed from start of run(t=0) in seconds
        DateTime = 0   %time elapsed from Year 0 in seconds
        iStep = 0      %current step number during run
        RunSteps = 0   %number of time steps after checking stability
        delta  = 0     %time step in seconds 
        StepTime       %time to be saved (seconds during run and converted to years)
        Grid           %model grid at timestep, t, struct with
                       % x,y,z co-ordinates
                       % t - timstep (years)
                       % Wz - array of column vectors for width at hw,mt,lw
                       % zdiff - difference over a timestep
                       % intidx - x-axis indices from mouth to tidal limit
                       % hrv - hydraulic depthof river (m)
                       % Wrv - width of river (m)
                       % xM = distance to mouth in the grid
        CLatT          %struct array for x and y coordinates of meander centreline at each time step
        zGrid          %z grids at each sampled time step  
        zWL            %alongchannel water levels at each sampled time step
        tPlan          %plan form, Whw,Wmt,Wlw, at each sampled time step
        dTrans         %incremental changes in timestep, t, includes:
                       % delX - unadjusted transgression distance
                       % estdX - adjusted transgression distance (inc sediment flux)
                       % cstdX - open coast transgression distance
                       % dSLR - increase in mean sea level at mouth
                       % Lt - distance to tidal limit (to record change - model uses CF_HydroData.xTidalLimit)
                       % vdiffx - volume difference for [0,delX/2,delX,3delX/2]
                       % sedVol - sediment flux (+ve=sediment import)
                       % FPA - flood plain area
        Trans          %struct used to store transgresson output at each time step                       
                       % dTrans struct variables at each sampled time step
                       %        delX,estdX,csrsX and sedVol values of 
                       %        dTrans at each timestep are converted
                       %        to cumulative values at end of run
        cns            %struct of default model constants
    end
    
    methods
        function obj = CF_TransModel()                   
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
            obj = CF_TransModel;             
            %now check that the input data has been entered
            %isValidModel checks the InputHandles defined in ChannelForm
            if ~isValidModel(mobj, metaclass(obj).Name)  
                warndlg('Use Setup to define model input parameters');
                return;
            end
            %select baseline models to be used and assign to Channel and
            %Valley propeties. RunParam are also copied from Channel to
            %new model instance
            [obj,isok] = selectInputModels(obj,mobj);
            if ~isok, return; end
%--------------------------------------------------------------------------
% Model code
%--------------------------------------------------------------------------
            %iniatialise model setup 
            obj = InitialiseModel(obj,mobj); %uses inputs assigned to RunParam
            if isempty(obj), return; end

            msg = sprintf('ChannelForm processing, please wait');
            hw = waitbar(0,msg);
            %run model
            for jt = 1:obj.RunSteps
                InitTimeStep(obj,mobj,jt)
                RunTimeStep(obj)
                PostTimeStep(obj);
                %if waitbar has been deleted then program has failed to
                %find a solution (propable in CSTmodel)
                if ~isvalid(hw), return; end
                %to report time step during run use the following
                msg = sprintf('ChannelForm processing, step %d of %d',...
                                                         jt,obj.RunSteps);
                waitbar(jt/obj.RunSteps,hw,msg);                
            end
            close(hw);
            %cumulative values for first 4 variables in Trans output table
            %order is: 'delX','estdX','cstdX','dSLR','Lt','FPA','sedVol','vdiffx'                                                
            obj.Trans(:,1:4) = varfun(@cumsum,obj.Trans(:,1:4)); 
            obj.Trans.waterVol = cumsum(obj.Trans.waterVol,1);
            obj.Trans.sedVol = cumsum(obj.Trans.sedVol,1);
            obj.Trans.vdiffx = cumsum(obj.Trans.vdiffx,1);
            %plots of run
            slr = obj.Trans.dSLR(end);
%             summaryGridPlots(obj,slr);
%             crossectionPlot(obj,slr);
            changePlot(obj,slr);         
%             thalwegPlot(obj,slr);

            %now assign results to object properties  
            mtime = years(obj.StepTime/obj.cns.y2s);
            dims = struct('x',obj.Grid.x,'y',obj.Grid.y,'t',mtime,...
                          'ishead',false,'xM',obj.Trans.cstdX,'cline',obj.CLatT);
%--------------------------------------------------------------------------
% Assign model output to dstable using the defined dsproperties meta-data
%--------------------------------------------------------------------------                   
            %assign metadata about model and save grid
            meta.source = obj.MetaData;
            meta.data = obj.Channel.MetaData;
            obj = setGrid(obj,{obj.zGrid},dims,meta);
            obj = setPlanProps(obj,obj.tPlan,meta);  %half width data           
            obj = setWLProps(obj,obj.zWL,meta); 
%--------------------------------------------------------------------------
% Add property dstables in function GDinterface.setFormProperties
%--------------------------------------------------------------------------  
            obj = setHyps_SP_GP(obj,meta);
            obj = setTransModel(obj,meta);
%--------------------------------------------------------------------------
% Save results
%--------------------------------------------------------------------------             
            setCase(mobj.Cases,obj,'trans_model');
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
%% ------------------------------------------------------------------------
% Methods called by timestep loop
%--------------------------------------------------------------------------
    methods (Access=private)
        function obj = InitialiseModel(obj,mobj)
            %initialise ChannelForm properties and run parameters
            rnpobj = obj.RunParam.RunProperties;
            %initialise time step paramters
            tstep = rnpobj.TimeStep; 
            y2s   = mobj.Constants.y2s;     %year to seconds conversion factor
            obj.delta = tstep*y2s;          %time step in seconds
            obj.DateTime = rnpobj.StartYear*y2s;%time elapsed from Year 0 in seconds
            obj.RunSteps = rnpobj.NumSteps; %not needed unless stability check added
            %model constants
            obj.cns = getConstantStruct(mobj.Constants); 
            
            %use selected channel and valley forms to create initial grid
            fgrid = getGrid(obj.Channel,1);
            %extract plan form width data
            Wz =obj.Channel.Data.Plan.DataTable(1,:);
            fgrid.Wz = table2cell(Wz);  %COULD change to use a table***
            obj.Grid = updateValley(obj,fgrid,false);
            obj.Grid.xM = 0;  %initialise mouth location at x=0
%             obj.CLatT = []; %initialise meander centreline
            %model selection options(wlflag,modeltype,etc)
            obj.Selection = obj.Channel.Selection;
            obj.CSTparams = obj.Channel.CSTparams;
            
            %initialise water levels for Time=0
            setTransHydroProps(obj.RunParam.CF_HydroData,mobj);
            %along channel water levels are transient properties. Set 
            %initial calues using Channel water levels          
            setChannelWaterLevels(obj);
            %get the distance to the tidal limit and index on the x-axis
            updateTidalLimit(obj); %sets Grid.intidx
            %get the dimensions of the river for the initial river discharge
            [obj.Grid.hrv,obj.Grid.Wrv] = riverDimensions(obj);
            
            %transgression ouput table for 
            % delX - unadjusted transgression distance
            % estdX - adjusted transgression distance (inc sediment flux)
            % cstdX - open coast transgression distance
            % dSLR - rate of sea level rise in time step
            % Lt - distance to tidal limit (to record change - model uses CF_HydroData.xTidalLimit)
            % FPA - flood plain area
            % waterVol - water volume change due to changes in slr and tidal range
            % sedVol - sediment flux (+ve=sediment import)            
            % vdiffx - volume difference for [0,delX/2,delX,3delX/2]
            varnames = {'delX','estdX','cstdX','dSLR','Lt','FPA',...
                                           'waterVol','sedVol','vdiffx'};                                                        
            Lt = obj.RunParam.CF_HydroData.xTidalLimit;
            FPA = floodPlainArea(obj);
            %sedflux = sedimentFlux(obj);
            vars = {0,0,0,0,Lt,FPA,0,0,[0,0,0,0]};                    
            obj.dTrans = table(vars{:},'VariableNames',varnames);
            %write data for initial time step (t=0)
            PostTimeStep(obj); 
        end
%%
        function InitTimeStep(obj,mobj,jt)
            %initialise model parameters for next time step
            obj.iStep = jt;
            obj.Time = jt*obj.delta;
            obj.DateTime = obj.DateTime+obj.delta;
            obj.Grid.t = obj.DateTime/obj.cns.y2s;
            %update form input parameters
            obj = updateInputParams(obj);
            %get the distance to the tidal limit and index on the x-axis
            updateTidalLimit(obj); 
            %update water levels at boundary
            newWaterLevels(obj.RunParam.CF_HydroData,mobj,obj);
            dt = obj.delta/obj.cns.y2s;
            obj.dTrans.dSLR = obj.RunParam.CF_HydroData.dslr*dt;
            %could also update river discharge if required
            % Qr = source of river discharge time series;
            % obj.RunParam.CF_HydroData.RiverDischarge = Qr;
        end
 %%
        function RunTimeStep(obj)
            %run model for the time step jt
            trnobj = obj.RunParam.CF_TransData;
            [obj.dTrans.sedVol,obj.dTrans.waterVol] = sedimentFlux(obj);
            %get transgression distance for the combined form
            % profile on
            obj = getTransDist(obj);
            % s = profile('info');
            % profile off
            %get the transgression distance for zero volume change
            delX = obj.dTrans.delX;
            vdiffx = obj.dTrans.vdiffx;
            delint = [0,delX/2,delX,delX+delX/2];
            %adjusted delX accounts for sediment input/output
            obj.dTrans.estdX = interp1(vdiffx,delint,obj.dTrans.sedVol,...
                                                        'linear','extrap');            
            %shoreline change on open coast
            hydobj = obj.RunParam.CF_HydroData;
            incslr = hydobj.dslr*obj.delta/obj.cns.y2s; %slr in time step
            obj.dTrans.cstdX = trnobj.BruunRatio*incslr;%coastal transgression 
            obj.Grid.xM = sum(obj.Trans.cstdX)+obj.dTrans.cstdX;
                    
            %update grid based on transgression
            fgrid = getNewChannel(obj);            %new channel form
            if isempty(fgrid.x), return; end
            fgrid = updateMeander(obj,fgrid);            
            fgrid = updateValley(obj,fgrid,true);  %new combined form            
            obj.Grid = updateShoreline(obj,fgrid); %shoreline adjusted 
            obj = cf_offset_wls(obj);              %translate wls if coast erodes
            obj.dTrans.FPA = floodPlainArea(obj);  %modified flood plain area
            %summary plot of change from initial form (used for checking)
            % summaryGridPlots(obj);
        end 
%%
        function obj = PostTimeStep(obj)
            %store the results for each time step
            jr = length(obj.StepTime)+1;
            obj.StepTime(jr,1) = obj.DateTime;   %numeric seconds (double)
            %add variables to transgression table
            obj.Trans = [obj.Trans;obj.dTrans];
            %add grid of new form (incl valley)
            sz = num2cell(size(obj.Grid.z));
            obj.zGrid = [obj.zGrid;reshape(obj.Grid.z,1,sz{:})];
            %add meander centreline coordinates
            obj.CLatT = [obj.CLatT,obj.Grid.cline];
            %add plan form description
            Wz = obj.Grid.Wz';
            tPlani = table(Wz{:},'VariableNames',{'Whw','Wmt','Wlw'});
            obj.tPlan = [obj.tPlan;tPlani];
            %add alongchannel water levels
            hydobj = obj.RunParam.CF_HydroData;
            zwli = table(hydobj.zhw',hydobj.zmt',hydobj.zlw',...
                                'VariableNames',{'zhw','zmt','zlw'});                
            obj.zWL = [obj.zWL;zwli];
        end       
%% ------------------------------------------------------------------------
% Methods called during runtime
%--------------------------------------------------------------------------
        function obj = getTransDist(obj)
            %find transgression distance for a given SLR increment
            %using fzero is slower but more accurate  
            F = obj.Grid;           
            Lt = obj.RunParam.CF_HydroData.xTidalLimit;
            slr = obj.RunParam.CF_HydroData.dslr*obj.delta/obj.cns.y2s;
            adist = Lt/abs(min(min(F.z,[],'omitnan')))*slr; %used as initial guess
            %reducing TolX to 1 reduces time by ~7% but increases vdiff@delX
            %using fminsearch is much slower
            % vdiff = @(delx) getZdiff(obj,delx); 
            % options = optimset('TolX',1e-3); 
            % [delX,~] = fzero(vdiff,adist,options);
            delX = interp_delX(obj,adist);  %quick and dirty option (comment-out fzero)
            [~,zdiff] = getZdiff(obj,delX);
            delint = [0,delX/2,delX,delX+delX/2];
            for i=1:length(delint)
                obj.dTrans.vdiffx(1,i) = getZdiff(obj,delint(i));
            end             
            obj.dTrans.delX = delX;
            obj.Grid.zdiff = zdiff;
        end
%%
        function [vdiff,zdiff] = getZdiff(obj,delX)
            %compute the difference between the surfaces for given distance 
            %of translation, delX. Uses current grid+dhw to define new grid             
            F = obj.Grid;
            dxdy = abs(F.x(2)-F.x(1))*abs(F.y(2)-F.y(1));
            %add change in high water over a time step
            z2 = F.z+obj.RunParam.CF_HydroData.dhw;
            x2 = F.x+delX;    
            [X,Y] = ndgrid(F.x,F.y);   
            z1 = F.z;
            z2 = griddata(x2,F.y,z2',X,Y);    
            zdiff = (z2-z1);
            
            %apply limits to domain and any vertical constraints
            zdiff = applyConstraints(obj,zdiff);
            
            %subsample grid to integration length and find volume change
            subx = obj.Grid.intidx;  %indices from mouth to tidal limit
            if rem(obj.Grid.xM,(F.x(2)-F.x(1)))>0
                %if coast erosion occupies part of a grid interval in x
                %only include complete grid cells in integration
                subx = subx(2:end);
            end
            dz = zdiff(subx,:);
            dz(isnan(dz)) = 0;
            vdiff = trapz(trapz(dz))*dxdy;
        end
%%
        function zdiff = applyConstraints(obj,zdiff)
            %apply limits to domain and any vertical constraints
            %control lateral spatial extent
            trnobj = obj.RunParam.CF_TransData;                        
            [Y,Yhw0,Yhw,~,zFP,dely] = constraintGrids(obj,obj.Grid);
            
            if trnobj.inclHWConstraint       %sea wall fixes HW line                
                zdiff(Y<Yhw0.left) = NaN;    %left side of initial hw
                zdiff(Y>Yhw0.right) = NaN;   %right side of initial hw
            elseif trnobj.inclFloodPlain     %open flood plain
                %flood plain surface that determines intersection with valley slope
                zdiff(obj.Grid.z>zFP) = NaN; %use combined channel+valley to HW limit        
            else
                %HW changes but flood plain excluded, sediment demand limited to near bank
                fact = 4; %needs to be consistent with value in updateValley
                %use channel HW+fact.dy as limit of volume calcs   
                zdiff(Y<(Yhw.left-fact*dely)) = NaN;
                zdiff(Y>(Yhw.right+fact*dely)) = NaN;
            end

            %control vertical constraints
            if trnobj.inclGeoConstraint && ~isempty(trnobj.StConstraints)
                ids = length(trnobj.StConstraints);
                for i=1:ids
                    x = obj.Grid.x;
                    idx = x>trnobj.StConstraints(i) & x<trnobj.NdConstraints(i) & zdiff<0;
                    zdiff(idx)=0;  %remove any erosion over constrained bed
                end
            end
        end
%%
        function delX = interp_delX(obj,adist)
            %interpolate end values to find delX
            dVpos = getZdiff(obj,0);
            dVneg = dVpos; 
            count = 2;           
            while dVneg>0
                %ensure the largest delX value gives a negative vdiffx
                newdel = count*adist;
                dVneg = getZdiff(obj,newdel);
                count = count+1;
            end
            delX = dVpos/((dVpos-dVneg)/newdel);
        end
%%
        function grid = getNewChannel(obj)
            %generate a new channel based on updated input parameters
            %with no offset (ie mouth is at x=0)
            grid = obj.Grid;
            z0 = grid.z;   %existing grid before updating (includes valley)
            option = obj.Channel.Selection.modeltype;
            switch option
                case 'Exponential'                    
                    %model selection assigned to obj.Selection struct
                    [grid.x,grid.y,grid.z,grid.Wz] = cf_exp_models(obj);
                case 'Power'
                    [grid.x,grid.y,grid.z,grid.Wz] = cf_pow_model(obj);
                case 'CKFA'
                    [grid.x,grid.y,grid.z,grid.Wz] = ckfa_form_model(obj);
            end    
            if isempty(grid.x)
                return;
            elseif isrow(grid.x)
                grid.x = grid.x'; 
            end  
            
            %if set apply the modifications defined in CF_ModsData
            if obj.Selection.incmods
                grid.z = squeeze(setMorphMods(obj.Channel,grid));
            end
            
            %apply any vertical constraints
            trnobj = obj.RunParam.CF_TransData; 
            if trnobj.inclGeoConstraint && ~isempty(trnobj.StConstraints)
                 %elevation of initial flood plain
                hydobj0 = obj.Channel.Data.WaterLevels; %source form water levels
                zFP0 = hydobj0.zhw+trnobj.FPoffset;     %initial flood plain levels
                zFP0 = repmat(zFP0',1,length(grid.y));
                ids = length(trnobj.StConstraints);
                for i=1:ids
                    %fixes channel only because only modifies below initial
                    %flood plain elevation, zFP0
                    x = grid.x;
                    idx = x>trnobj.StConstraints(i) & ...
                          x<trnobj.NdConstraints(i) & grid.z<zFP0;
                    grid.z(idx)=z0(idx);  %remove any erosion over constrained bed
                end
            end            
            grid.metadata = obj.Channel.Data.Form.MetaData;
            %adjust flood plain in channel model for any constraints
            grid = updateFloodPlain(obj,grid);            
        end
%%
        function grid = updateFloodPlain(obj,grid)
            %apply any constraints at high water, and/or exclude flood plain
            %note: geological constraints are applied when creating fgrid
            trnobj = obj.RunParam.CF_TransData;
            [Y,Yhw0,Yhw,zFP0,~,dely] = constraintGrids(obj,grid);

            if ~trnobj.inclFloodPlain && ~trnobj.inclHWConstraint
                %exclude flood plain with no constraint at HW so reset levels
                %to initial flood plain with a fact*dely border along bank
                fact = 4; %needs to be consistent with value in applyConstraints
                idzleft = Y<(Yhw.left-fact*dely);    %current hw widths
                idzright = Y>(Yhw.right+fact*dely);
                idz = idzleft | idzright;
                grid.z(idz) = zFP0(idz);
            end
            
            if trnobj.inclHWConstraint
                %include constraint at HW so reset levels to intial flood
                %plain and add a sea wall 
                idzleft = Y<(Yhw0.left-dely);        %initial hw widths
                idzright = Y>(Yhw0.right+dely);
                idz = idzleft | idzright;
                grid.z(idz) = zFP0(idz);    
                %
                idzleft =  Y>(Yhw0.left-dely) & Y<Yhw0.left;   %wall on left bank
                idzright = Y>Yhw0.right & Y<(Yhw0.right+dely); %wall on right bank
                idz = idzleft | idzright;
                grid.z(idz) = grid.z(idz)+2; %offset to define a wall crest level
            end            
        end
%%
        function [Y,Yhw0,Yhw,zFP0,zFP,dely] = constraintGrids(obj,grid)
            %grid arrays for the y-dimension, alongchannel high water
            %widths for left and right bank flood plain elevations. Returns
            %values for initial condition and current grid
            
            %grid for Y dimension
            grdobj = obj.RunParam.GD_GridProps;
            [x,y,~,dely] = getGridDimensions(grdobj);
            Y = repmat(y',length(x),1);

            %y-coordinates of centre-line (channel thalweg)
            yCL = y(1)+(y(end)-y(1))/2;   %scalar constant, assume symetric form about x-axis

            %grid of initial high water y co-ordinates
            xhw0 = obj.Channel.Data.Plan.Whw/2;     %initial widths
            Yhwo = repmat(xhw0',1,length(y));
            Yhw0.left = yCL-Yhwo;    
            Yhw0.right = yCL+Yhwo;

            %grid of high water y co-ordinates
            xhw = grid.Wz{1}/2;                     %current state
            Yhwi = repmat(xhw',1,length(y)); 
            Yhw.left = yCL-Yhwi;    
            Yhw.right = yCL+Yhwi;

            trnobj = obj.RunParam.CF_TransData;
            %elevation of initial flood plain
            hydobj0 = obj.Channel.Data.WaterLevels; %source form water levels
            zFP0 = hydobj0.zhw+trnobj.FPoffset;     %initial flood plain levels
            zFP0 = repmat(zFP0',1,length(y));

            %elevation of current flood plain based on current zhw
            hydobj = obj.RunParam.CF_HydroData;
            zFP = hydobj.zhw+trnobj.FPoffset;       %flood plain levels
            zFP = repmat(zFP',1,length(y));
        end
%%
        function grid = updateValley(obj,fgrid,isoffset)
            %combine channel form (fgrid) with valley form (obj.Valley.Data.Form)
            % fgrid - new channel grid to be added to valley grid
            % isoffset - flag to include offset if xM>0, model uses true,
            %            summary plots use false to generate initial form
            grid = fgrid;                  %copy channel grid to retain struct data
            vgrid = getGrid(obj.Valley,1); %valley grid (does not move)
            
            if isoffset && ~isempty(obj.Grid) && obj.Grid.xM>0
                %fgrid.x = fgrid.x+obj.Grid.xM; %moves grid by erosion distance
                %subgrid scale movement makes control of property estimates
                %complex. Only make coast when it has moved integer delx
                ixM = floor(grid.xM/(grid.x(2)-grid.x(1)))+1;
                fgrid.x = fgrid.x+fgrid.x(ixM);            
            end
            %combine channel and valley grids
            [X,Y] = ndgrid(vgrid.x,vgrid.y);
            zf = griddata(fgrid.x,fgrid.y,fgrid.z',X,Y);%ensure both use the same grid
            grid.x = vgrid.x;
            grid.y = vgrid.y;
            grid.z = max(zf,vgrid.z);                   %combined forms
        end
%%
        function grid = updateShoreline(obj,grid)
            %modify the shoreline if coast has retreated
            if ~isempty(obj.Grid) && obj.Grid.xM>0
                ixM = floor(grid.xM/(grid.x(2)-grid.x(1)));
                if ixM>0
                    subgrid = grid;
                    subgrid.x = grid.x(ixM+1:end);
                    subgrid.z = grid.z(ixM+1:end,:);
                    shore = setShoreline(obj.RunParam.CF_ShoreData,subgrid,false); %shore strip grid
                    grid.z(1:ixM,:) = shore.z(end-(ixM-1):end,:);
                end
                % %defined slope shoreline  
                % hydobj = obj.RunParam.CF_HydroData;
                % trnobj = obj.RunParam.CF_TransData;
                % zmx = hydobj.zhw(1)+trnobj.FPoffset;
                % zmn = min(obj.Valley.Data.Form.Z,[],'all');
                % zM = grid.z(ixM,:); %elevation at mouth
                % %add slope that varies from 1:20 at +zhw to 1:200 at the
                % %invert of the valley at x=0 ie zmn
                % slope = (zM-zmn)*(0.05-0.005)/(zmx-zmn)+0.005;
                % delz = zM-(grid.x(ixM)-grid.x(1:ixM))*slope;
                % grid.z(1:ixM,:) = min(grid.z(ixM+1,:),delz);
            end
        end
%%      
        function grid = updateMeander(obj,grid)
            %modify the grid if a meander is to be included
            ismeander = obj.RunParam.CF_TransData.isMeander;
            if ~isnan(ismeander) && ~isempty(grid.cline)
                if ismeander>0
                    %meander migrates with the coast
                    cline = grid.cline;
                    ixM = floor(grid.xM/(grid.x(2)-grid.x(1)));
                    if ixM>0
                        cline.y = [ones(ixM,1)*cline.y(1);cline.y(1:end-ixM)];
                        grid.cline = cline;
                    end
                    grid = gd_xy2sn(grid,cline,true,false);
                else      
                    %meander is fixed as initially defined
                    cline = grid.cline;
                    grid = gd_xy2sn(grid,cline,true,false);
                end
            end
        end
%%
        function obj = updateTidalLimit(obj)
            %find the new tidal limit and update integration distance indices
    %NB: water levels corrected for shoreline offset in cf_set_hydroprops
    %    and Lt is distance from x=0
            obj.dTrans.Lt = tidalLimit(obj);
            obj.RunParam.CF_HydroData.xTidalLimit = obj.dTrans.Lt; %used in models
            
            %upstream limit (tidal limit)
            intdist = obj.RunParam.CF_TransData.IntDist;
            if intdist>0           %user has defined integration distance
                subxu = find(obj.Grid.x>intdist,1,'first');
            else                   %use dynamic tidal limit
                subxu = find(obj.Grid.x>obj.dTrans.Lt,1,'first');
            end
            %downstream limit (mouth)
            subxl = 1;
            if ~isempty(obj.Grid) && obj.Grid.xM>0
                %record change as distance from x=0 and not channel mouth
                dx = obj.Grid.x(2)-obj.Grid.x(1);
                subxl = floor(obj.Grid.xM/dx)+1;
                if intdist>0
                    %increment upper to maintain specified integration distance
                    subxu = subxu+subxl-1; 
                end
            end
            obj.Grid.intidx = subxl:subxu;
        end
%%
        function Lt = tidalLimit(obj)
           %find tidal limit, where amplitude<0.05 or river bed level = HW level
           hydobj = obj.RunParam.CF_HydroData; %water levels at time t
           amp = hydobj.zhw-hydobj.zlw;
           idx = find(amp<0.05,1,'first');       %index of negligible tidal amplitude 
           thalweg = obj.Grid.z(:,ceil(size(obj.Grid.z,2)/2));
           hwl = hydobj.zhw;
           
           %with a river channel and shallow valley slope the HW may not
           %intersect the valley slope
           %figure;  plot(obj.Grid.x,thalweg,obj.Grid.x,hwl);
           Pxz =  InterX([obj.Grid.x';thalweg'],[obj.Grid.x';hwl']);
           if ~isempty(Pxz) && isempty(idx)      
               Lt = Pxz(1);                      %does not converge but intersects valley
           elseif isempty(Pxz) && ~isempty(idx)
               Lt = obj.Grid.x(idx);             %converges but does not intersect valley
           elseif ~isempty(Pxz) && ~isempty(idx)
               Lt = min(obj.Grid.x(idx),Pxz(1)); %converges and intersects valley
           else
               Lt = [];   %no solution found
           end
        end
%%
        function [hrv,Wrv] = riverDimensions(obj)
            %get the dimensions of the river
            hydobj = obj.RunParam.CF_HydroData;
            sedobj = obj.RunParam.CF_SediData;
            amp = (hydobj.zhw(1)-hydobj.zlw(1))/2;
            d50riv = sedobj.d50river;
            tauriv = sedobj.tauriver;
            rhos = obj.cns.rhos;
            rhow = obj.cns.rhow;
            Qr = hydobj.RiverDischarge;
            Le = hydobj.xTidalLimit;   %distance to tidal limit 
            Sr  = 2*amp/Le;

            [hrv,Wrv,~] = river_regime(Qr,Sr,d50riv,tauriv,rhos,rhow);
        end
%%
        function fparea = floodPlainArea(obj)
            %get the flood plain area based on elevations at HW+offset           
            fp_offset = obj.RunParam.CF_TransData.FPoffset;
            %model water level surfaces
            hydobj = obj.RunParam.CF_HydroData;
            zfp = hydobj.zhw+fp_offset;
            Zfp = repmat(zfp,1,length(obj.Grid.y));

            zf = obj.Grid.z;
            zf(zf<Zfp-0.05 | zf>Zfp+0.05) = NaN;
            ndx = sum(sum(~isnan(zf)));
            [~,~,delx,dely] = getGridDimensions(obj.RunParam.GD_GridProps);
            fparea = delx*dely*ndx;
        end
%%
        function [sedvol,watervol,conc] = sedimentFlux(obj)
            %get the sediment flux in a time step
            % conc is an optional output but is not currently used
            trnobj = obj.RunParam.CF_TransData;
            dt = obj.delta/obj.cns.y2s; 
            hydobj = obj.RunParam.CF_HydroData;
            slr = hydobj.dslr;     %slr in time step
            if isempty(obj.CSTparams) || ~isfield(obj.CSTparams,'Vhw')
                obj = setCSTparams(obj);
            end
            
            if trnobj.SedFlux>0
                %user defined import/export flux in m3/yr              
                sedvol = trnobj.SedFlux*dt;                %+ve volume for import
                watervol = obj.CSTparams.Shw*slr*dt;
            else
                %compute flux from form and hydraulics in m3/yr
                %parameters assigned in updateInputParameters
                inp.Volume = obj.CSTparams.Vhw;            %element volume at start of timestep (m^3)
                inp.SurfaceArea = obj.CSTparams.Shw;       %element surface area (m^2)
                inp.Prism = obj.CSTparams.Pr;              %tidal prism of channel (m^3)
                inp.RiverDischarge = hydobj.RiverDischarge;%river discharge (m^3/s) +ve downstream
                %sediment input data
                sedobj = obj.RunParam.CF_SediData;
                inp.TransportCoeff = sedobj.TransportCoeff;%transport coefficient n (3-5)
                inp.EqConc = sedobj.EqConcentration;       %equilibrium concentration (-)
                inp.RiverConc = sedobj.RiverConcentration; %river load imported by advection (-)
                inp.BedConc = sedobj.BedConcentration;     %concentration of bed (-)
                inp.y2s = obj.cns.y2s;                     %factor to convert from years to seconds                
                
                %equilibrium volume definition
                if sedobj.EqScaleCoeff==0
                    V0 = obj.Channel.Data.GrossProps.Vhw(1);
                    Pr0 = obj.Channel.Data.GrossProps.Pr(1);
                    inp.EqScaleCoeff = V0/Pr0; 
                    % inp.EqScaleCoeff = inp.Volume/inp.Prism; %used in v1 code  - makes no difference in linear case
                    inp.EqShapeCoeff = 1;
                else
                    inp.EqScaleCoeff = sedobj.EqScaleCoeff;
                    inp.EqShapeCoeff = sedobj.EqShapeCoeff;
                end     
                
                %vertical and horizontal rates of sediment exchange
                inp = sedExchanges(obj,inp); 
                
                %sediment flux with external environment
                [sedvol,watervol,conc] = get_sed_flux(inp,slr);
                %report sedvol as +ve volume for sediment import       
                sedvol = -sedvol*dt; watervol = watervol*dt;
            end
        end
%%
        function inp = sedExchanges(obj,inp)
            %set the vertical and horizontal reate sediment exchange
            hydobj = obj.RunParam.CF_HydroData;
            sedobj = obj.RunParam.CF_SediData;
            cn = obj.cns;                              %model constants 
            %vertical exchange   
            ws = settling_velocity(sedobj.SedimentSize,cn.g,cn.rhow,...
                                    cn.rhos,cn.visc,sedobj.EqDensity);                  
            inp.VerticalExchange = ws;                 %vertical exchange (m/s)

            %horizontal exchange
            if obj.Selection.wlflag==0 && ~isempty(hydobj.cstres)  
                %hydraulic results from CST model
                u = hydobj.cstres.U(1);                %velocity at mouth
                H = mean(hydobj.cstres.d);             %average hydraulic depth
            else
                %no hydraulic data
                u = 1;                                 %assumed velocity at mouth
                H = obj.CSTparams.Vhw/obj.CSTparams.Shw;%estimated hydraulic depth
            end
            A = obj.CSTparams.Am;                      %CSA at mouth
            D = u^2*H/ws;

            delx = u*hydobj.tidalperiod/4;             %tidal excursion length
            inp.HorizontalExchange = D*A/delx;         %horizontal exchange (m/s)
        end
%%
        function obj = setCSTparams(obj)
            %extract the CSTparams from the model GrossProps table
            %only called when initialising model
            gp = obj.Channel.Data.GrossProps(1,:);
            %used in CSTmodel
            obj.CSTparams.Wm = gp.Wm;     %mean tide width at mouth (m)
            obj.CSTparams.Lw = gp.Lw;     %width convergence length at MT (m)
            obj.CSTparams.Am = gp.Am;     %mean tide csa at mouth (m^2)
            obj.CSTparams.La = gp.La;     %csa convergenve length at MT (m)
            %used in get_sed_flux (when called in CF_TransModel)
            obj.CSTparams.Vhw = gp.Vhw;   %element volume (m^3)
            obj.CSTparams.Shw = gp.Shw;   %element surface area (m^2)
            obj.CSTparams.Pr = gp.Pr;     %tidal prism of channel (m^3) 
        end
%%
        function obj = updateInputParams(obj)
            %update the input properties based on the new offset
            % edX
            hydobj = obj.RunParam.CF_HydroData; 
            grdobj = obj.RunParam.GD_GridProps;
            grid = obj.Grid;
            
            %net change in position of "initial" mouth due transgression 
            %of coast and estuary NB - relative to translating channel form
            %and not x=0
            edX = obj.dTrans.estdX-obj.dTrans.cstdX; 
            
            %update parameters used in call to newWaterLevels/CSTmodel
            %should exist already for ckfa model but not for others  
            hyps = gd_channel_hypsometry(grid,hydobj,grdobj.histint,0);
            [w,csa,~] = gd_section_properties(grid,hydobj);
            gp = gd_gross_properties(grid,hydobj,hyps,w{2},csa{2});
            %write width and csa to command window
            time = obj.StepTime(end)/obj.cns.y2s;
            fprintf('T = %g yrs: Width: Wm %g m, CSA: Am %g m^2\n',time,gp.Wm,gp.Am)
            
            obj.CSTparams.Wm = gp.Wm;   %mean tide width at mouth (m)
            obj.CSTparams.Lw = gp.Lw;   %width convergence length at MT (m)
            obj.CSTparams.Am = gp.Am;   %mean tide csa at mouth (m^2)
            obj.CSTparams.La = gp.La;   %csa convergence length at MT (m)  
            %update parameters used in call to get_sed_flux
            obj.CSTparams.Vhw = gp.Vhw; %element volume (m^3)
            obj.CSTparams.Shw = gp.Shw; %element surface area (m^2)
            obj.CSTparams.Pr = gp.Pr;  %tidal prism of channel (m^3) 

            option = obj.Channel.Selection.modeltype;
            if strcmp(option,'Exponential') || strcmp(option,'Power')
                %movement of open coast and x index for offset from x=0
                thalweg = grid.z(:,ceil(size(grid.z,2)/2));
                zm = interp1(grid.x,thalweg,obj.Grid.xM,'linear','extrap');
                modsel = obj.Channel.Selection;
                if isfield(modsel,'channelform') && ...
                                    strcmp(modsel.channelform,'Rectangular')
                    %adjust hydraulic depth to effective thalweg depth
                    %zm=zlw-(nc+1)/nc)*hbar
                    nc = obj.Channel.RunParam.CF_FormData.ChannelShapeParam;
                    zlw = hydobj.zlw(1);        %low water at mouth
                    %hd = mc.Wlw/2/(nc+1), or %hd = nc.hc/(nc+1),
                    hc = zlw-zm;   %depth of channel at centre-line
                    zm = zlw-nc*hc/(nc+1); 
                end   

                %for the ckfa model the tidal limit and tidal amplitude 
                %are the input parameters that change and these are updated 
                %when CF_HydroData is updated in cf_set_hydroprops and
                %updateTidalLimit             
                classname = 'CF_FormData';
                Wu = obj.RunParam.(classname).HWmouthWidth;
                Wl = obj.RunParam.(classname).LWmouthWidth;
                
                ixM = floor(obj.Grid.xM/(grid.x(2)-grid.x(1)))+2; %+2 to avoid mouth
                Lwu = -getconvergencelength(grid.x(ixM:end),grid.Wz{1}(ixM:end));
                Lwl = -getconvergencelength(grid.x(ixM:end),grid.Wz{3}(ixM:end));
                
                %increase in mouth width due to landward transgression
                %new widths based on increment, edX, from previous time step
                Wu = (Wu-grid.Wrv)*exp(edX/Lwu)+grid.Wrv; 
                Wl = (Wl-grid.Wrv)*exp(edX/Lwl)+grid.Wrv;
                
                %update model run time parameters
                obj.RunParam.(classname).HWmouthWidth = Wu;
                obj.RunParam.(classname).LWmouthWidth = Wl;
                obj.RunParam.(classname).zMouthInvert = zm;
                
                if strcmp(option,'Exponential')
                    %update the convergence lengths
                    obj.RunParam.(classname).HWwidthELength = Lwu;
                    obj.RunParam.(classname).HWwidthELength = Lwl;
                end                
            end                
        end
        %%
        function setChannelWaterLevels(obj)
            wls = obj.Channel.Data.WaterLevels(1,:);
            obj.RunParam.CF_HydroData.zhw = wls.zhw'; %high water
            obj.RunParam.CF_HydroData.zmt = wls.zmt'; %mean tide level
            obj.RunParam.CF_HydroData.zlw = wls.zlw'; %low water
        end
%%
        function obj = setTransModel(obj,meta)
            %add the transgression results table to the model output
            tstep = obj.Data.Form.RowNames;
            tdsp = modelDSproperties(obj);
            tdst = dstable(obj.Trans,'RowNames',tstep,'DSproperties',tdsp); 
            tdst.Source = meta.source;
            tdst.MetaData = meta.data;
            tdst.Dimensions.delx = 1:4;  %dummy dimensions for 4 cases at increments of delX
            obj.Data.Transgression = tdst;
        end
 %% ------------------------------------------------------------------------
% Model input selection and initialisation
%--------------------------------------------------------------------------
        function [obj,isok] = selectInputModels(obj,mobj)
            %select Form and Valley models
            isok = true;
            muicat = mobj.Cases;
            ftxt = 'Select Form Model to use:';
            obj.Channel = selectCaseObj(muicat,[],{'CF_FormModel'},ftxt);
            if isempty(obj.Channel), isok = false; return; end
            vtxt = 'Select Valley Model to use:';
            obj.Valley = selectCaseObj(muicat,[],{'CF_ValleyModel','GD_GridImport'},vtxt);
            if isempty(obj.Valley), isok = false; end
            
            %Metadata of selection
            obj.MetaData = sprintf('Transgression model using "%s" for channel and "%s" for valley',...
                           obj.Channel.Data.Form.Description,obj.Valley.Data.Form.Description);
            
            %assign Channel RunParam to the new model
            fnames = fieldnames(obj.Channel.RunParam);
            for i=1:length(fnames)
                obj.RunParam.(fnames{i}) = copy(obj.Channel.RunParam.(fnames{i}));
            end
            %add the additional input data needed for transgression model
            obj.RunParam.RunProperties = copy(mobj.Inputs.RunProperties);
            obj.RunParam.CF_TransData = copy(mobj.Inputs.CF_TransData);
            if isfield(mobj.Inputs,'CF_ShoreData') && ~isfield(obj.RunParam,'CF_ShoreData')
                obj.RunParam.CF_ShoreData = copy(mobj.Inputs.CF_ShoreData);
            end
            %prompt user to modify
            obj = updateRunParameters(obj);
        end
%%
        function obj = updateRunParameters(obj)
            %allow user to modify parameters that do not influence initial 
            %form including: sea level rise, tidal range cycles, and the 
            %coefficients used in the sediment flux calculations
            titletxt = 'Update Parameters';
            promptxt = {'Rate of sea level rise (m/year)', ...
                        'Amplitude of Cycles (m)',...
                        'Equilibrium sediment density (kg/m^3)', ...
                        'Sediment load in river (kg/m^3)',...
                        'Transport coefficient (+/-n)',...                         
                        'Equilibrium scale coefficient (0=scale to initial)',...                          
                        'Equilibrium shape coefficient (-)'};
            defaults = {num2str(obj.RunParam.WaterLevels.SLRrate),...
                        num2str(obj.RunParam.WaterLevels.CycleAmp),...
                        num2str(obj.RunParam.CF_SediData.EqDensity),...
                        num2str(obj.RunParam.CF_SediData.RiverDensity),...
                        num2str(obj.RunParam.CF_SediData.TransportCoeff),...
                        num2str(obj.RunParam.CF_SediData.EqScaleCoeff),...
                        num2str(obj.RunParam.CF_SediData.EqShapeCoeff)};
            answer = inputdlg(promptxt,titletxt,1,defaults);
            if ~isempty(answer)
                obj.RunParam.WaterLevels.SLRrate = str2num(answer{1}); %#ok<ST2NM>
                obj.RunParam.WaterLevels.CycleAmp = str2num(answer{2}); %#ok<ST2NM>
                obj.RunParam.CF_SediData.EqDensity = str2double(answer{3});
                obj.RunParam.CF_SediData.RiverDensity = str2double(answer{4});
                obj.RunParam.CF_SediData.TransportCoeff = str2double(answer{5});
                obj.RunParam.CF_SediData.EqScaleCoeff = str2double(answer{6});
                obj.RunParam.CF_SediData.EqShapeCoeff = str2double(answer{7});
            end                       
        end
%% ------------------------------------------------------------------------
% Plotting methods
%--------------------------------------------------------------------------
        function changePlot(obj,slr)
            %generate a plot of change
            hf = figure('Name','Difference Plot','Units','normalized',...
                                                    'Tag','PlotFig');                                    
            hf.Position(1) = 0.1;
            hf.Position(3) = hf.Position(3)*2;
            p = uipanel('Parent',hf,'BorderType','none'); 
            p.Title = sprintf('Change plot, slr=%0.1g',slr);
            p.TitlePosition = 'centertop'; 
            p.FontSize = 12;
            p.FontWeight = 'bold';
            
            ax1 = subplot(1,2,1,'Parent',p);
            netchangePlot(obj,ax1,slr);
            
            ax2 = subplot(1,2,2,'Parent',p);
            rateofchangePlot(obj,ax2);

            %to replace the default background with a white background uncomment
            %used for figure plots for papers
            %hf.Color = [1,1,1];
            %ax1.Color = [1,1,1];
            %ax2.Color = [1,1,1];            
        end
%%
        function netchangePlot(obj,ax,slr)
            %plot the difference between the original and translated channel
            %origin of the plot is at the offset distance (amount translated) 
            % NB uses transient properties so only works in current session
            subx = obj.Grid.intidx; %index to subsample domain for integration
            xs = obj.Grid.x(subx);                
            zs = obj.Grid.z(subx,:);%current grid
            ys = obj.Grid.y;  
            initgrid = getGrid(obj.Channel);
            initgrid = updateValley(obj,initgrid,false);
            z0 = initgrid.z(subx,:);
            zdiff = zs-z0;
            
            yzhw = obj.Grid.Wz{1}(subx)/2;
            yzlw = obj.Grid.Wz{3}(subx)/2;
            zyz = ones(size(xs))*max(max(zdiff));
            adX = obj.Trans.estdX(end);  %cumulative estuary trangression
            cdX = obj.Trans.cstdX(end);  %cumulative coastal trangression

            %contourf(ax,xi,yi2,zgrd','LineStyle', 'none');
            surf(ax,xs,ys,zdiff','FaceColor','interp','EdgeColor', 'none');
            colormap(cmap_selection(20));
            caxis([-2*slr,2*slr])      %only show +/-2*slr change
            view(2);
            hold on
            plot3(xs,yzhw,zyz,'--k');  %high water line
            plot3(xs,-yzhw,zyz,'--k');
            plot3(xs,yzlw,zyz,'-.k');  %low water line
            plot3(xs,-yzlw,zyz,'-.k');
            hold off
            xlabel('Distance along channel (m)');
            ylabel('Width (m)');
            h_c = colorbar(ax);   
            h_c.Label.String = 'Change in elevation (m) for +/-2.slr';
            timetxt = num2str(obj.Grid.t);
            infotxt = sprintf('Dashed lines are HW/LW at T=%s yrs',timetxt);
            tltxt = sprintf('dV=0 for: SLR = %0.1g m, Transgression = %d m, Coast erosion = %d m\n%s',...
                                        slr,round(adX),round(cdX),infotxt);
            title(tltxt,'FontSize',10);
        end  
%%
        function rateofchangePlot(obj,ax)
            %plot the volume difference for different distanes of transgression
            % NB uses transient properties so only works in current session
            delX = obj.Trans.delX(end);
            dx = [0,delX/2,delX,delX+delX/2];  %x intervals used            
            dV = obj.Trans.vdiffx(end,:);      %volume change for dx
            trnobj = obj.RunParam.CF_TransData;
            sedvol = obj.Trans.sedVol(end);
            estdX = obj.Trans.estdX(end);
            
            plot(ax,dx,dV,'-k','LineWidth',1)
            ax.YGrid = 'on';
            xlabel('Landward transgression (m)');
            ylabel('Volume change (m^3)');
            %add explanatory text
            mxX = ax.XLim(2); mxY = ax.YLim(2);
            hold on
            plot([0,mxX],[0,0],'--r');  %'Color',[0.96,0.96,0.96]
            plot([0,estdX],[sedvol,sedvol],'--b'); %transgression point
            plot([estdX,estdX],[0,sedvol],'--b'); %transgression point
            hold off
            txt1 = 'Channel only';
            if trnobj.inclFloodPlain
                txt1 = 'Include flood plain';
            end
            txt2 = 'No bed constraint';
            if trnobj.inclGeoConstraint
                txt2 = 'Including bed constraints';
            end
            txt3 = 'No HW constraint';
            if trnobj.inclHWConstraint
                txt3 = 'Constrained at HW';
            end
            title(sprintf('%s; %s; %s\n',txt1,txt2,txt3),'FontSize',10);
            text(mxX,mxY/10,'Sediment Import  ','HorizontalAlignment','right');
            text(mxX,-mxY/10,'Sediment Export  ','HorizontalAlignment','right');
        end
%%
        function crossectionPlot(obj,slr)
            %plot cross-sections along length of channel at start and end of run
            % NB uses transient properties so only works in current session
            xi = obj.Grid.x;  
            yi = obj.Grid.y;
            zi = obj.Grid.z;    %current grid
            initgrid = getGrid(obj.Channel);
            initgrid = updateValley(obj,initgrid,false);
            z0 = initgrid.z;
            if obj.Grid.ishead  %orientation of x-axis, x=0 is nearest the mouth if ishead=false
                zi = flipud(zi);                          
                z0 = flipud(z0);  
            end

            wl0 = obj.Channel.Data.WaterLevels.zmt(1);
            wli = obj.RunParam.CF_HydroData.zmt(1);

            figure('Name','XS Plot','Units','normalized','Tag','PlotFig');                                                    
            ax = axes('Tag','PropertyPlot');
            noxi=length(xi);
            ix0 = find(xi>=obj.Grid.xM-eps,1,'first'); 
            noxi = noxi-ix0;
            nx1=ix0;  nx2=ceil(ix0+0.1*noxi);  nx3=ceil(ix0+0.2*noxi);
 
            green = mcolor('green');
            plot(ax,yi,z0(nx1,:),'-r','LineWidth',0.6);
            hold on
            plot(ax,yi,zi(nx1,:),'--r','LineWidth',0.6);
            
            plot(ax,yi,z0(nx2,:),'-b','LineWidth',0.58);
            plot(ax,yi,zi(nx2,:),'--b','LineWidth',0.58);
            
            plot(ax,yi,z0(nx3,:),'-','Color',green,'LineWidth',0.56);
            plot(ax,yi,zi(nx3,:),'--','Color',green,'LineWidth',0.56);
            
            %add water levels at mouth
            plot(xlim, wl0*[1 1],'-','Color',[0.7,0.7,0.7]);
            plot(xlim, wli*[1 1],'--','Color',[0.7,0.7,0.7]);
            
            hold off
            xlabel('Width (m)'); 
            ylabel('Change in level (m)');
            hL=legend('0pre','0post','0.1pre','0.1post','0.2pre','0.2post','Location','SouthEast');
            set(hL, 'Color', 'none');
            casedesc = sprintf('Cross-sections difference plot, slr=%0.1g',slr);
            title(casedesc,'FontWeight','normal','FontSize',10);            
        end
%%
        function thalwegPlot(obj,slr)
            %plot the thalweg at the start and end of run
            % NB uses transient properties so only works in current session
            xi = obj.Grid.x;  
            zi = obj.Grid.z(:,ceil(size(obj.Grid.z,2)/2));
            
            initgrid = getGrid(obj.Channel);
            initgrid = updateValley(obj,initgrid,false);
            z0 = initgrid.z(:,ceil(size(initgrid.z,2)/2));
            if obj.Grid.ishead  %orientation of x-axis, x=0 is nearest the mouth if ishead=false
                zi = flipud(zi);                       
                z0 = flipud(z0);  
            end
            
            wl0 = obj.Channel.Data.WaterLevels.zmt(1);
            wli = obj.RunParam.CF_HydroData.zmt(1);
            
            figure('Name','CL Plot','Units','normalized','Tag','PlotFig');       
            ax = axes('Tag','PropertyPlot');
            plot(ax,xi,z0(:,1),'-k');            
            hold on
            plot(ax,xi,zi(:,1),'--k');
            
            %add water levels at mouth
            plot(xlim, wl0*[1 1],'-','Color',[0.7,0.7,0.7]);
            plot(xlim, wli*[1 1],'--','Color',[0.7,0.7,0.7]);  
            
            hold off
            xlabel('<= mouth       Distance along channel (m)       head =>'); 
            ylabel('Elevation (mAD)');
            hL=legend('Initial','Final','Location','SouthEast');
            set(hL, 'Color', 'none');
            casedesc = sprintf('Centre-line difference plot, slr=%0.1g',slr);
            title(casedesc,'FontWeight','normal','FontSize',10); 
        end
%%
        function summaryGridPlots(obj,slr)             
            %plot initial grid, new grid and difference
            % NB uses transient properties so only works in current session
            intigrid = getGrid(obj.Channel,1);            
            F0 = updateValley(obj,intigrid,false);
            F = obj.Grid;
            zdiff = F.z-F0.z;
            start = obj.RunParam.RunProperties.StartYear;
            
            figure('Name','Grid Plot','Units','normalized','Tag','PlotFig'); 
            subplot(3,1,1)
            surf(F0.x,F0.y,F0.z','FaceColor','interp','EdgeColor', 'none');    
            
            view(2); 
            h1 = colorbar;  
            h1.Label.String = 'Elevation (mAD)';
            title(sprintf('Intiial form, T = %s',num2str(start)))
            subplot(3,1,2)
            surf(F.x,F.y,F.z','FaceColor','interp','EdgeColor', 'none');
            view(2); 
            h2 = colorbar; 
            h2.Label.String = 'Elevation (mAD)';
            title(sprintf('Form at time = %s',num2str(F.t)))
            s3 = subplot(3,1,3);
            surf(F.x,F.y,zdiff','FaceColor','interp','EdgeColor', 'none');            
            view(2); 
            colormap(s3,cmap_selection(20));
            caxis([-2*slr,2*slr])      %only show +/-2*slr change
            h3 = colorbar; 
            h3.Label.String = 'Change in elevation (m)';
            title('Difference plot (showing +/-2.slr)')
            casedesc = sprintf('Summary of change, slr=%0.1g',slr);
            sgtitle(casedesc,'FontWeight','normal','FontSize',10); 
        end
%% ------------------------------------------------------------------------
% Ouput property definition for the Transgression specific data
%--------------------------------------------------------------------------
        function dsp = modelDSproperties(~)
            %define a dsproperties struct and add the model metadata
            dsp = struct('Variables',[],'Row',[],'Dimensions',[]);
            %define each variable to be included in the data table and any
            %information about the dimensions. dstable Row and Dimensions can
            %accept most data types but the values in each vector must be unique

            %struct entries are cell arrays and can be column or row vectors
            %static ouput (mean tide values)
            % delX - unadjusted transgression distance
            % estdX - adjusted transgression distance (inc sediment flux)
            % cstdX - open coast transgression distance
            % Lt - distance to tidal limit (to record change - model uses CF_HydroData.xTidalLimit)
            % FPA - flood plain area
            % waterVol - water volume change due to changes in slr and tidal range
            % sedVol - sediment flux (+ve=sediment import)            
            % vdiffx - volume difference for [0,delX/2,delX,3delX/2]
            
            dsp.Variables = struct(...
                        'Name',{'delX','estdX','cstdX','dSLR','Lt','FPA',...
                                            'waterVol','sedVol','vdiffx'},...
                        'Description',{'Unadjusted transgression distance',...
                                       'Adjusted transgression distance',...
                                       'Open coast transgression distance',...
                                       'Cumulative sea level rise',...
                                       'Distance to tidal limit',...
                                       'Flood plain area',...
                                       'Water volume change',...
                                       'Sediment flux (+ve=sediment import)',...
                                       'Volume change for [0,dx/2,dx,3dx/2]'},...
                        'Unit',{'m','m','m','m','m','m^2','m^3','m^3','m^3'},...
                        'Label',{'Distance (m)','Distance (m)','Distance (m)',...
                                 'Sea level rise (m)','Distance (m)',...
                                 'Area (m^2)','Volume (m^3)',...
                                 'Volume (m^3)','Volume (m^3)'},...                                 
                        'QCflag',repmat({'model'},1,9));
            dsp.Row = struct(...
                        'Name',{'Time'},...
                        'Description',{'Time'},...
                        'Unit',{'y'},...
                        'Label',{'Time (yr)'},...
                        'Format',{'y'});       
           dsp.Dimensions = struct('Name',{'delx'},'Description',{'Increments of delX'},...
                               'Unit',{'-'},'Label',{'Intervals for [0,dx/2,dx,3dx/2]'},'Format',{''});    
        end
    end
end

