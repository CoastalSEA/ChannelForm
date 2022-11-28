classdef CF_TransModel < FGDinterface  
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
        Selection %struct for plan, channel and intertidal form selection 
                  %along-channel water level model and modifications flags
                  %NB: needed for calls to models (e.g. cf_exp_models)
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
                       % t - timestep (years)
                       % ishead - %orientation of x-axis true if head
                       % xM - distance to mouth from grid origin
                       % cline - %x,y coordinates of meander centre-line  
                       % note: irow, desc and metadata of grid stuct not used
        TranProp       %additional grid properties at timestep, t, struct with
                       % Wz - table of column vectors for width at hw,mt,lw
                       % zdiff - difference over a timestep
                       % intidx - x-axis indices from mouth to tidal limit  
                       % Rv - struct of river regime properties Hr, Wr, Ar
        CLatT          %struct array for x and y coordinates of meander centreline at each time step
        zGrid          %z grids at each sampled time step  
        zWL            %alongchannel water levels at each sampled time step
        tPlan          %plan form, Whw,Wmt,Wlw, at each sampled time step
        dTrans         %incremental changes in timestep, t, includes:
                       % delX - unadjusted transgression distance
                       % estdX - adjusted transgression distance (inc sediment flux)
                       % cstdX - open coast transgression distance
                       % SLR - increase in mean sea level at mouth (m)
                       % Lt - distance to tidal limit (to record change - model uses CF_HydroData.xTidalLimit)
                       % vdiffx - volume difference for [0,delX/2,delX,3delX/2]
                       % sedVol - sediment flux (+ve=sediment import)
                       % FPA - flood plain area
        Trans          %struct used to store transgression output at each time step                       
                       % dTrans struct variables at each sampled time step
                       %        delX,estdX,cstdX and sedVol values of 
                       %        dTrans at each timestep are converted
                       %        to cumulative values at end of run
        cns            %struct of default model constants
    end
    
    properties (Dependent, SetAccess=private)
        zMouthInvert        %zm - thalweg bed level at mouth to zero datum (m)
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
            msg = 'Solution not found - run terminated';
            ok = InitialiseModel(obj,mobj); %uses inputs assigned to RunParam
            if ok<1, warndlg(msg); return; end
            
            %run model
            if isnan(obj.RunParam.CF_TransData.SedFlux)
                ok = kinematic_model(obj);
            else
                ok = dynamic_model(obj);
            end
            if ok<1, warndlg(msg); return; end
            %write end of run 'input parameters' to command window
            timestepInput(obj); 
            
            %cumulative values for first 4 variables in Trans output table
            %order is: 'delX','estdX','cstdX','SLR','Lt','FPA','sedVol','vdiffx'                                                
            obj.Trans(:,1:4) = varfun(@cumsum,obj.Trans(:,1:4)); 
            obj.Trans.waterVol = cumsum(obj.Trans.waterVol,1);
            obj.Trans.sedVol = cumsum(obj.Trans.sedVol,1);
            obj.Trans.vdiffx = cumsum(obj.Trans.vdiffx,1);
            %plots of run
            [~,slr,~] = netChangeWL(obj.RunParam.WaterLevels,obj);
            summaryGridPlots(obj,slr);
            crossectionPlot(obj,slr);
            changePlot(obj,slr);         
            thalwegPlot(obj,slr);

            %now assign results to object properties  
            mtime = years(obj.StepTime/obj.cns.y2s);
            xM0 = obj.Valley.Data.Grid.UserData.xM;
            dims = struct('x',obj.Grid.x,'y',obj.Grid.y,'t',mtime,...
                          'ishead',false,'xM',xM0+obj.Trans.cstdX,...
                          'Lt',obj.Trans.Lt,'Rv',obj.TranProp.Rv,...
                          'cline',obj.CLatT);
%--------------------------------------------------------------------------
% Assign model output to dstable using the defined dsproperties meta-data
%--------------------------------------------------------------------------                   
            %assign metadata about model and save grid
            meta.source = obj.MetaData;
            meta.data = obj.Channel.MetaData;
            obj = setGrid(obj,{obj.zGrid},dims,meta);
%--------------------------------------------------------------------------
% Add property dstables in function GDinterface.setFormProperties
%--------------------------------------------------------------------------  
            %zwl and Wz are empty and resolved in setProperties
            %limits=0 to use grid to determine hypsometry limits        
            %histint = obj.RunParam.GD_GridProps.histint;
            obj = setModelFormProps(obj);  
            setTransModel(obj,meta);
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
%% ------------------------------------------------------------------------
% Methods to initialise models and run kinematic model
%--------------------------------------------------------------------------
    methods (Access=private)
        function ok = InitialiseModel(obj,mobj)
            %initialise ChannelForm properties and run parameters
            rnpobj = obj.RunParam.RunProperties;
            %model constants
            obj.cns = getConstantStruct(mobj.Constants);
            %initialise time step paramters
            tstep = rnpobj.TimeStep; 
            obj.delta = tstep*obj.cns.y2s;              %time step in seconds
            obj.DateTime = rnpobj.StartYear*obj.cns.y2s;%time elapsed from Year 0 in seconds
            obj.RunSteps = rnpobj.NumSteps; %not needed unless stability check added

            %transgression time step table for 
            % delX - unadjusted transgression distance
            % estdX - adjusted transgression distance (inc sediment flux)
            % cstdX - open coast transgression distance
            % SLR - sea level rise in time step
            % Lt - distance to tidal limit (to record change - model uses CF_HydroData.xTidalLimit)
            % FPA - flood plain area
            % waterVol - water volume change due to changes in slr and tidal range
            % sedVol - sediment flux (+ve=sediment import)            
            % vdiffx - volume difference for [0,delX/2,delX,3delX/2]
            varnames = {'delX','estdX','cstdX','SLR','Lt','FPA',...
                                           'waterVol','sedVol','vdiffx'};                                                        
            vars = {0,0,0,0,0,0,0,0,[0,0,0,0]};                    
            obj.dTrans = table(vars{:},'VariableNames',varnames);
            
            %transgression summary ouput table - same variables as dTrans
            obj.Trans = array2table(zeros(0,9));
            obj.Trans.Properties.VariableNames = varnames; 

            %initialise model grid at timestep, t=0
            obj.Grid = getGrid(obj.Channel,1);
            %in case shoreline present reset xM=0 (shore is reassigned in
            %updateModelGrid).
            obj.Grid.xM = 0;
            %model selection options(wlflag,modeltype,etc)                 
            obj.Selection = obj.Channel.Selection;
            
            %use current settings in WaterLevels and RunProperties to
            %initialise hydraulic parameters (single values). 
            setTransHydroProps(obj.RunParam.CF_HydroData,obj); 
            
            %use selected channel and valley forms to create initial grid
            obj.TranProp = struct('Wz',[],'zdiff',[],'intidx',[],'Rv',[]);
            %initialise the dimensions of the river for the initial river discharge
            %assumed constant as the channel migrates
            obj.TranProp.Rv = obj.Grid.Rv;

            %set parameters used in CSTmodel and sed_flux model
            setCSTparams(obj);

            %to ensure water levels are correctly assigned, recreate the
            %form based on the selected model. Ensures consistency and
            %initialises transient properties in obj.RunParam.CF_HydroData 
            ok = updateModelGrid(obj);
            if ok<1, return; end

            %write data for initial time step (t=0)
            PostTimeStep(obj); 
        end
%%
        function ok = kinematic_model(obj)
            %compute transgression without time stepping
            %write initial conditions to command window
            timestepInput(obj)
            
            obj.Time = obj.RunSteps*obj.delta; %time period in seconds
            obj.DateTime = obj.DateTime+obj.Time;
            obj.Grid.t = obj.DateTime/obj.cns.y2s;
            %update water levels at boundary
            newWaterLevels(obj.RunParam.CF_HydroData,obj);
            %net change in MSL over the simulation period
            [~,obj.dTrans.SLR] = netChangeWL(obj.RunParam.WaterLevels,obj);
            %change in water volume over run duration
            obj.dTrans.waterVol =  obj.CSTparams.Shw*obj.dTrans.SLR;
            
            %get transgression distance for the combined form
            hd = setdialog('Processing transgression');
            pause(0.1)
            obj = getTransDist(obj);
            obj.dTrans.estdX = obj.dTrans.delX;    %ignore sediment flux
            %coastal transgression 
            trnobj = obj.RunParam.CF_TransData;
            obj.dTrans.cstdX = trnobj.BruunRatio*obj.dTrans.SLR;
            
            %update form input parameters
            updateFormParams(obj);
            
            %update grid based on transgression
            ok = updateModelGrid(obj);
            close(hd);
            if ok<1, return; end
            %save results
            PostTimeStep(obj);
        end
%% ------------------------------------------------------------------------
% Methods for dynamic model with a time stepping loop
%--------------------------------------------------------------------------
        function ok = dynamic_model(obj)
            %run time stepping for the quasi-dynamic transgression model
            ok = 1;
            msg = sprintf('ChannelForm processing, please wait');
            hw = waitbar(0,msg);
            for jt = 1:obj.RunSteps
                InitTimeStep(obj,jt)
                ok = RunTimeStep(obj);
                if ok <1  || ~isvalid(hw)
                    close(hw);
                    return; 
                end
                PostTimeStep(obj);
                %to report time step during run use the following
                msg = sprintf('ChannelForm processing, step %d of %d',...
                                                         jt,obj.RunSteps);
                waitbar(jt/obj.RunSteps,hw,msg);                
            end
            close(hw);  
        end
%%
        function InitTimeStep(obj,jt)
            %initialise model parameters for next time step
            obj.iStep = jt;
            obj.Time = jt*obj.delta;
            obj.DateTime = obj.DateTime+obj.delta;
            obj.Grid.t = obj.DateTime/obj.cns.y2s;

            %update flux model parameters (call be fore updating waterlevels)
            updateFluxParams(obj);
            
            %update water levels at boundary (mouth) - NB: this method
            %is in CF_HydroData and NOT WaterLevels
            newWaterLevels(obj.RunParam.CF_HydroData,obj);
            %set slr increment for timestep
%             dt = obj.delta/obj.cns.y2s;
            obj.dTrans.SLR = obj.RunParam.CF_HydroData.slr; %change in msl over a time step
            %could also update river discharge if required
            % Qr = source of river discharge time series;
            % obj.RunParam.CF_HydroData.RiverDischarge = Qr;

            %sediment flux for new time step, jt, based on form at jt-1
            %needs slr for current time step, hence after waterlevel update
            [obj.dTrans.sedVol,obj.dTrans.waterVol] = sedimentFlux(obj); 

            %write selected parameters to command window
            timestepInput(obj)
        end
%%
        function ok = RunTimeStep(obj)
            %run model for the time step jt
            trnobj = obj.RunParam.CF_TransData;            
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
            %shoreline change on open coast, ie coastal transgression 
            obj.dTrans.cstdX = trnobj.BruunRatio*obj.dTrans.SLR;

            %update form input parameters
            updateFormParams(obj);
                    
            %update grid based on transgression
            ok = updateModelGrid(obj);

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
            % figure; mesh(obj.Grid.z);
            obj.zGrid = [obj.zGrid;reshape(obj.Grid.z,1,sz{:})];
            %add meander centreline coordinates
            obj.CLatT = [obj.CLatT,obj.Grid.cline];
            %add plan form description
            obj.tPlan = [obj.tPlan;obj.TranProp.Wz];
            %add alongchannel water levels
            hydobj = obj.RunParam.CF_HydroData;
            zwli = table(hydobj.zhw,hydobj.zmt,hydobj.zlw,...
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
%             slr = obj.RunParam.CF_HydroData.dslr*obj.delta/obj.cns.y2s;
            hc = obj.RunParam.CF_FormData.MTmouthDepth; %depth at mouth - *****DOES NOT work for CKFA model
            adist = Lt/hc*obj.dTrans.SLR;               %used as initial guess
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
            obj.TranProp.zdiff = zdiff;
        end
%%
        function [vdiff,zdiff] = getZdiff(obj,delX)
            %compute the difference between the surfaces for given distance 
            %of translation, delX. Uses current grid+dhw to define new grid             
            F = obj.Grid;
            [~,ixM] = gd_basin_indices(obj.Grid); %nearest grid point
            %remove shore whilst gridding differneces
            F.z = F.z(ixM:end,:);
            F.x = F.x(ixM:end);
            dx = abs(F.x(2)-F.x(1));
            dy = abs(F.y(2)-F.y(1));
            %add change in high water over a time step
            z2 = F.z+obj.RunParam.CF_HydroData.dhw;
            x2 = F.x+delX;    
            [X,Y] = ndgrid(F.x,F.y);   
            z1 = F.z;
            z2 = griddata(x2,F.y,z2',X,Y);    
            zdiff = (z2-z1);
            
            %restore shore with NaN values
            zdiff = [NaN(ixM-1,length(F.y));zdiff];
            %apply limits to domain and any vertical constraints
            zdiff = applyConstraints(obj,zdiff);
            
            %subsample grid to integration length and find volume change
            subx = obj.TranProp.intidx;  %indices from mouth to tidal limit
            %as coast erodes ixM=subx(1) is rounded to the nearest grid
            %point. When xM<dx/2 increment the index to avoid grid cells 
            %outside channel
            if rem(obj.Grid.xM,dx)>0 && rem(obj.Grid.xM,dx)<dx/2
                subx = subx(2:end);
            end
            dz = zdiff(subx,:);
            dz(isnan(dz)) = 0;
            vdiff = trapz(trapz(dz))*dx*dy;
        end
%%
        function zdiff = applyConstraints(obj,zdiff)
            %apply limits to domain and any vertical constraints
            %control lateral spatial extent
            trnobj = obj.RunParam.CF_TransData;                        
            [Y,Yhw0,Yhw,~,idV,dely] = constraintGrids(obj,obj.Grid);
            
            if trnobj.inclHWConstraint       %sea wall fixes HW line                
                zdiff(Y<Yhw0.left) = NaN;    %left side of initial hw
                zdiff(Y>Yhw0.right) = NaN;   %right side of initial hw
            elseif trnobj.inclFloodPlain     %open flood plain
                %flood plain surface that determines intersection with valley slope
                zdiff(idV) = NaN; %use combined channel+valley to HW limit        
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
                xM = obj.Valley.Data.Grid.UserData.xM;
                for i=1:ids
                    x = obj.Grid.x;
                    idx = x>xM+trnobj.StConstraints(i) & x<xM+trnobj.NdConstraints(i) & zdiff<0;
                    zdiff(idx)=0;  %remove any erosion over constrained bed
                end
            end
        end
%%
        function delX = interp_delX(obj,adist)
            %interpolate end values to find delX
            % zf = squeeze(obj.Channel.Data.Grid.Z);
            % zv = squeeze(obj.Valley.Data.Grid.Z);
            % figure; 
            % subplot(2,1,1)
            % contourf(zf')
            % subplot(2,1,2)
            % contourf(zv');
            
%             dVpos = getZdiff(obj,0);
%             sgn = 1;
%             if dVpos<0
%                 sgn = -1;
%             end
%             dVneg = dVpos; 
%             count = 2;           
%             while sgn*dVneg>0
%                 %ensure the largest delX value gives a negative vdiffx
%                 newdel = sgn*count*adist;
%                 dVneg = getZdiff(obj,newdel);
%                 count = count+1;
%             end
%             delX = dVpos/((dVpos-dVneg)/newdel);

            dV0 = getZdiff(obj,0);      %diffference with zero translation
            dV1 = getZdiff(obj,adist);  %guestimate with dx=adist
            sgnx = +1;
            if (dV0>0 && dV0<dV1) || (dV0<0 && dV0>dV1) 
                sgnx = -sgnx;           %going in wrong direction
            end
            count = 2;
            
            if (dV0>0 && dV1<0) || (dV0<0 && dV1>0)
                %opposited signs so solution found for +ve dx
            elseif dV0>0                %increase in volume with dx=0
                while dV1>0
                    newdel = sgnx*count*adist;
                    dV1 = getZdiff(obj,newdel);
                    if dV0<dV1, sgnx = -sgnx; end %direction should be established but check
                    count = count+1;
                end
            elseif dV0<0                %decrease in volume with dx=0
                while dV1<0
                    newdel = sgnx*count*adist;
                    dV1 = getZdiff(obj,newdel);
                    if dV0>dV1, sgnx = -sgnx; end %direction should be established but check
                    count = count+1;
                end
            else                        %dV0=0
                dV1 = 0;                %implies no change in sea level
                newdel = 1;             %hence no transgression
            end

            delX = dV0/((dV0-dV1)/newdel);    
        end
%%
        function ok = updateModelGrid(obj)
            %update grid based on transgression
            chgrid = getNewChannel(obj);             %new channel form and updates
            if isempty(chgrid.x), ok = 0; return; end%water levels in obj.RunParam.CF_HydroData
            chgrid = updateMeander(obj,chgrid);      %apply meander if included    
            chgrid = updateGeoConstraints(obj,chgrid);%adjust channel for any non-erodible geology
            chgrid = updateFloodPlain(obj,chgrid);   %adjust flood plain in channel model for any constraints
            newgrid = updateValley(obj,chgrid,true); %new combined channel+valley form            
            obj.Grid = updateShoreline(obj,newgrid); %shoreline adjusted 
            obj = cf_offset_wls(obj,true);           %translate wls, ??? length of wl vectors
            obj.dTrans.FPA = floodPlainArea(obj);    %modified flood plain area
            %get the distance to tidal limit and indices to integrate over
            updateTidalLimit(obj); 
            ok = 1;
        end
%%
        function grid = getNewChannel(obj)
            %generate a new channel based on updated input parameters
            %with no offset (ie mouth is at x=0)
            grid = obj.Grid;            
            option = obj.Selection.modeltype;
            %new channel grid and water levels in obj.RunParam.CF_HydroData
            %is also updated in the function called
            switch option
                case 'Exponential'                    
                    %model selection assigned to obj.Selection struct
                    [grid.x,grid.y,grid.z,Wz] = cf_exp_models(obj);
                case 'Power'
                    [grid.x,grid.y,grid.z,Wz] = cf_pow_model(obj);
                case 'CKFA'
                    [grid.x,grid.y,grid.z,Wz] = ckfa_form_model(obj);
            end   
            
            if isempty(grid.x)
                return;
            elseif isrow(grid.x)
                grid.x = grid.x'; 
            end  
            
            %update plan widths in TranProp
            obj.TranProp.Wz = Wz;
            
            %if set apply the modifications defined in CF_ModsData
            if obj.Selection.incmods
                modobj = obj.RunParam.CF_ModsData;
                grid.z = squeeze(setMorphMods(modobj,grid));
            end      
            grid.metadata = obj.Channel.Data.Grid.MetaData;       
        end
%%
        function grid = updateGeoConstraints(obj,grid)
            %apply any vertical constraints
            trnobj = obj.RunParam.CF_TransData; 
            if trnobj.inclGeoConstraint && ~isempty(trnobj.StConstraints)
                %find areas of erosion relative to initial grid
                ids = length(trnobj.StConstraints);
                z0 = squeeze(obj.Channel.Data.Grid.Z);  %source channel grid
                for i=1:ids
                    %constrains erosion relative to initial channel
                    %form over the constraint distances defined
                    x = grid.x;                    
                    %index is for channel grid
                    idx = x>trnobj.StConstraints(i) & ...
                          x<trnobj.NdConstraints(i) & grid.z<z0;
                    grid.z(idx)=z0(idx);  %remove any erosion over constrained bed
                end
            end
        end
%%
        function grid = updateFloodPlain(obj,grid)
            %apply any constraints at high water, and/or exclude flood plain
            %note: geological constraints are applied when creating fgrid
            trnobj = obj.RunParam.CF_TransData;
            [Y,Yhw0,Yhw,zFP0,idV,dely] = constraintGrids(obj,grid);

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
            
            if obj.RunParam.CF_TransData.MaxFPwidth>0
                %impose a maximum width on the valley flood plain
                %zFP is the flood plain elevation along the channel and the
                %valley dtm outside of the channel
                %idv = grid.z>zFP;       %model grid above flood plain
                grid.z(idV) = NaN; %apply flood plain/valley values
            end
        end
%%
        function [Y,Yhw0,Yhw,zFP0,idV,dely] = constraintGrids(obj,grid)
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
            xhw = obj.TranProp.Wz.Whw/2;            %current state
            Yhwi = repmat(xhw',1,length(y)); 
            Yhw.left = yCL-Yhwi;    
            Yhw.right = yCL+Yhwi;

            trnobj = obj.RunParam.CF_TransData;
            %elevation of initial flood plain
            hydobj0 = obj.Channel.Data.WaterLevels; %source form water levels
            zFP0 = hydobj0.zhw+trnobj.FPoffset;     %initial flood plain levels
            zFP0 = repmat(zFP0',1,length(y));       %returns array with size(z)

            %elevation of current flood plain based on current zhw
            hydobj = obj.RunParam.CF_HydroData;
            zFP = hydobj.zhw+trnobj.FPoffset;       %flood plain levels
            if isscalar(zFP)
                zFP = repmat(zFP,length(x),length(y));
            else
                zFP = repmat(zFP',1,length(y));     %returns array with size(z)
            end
            
            idV = grid.z>zFP;                       %valley above flood plain  
            if obj.RunParam.CF_TransData.MaxFPwidth>0
                %impose a maximum width on the valley flood plain
                %used when modelling Holocene and starting offshore of
                %present day coastline/valley formation to limit extent of
                %sediment demand on the flood plain
                zV = squeeze(obj.Valley.Data.Grid.Z);
                vb = obj.RunParam.CF_TransData.MaxFPwidth/2;  %limit valley;
                idv = ~(Y>=yCL-vb & Y<=yCL+vb); %indices for valley outwith limit
                idV = idv | idV;
                %zFP0(idV) = zV(idV);            %restore valley levels outside limit
            end
        end
%%
        function newgrid = updateValley(obj,fgrid,isoffset)
            %combine channel form (fgrid) with valley form (obj.Valley.Data.Grid)
            % fgrid - new channel grid to be added to valley grid
            % isoffset - flag to include offset if xM>0, model uses true,
            %            summary plots use false to generate initial form
            
            %apply an offset equivalent to any coastal erosion/accretion
            fgrid.xM = sum(obj.Trans.cstdX)+obj.dTrans.cstdX;
            %valley grid (does not move) but may be offset if shore present
            vgrid = getGrid(obj.Valley,1);     
            %copy channel grid to retain struct data
            newgrid = fgrid;               
            newgrid.xM = vgrid.xM+fgrid.xM;  %combined offset of shore & erosion/accretion

            if isoffset && (vgrid.xM>0 || fgrid.xM>0)
                %intial grid may include a shore so find the offset
                %relative to the initial coast
                [~,~,delx] = getGridDimensions(obj.RunParam.GD_GridProps);
                idM = round((vgrid.xM+fgrid.xM)/delx)+1;
                fgrid.x = fgrid.x+fgrid.x(idM);
                Wz = obj.TranProp.Wz;
                %pad vectors for shore and remove any landward point
                %outside grid
                varfunc = @(x) [NaN(1,idM-1),x(1:length(vgrid.x)-idM+1)];
                obj.TranProp.Wz = varfun(varfunc,Wz);
                obj.TranProp.Wz.Properties.VariableNames = Wz.Properties.VariableNames;
            end
            
            %combine channel and valley grids
            [X,Y] = ndgrid(vgrid.x,vgrid.y);
            zf = griddata(fgrid.x,fgrid.y,fgrid.z',X,Y);%ensure both use the same grid
            newgrid.x = vgrid.x;
            newgrid.y = vgrid.y;
            newgrid.z = max(zf,vgrid.z);                   %combined forms
        end
%%
        function grid = updateShoreline(obj,grid)
            %modify the shoreline if coast has retreated
            if grid.xM~=0
                [~,ixM] = gd_basin_indices(grid); %nearest grid point
                if ixM>1
                    subgrid = grid;
                    subgrid.x = grid.x(ixM:end);
                    subgrid.z = grid.z(ixM:end,:);
                    obj.RunParam.CF_ShoreData.ShoreWidth = grid.xM;
                    z0 = obj.RunParam.CF_HydroData.zmt(ixM); %msl at the mouth
                    shore = setShoreline(obj.RunParam.CF_ShoreData,grid,z0,false); %false returns a shore strip grid
                    grid.z(1:ixM-1,:) = shore.z(end-(ixM-2):end,:);                    
                elseif ixM<1
                    warndlg('Seaward migration not handled in CF_TransModel.updateShoreline')
                end
                % %defined slope shoreline  
                % hydobj = obj.RunParam.CF_HydroData;
                % trnobj = obj.RunParam.CF_TransData;
                % zmx = hydobj.zhw(1)+trnobj.FPoffset;
                % zmn = min(obj.Valley.Data.Grid.Z,[],'all');
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
                    [~,ixM] = gd_basin_indices(grid);
                    if ixM>0
                        %offset y-dimensions of meander by the amount of
                        %coast erosion
                        xM = sum(obj.Trans.cstdX)+obj.dTrans.cstdX;
                        if grid.xM~=xM
                            %grid includes a shore so find coastal erosion
                            delx = abs(grid.x(2)-grid.x(1));
                            ixM = round(xM/delx);
                        end
                        cline.y = [ones(ixM,1)*cline.y(1);cline.y(1:end-ixM)];
                        grid.cline = cline;
                    end
                else      
                    %meander is fixed as initially defined
                    cline = grid.cline;
                end
                grid = gd_xy2sn(grid,cline,true,false);
            end
        end
%%
        function obj = updateTidalLimit(obj)
            %find the new tidal limit and update integration distance indices
            %NB: water levels corrected for shoreline offset in cf_set_hydroprops
            %and Lt is distance from grid origin
            obj.dTrans.Lt = tidalLimit(obj);
            obj.RunParam.CF_HydroData.xTidalLimit = obj.dTrans.Lt; %used in non-cst form models
            
            %user specified limit for integration
            intdist = obj.RunParam.CF_TransData.IntDist;
            [ich,subxl] = gd_basin_indices(obj.Grid,[],obj.dTrans.Lt); %x-indices of channel
            if intdist>0  %user has defined integration distance from mouth
                subxu = find(obj.Grid.x>intdist,1,'first');
                %increment upper to maintain specified integration distance
                subxu = subxu+subxl-1; 
                ich = subxl:subxu;
            end
            obj.TranProp.intidx = ich;
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

           Pxz =  InterX([obj.Grid.x';thalweg'],[obj.Grid.x';hwl]);
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
        function fparea = floodPlainArea(obj)
            %get the flood plain area based on elevations at HW+offset           
            fp_offset = obj.RunParam.CF_TransData.FPoffset;
            %model water level surfaces
            hydobj = obj.RunParam.CF_HydroData;
            ich = gd_basin_indices(obj.Grid);  %account for offset to mouth
            zfp = hydobj.zhw(ich)+fp_offset;
            Zfp = repmat(zfp',1,length(obj.Grid.y));

            zf = obj.Grid.z(ich,:);
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
            hydobj = obj.RunParam.CF_HydroData;
            dt = obj.delta/obj.cns.y2s; %time step in years            
            slr = obj.dTrans.SLR;       %slr in time step (m)

            if trnobj.SedFlux>0
                %user defined import/export flux in m3/yr              
                sedvol = trnobj.SedFlux*dt;                %+ve volume for import
                watervol = obj.CSTparams.Shw*slr;
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
                    V0 = obj.Channel.Data.GrossProps.Vhw;
                    Pr0 = obj.Channel.Data.GrossProps.PrA;
                    inp.EqScaleCoeff = V0/Pr0; 
                    %used in v1 code - but updates each time step??
                    % inp.EqScaleCoeff = inp.Volume/inp.Prism; 
                    inp.EqShapeCoeff = 1;
                else
                    inp.EqScaleCoeff = sedobj.EqScaleCoeff;
                    inp.EqShapeCoeff = sedobj.EqShapeCoeff;
                end     
                
                %vertical and horizontal rates of sediment exchange
                inp = sedExchanges(obj,inp); 
                
                %sediment flux with external environment
                [sedvol,watervol,conc] = get_sed_flux(inp,slr/dt);
                %report sedvol as +ve volume for sediment import       
                sedvol = -sign(inp.TransportCoeff)*sedvol*dt; 
                watervol = watervol*dt;
            end
        end
%%
        function inp = sedExchanges(obj,inp)
            %set the vertical and horizontal rate sediment exchange
            hydobj = obj.RunParam.CF_HydroData;
            sedobj = obj.RunParam.CF_SediData;
            cn = obj.cns;                              %model constants 
            %vertical exchange   
            ws = settling_velocity(sedobj.SedimentSize,cn.g,cn.rhow,...
                                    cn.rhos,cn.visc,sedobj.EqDensity);                  
            inp.VerticalExchange = ws;                 %vertical exchange (m/s)

            %horizontal exchange
            if obj.Selection.wlflag==0 && ~isempty(hydobj.cstres)  
                %hydraulic results from CST model (cstres are results from
                %mouth interpolated onto the static form grid)
                u = hydobj.cstres.U(1);                %velocity at mouth
                L = min(3*obj.CSTparams.La,obj.dTrans.Lt);
                ich = gd_basin_indices(obj.Grid,[],L); %indices from mouth to tidal limit
                H = mean(hydobj.cstres.d(ich));        %average hydraulic depth
                %H = mean(hydobj.cstres.d);
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
        function obj = updateFormParams(obj)
            %update the form model input properties based on the new offset
            hydobj = obj.RunParam.CF_HydroData; 
            grid = obj.Grid;
            tp = obj.TranProp;
            %net change in position in a timestep of "initial" mouth 
            %due transgression of coast and estuary NB - relative to 
            %translating channel form and not x=0
            edX = obj.dTrans.estdX-obj.dTrans.cstdX; 

            %for the ckfa model the tidal limit and tidal amplitude 
            %are the input parameters that change and these are updated 
            %when CF_HydroData is updated in cf_set_hydroprops and
            %updateTidalLimit  
            
            %for the other form models, using CSTmodel, the mouth parameters
            %need updating. Note that the convergence lengths are set in 
            %InitiailiseModel>setCSTparams and are assumed to remain 
            %constant in these models (this is a key difference to the
            %ckfa model, which uses La as a dependent parameter.
            Wm0 = obj.CSTparams.Wm;
            Lw = obj.CSTparams.Lw; 
            obj.CSTparams.Wm = (Wm0-tp.Rv.Wr)*exp(edX/Lw)+tp.Rv.Wr;
            
            Am0 = obj.CSTparams.Am;
            La = obj.CSTparams.La;
            obj.CSTparams.Am = (Am0-tp.Rv.Ar)*exp(edX/La)+tp.Rv.Ar;

            option = obj.Selection.modeltype;
            if strcmp(option,'Exponential') || strcmp(option,'Power')
                %location of mouth in grid
                [~,ixM] = gd_basin_indices(grid);
                %depth at mouth to mtl
                subgrid.x = grid.x(ixM:end);    %remove shore
                subgrid.z = grid.z(ixM:end,:);  
                thalweg = subgrid.z(:,ceil(size(subgrid.z,2)/2)); %centreline
                %invert x and add offset of mouth to max x to interpolate
                %seaward of the exisitng form if edX is +ve
                zm = interp1(flipud(subgrid.x),thalweg,max(subgrid.x)+edX,'linear','extrap');
                hc = hydobj.zmt(1)-zm;  %zmt is the new boundary value set in InitTimeStep and is scalar
                           
                %widths before transgression
                Wu = obj.RunParam.CF_FormData.HWmouthWidth;
                Wl = obj.RunParam.CF_FormData.LWmouthWidth;                
                
                %option 1 - update the convergence lengths and mouth widths
                %convergence length from mouth accounting for any offset
                % [ich,~] = gd_basin_indices(grid);
                % ewidth = tp.Wz.Whw-tp.Rv.Wr; ewidth(ewidth<0) = 0;
                % Lwu = -getconvergencelength(grid.x(ich),ewidth(ich)); %width convergence length at mean tide level
                % ewidth = tp.Wz.Wlw-tp.Rv.Wr; ewidth(ewidth<0) = 0;
                % Lwl = -getconvergencelength(grid.x(ich),ewidth(ich)); %csa convergence length at mean tide level
                % %increase in mouth width due to landward transgression
                % %new widths based on increment, edX, from previous time step
                % Wu = (Wu-tp.Rv.Wr)*exp(edX/Lwu)+tp.Rv.Wr; 
                % Wl = (Wl-tp.Rv.Wr)*exp(edX/Lwl)+tp.Rv.Wr;
                % if strcmp(option,'Exponential')
                %     %update the convergence lengths
                %     obj.RunParam.(classname).HWwidthELength = Lwu;
                %     obj.RunParam.(classname).LWwidthELength = Lwl;                        
                % elseif strcmp(option,'Power')    
                %     warndlg('Not yet coded in CF_TransModel.updateFormParams')
                %     %means that Wu and Wl do not change
                % end 
                %NB - should really update Lw and La to be consitent it
                %this option is used.
                
                %option 2 - update mouth widths only
                if strcmp(option,'Exponential')
                    Lwu = obj.RunParam.CF_FormData.HWwidthELength;
                    Lwl = obj.RunParam.CF_FormData.LWwidthELength;
                    %increase in mouth width due to landward transgression
                    %new widths based on increment, edX, from previous time step
                    Wu = (Wu-tp.Rv.Wr)*exp(edX/Lwu)+tp.Rv.Wr; 
                    Wl = (Wl-tp.Rv.Wr)*exp(edX/Lwl)+tp.Rv.Wr;    
                elseif strcmp(option,'Power')    
                    warndlg('Not yet coded in CF_TransModel.updateFormParams')
                    %means that Wu and Wl do not change
                end                 

                %update model run time parameters   
                obj.RunParam.CF_FormData.HWmouthWidth = Wu;
                obj.RunParam.CF_FormData.LWmouthWidth = Wl;
                obj.RunParam.CF_FormData.MTmouthDepth = hc; 
            elseif strcmp(option,'CKFA')
                %No additional parameters needed for CKFA model
            else
                warndlg('Unkonwn model type in CF_TransModel.updateInputParameters')
            end    
        end
%%
        function updateFluxParams(obj)
            %update the sediment flux input parameters based on new form
            hydobj = obj.RunParam.CF_HydroData; 
            grdobj = obj.RunParam.GD_GridProps;
            grid = obj.Grid;
            [~,hypdst] = gd_basin_hypsometry(grid,hydobj,grdobj.histint,0,true);
            spdst = gd_section_properties(grid,hydobj,hypdst);
            gpdst = gd_gross_properties(grid,hydobj,spdst);
            %update parameters used in call to get_sed_flux
            obj.CSTparams.Vhw = gpdst.Vhw; %element volume (m^3)
            obj.CSTparams.Shw = gpdst.Shw; %element surface area (m^2)
            obj.CSTparams.Pr = gpdst.PrA;  %tidal prism of channel (m^3) 
        end
%%
        function setChannelWaterLevels(obj)
            %use the Channel model water levels to initialise CF_HydroData
            hydobj = obj.RunParam.CF_HydroData;        %input hydro data
            wls = obj.Channel.Data.WaterLevels(1,:);   %Channel water levels
            
            hydobj.zhw = wls.zhw;          %high water
            hydobj.zmt = wls.zmt;          %mean tide level
            hydobj.zlw = wls.zlw;          %low water             
            obj.RunParam.CF_HydroData = hydobj;
        end
%%
        function setCSTparams(obj)
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
            obj.CSTparams.Pr = gp.PrA;    %tidal prism of channel (m^3) 
        end        
%%
        function setTransModel(obj,meta)
            %add the transgression results table to the model output
            tstep = obj.Data.Grid.RowNames;
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
            muicat = mobj.Cases;
            ftxt = 'Select Form Model to use:';
            obj.Channel = selectCaseObj(muicat,[],{'CF_FormModel'},ftxt);
            if isempty(obj.Channel), isok = false; return; end
            vtxt = 'Select Valley Model to use:';
            obj.Valley = selectCaseObj(muicat,[],{'CF_ValleyModel','GD_ImportData'},vtxt);
            if isempty(obj.Valley), isok = false; return; end
            
            %Metadata of selection
            obj.MetaData = sprintf('Transgression model using "%s" for channel and "%s" for valley',...
                           obj.Channel.Data.Grid.Description,obj.Valley.Data.Grid.Description);
            
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
            [obj,isok] = updateRunParameters(obj);
        end
%%
        function [obj,isok] = updateRunParameters(obj)
            %allow user to modify parameters that do not influence initial 
            %form including: sea level rise, tidal range cycles, and the 
            %coefficients used in the sediment flux calculations
            titletxt = 'Update Parameters';
            promptxt = {'Rate of sea level rise (m/year)', ...
                        'Tidal ampltude (m)',...
                        'Amplitude of Cycles (m)',...
                        'Equilibrium sediment density (kg/m^3)', ...
                        'Sediment load in river (kg/m^3)',...
                        'Transport coefficient (+/-n)',...                         
                        'Equilibrium scale coefficient (0=scale to initial)',...                          
                        'Equilibrium shape coefficient (-)',...
                        'Constant sediment flux (m^3/yr): 0,NaN,or value'};
            defaults = {num2str(obj.RunParam.WaterLevels.SLRrate),...
                        num2str(obj.RunParam.WaterLevels.TidalAmp),...
                        num2str(obj.RunParam.WaterLevels.CycleAmp),...
                        num2str(obj.RunParam.CF_SediData.EqDensity),...
                        num2str(obj.RunParam.CF_SediData.RiverDensity),...
                        num2str(obj.RunParam.CF_SediData.TransportCoeff),...
                        num2str(obj.RunParam.CF_SediData.EqScaleCoeff),...
                        num2str(obj.RunParam.CF_SediData.EqShapeCoeff),...
                        num2str(obj.RunParam.CF_TransData.SedFlux)};
            answer = inputdlg(promptxt,titletxt,1,defaults);
            if ~isempty(answer)
                isok = true;
                obj.RunParam.WaterLevels.SLRrate = str2num(answer{1}); %#ok<ST2NM>
                obj.RunParam.WaterLevels.TidalAmp = str2double(answer{2}); 
                obj.RunParam.WaterLevels.CycleAmp = str2num(answer{3}); %#ok<ST2NM>
                obj.RunParam.CF_SediData.EqDensity = str2double(answer{4});
                obj.RunParam.CF_SediData.RiverDensity = str2double(answer{5});
                obj.RunParam.CF_SediData.TransportCoeff = str2double(answer{6});
                obj.RunParam.CF_SediData.EqScaleCoeff = str2double(answer{7});
                obj.RunParam.CF_SediData.EqShapeCoeff = str2double(answer{8});
                obj.RunParam.CF_TransData.SedFlux = str2double(answer{9});
            else
                isok = false;
            end                       
        end
%%
        function timestepInput(obj)
            %write the inputs used for each time step to the command window
            time = obj.StepTime(end)/obj.cns.y2s;
            ixM = obj.TranProp.intidx(1);
            edX = obj.dTrans.estdX-obj.dTrans.cstdX; 
            gp = obj.CSTparams;
            hc = obj.RunParam.CF_FormData.MTmouthDepth; %*** only for exp/power models
            fprintf('T = %g yrs: ixM=%g; xM=%g m; edX=%g; Width Wm=%g m; CSA Am=%g m^2; Depth hc=%g m\n',...
                                   time,ixM,obj.Grid.xM,edX,gp.Wm,gp.Am,hc)
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
            p.Title = sprintf('Change plot, slr=%0.2g m',slr);
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
            
            % subx = obj.TranProp.intidx; %index to subsample domain for integration
            xs = obj.Grid.x;                
            zs = obj.Grid.z;%current grid
            ys = obj.Grid.y;  
            z0 = squeeze(obj.zGrid(1,:,:));
            zdiff = zs-z0;

            yzhw = obj.TranProp.Wz.Whw/2;
            yzlw = obj.TranProp.Wz.Wlw/2;
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
            tltxt = sprintf('dV=0 for: SLR = %0.2g m, Transgression = %d m, Coast erosion = %d m\n%s',...
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
            edX = obj.Trans.estdX(end);
            
            plot(ax,dx,dV,'-k','LineWidth',1)
            ax.YGrid = 'on';
            xlabel('Landward transgression (m)');
            ylabel('Volume change (m^3)');
            %add explanatory text
            mxX = ax.XLim(2); mxY = ax.YLim(2);
            hold on
            plot([0,mxX],[0,0],'--r');  %'Color',[0.96,0.96,0.96]
            plot([0,edX],[sedvol,sedvol],'--b'); %transgression point
            plot([edX,edX],[0,sedvol],'--b'); %transgression point
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
            z0 = squeeze(obj.zGrid(1,:,:));
            gd_dir = gd_ax_dir(obj.Grid);
            if gd_dir.x==1 || gd_dir.x==4                 
                %orientation of x-axis, x=0 is nearest the mouth if ishead=false
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
            %reverse x-axis - sections are oriented looking from mouth up-estuary
            ax.XDir = 'reverse';
            %add meta-data
            xlabel('Width (m)'); 
            ylabel('Change in level (m)');
            hL=legend('0pre','0post','0.1pre','0.1post','0.2pre','0.2post','Location','SouthEast');
            set(hL, 'Color', 'none');
            casedesc = sprintf('Cross-sections difference plot, slr=%0.2g m',slr);
            title(casedesc,'FontWeight','normal','FontSize',10);            
        end
%%
        function thalwegPlot(obj,slr)
            %plot the thalweg at the start and end of run
            % NB uses transient properties so only works in current session
            xi = obj.Grid.x;  
            zi = obj.Grid.z(:,ceil(size(obj.Grid.z,2)/2));

            z0 = squeeze(obj.zGrid(1,:,:));
            z0 = z0(:,ceil(size(z0,2)/2));
            gd_dir = gd_ax_dir(obj.Grid);
            if gd_dir.x==1 || gd_dir.x==4                 
                %orientation of x-axis, x=0 is nearest the mouth if ishead=false
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
            casedesc = sprintf('Centre-line difference plot, slr=%0.2g m',slr);
            title(casedesc,'FontWeight','normal','FontSize',10); 
        end
%%
        function summaryGridPlots(obj,slr)             
            %plot initial grid, new grid and difference
            % NB uses transient properties so only works in current session
            z0 = squeeze(obj.zGrid(1,:,:));
            F = obj.Grid;
            zdiff = F.z-z0;
            start = obj.RunParam.RunProperties.StartYear;
            
            figure('Name','Grid Plot','Units','normalized','Tag','PlotFig'); 
            subplot(3,1,1)
            surf(F.x,F.y,z0','FaceColor','interp','EdgeColor', 'none');    
            
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
            casedesc = sprintf('Summary of change, slr=%0.2g',slr);
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
            % SLR - cumulative sea level rise
            % Lt - distance to tidal limit (to record change - model uses CF_HydroData.xTidalLimit)
            % FPA - flood plain area
            % waterVol - water volume change due to changes in slr and tidal range
            % sedVol - sediment flux (+ve=sediment import)            
            % vdiffx - volume difference for [0,delX/2,delX,3delX/2]
            
            dsp.Variables = struct(...
                        'Name',{'delX','estdX','cstdX','SLR','Lt','FPA',...
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

