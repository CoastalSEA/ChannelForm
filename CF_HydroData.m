classdef CF_HydroData < muiPropertyUI                
%
%-------class help---------------------------------------------------------
% NAME
%   CF_HydroData.m
% PURPOSE
%  Class for hydraulic parameters used by CSTmodel in ChannelForm model
% USAGE
%   obj = CF_HydroData.setInput(mobj); %mobj is a handle to Main UI
% NOTES
%   Co-ordinate convention is that the x-axis is landward from the mouth to
%   the tidal limit. The y-axis is the cross-channel axis, mirrored about the
%   centre-line.
% SEE ALSO
%   inherits muiPropertyUI. Very similar to CSThydraulics in Asmita
%   used in CF_FormModel and CT_Transgression in the ChannelForm model
%
% Author: Ian Townend
% CoastalSEA (c) Jan 2022
%--------------------------------------------------------------------------
%      
    properties (Hidden)
        %abstract properties in muiPropertyUI to define input parameters
        PropertyLabels = {'River discharge (m3/s)',...
                          'Wind speed at 10m (m/s)',...
                          'Distance to estuary/river switch (m)',...
                          'Distance to tidal limit (m)',...
                          'Manning friction coefficient [mouth switch head]',...
                          'Intertidal storage ratio [mouth switch head]',...
                          'CSTmodel distance increment (m)',...
                          'River discharge range [Q1 Q2 ...Qn]'};
        %abstract properties in muiPropertyUI for tab display
        TabDisplay   %structure defines how the property table is displayed 
    end
    
    properties
        RiverDischarge = 0    %river discharge (m^3/s) +ve downstream
        WindSpeed = 0         %wind speed at 10m (m/s)
        xTideRiver            %distance from mouth to estuary/river switch (m)
        xTidalLimit           %distance from mouth to tidal limit (m)
        Manning               %Manning friction coefficient [estuary river]
        StorageRatio          %intertidal storage ratio [estuary river]
        DistInt = 5000        %distance increment along estuary (m) 
        Qrange                %vector of input river discharges (m^3/s)
    end     
    
    properties (Transient)
        tidalperiod     %tidal period in seconds (s)
        zhw = 0         %high and low water level from amplitude+mtl 
        zmt = 0         %or CST model if varying along channel and/or
        zlw = 0         %time dependent (mOD)
        dhw = 0         %change in high water during time step
        dslr = 0        %rate of sea level rise at time t (only changes if SLRrate is not constant)
        Qr = 0          %time dependent river discharge
        FormModel       %definition of morphological form
        cstres          %struct with velocity and hyd.depth from cst model
        cstmsg = 'CSTmodel not found. Check that App is installed';
    end

    properties (Hidden)
        CSTmodel         %dstable from CSTrunmodel with following fields
        % MeanTideLevel  %z - mean water suface elevation along estuary
        % TidalElevAmp   %a - elevation amplitude along estuary
        % TidalVelAmp    %U - tidal velocity amplitude along estuary
        % RiverVel       %v - river flow velocity along estuary
        % HydDepth       %d - along estuary hydraulic depth 
        %used for Hydraulics display utility (NOT used in models)
    end  
%%   
    methods (Access=protected)
        function obj = CF_HydroData(mobj)       
            %constructor code:            
            %TabDisplay values defined in UI function setTabProperties used to assign
            %the tabname and position on tab for the data to be displayed
            obj = setTabProps(obj,mobj);  %muiPropertyUI function            
            ok = initialise_mui_app('CSTmodel',obj.cstmsg,'CSTfunctions');
            if ok<1, return; end
        end 
    end
%%  
    methods (Static)  
        function obj = setInput(mobj,editflag)
            %gui for user to set Parameter Input values
            classname = 'CF_HydroData';           
            obj = getClassObj(mobj,'Inputs',classname);
            if isempty(obj)
                obj = CF_HydroData(mobj);      
            end
            
            %use muiPropertyUI function to generate UI
            if nargin<2 || editflag
                %add nrec to limit length of props UI (default=12)
                obj = editProperties(obj);  
                %add any additional manipulation of the input here
            end
            setClassObj(mobj,'Inputs',classname,obj);
        end  
%%
        function displayRiverDims(mobj)
            %display river dimensions for current property settings
            msgtxt = 'Hydraulic Parameters have not been defined';
            obj = getClassObj(mobj,'Inputs','CF_HydroData',msgtxt);
            if isempty(obj), return; end
            msgtxt = 'Sediment Parameters have not been defined';
            sedobj = getClassObj(mobj,'Inputs','CF_SediData',msgtxt);
            if isempty(sedobj), return; end
            msgtxt = 'Water level Parameters have not been defined';
            wlvobj = getClassObj(mobj,'Inputs','WaterLevels',msgtxt);
            if isempty(wlvobj), return; end

            am = wlvobj.TidalAmp;
            d50riv = sedobj.d50river;
            tauriv = sedobj.tauriver;
            rhos = mobj.Constants.SedimentDensity;
            rhow = mobj.Constants.WaterDensity;
            Qr = obj.RiverDischarge;
            Le = obj.xTidalLimit;   %distance to tidal limit 
            Sr  = 2*am/Le;

            [hrv,Wrv,~] = river_regime(Qr,Sr,d50riv,tauriv,rhos,rhow);
            msgbox(sprintf('River input is %gm^3/s, slope is %g\nWidth is %0.1fm and hydraulic depth is %0.2fm',...
                      Qr,Sr,Wrv,hrv));
        end       
    end
%%
    methods
        function runModel(obj,mobj)
            %compile input data and run model. Uses definitions based on
            %CF_CKFAdata properties when run from Utilities>Hydraulic Model.
            %Whereas runModelatT (see below) which uses the form model currently 
            %being manipulated (water surface added to morphological form)            
            ok = initialise_mui_app('CSTmodel',obj.cstmsg,'CSTfunctions');
            if ok<1, return; end
            
            if isempty(obj.xTideRiver)  %check first property has been set
                warndlg('Hydraulic properties have not been defined')
                return;
            end
            dsp = modelDSproperties(obj);

            %now check that the input data has been entered
            %isValidModel checks the InputHandles defined in ModelUI
            if ~isValidModel(mobj, metaclass(obj).Name)
                warndlg('Use Setup to define Form and Forcing parameters');
                return;
            end
            %--------------------------------------------------------------------------
            % Model code
            %--------------------------------------------------------------------------
            %input parameters for model 
            setTransHydroProps(obj,mobj);
            obj.FormModel = selectFormModel(obj,mobj);
            if isempty(obj.FormModel), return; end  %no form model retrieved
            
            %assign the run parameters to the model instance          
            [inp,rnp] = getHydroModelParams(obj,false);            
            est = [];  %observed values of estuary form so can be empty

            %run model iteratively for a range of river discharges
            if isempty(obj.Qrange)
                obj.Qrange = obj.Qr;
            end
            nrow = length(obj.Qrange);
            resX{nrow,5} = [];
            for i=1:nrow
                inp = updateModelParameters(obj,inp,i);
                try
                    [res,xy,~,~] = cst_model(inp,rnp,est);
                    resX(i,:) = res;
                catch
                    %remove the waitbar if program did not complete
                    hw = findall(0,'type','figure','tag','TMWWaitbar');
                    delete(hw);
                    inpQr = inp.RiverDischarge;
                    msg = sprintf('Failed to find solution in cst_model for Qr=%d',inpQr);
                    warndlg(msg);
                    return;
                end
            end
            resXQ = cell(1,5);
            for col = 1:5
               resXQ{col} = vertcat(resX{:,col});
            end
            %now assign results to object properties
            %--------------------------------------------------------------------------
            % Assign model output to a dstable using the defined dsproperties meta-data
            %--------------------------------------------------------------------------
            dst = dstable(resXQ{:},'RowNames',inp.Qrange','DSproperties',dsp);
            dst.Dimensions.X = xy{:,1};     %grid x-coordinate
            %--------------------------------------------------------------------------
            % Save results
            %--------------------------------------------------------------------------
            %assign metadata about model
            dst.Source = metaclass(obj).Name;
            %save results
            obj.CSTmodel = dst;
            setClassObj(mobj,'Inputs','CSThydraulics',obj);
            getdialog('Run complete');            
        end
%%
        function [resX,xyz,resXT,time] = runHydroModel(obj,estobj)
            %run model when updating models eg adding surface to form model or in
            %transgression model (i) no checks made; (ii) uses current transient
            %properties for water levels; and (iii) specified form model. 
            %See cst_model help for definitions of output.
            resX = []; resXT = []; time = []; xyz = [];
            ok = initialise_mui_app('CSTmodel',obj.cstmsg,'CSTfunctions');
            if ok<1, return; end
            
            obj.FormModel = estobj;
            [inp,rnp] = getHydroModelParams(obj,true);
            
            try
                [resX,xyz,resXT,time] = cst_model(inp,rnp,[]);
            catch ME
                %remove the waitbar if program did not complete
                model_catch(ME,'cst_model','Qr',inp.RiverDischarge);                
            end
        end
%%
        function tabPlot(obj,src,~)
            %create a summary of the hydraulic model results
            if isempty(obj.CSTmodel)
                warndlg('No model results to display')
                return;
            end
            
            if nargin<2
                src = figure('Name','Hydraulic plot','Tag','PlotFig');
            else
                ax = findobj(src,'Tag','PlotFigAxes');
                delete(ax)
            end
            ax = axes('Parent',src,'Tag','PlotFigAxes');
            ax.Position = [0.16,0.18,0.65,0.75]; %make space for slider bar
            setYaxisLimits(obj,ax);  %wakes dynamic properties if necessary
            
            Q = obj.Qrange;
            setSlideControl(obj,src,Q(1),Q(end));
            
            cstPlot(obj,ax,Q(1))
        end
%%
        function obj = setTransHydroProps(obj,mobj)
            %initialise the transient properties used in the models
            wlvobj = getClassObj(mobj,'Inputs','WaterLevels');
            [obj.zhw,obj.zmt,obj.zlw] = newWaterLevels(wlvobj,0,0);
            obj.Qr = obj.RiverDischarge;  %initialise transient river discharge
            obj.tidalperiod = wlvobj.TidalPeriod*3600; %tidal period in seconds
        end 
%%
        function newWaterLevels(obj,mobj,robj)
            %update water levels when running transgression model
            % robj is the run time object that defines the time step
            % ie CF_Transgression in the ChannelForm model
            % WaterLevels provides water levels at the mouth. To get 
            % alongchannel values call cf_set_hydroprops after calling
            % newWaterLevels.
            WaterLevels.setWaterLevels(mobj,robj);
            wlvobj = getClassObj(mobj,'Inputs','WaterLevels');
            obj.zhw = wlvobj.HWaterLevel;
            obj.zmt = wlvobj.MeanSeaLevel;
            obj.zlw = wlvobj.LWaterLevel;
            obj.dhw = wlvobj.dHWchange;  %change in high water over a time step
            obj.dslr = wlvobj.dslr;
        end        
    end
%%    
    methods (Access = private)
        function cst = interpCSTproperties(obj,Q)
            %interpolate the CSTmodel output for given value of river discharge
            dst = obj.CSTmodel;
            Qrows = dst.RowNames;
            cst.x = dst.Dimensions.X; 
            if ~isprop(dst,'MeanTideLevel') 
                %check if dynamic properties are active
                dst = activatedynamicprops(dst);
            end
            cst.z = interp1(Qrows,dst.MeanTideLevel,Q,'linear'); %mean tide level
            cst.a = interp1(Qrows,dst.TidalElevAmp,Q,'linear');  %tidal amplitude
			cst.U = interp1(Qrows,dst.TidalVelAmp,Q,'linear');   %tidal velocity amplitude
			cst.v = interp1(Qrows,dst.RiverVel,Q,'linear');      %river velocity 
			cst.d = interp1(Qrows,dst.HydDepth,Q,'linear');      %hydraulic depth
        end
%%
        function cstPlot(obj,ax,Q)
            %plot the along channel variables from the CSTmodel (as per tab)
            cst = interpCSTproperties(obj,Q);
            green = mcolor('green');
            orange = mcolor('orange');
			yyaxis(ax,'left')
            cla                                  %clear any existing plot lines
            plot(ax,cst.x,cst.z,'-r','LineWidth',1.0);           %plot time v elevation
            hold on
            plot(ax,cst.x,(cst.z+cst.a),'-.b','LineWidth',0.8)   %plot high water level
            plot(ax,cst.x,(cst.z-cst.a),'-.b','LineWidth',0.8)   %plot low water level
			plot(ax,cst.x,(cst.z-cst.d),'-k','LineWidth',0.6);   %hydraulic depth below mean tide level
            ylabel('Elevation (mOD)'); 
			yyaxis(ax,'right')
            cla                                  %clear any existing plot lines
			plot(ax,cst.x,cst.U,'--','Color',orange,'LineWidth',0.6)%plot tidal velocity
			plot(ax,cst.x,cst.v,'--','Color',green,'LineWidth',0.6) %plot river velocity
            hold off
            xlabel('Distance from mouth (m)'); 
            ylabel('Velocity (m/s)'); 
			legend('MTL','HWL','LWL','Hydraulic depth',...
                'Tidal velocity','River velocity','Location','east');			
            title('Along channel variation');
        end

%%
        function hm = setSlideControl(obj,hfig,qmin,qmax)
            %intialise slider to set different Q values   
            invar = struct('sval',[],'smin',[],'smax',[],'size', [],...
                           'callback','','userdata',[],'position',[],...
                           'stxext','','butxt','','butcback','');            
            invar.sval = qmin;     %initial value for slider 
            invar.smin = qmin;     %minimum slider value
            invar.smax = qmax;     %maximum slider value
            invar.callback = @(src,evt)updateCSTplot(obj,src,evt); %callback function for slider to use
            invar.position = [0.15,0.005,0.45,0.04]; %position of slider
            invar.stext = 'River discharge = ';   %text to display with slider value, if included          
            hm = setfigslider(hfig,invar);   
        end   
%%
        function updateCSTplot(obj,src,~)
            %use the updated slider value to adjust the CST plot
            stxt = findobj(src.Parent,'Tag','figsliderval');
            Q = round(src.Value);
            stxt.String = num2str(Q);     %update slider text
            %figure axes and update plot
            figax = findobj(src.Parent,'Tag','PlotFigAxes'); 
            cstPlot(obj,figax,Q)
        end
%%
        function setYaxisLimits(obj,ax)
            %set the Y axis limits so they do not change when plot updated
            dst = obj.CSTmodel;
            if ~isprop(dst,'MeanTideLevel')
                dst = activatedynamicprops(dst);
            end
            z = dst.MeanTideLevel;
            a = dst.TidalVelAmp;
            d = dst.HydDepth;
            U = dst.TidalVelAmp;
            v = dst.RiverVel;
            
            lim1 = floor(min(z-d,[],'All'));
            lim2 = ceil(max(z+a,[],'All'));
            lim3 = floor(min(v,[],'All'));
            lim4 = ceil(max(U,[],'All'));
            yyaxis left                          %fix y-axis limits
            ax.YLim = [lim1,lim2];
            yyaxis right
            ax.YLim = [lim3,lim4];
        end        
%%
        function [inp,rnp] = getHydroModelParams(obj,incriver)
            %extract the additional parameters needed to run the CSTmodel
            %from the CF_HydroData and CF_FormModel classes
            inp = getPropertiesStruct(obj);         
            
            %inp parameters requires the following form parameters
            %in the CSTmodel estuary length is the length of model domain
            %and inp.xTidalLimit is not used
            xlength = obj.FormModel.RunParam.GD_GridProps.XaxisLimits(2);
            inp.EstuaryLength = xlength;%total length of channel (m)
            
            %width, CSA and convergence length are form model dependent
            inp = getHydroFormProps(obj,inp);

            %get time dependent water level properties (due to slr and ntc)          
            amp = (obj.zhw(1)-obj.zlw(1))/2;     %tidal amplitude at mouth
            %inp parameters requires the following water level parameters
            inp.MTLatMouth = obj.zmt(1);         %mean tide level at mouth (mOD)
            inp.TidalAmplitude = amp;            %tidal amplitude (m)
            inp.TidalPeriod = obj.tidalperiod/3600 ;%tidal period (hr)
            
            if incriver
                inp.RiverDischarge = obj.Qr;     %transient river discharge
                [~,Wrv,Arv] = get_river_regime(obj.FormModel,2*amp,obj.Qr);                
                inp.RiverWidth = Wrv;            %upstream river width (m) 
                inp.RiverCSA = Arv;              %upstream river cross-sectional area (m^2)
            else
                inp.RiverDischarge = 0;
                inp.RiverWidth = 0;
                inp.RiverCSA = 0;
            end

            %rnp parameters requires the following
            rnp.TimeInt = 0;           %time increment in analytical model (hrs) - only needed if tidal cycle output required
            rnp.DistInt = obj.DistInt; %distance increment along estuary (m)
            rnp.useObs = false;        %flag to indicate whether to use observations
        end
%%
        function inp = getHydroFormProps(obj,inp)
            %initialise the properties required for the CST model using the
            %Form model which was used to call CF_HydroData and assigned to
            %obj.FormModel
            if isempty(obj.FormModel.Data)
                frm = obj.FormModel.CSTparams;                   
            else 
                frm = obj.FormModel.Data.GrossProps; 
            end
            inp.MouthWidth = frm.Wm;    %width of mouth at high water(m)
            inp.WidthELength = frm.Lw;  %width convergence length at high water (m)
            inp.MouthCSA = frm.Am;      %CSA of mouth at high water(m2)
            inp.AreaELength = frm.La;   %area
        end
%%
        function [formobj,caserec] = selectFormModel(~,mobj)    
            %prompt user to select a form model from existing cases
            muicat = mobj.Cases;
            promptxt = 'Select Form Model to use:';            
            [caserec,ok] = selectRecord(muicat,'PromptText',promptxt,...
                              'CaseType',{'form_model'},'CaseClass',[],...
                              'SelectionMode','single','ListSize',[250,200]);
            if ok<1, formobj = []; return; end
            formobj = getCase(muicat,caserec);             
        end
%%
        function inp = updateModelParameters(obj,inp,idx)
            %update additional properties that change with river discharge
            inp.RiverDischarge = obj.Qrange(idx);%value to use in cst_model
            [~,Wrv,Arv] = get_river_regime(obj.FormModel,obj.Qrange(idx));
            inp.RiverWidth = Wrv;                %upstream river width (m) 
            inp.RiverCSA = Arv;                  %upstream river cross-sectional area (m^2)
        end
%%
        function dsp = modelDSproperties(~)
            %define a dsproperties struct and add the model metadata
            dsp = struct('Variables',[],'Row',[],'Dimensions',[]);
            %define each variable to be included in the data table and any
            %information about the dimensions. dstable Row and Dimensions can
            %accept most data types but the values in each vector must be unique

            %struct entries are cell arrays and can be column or row vectors
            %static ouput (mean tide values)
            dsp.Variables = struct(...
                'Name',{'MeanTideLevel','TidalElevAmp','TidalVelAmp',...
                        'RiverVel','HydDepth'},...
                'Description',{'Mean water level',...
                               'Tidal elevation amplitude',...
                               'Tidal velocity amplitude',...
                               'River flow velocity',...
                               'Hydraulic depth'},...
                'Unit',{'m','m','m/s','m/s','m'},...
                'Label',{'Mean water level (m)',...
                         'Elevation amplitude (m)',...
                         'Velocity amplitude (m/s)',...
                         'Velocity (m/s)','Depth (m)'},...
                'QCflag',repmat({'model'},1,5));
            dsp.Row = struct(...
                'Name',{''},...
                'Description',{''},...
                'Unit',{''},...
                'Label',{''},...
                'Format',{''});
            dsp.Dimensions = struct(...
                'Name',{'X'},...
                'Description',{'Chainage'},...
                'Unit',{'m'},...
                'Label',{'Distance from mouth (m)'},...
                'Format',{'-'});
        end
    end   
end