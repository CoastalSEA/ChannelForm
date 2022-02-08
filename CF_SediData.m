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
                          'Erosion rate constant (kg/N/s)',...
                          'Equilibrium sediment density (kg/m^3)', ...                          
                          'Bulk density of bed (kg/m^3)',...
                          'Transport coefficient (+/-n)',...
                          'Equilibrium scale coefficient (0=scale to initial)',...
                          'Equilibrium shape coefficient (-)',...
                          'Average depth over marsh (m)',...
                          'Maximum marsh depth (m), 0=no marsh',...
                          'Sediment load in river (kg/m^3)',...
                          'Sediment grain size in river, d50 (m)',...
                          'Critical bed shear stress in river (Pa)'};
        %abstract properties in muiPropertyUI for tab display
        TabDisplay   %structure defines how the property table is displayed 
    end
    
    properties
        SedimentSize         %sediment grain size, D50 (m)
        CritBedShear         %critical bed shear stress (Pa)
        ErosionRate = 0.002  %erosion rate constant (kg/N/s)
        EqDensity = 0        %equilibrium concentration density (kg/m^3)        
        BedDensity = 0       %bed density (kg/m^3)
        TransportCoeff = 3   %transport coefficient n (3-5) 
        EqScaleCoeff = 0.84  %equilibrium scale coefficient, alpha, UK default - set to zero to scale to initial volume
        EqShapeCoeff = 1     %equilibrium shape coefficient, beta, UK default      
        AvMarshDepth = 0     %average depth of marsh surface (m)
        MaxMarshDepth = 0    %maximum depth of saltmarsh (m)
        RiverDensity = 0     %river load imported by advection (kg/m^3)
        d50river             %sediment grain size in river, D50 (m)
        tauriver             %critical bed shear stress in river (Pa)
    end    
    
    properties (Dependent)
        EqConcentration      %equilibrium concentration (-)
        BedConcentration     %concentration of bed (-)
        RiverConcentration   %concentration of river load (-)
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
                obj = editProperties(obj,13);  
                %add any additional manipulation of the input here
            end
            setClassObj(mobj,'Inputs',classname,obj);
        end   
%%
        function displayMorphTime(mobj)
            %display the morphological time scale for a single element using
            %the model parameters currently defined
            obj = getClassObj(mobj,'Inputs','CF_SediData');
            if isempty(obj)
                warndlg('Model properties, including sediment properties, need to be defined')
                return;
            end
            
            %prompt user to select an existing form model to modify
            promptxt = 'Select an existing Form Model?';
            formobj = selectCaseObj(mobj.Cases,{'form_model'},[],promptxt);
            if isempty(formobj), return; end
           
            sedinp = getSedFluxInputs(obj,mobj,formobj);
            V = sedinp.Volume;            %element volume at start of run (m^3)
            S = sedinp.SurfaceArea;       %element surface area (m^2)
            n = sedinp.TransportCoeff;    %transport coefficient n (3-5)
            cE = sedinp.EqConc;           %equilibrium concentration (-)
            w = sedinp.VerticalExchange;  %vertical exchange (m/s)
            d = sedinp.HorizontalExchange;%horizontal exchange (m/s)
            
            tau = 1/(n*cE)*(V/(w*S)+V/d)/sedinp.y2s;
  
            answer = inputdlg('Rate of slr (m/yr)','SLR rate',1,{'0.002'});
            slr = str2double(answer);

            %get_sed_flux_returns change in morphological volume/yr, dvol 
            %hence negative is infilling and import of sediment
            %delV is the water volume change (S x slr)
            [dvol,delV] = get_sed_flux(sedinp,slr);
            msg1 = sprintf('SEM morphological timescale = %0.3f years',tau);
            msg2 = sprintf('Volume change due to SLR (m^3/yr) = %0.3e',delV);
            msg3 = sprintf('Volume imported/exported (m^3/yr) = %0.3e',-dvol);
            msgtxt = sprintf('%s\n%s\n%s',msg1,msg2,msg3);
            msgbox(msgtxt,'Morphological Time');
        end
    end
%%        
    methods
        function eqConc = get.EqConcentration(obj)
            %dependent property derived from EqDensity
            cn = muiConstants.Evoke;
            eqConc = obj.EqDensity/cn.SedimentDensity;
        end
%%
        function eqConc = get.BedConcentration(obj)
            %dependent property derived from BedDensity
            cn = muiConstants.Evoke;
            nomi = obj.BedDensity-cn.WaterDensity;
            dnom = cn.SedimentDensity-cn.WaterDensity;
            eqConc = nomi/dnom;
        end
%%
        function eqConc = get.RiverConcentration(obj)
            %dependent property derived from RiverDensity
            cn = muiConstants.Evoke;
            eqConc = obj.RiverDensity/cn.SedimentDensity;
        end
    end
%%    
    methods (Access=private)
        function sedinp = getSedFluxInputs(obj,mobj,formobj)
            %use model definition to construct sedflux input parameters
            answer = questdlg('Use saved or current sediment properties?',...
                                     'Sed Props','Saved','Current','Saved');
            if strcmp(answer,'Current')     %sediment properties
                sed = obj;
            else
                sed = formobj.RunParam.CF_SediData;                
            end

            cn = getConstantStruct(mobj.Constants); %model constants 
            sedinp.y2s = cn.y2s;
            ws = settling_velocity(sed.SedimentSize,cn.g,cn.rhow,...
                                          cn.rhos,cn.visc,sed.EqDensity);                                             
            sedinp.VerticalExchange = ws;             %vertical exchange (m/s)
            sedinp.BedConc = sed.BedConcentration;    %concentration of bed (-)
            sedinp.EqConc = sed.EqConcentration;      %equilibrium concentration (-)
            sedinp.RiverConc = sed.RiverConcentration;%river load imported by advection (-)

            gprop = formobj.Data.GrossProps(1,:);
            sedinp.Volume = gprop.Vhw;              %element volume at start of run (m^3)
            sedinp.SurfaceArea = gprop.Shw;         %element surface area (m^2)
            sedinp.Prism = gprop.Pr;                %tidal prism of channel (m^3)
            
            hydobj = formobj.RunParam.CF_HydroData;
            sedinp.RiverDischarge = hydobj.RiverDischarge;  %river discharge (m^3/s)

            if ~isempty(hydobj.cstres)  %hydraulic results from CST model
                u = hydobj.cstres.U(1);             %velocity at mouth
                H = mean(hydobj.cstres.d);               %average hydraulic depth
            else                        %no hydraulic data
                u = 1;                              %assumed velocity at mouth
                H = gprop.Vhw/gprop.Shw;            %estimated hydraulic depth
            end
            A = gprop.Am;                    %CSA at mouth
            D = u^2*H/ws;
            wlvobj = mobj.Inputs.WaterLevels;
            %delx = cfm.Le/2;                       %half estuary length
            delx = u*wlvobj.TidalPeriod/4*3600;     %tidal excursion length
            sedinp.HorizontalExchange = D*A/delx;   %horizontal exchange (m/s)
            sedinp.TransportCoeff = sed.TransportCoeff;%transport coefficient n (3-5)
            %equilibrium volume definition
            if sed.EqScaleCoeff==0
                sedinp.EqScaleCoeff = gprop.Vhw/gprop.Pr; 
                sedinp.EqShapeCoeff = 1;
            else
                sedinp.EqScaleCoeff = sed.EqScaleCoeff;
                sedinp.EqShapeCoeff = sed.EqShapeCoeff;
            end
        end        
    end
end