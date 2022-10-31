function ckfa_dimensions(mobj)
%
%-------function help------------------------------------------------------
% NAME
%   ckfa_dimensions.m
% PURPOSE
%   utility to get the CKFA channel form properties: hm,Lw,Ucr
% USAGE
%   ckfa_dimensions(mobj)
% INPUTS
%   mobj - instance of ChannelForm model UI
% OUTPUTS
%   table figure of output
% NOTES
%   CKFA cross-section comprises a channel using Cao&Knight section and an
%   intertidal using the profile proposed by Friedrichs & Aubrey (tide only)
% SEE ALSO
%   used in CF_FormModel as part of ChannelForm model
%
% Author: Ian Townend
% CoastalSEA (c) Jan 2022
%--------------------------------------------------------------------------
%        
    hydobj = mobj.Inputs.CF_HydroData;
    sedobj = mobj.Inputs.CF_SediData;
    wlvobj = mobj.Inputs.WaterLevels;
    if ~isValidModel(mobj, 'CF_HydroData') 
        warndlg('Need hydraulic model, sediment and water level data')
        return;
    end
    cn = getConstantStruct(mobj.Constants);

    Le = hydobj.xTidalLimit;         %distance to tidal limit
    Lx = hydobj.xTideRiver;          %distance to tide/river switch

    %channel sediment properties
    d50 = sedobj.SedimentSize;       %sediment grain size, D50 (m)
    rhoc = sedobj.EqDensity;         %suspended sediment concentration (kg/m3)
    taucr = sedobj.CritBedShear;     %critical bed shear stress (Pa)

    %river sediment properties
    d50riv = sedobj.d50river;        %sediment grain size, D50 (m)
    tauriv = sedobj.tauriver;        %critical bed shear stress (Pa)

    %hydraulic properties
    am = wlvobj.TidalAmp;            %tidal amplitude (m)
    tp = wlvobj.TidalPeriod*3600;    %tidal period (s)
    
    %river properties
    Qr = hydobj.RiverDischarge;      %river discharge (m^3/s)
    Sr  = 2*am/Le;   %energy slope at tidal limit (-); **estimate**

    % calc fall velocity.  Mud Manual, eqn 5.7 including floculation
    ws = settling_velocity(d50,cn.g,cn.rhow,cn.rhos,cn.visc,rhoc);   

    % River width and depth
    [hrv,Wrv,~] = river_regime(Qr,Sr,d50riv,tauriv,cn.rhos,cn.rhow);
    if hrv==0
        warndlg(sprintf('Zero river depth for %dm^3/s discharge in channel_form.m',Qr))
        return;
    end

    % Channel properties
    Arv = hrv*Wrv;
    params = struct('am',am,...                  %tidal amplitude at mouth (m)
                'tp',tp,...                      %tidal period (s)
                'Le',hydobj.xTidalLimit,...      %channel length (m)
                'Uw',hydobj.WindSpeed,...        %wind speed at 10m (m/s)
                'Qr',hydobj.Qr,...               %river discharge (m3/s)
                'hrv',hrv,'Wrv',Wrv,'Arv',Arv,...%from river_regime
                'g',cn.g,'rhow',cn.rhow,'rhos',cn.rhos,'rhoc',rhoc,...
                'taucr',taucr,'d50',d50,'ws',ws,...%see above   
                'tauriv',tauriv,'d50riv',d50riv,...%see above
                'me',sedobj.ErosionRate,...      %erosion rate (kg/N/s)
                'Dsm',sedobj.AvMarshDepth,'Dmx',sedobj.MaxMarshDepth); 
   
    initdepth = 1;  %initial guess of hydraulic depth
    formdims1 = ckfa_form_solver(initdepth,params);
    if isempty(formdims1)
        warndlg('Unable to find solution in channel_form.m')
        return;
    end
    tablefigure('CKFA model','CKFA channel form properties',formdims1);
    % formdims2 = ckfa_form_solver(initdepth,params);
    % if isempty(formdims2)
    %     warndlg('Unable to find solution in channel_form.m')
    %     return;
    % end
    % 
    % gp_form = vertcat(formdims1,formdims2);
    % tablefigure('CKFA model','CKFA channel form properties using two solvers',gp_form);
end