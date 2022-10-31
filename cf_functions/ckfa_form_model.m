function [xi,yi,zgrd,Wz,Rv] = ckfa_form_model(obj,isfull)
%
%-------function help------------------------------------------------------
% NAME
%   ckfa_form_model.m
% PURPOSE
%   construct idealised channel form using 3D CKFA model
% USAGE
%   [xi,yi,zgrd,yz] = ckfa_form_model(obj,isfull)
% INPUTS
%   obj - CF_FormModel class instance
%         obj.Selection.wlflag - indicates type of water surface to use
%            0=CSTmodel used to define water levels
%            1=constant HW tapering LW 
%            2=constant HW & LW
%   isfull - true returns full grid, false half-grid
% OUTPUTS
%   xi - x co-ordinate (m) origin at head/river
%   yi - y co-ordinate (m) origin on centre-line
%   zgrd - bed elevation grid (m)
%   Wz -  table of Whw,Wmt,Wlw for width at hw,mt,lw (m)
%   Rv - struct of river regime properties Hr, Wr, Ar
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
    xi = []; yi = []; zgrd = []; Wz = [];
    if nargin<3
        isfull = true;
    end

    % Channel properties
    params = ckfa_properties(obj);

    %set the water level variations along the estuary
    obj.CSTparams = params.form;
    [obj,ok] = cf_set_hydroprops(obj);
    if ok<1, return; end
    
    Rv.Hr = params.input.hav; 
    Rv.Wr = params.input.Wrv;
    Rv.Ar = params.input.WArv;
    %generate 3D surface of CKFA model
    [xi,yi,zi] = ckfa_3D_form(obj,params);
    if isempty(xi),return; end
    
    %model x-axis is from mouth. no need to reverse data for use in ChannelForm
    [yu,yo,yl] = expPlan(params,xi);
    Wz = [yu,yo,yl]*2;      %full width  
    yzcell = num2cell(Wz',2)';  %formatted to load into table
    Wz = table(yzcell{:},'VariableNames',{'Whw','Wmt','Wlw'});
    %generate complete 3D channel form by mirroring half section
    if isfull                              %return full grid
        zgrd = cat(2,fliplr(zi(:,2:end)),zi);
        yi  = [-flipud(yi(2:end)); yi];
    else                                   %return half grid
        zgrd = zi;
    end
%     zgrd = flipud(zgrd); %make orientation consistent with other models
end
%%
function params = ckfa_properties(obj)
    %get the input parameters and then call the ckfa model components for
    %form, flow and tide+wave properties
    params.input = ckfa_parameters(obj);          
    %call solver function 
    initdepth = 1;  %initial guess of hydraulic depth
    params.form = ckfa_form_solver(initdepth,params.input);
    % flow only gross properties
    params.flow = ckfa_flowprops(params);
    % tide+wave gross properties (incl saltmarsh if included)
    params.wave = ckfa_waveprops(params);
end
%%
function params = ckfa_parameters(obj)
    %intialise the properties required for the CKFA_solver function
    
    %default model constants (NB: not as held and modifiable in UI)
    cn = getConstantStruct(muiConstants.Evoke);    
    
    sedobj = obj.RunParam.CF_SediData;    
    %channel sediment properties
    d50 = sedobj.SedimentSize;    %sediment grain size, D50 (m)
    rhoc = sedobj.EqDensity;      %suspended sediment concentration (kg/m3)
    taucr = sedobj.CritBedShear;  %critical bed shear stress (Pa)
    d50riv = sedobj.d50river;%sediment grain size in river (m)
    tauriv = sedobj.tauriver;%critical bed shear stress (Pa)
    
    % calc fall velocity.  Mud Manual, eqn 5.7 including floculation
    ws = settling_velocity(d50,cn.g,cn.rhow,cn.rhos,cn.visc,rhoc);    

    %parameters required by solver
    hydobj = obj.RunParam.CF_HydroData;
    am0 = (hydobj.zhw(1)-hydobj.zlw(1))/2;    
    [hrv,Wrv,Arv] = get_river_regime(obj,2*am0); %Cao & Knight, 1996, uses Qr in RunParam.CF_HydroData
    
    params = struct('am',am0,...                     %tidal amplitude at mouth (m)
                    'tp',hydobj.tidalperiod,...      %tidal period (s)
                    'Le',hydobj.xTidalLimit,...      %channel length (m)
                    'Uw',hydobj.WindSpeed,...        %wind speed at 10m (m/s)
                    'Qr',hydobj.Qr,...               %river discharge (m3/s)
                    'hrv',hrv,'Wrv',Wrv,'Arv',Arv,...%from river_regime
                    'g',cn.g,'rhow',cn.rhow,'rhos',cn.rhos,'rhoc',rhoc,...
                    'taucr',taucr,'d50',d50,'ws',ws,...%see above  
                    'tauriv',tauriv,'d50riv',d50riv,...%see above
                    'me',sedobj.ErosionRate,...      %erosion rate (kg/N/s)
                    'Dsm',sedobj.AvMarshDepth,'Dmx',sedobj.MaxMarshDepth);        
end
%%
function output = ckfa_flowprops(params)
    %flow only - gross properties using the CKFA model
    nbk = 2; %number of banks is assumed to be 2 in the ChannelForm model
    
    %CKFA form properties from ckfa_form_solver    
    hm = params.form.hm;           %MTL hydraulic depth at mouth (m)
    LW = params.form.Lw;           %e-folding length for width (m)
    Wm = params.form.Wm;           %MTL Width at mouth (m)
    LA = params.form.La;           %e-folding length for CSA (m2)
    Am = params.form.Am;           %MTL CSA at mouth (m2)
    Ucr = params.form.Ucr;         %peak tidal amplitude (m/s)    
    
    %Input parameters
    am = params.input.am;          %tidal amplitude (m)
    tp = params.input.tp;          %tidal period (s)
    Le = params.input.Le;          %estuary length (m)
    hrv = params.input.hrv;        %hydraulic depth of regime river channel (m)
    Wrv = params.input.Wrv;        %width of regime river channel (m)

    omega = 2*pi/tp;               %tidal frequency (1/s)
    g = muiConstants.Evoke.Gravity;
    %
    % Lower intertidal slope
    phi = Wm/2/Le*(1-exp(-Le/LW)); %average angle of bank to direction of flow
    Lst = Ucr/omega*tan(phi);      %lower intertidal width
    meq = Lst/am;                  %lower intertidal slope
    mSeq= meq*Le;                  %lower intertidal slope-area
    %
    % Tidal prism
    etaA = omega*LA/sqrt(g*hm);    %wave number ratio
    etah = am/hm;                  %amplitude-depth ratio
    Vph = Ucr*Wm*hm/2/omega*(4-pi*etah*sin(etaA));  %prism
    %
    % Gross properties for flow only
    Soc = Wm*LW*(1-exp(-Le/LW));   %surface area at mtl
    Voc = Am*LA*(1-exp(-Le/LA));   %volume at mtl
    Slw = Soc - nbk*am*mSeq;       %surface area at low water
    if Slw<0, Slw=Wrv*Le; end      %minimum area based on river section
    Vlw = Voc-am*(Soc+Slw)/2;      %volume at low water
    if Vlw<0, Vlw=Wrv*hrv*Le; end  %minimum volume based on river section
    Sfl = (1+pi/2)*nbk*am*mSeq;    %surface area of tidal flat
    Vfl = (3+pi)/2*nbk*am^2*mSeq;  %volume of tidal flat %corrected 9/12/21, was (1+pi)/2
    Vp  = 2*am*Slw+Vfl;            %form based prism
    Arv = Wrv*hrv;                 %CSA of the regime river 
    %
    % Set up output array for flow only
    output = table(Slw,Vlw,Soc,Voc,Sfl,Vfl,Vp,Wm,Am,Wrv,Arv,Lst,Vph);
end
%%
function output = ckfa_waveprops(params)
    %flow+waves - gross properties using the CKFA model
    nbk = 2; %number of banks is assumed to be 2 in the ChannelForm model
    type= 4; %estuary classification type (1-7) - assumed fixed in ChannelForm
    
    % CKFA form properties from ckfa_form_solver    
    hm = params.form.hm;         %MTL hydraulic depth at mouth (m) 
    LA  = params.form.La;        %CSA convergence length (m)
    LW  = params.form.Lw;        %width convergence length (m)
    
    % Flow only properties
    Slw = params.flow.Slw;       %surface area at low water (m2)
    Vlw = params.flow.Vlw;       %volume at low water (m3)
    Sfl = params.flow.Sfl;       %surface area of tidal flats (m2)
    Vfl = params.flow.Vfl;       %volume of tidal flats (m3)
    Vp  = params.flow.Vp;        %tidal prism volume (m3)    
    Lst = params.flow.Lst;       %width of lower intertidal - LW to MT (m)
    
    % Input parameters
    am = params.input.am;        %tidal amplitude (m)
    tp = params.input.tp;        %tidal period (s)    
    Le = params.input.Le;        %estuary length (m)
    Uw = params.input.Uw;        %wind speed at 10m (m/s)
    d50 = params.input.d50;      %sediment grain size, D50 (m)
    tau = params.input.taucr;    %critical bed shear stress (Pa)
    me = params.input.me;        %erosion rate (kg/N/s)
    ws = params.input.ws;        %sediment fall velocity (m/s)
    rhoc = params.input.rhoc;    %suspended sediment concentration (kg/m3)
    Dsm = params.input.Dsm;      %average depth over saltmarsh (m)
    Dmx = params.input.Dmx;      %maximum depth of salt marsh (m)
    d50riv = params.input.d50riv;%sediment grain size in river (m)
    tauriv = params.input.tauriv;%critical bed shear stress (Pa)
    
    % Constant properties
    g = params.input.g;          %acceleration due to gravity (m/s2)
    rhow = params.input.rhow;    %density of water (default = 1025 kg/m^3)
    rhos = params.input.rhos;    %density of sediment (default = 2650 kg

    % Initialise model settings
    if Dmx==0, ism = false; else, ism = true; end %include saltmarsh if Dmx>0
    conc = rhoc/rhos;            %volume concentration (-)   
    meq = Lst/am;                %lower intertidal slope
    mSeq= meq*Le;                %lower intertidal slope-area
    
    % Obtain depth and width of wave formed profile at high and low water
    Wlw = Slw/Le;
    Whw = (Slw+Sfl)/Le;           %width at high water
    [dhw,yhw,dlw,ylw] = ckfa_wave_form(am,Uw,d50,tau,me,ws,conc,hm,Wlw,Whw);
    Whw = (Slw+Sfl)/Le;          %width at high water
    %
    % Area and volume adjustments due to wave effects
    Lc = Slw/Le/nbk;             %flow only half-width of channel
    mu = Vlw*Le*nbk*3/Slw^2;     %low water channel aspect ratio, mu
    dc = mu*Lc/2;                %centre-line depth of idealised low water channel
    % Low water adjustments
    y0lw = 0;
    %
    if dlw<dc, y0lw = Lc*(1+2*(-dlw)/mu/Lc)^0.5; end  %distance from c.l. to dlw
    dLlw = y0lw+ylw-Lc;          %increase in channel width due to waves
    %
    if dhw>2*am %handle case where hw wave profile is greater than tidal range
        y0hwc= Lc*(1+2*(-(dhw-2*am))/mu/Lc)^0.5; %distance from cl to dhw
        if y0hwc<0, y0hwc=0; end
        ylwc = yhw*(1-(2*am/dhw)^1.5); %distance from dhw to 2a intersection on wave profle
        dLlw = y0hwc+ylwc-Lc;    %increase in channel width due to waves
    end 
    %
    if dLlw<0, dLlw=0; end
    Slww = Slw;
    Vlww = Vlw;
    %
    if dLlw>0, Slww = Slww+nbk*dLlw*Le; end  %adjusted channel area
    %
    Awlw = nbk*3/5*ylw*dlw;                 %csa of wave profile at lw
    Aclw = nbk*dc/3*(2*Lc-y0lw*((3*Lc^2-y0lw^2)/Lc^2)); %area of flow only profile from y0lw to Lc
    %
    if (Awlw-Aclw)>0, Vlww = Vlw+(Awlw-Aclw)*Le; end     %adjusted channel volume
    if dhw>2*am %handle case where hw wave profile is greater than tidal range
        %Aclc = nbk*dc/3*(2*Lc-y0hwc*((3*Lc^2-y0hwc^2)/Lc^2)); %area of lw profile from y0hwc to Lc
        Awlc = nbk*3/5*yhw*dhw*(1-((yhw-ylwc)/yhw)^(5/3)); 
        Awlc = Awlc-nbk*2*am*ylwc; %area of hw profile from y0hwc to lw
        ylww = ylw*(1-((dhw-2*am)/dlw)^1.5); %distance from dlw to r of (dhw-2a) wiht lw wave profile
        Alww = nbk*3/5*ylw*dlw*(1-((ylw-ylww)/ylw)^(5/3));
        Alww = Alww-nbk*(dhw-2*am)*ylww; %area of lw profile from y0lw to intersection of (dhw-2a)
        dAlw = nbk*(y0hwc-y0lw)*(dhw-2*am); %area of box between y0lw and y0hwc
        Vlww = Vlw+(Awlc+dAlw+Alww-Aclw)*Le;
    end
    % High water adjustments
    y0hwu= Lst*(asin((am-dhw)/am)+1);  %distance from lw to dhw if dhw<a
    y0hwl= Lst*((am-dhw)/am+1);        %distance from lw to dhw if dhw>a
    %
    if dhw>am, y0hw = y0hwl; else y0hw = y0hwu; end
    Whww = nbk*(y0hw+yhw+(y0lw+ylw)); %width at high water
    if dhw>2*am,  Whww = nbk*(y0hwc+yhw); end
    %
    if Whww<Whw, Whww = Whw; end
    Sflw = Whww*Le-Slww;               %adjusted flat area
    Awhw = nbk*3/5*yhw*dhw;            %csa of wave profile at hw
    ycu = (1+pi/2)*Lst-y0hwu;   %width of flow profile replaced by wave profile
    ycl = Lst-y0hwl;   %width of flow profile below z=0 replaced by wave profile
    %
    if dhw>am, yc = ycl; else yc = ycu; end
    Acu = nbk*am*(yc-Lst*(0.54*cos(y0hw/Lst)+0.841*sin(y0hw/Lst)));%area of flow only profile from y0hw to Li if dhw<a
    Acl = nbk*((pi/2-1)*am*Lst+((am+dhw)/2)*yc);%area of flow only profile from y0hw to Li if dhw>a
    %
    if dhw>am Achw = Acl; else Achw = Acu; end
    Vflw = Vfl;
    %
    if (Awhw-Achw)>0, Vflw = Vflw+(Awhw-Achw)*Le; end  %adjusted flat volume
    if dhw>2*am 
        Vflw = (Awhw-Awlc-nbk*(Lc-y0hwc)*2*am)*Le; 
    end

    % Calculate mtl area
    Sow = Slww+nbk*am*mSeq;   %adjusted mtl area when dhw<a
    %Wmw = Sow/Le;
    if dhw>2*am %hw wave profile is greater than tidal range
               %y0hwc is distance from cl to dhw, when dhw>2a (see above)
        ymtl = yhw*(1-(am/dhw)^1.5);%distance from dhw to a intersection on wave profle
        ymtl = y0hwc+ymtl;          %mtl half-width inc influence of waves
        Sow  = nbk*Le*ymtl;         %mtl area inc influence of waves
    elseif dhw>am %hw wave profile is greater than tidal amplitude
        y0ow = (dLlw+Lc)+(2*am-dhw)*meq;   %distance from cl to dhw, when dhw>a
        ymtl = yhw*(1-(am/dhw)^1.5);%distance from dhw to a intersection on wave profle
        ymtl = y0ow+ymtl;           %mtl half-width inc influence of waves
        Sow  = nbk*Le*ymtl;         %mtl area inc influence of waves
    end
    %
    Vpw = 2*am*Slww+Vflw;           %form based estimate of prism
    % Correct Low water volume for effect of increased prism (pro-rata)
    % [should be able to do this using prism to estimate hm and whence h(lw)]
    Qp = pi*Vpw/tp/2;
    if type==2 || type==3 || type==7 %fjards, rias and inlets
        Qp = Le/LA*Qp;
    end
    if Vpw>Vp
        Slope = 4*am/tp/sqrt(g*hm);
        [hom,Wom,~] = river_regime(Qp,Slope,d50riv,tauriv,rhos,rhow);
        Aav  = (hom*Wom)*LA*(1-exp(-Le/LA))/Le;
        Vlww = Aav*Le;
        mu=3*nbk*Vlw*Le/Slw^2;
        Slww = sqrt(3*nbk*Vlww*Le/mu);
    end

    if dLlw>0 && type==7, Slww = Slww+nbk*dLlw*Le; end  %adjusted channel area
    Vow = Vlww+am*(Slww+Sow)/2;                         %adjusted mtl volume
    Vpw = 2*am*Slww+Vflw; 
    %
    Wmw = Le*Sow/Le/LW/(1-exp(-Le/LW));                 %width at mouth
    Amw = Le*Vow/Le/LA/(1-exp(-Le/LA));                 %CSA at mouth   
    %
    % Correct tidal flat area and volume if a saltmarsh is present and
    % calculate area and volume of the profile occupied by the marsh
    %
    if ism
        %find values for flow only form (code as per bio_init.m)
        usf = Sfl/Le;       %area per unit length
        uvf = Vfl/Le;       %volume per unit length
        %
        usm = nbk.*am.*meq.*(3.142/2-asin((am-Dmx)/am)); %marsh area per unit length
        sfm = (usf-usm)/nbk/am/meq;
        uvm = am*(usm-nbk*am*meq*(0.54*cos(sfm)+0.841*sin(sfm))); % unvegetated marsh volume per unit length
        if Dmx<=dhw %if marsh is within wave profile
            %find values with inclusion of waves
            ysm  = yhw*(Dmx/dhw)^1.5;   %distance from Dmx to hw
            Ssmw = nbk*Le*ysm;          %area of wave profile occupied by marsh
            Vsw  = 3/5*nbk*dhw*yhw*Le*(ysm/yhw)^(5/3); %volume of unvegeated wave profile occupied by marsh
            Vsmw = Ssmw*Dsm;            %wet volume of marsh
            Sflw = Sflw - Ssmw;         %flat area with marsh area removed
            Vflw = Vflw - Vsw;          %flat volume with marsh volume removed
        else %if marsh extends beyond wave profile revert to flow only form
            Ssmw = usm*Le;              %area of marsh on flow only profile
            Vsmw = Ssmw*Dsm;            %wet volume of marsh
            Sflw = Sfl - Ssmw;          %flat area with marsh area removed
            Vflw = Vfl - uvm*Le;        %flat volume with unvegetated marsh volume removed
            %uses Sfl and Vfl rather than Sflw and Vflw because assume that wave
            %profile is infilled by marsh sediment
        end
    else
        Ssmw = 0; Vsmw = 0;
    end
    %
    %alternative output V,S for the wave profile instead of low water channel
    %propc = [Ssmw,Vsmw,Soc,Voc,Sfl,Vfl,Vp];  
    output = table(Slww,Vlww,Sow,Vow,Sflw,Vflw,Ssmw,Vsmw,Vpw,Wmw,Amw);
end
%%
function [yu,yo,yl] = expPlan(params,xi)
    %compute controlling dimensions (y) for an exponential plan form                                        
    Wm = params.flow.Wm;
    Wmw = params.wave.Wmw;
    Lst = params.flow.Lst;  %width of lower intertidal - LW to MT (m)
    Wm = max(Wm,Wmw);
    bl = (Wm-2*Lst)/2;
    bu = (Wm+2*1.57*Lst)/2;
    bh = params.flow.Wrv/2;
    Lw = params.form.Lw;

    Ls = (bu-bl)/2.57;            %Lstar in F&A, lower intertial width (lw to mtl)(m) 
    bo = bl+Ls;                   %half-width at mtl(m)
    yu = (bu-bh)*exp(-xi/Lw)+bh;  %distance from centre line to hw
    yo = (bo-bh)*exp(-xi/Lw)+bh;  %distance from centre line to mtl
    yl = (bl-bh)*exp(-xi/Lw)+bh;  %distance from centre line to lw       
end