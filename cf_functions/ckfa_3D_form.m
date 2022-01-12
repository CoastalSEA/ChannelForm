function [xi,yi,zi] = ckfa_3D_form(obj,params)
%
%-------function help------------------------------------------------------
% NAME
%   ckfa_3D_form.m
% PURPOSE
%   construct idealised channel form using 3D CKFA model
% USAGE
%   [xi,yi,zgrd,yz] = ckfa_3D_form(obj)
% INPUTS
%   obj - CF_FormModel class instance
%   params - ckfa model parameters defined in ckfa_form_model.ckfa_parameters
% OUTPUTS
%   xi - x co-ordinate (m)
%   yi - y co-ordinate (m)
%   zi - bed elevation grid (m)
% NOTES
%   CKFA cross-section comprises a channel using Cao&Knight section and an
%   intertidal using the profile proposed by Friedrichs & Aubrey (tide only)
% SEE ALSO
%   used in ckfa_form_model which is one of the models called by 
%   CF_FormModel as part of the ChannelForm model 
%
% Author: Ian Townend
% CoastalSEA (c) Jan 2022
%--------------------------------------------------------------------------
%
    nbk = 2; %number of banks is assumed to be 2 in the ChannelForm model
    %get the required input parameter classes
    grdobj = obj.RunParam.GD_GridProps;
    hydobj = obj.RunParam.CF_HydroData;
    
    %model run parameters
    offset = 2*grdobj.histint;       %level above hw for land outside section

    % CKFA form properties from ckfa_form_solver    
    hm = obj.CKFAform.form.hm;   %MTL hydraulic depth at mouth (m) 
    LA  = obj.CKFAform.form.La;  %CSA convergence length (m)
    LW  = obj.CKFAform.form.Lw;  %width convergence length (m)
    Ucr = obj.CKFAform.form.Ucr; %peak tidal velocity

    % Flow only properties
    Slw = obj.CKFAform.flow.Slw; %surface area at low water (m2)
    Sfl = obj.CKFAform.flow.Sfl; %surface area of tidal flats (m2)
    Wm = obj.CKFAform.flow.Wm;   %Width at mouth (m)    
    Am = obj.CKFAform.flow.Am;   %CSA at mouth (m2)
    Wrv = obj.CKFAform.flow.Wrv; %width of regime river (m)
    Arv = obj.CKFAform.flow.Arv; %CSA of regime river (m2)
    % tide+wave properties (used for Qp - commented out)
    % Vpw = obj.CKFAform.wave.Vpw;  %tidal prism including wave effects
    
    % Input parameters
    am = params.am;              %tidal amplitude (m)
    tp = params.tp;              %tidal period (s)    
    Le = params.Le;              %estuary length (m)
    Uw = params.Uw;              %wind speed at 10m (m/s)
    d50 = params.d50;            %sediment grain size, D50 (m)
    tau = params.taucr;          %critical bed shear stress (Pa)
    me = params.me;              %erosion rate (kg/N/s)
    ws = params.ws;              %sediment fall velocity (m/s)
    rhoc = params.rhoc;          %suspended sediment concentration (kg/m3)
    Dsm = params.Dsm;            %average depth over saltmarsh (m)
    Dmx = params.Dmx;            %maximum depth of salt marsh (m)          
    
    % Constant properties
    g = params.g;                %acceleration due to gravity (m/s2)
    rhos = params.rhos;          %density of sediment (default = 2650 kg
    conc = rhoc/rhos;            %volume concentration (-)  
    omega = 2*pi/tp;             %tidal frequency (1/s)
    nc  = 2.0;                   %channel shape factor
    cflag = 0;                   %constraint based on basin area at 
                                 %high water (1=constraint included)    
    % Qp = pi*Vpw/tp/2;              
    % type= 4;          estuary classification type (1-7)
    % if type==2 || type==3 || type==7 %fjards, rias and inlets
    %     Qp = Le/LA*Qp;
    % end
    %set-up co-ordinate system
    [xi,yi] = getGridDimensions(grdobj);
    yi = yi(yi>=0);  %half the grid    
    zi = zeros(length(xi),length(yi));

    %water level properties based on amplitude+mtl or CST model (mAD)
    %ckfa model has origin at mouth, x positive upstream
    zHWxi = hydobj.zhw;           %high water level(mAD)
    zLWxi = hydobj.zlw;           %low water level(mAD) 
    am0 = (zHWxi(1)-zLWxi(1))/2;  %tidal amplitude at mouth

    Slope = 4*am0/tp/sqrt(g*hm);
%     [hlw,wlw,~] = river_regime(Qp,Slope,d50river,tauriver,rhos,rhow);
%     wlwm=wlw; alwm=wlw*hlw;
    eqtr = Slope*Le;
    [~,wlwm,alwm] = get_river_regime(obj,eqtr);
    

    Wbs = 1.2*(Slw+Sfl)/LW/(1-exp(-Le/LW)); %high water width to generate basin area
    Lsll= Ucr/omega;       %lower intertidal width along system axis
    
    for j = 1:length(xi)
        zhw = zHWxi(j); 
        zlw = zLWxi(j);
        ax = (zhw-zlw)/2;               %tidal amplitude(m) at x(j)
        mtl = zhw-ax;                   %mean tide level(mOD)
        phi  = Wm/2/LW*exp(-xi(j)/LW);  %angle to bank at xj
        Lst  = Lsll*tan(phi);           %width of lower intertidal at xj
        meq  = Lst/ax;                  %equilibrium lower intertidal slope at xj
        if isinf(meq), meq = Lst/0.1; end
        
        Wx   = Wm*exp(-xi(j)/LW);       %mtl width at xj
        Wlw  = Wx - nbk*ax*meq;         %flow only width at low water
        Whw  = Wx + nbk*pi/2*ax*meq;    %flow only width at high water
        Wbx  = Wbs*exp(-xi(j)/LW);      %max width at xj
        %
        Ax  = Am*exp(-xi(j)/LA);        %mtl csa at xj 
        Alw = Ax-ax*(Wx+Wlw)/2;         %lw csa at xj
%         hx  = Ax/Wx;                    %hydraulic depth (hAb) at xj
        %
        % Obtain depth and width of wave formed profile at high and low water
        [dhw,yhw,dlw,ylw] = ckfa_wave_form(am,Uw,d50,tau,me,ws,conc,hm,Wlw,Whw);
        %
        % correction for enlarged prism due to influence of waves
        %
        Alx  = alwm*exp(-xi(j)/LA);     %low water area at xj
        Wlx  = wlwm*exp(-xi(j)/LW);     %low water width at xj
        if Wlx>Wlw                      %correct for wave enlarged prism
            Wlw = Wlx;
        end  
        if Wlw<Wrv, Wlw=Wrv; end         %make river regime channel minimum section
        Lc  = Wlw/nbk;                  %half width of channel
        %
        if Alx>Alw                      %correct for wave enlarged prism
            Alw = Alx;
        end 

        if Alw<Arv, Alw=Arv; end         %make river regime channel minimum section
        %
        % calculate channel dimensions
        mu  = Alw*nbk*3/Wlw^nc;        %friction angle for low water channel
        dc  = mu*Lc/2;                 %depth at centre-line of low water channel
        %
        if dlw<dc
            y0lw = Lc*(1+2*(-dlw)/mu/Lc)^(1/nc); %distance from c.l. to dlw
        else
            y0lw = 0;
        end
        Lcw = max(Lc,y0lw+ylw);  %half width at lw including influence of waves
        if dhw<ax
            y0hw = Lcw + Lst*(asin((ax-dhw)/ax)+1);     %distance from cl to dhw
        elseif dhw>ax && dhw<2*ax
            y0hw = Lcw + Lst*((ax-dhw)/ax+1);
        elseif dhw>2*ax && dhw<2*ax+dlw
            y0hw = y0lw + ylw*(1-((dhw-2*ax)/dlw)^1.5);
        elseif dhw>2*ax+dlw && dhw<2*ax+dc
            y0hw = Lc*(1-2*(dhw-2*ax)/mu/Lc)^(1/nc);
        else
            y0hw=0;
        end
        Lhw  = max(Lcw+Lst*(1+pi/2),y0hw+yhw);     %half width at hw incl waves
        %
        %if dhw>Dmx  %wave profile extends beyond initiation depth of saltmarsh
        %this causes switching in the upper reaches where intertidal is very
        %narrow.
        if dhw>0   %wave profile exists
            ysm = yhw*(Dmx/dhw)^1.5;   %distance from Dmx to hw on wave profile
        elseif ax<=Dmx
            ysm = Lst*pi/2;            %distance from mtl to hw on flow profile
        else
            ysm = Lst*(pi/2-asin((ax-Dmx)/ax)); %distance from Dmx to hw on flow profile
        end
        %
        y0sm = (y0hw+yhw)-ysm;      %distance from c.l. to Dmx (saltmarsh)
        %
        %set up profile and calculate z values
        for k = 1:length(yi)
            %
            if  yi(k)<=y0lw && yi(k)<=y0hw
                zi(j,k) = mtl-(ax+mu*Lc/2*(1-(yi(k)/Lc)^nc));    %lw channel
            elseif yi(k)<y0lw && yi(k)>y0hw && yi(k)<y0hw+yhw
                zch    = -(ax+mu*Lc/2*(1-(yi(k)/Lc)^nc));       %lw channel
                zhw    = ax-dhw*(1-(yi(k)-y0hw)/yhw)^(2/3);      %hw wave profile
                zi(j,k) = mtl+min(zch,zhw);
            elseif yi(k)>y0lw && yi(k)<y0lw+ylw && yi(k)<y0hw
                zi(j,k) = mtl-(ax+dlw*(1-(yi(k)-y0lw)/ylw)^(2/3));%lw wave profile
            elseif yi(k)>y0lw && yi(k)>y0lw+ylw && yi(k)<Lcw && yi(k)<y0hw
                zi(j,k) = mtl-(ax+mu*Lc/2*(1-(yi(k)/Lc)^nc));    %lw channel profile above lw wave profile    
            elseif yi(k)>y0lw && yi(k)<Lcw && yi(k)>y0hw && yi(k)<y0hw+yhw && y0hw+yhw>Lcw
                zhw    = ax-dhw*(1-(yi(k)-y0hw)/yhw)^(2/3);      %hw wave profile
                if yi(k)<y0lw+ylw
                zlw    = -(ax+dlw*(1-(yi(k)-y0lw)/ylw)^(2/3));   %lw wave profile
                else
                    zlw=zhw;
                end
                zi(j,k) = mtl+min(zlw,zhw);
            elseif yi(k)>y0lw && yi(k)<Lcw
                zi(j,k) = mtl-(ax+mu*Lc/2*(1-(yi(k)/Lc)^nc));      %lw channel
            elseif yi(k)>Lcw && yi(k)<Lcw+Lst && (yi(k)<y0hw || yi(k)>y0hw+yhw)
                zi(j,k) = mtl+ax*((yi(k)-Lcw)/Lst-1);             %lower tidal flat
            elseif yi(k)>Lcw && yi(k)<Lcw+Lst && yi(k)>y0hw && yi(k)<y0hw+yhw
                zi(j,k) = mtl+ax-dhw*(1-(yi(k)-y0hw)/yhw)^(2/3);  %hw wave profile
            elseif yi(k)>Lcw+Lst && yi(k)<y0hw && yi(k)<y0sm
                zi(j,k) = mtl+ax*sin((yi(k)-Lcw)/Lst-1);          %upper tidal flat
            elseif yi(k)>Lcw+Lst && yi(k)>y0hw && yi(k)<y0hw+yhw && yi(k)<y0sm
                zi(j,k) = mtl+ax-dhw*(1-(yi(k)-y0hw)/yhw)^(2/3);  %hw wave profile
            elseif yi(k)>y0hw+yhw && yi(k)<Lhw && yi(k)<y0sm
                zi(j,k) = mtl+ax*sin((yi(k)-Lcw)/Lst-1);          %upper tidal flat
            elseif yi(k)>y0sm && yi(k)<Lhw && Dmx>0
                zi(j,k) = mtl+ax-Dsm;                          %horizontal marsh surface
            else
                zi(j,k) = mtl+ax+offset;                       %values outside section
            end
            if ~isreal(zi(j,k))
                zi(j,k) = 9.999;
            end
            %limit profile to size of basin
            if yi(k)>Wbx/2 && cflag==1
                zi(j,k) = mtl+ax+offset;
            end
        end
    end   
end
    
    
    