function [x,y,z] = ckfa3Dform(inputs,Le,me,ws,conc,gp_form,gp_flow,gp_wave,cstWL,nbk)
%generate surface for 3D form based on CKFA model (incl waves and marsh)
rnp = inputs{2};
hyd = inputs{3};
sed = inputs{4};
cns = inputs{5};

%model run parameters
Lt = rnp.ModelLength; %length of model domain (m)
nintx=rnp.intx;       %no of intervals in the x direction
bt = rnp.ModelWidth/2;%half width of model domain (m)
ninty=rnp.inty;       %no of intervals in the y direction

zo  = hyd.MTLatMouth;           %mean tide level at mouth (mOD)
am  = hyd.TidalAmplitude;       %tidal amplitude (m)
tp  = hyd.TidalPeriod*3600;     %tidal period (s) 
d50river = sed.d50river;         %sediment grain size, D50 (m)
tauriver = sed.tauriver;         %critical bed shear stress (Pa)
Dsm = sed.AvMarshDepth;         %average depth over saltmarsh (m)
Dmx = sed.MaxMarshDepth;        %maximum depth of salt marsh (m)
type= 4;          %estuary classification type (1-7)********************
g = cns.Gravity;              %acceleration due to gravity (m/s2)
rhow = cns.WaterDensity;      %density of water (default = 1025 kg/m^3)
rhos = cns.SedimentDensity;   %density of sediment (default = 2650 kg/m^3)
Slw = gp_flow.Slw;
Sfl = gp_flow.Sfl;
Vpw = gp_wave.Vpw;  %tidal prism including wave effects
Ucr = gp_form.Ucr;  %peak tidal velocity
hm  = gp_form.hm;   %average hydraulic depth
Wrv = gp_form.Wrv;  %regime width of river
hrv = gp_form.hrv;  %regime hydraulic depth of river
Wm  = gp_flow.Wm;   %width at mouth  ????????????????why not Wmw?
Am  = gp_flow.Am;   %area at mouth
LA  = gp_flow.LA;
LW  = gp_flow.LW;
Lst = gp_flow.Lst;

omega = 2*pi/tp;    %tidal frequency (1/s)
% profile plotting controls
csf  = 2.0;         %channel shape factor 
offset = 2*rnp.histint; %level above am for land outside section
% constraint based on basin area at high water (1=constraint included)
cflag = 0; 
Qp = pi*Vpw/tp/2;
if type==2 || type==3 || type==7 %fjards, rias and inlets
    Qp = Le/LA*Qp;
end
Slope = 4*am/tp/sqrt(g*hm);
% Slope = 2*am/Le;
[hlw,wlw,~] = river_regime(Qp,Slope,d50river,tauriver,rhos,rhow);
wlwm=wlw; alwm=wlw*hlw;

Wbs = 1.2*(Slw+Sfl)/LW/(1-exp(-Le/LW)); %high water width to generate basin area
%
Lsll= Ucr/omega;                %lower intertidal width along system axis
%
dx = Lt/nintx;                    %along estuary x increment
xj = 0:dx:Lt;                     %x-intervals along estuary
for j = 1:length(xj)
    zhw = cstWL(j,1); 
    zlw = cstWL(j,2);
    ax = (zhw-zlw)/2;               %tidal amplitude(m) at x(j)
    mtl = zhw-ax+zo;                %mean tide level(mOD)
    phi  = Wm/2/LW*exp(-xj(j)/LW);  %angle to bank at xj
    Lst  = Lsll*tan(phi);           %width of lower intertidal at xj
    meq  = Lst/ax;                  %equilibrium lower intertidal slope at xj

    Wx   = Wm*exp(-xj(j)/LW);       %mtl width at xj
    Wlw  = Wx - nbk*ax*meq;         %flow only width at low water
    Whw  = Wx + nbk*pi/2*ax*meq;    %flow only width at high water
    Wbx  = Wbs*exp(-xj(j)/LW);      %max width at xj
    %
    Ax  = Am*exp(-xj(j)/LA);        %mtl csa at xj 
    Alw = Ax-ax*(Wx+Wlw)/2;         %lw csa at xj
    hx  = Ax/Wx;                    %hydraulic depth (hAb) at xj
    %
    % Obtain depth and width of wave formed profile at high and low water
    [dhw,yhw,dlw,ylw] = waveXdims(hyd,sed,ax,hx,Wlw,Whw,me,ws,conc);
    %
    % correction for enlarged prism due to influence of waves
    %
    Alx  = alwm*exp(-xj(j)/LA);     %low water area at xj
    Wlx  = wlwm*exp(-xj(j)/LW);     %low water width at xj
    if Wlx>Wlw                      %correct for wave enlarged prism
        Wlw = Wlx;
    end  
    if Wlw<Wrv, Wlw=Wrv; end         %make river regime channel minimum section
    Lc  = Wlw/nbk;                  %half width of channel
    %
    if Alx>Alw                      %correct for wave enlarged prism
        Alw = Alx;
    end 
    Arv = Wrv*hrv;
    if Alw<Arv, Alw=Arv; end         %make river regime channel minimum section
    %
    % calculate channel dimensions
    mu  = Alw*nbk*3/Wlw^csf;        %friction angle for low water channel
    dc  = mu*Lc/2;                 %depth at centre-line of low water channel
    %
    if dlw<dc
        y0lw = Lc*(1+2*(-dlw)/mu/Lc)^(1/csf); %distance from c.l. to dlw
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
        y0hw = Lc*(1-2*(dhw-2*ax)/mu/Lc)^(1/csf);
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
    dy = bt/ninty;                          %horizontal y-interval
    yk = 0:dy:bt;   yk(1)=0.01;             %y-intervals over half-section
    for k = 1:length(yk)
        %
        if  yk(k)<=y0lw && yk(k)<=y0hw
            z(j,k) = mtl-(ax+mu*Lc/2*(1-(yk(k)/Lc)^csf));    %lw channel
        elseif yk(k)<y0lw && yk(k)>y0hw && yk(k)<y0hw+yhw
            zch    = -(ax+mu*Lc/2*(1-(yk(k)/Lc)^csf));       %lw channel
            zhw    = ax-dhw*(1-(yk(k)-y0hw)/yhw)^(2/3);      %hw wave profile
            z(j,k) = mtl+min(zch,zhw);
        elseif yk(k)>y0lw && yk(k)<y0lw+ylw && yk(k)<y0hw
            z(j,k) = mtl-(ax+dlw*(1-(yk(k)-y0lw)/ylw)^(2/3));%lw wave profile
        elseif yk(k)>y0lw && yk(k)>y0lw+ylw && yk(k)<Lcw && yk(k)<y0hw
            z(j,k) = mtl-(ax+mu*Lc/2*(1-(yk(k)/Lc)^csf));    %lw channel profile above lw wave profile    
        elseif yk(k)>y0lw && yk(k)<Lcw && yk(k)>y0hw && yk(k)<y0hw+yhw && y0hw+yhw>Lcw
            zhw    = ax-dhw*(1-(yk(k)-y0hw)/yhw)^(2/3);      %hw wave profile
            if yk(k)<y0lw+ylw
            zlw    = -(ax+dlw*(1-(yk(k)-y0lw)/ylw)^(2/3));   %lw wave profile
            else
                zlw=zhw;
            end
            z(j,k) = mtl+min(zlw,zhw);
        elseif yk(k)>y0lw && yk(k)<Lcw
            z(j,k) = mtl-(ax+mu*Lc/2*(1-(yk(k)/Lc)^csf));      %lw channel
        elseif yk(k)>Lcw && yk(k)<Lcw+Lst && (yk(k)<y0hw || yk(k)>y0hw+yhw)
            z(j,k) = mtl+ax*((yk(k)-Lcw)/Lst-1);             %lower tidal flat
        elseif yk(k)>Lcw && yk(k)<Lcw+Lst && yk(k)>y0hw && yk(k)<y0hw+yhw
            z(j,k) = mtl+ax-dhw*(1-(yk(k)-y0hw)/yhw)^(2/3);  %hw wave profile
        elseif yk(k)>Lcw+Lst && yk(k)<y0hw && yk(k)<y0sm
            z(j,k) = mtl+ax*sin((yk(k)-Lcw)/Lst-1);          %upper tidal flat
        elseif yk(k)>Lcw+Lst && yk(k)>y0hw && yk(k)<y0hw+yhw && yk(k)<y0sm
            z(j,k) = mtl+ax-dhw*(1-(yk(k)-y0hw)/yhw)^(2/3);  %hw wave profile
        elseif yk(k)>y0hw+yhw && yk(k)<Lhw && yk(k)<y0sm
            z(j,k) = mtl+ax*sin((yk(k)-Lcw)/Lst-1);          %upper tidal flat
        elseif yk(k)>y0sm && yk(k)<Lhw && Dmx>0
            z(j,k) = mtl+ax-Dsm;                          %horizontal marsh surface
        else
            z(j,k) = mtl+ax+offset;                       %values outside section
        end
        if ~isreal(z(j,k))
            z(j,k) = 9.999;
        end
        %limit profile to size of basin
        if yk(k)>Wbx/2 && cflag==1
            z(j,k) = mtl+ax+offset;
        end
    end
end

%generate histogram and 3D plots
%plotCKFAhist(z,am)
%
x = xj;
y  = [-fliplr(yk), yk];
z = flipud(cat(2,fliplr(z),z));
%plot3Dform(x,y,z)
end
%--------------------------------------------------------------------------
function [dhw,yhw,dlw,ylw] = waveXdims(hyd,sed,ax,hx,Wlw,Whw,me,ws,conc)
    % Obtain depth and width of wave formed profile at high and low water
    d50 = sed.SedimentSize;   %sediment grain size, D50 (m)
    tau = sed.CritBedShear; %critical bed shear stress (Pa)
    Uw = hyd.WindSpeed;       %wind speed at 10m (m/s)
    zw   = 10;                %elevation of wind speed (m)
    %
    Flw  = sqrt(2)*Wlw;            %fetch at low water
    hlw  = hx-ax;                  %hydraulic depth at low water
    if hlw>0
        cnc = conc*hx/hlw;         %modified concentration at low water
        %depth and width of wave profile
        [dlw, ylw] = wave_profile(Uw,zw,Flw,cnc,hlw,hlw,d50,tau,me,ws);
    else
        dlw=0; ylw=0;
    end
    %
    Fhw  = sqrt(2)*Whw;            %fetch at high water
    hhw  = hx+ax;                  %hydraulic depth at high water
    %conc = 5*0.0001*rhow/rhos*Ucr;  %SSC based on eqn of Prandle etal, GRL, 2005
    cnc = conc*hx/hhw;             %modified concentration at high water
    %depth and width of wave profile
    [dhw, yhw] = wave_profile(Uw,zw,Fhw,cnc,hhw,(dlw+2*ax),d50,tau,me,ws);       
end
%--------------------------------------------------------------------------
function plot3Dform(x,y,z)
    %surface plot of the resultant 3D form
    hf = figure('Name','Results Plot','Tag','PlotFig','Units','normalized'); 
    hf.Position(1) = 0.5;                                    
    surf(y,x,z, 'FaceColor','interp','EdgeColor', 'none');
    ylabel('Distance from mouth (m)');
    xlabel('Width (m)');
    zlabel('Elevation (mOD)')
end
%--------------------------------------------------------------------------
function plotCKFAhist(z,am)
    %histogram of the resultant hypsometry
    hf = figure('Name','Results Plot','Tag','PlotFig','Units','normalized'); 
    hf.Position(1) = 0.5;   
    % range for histogram data
    lowlim = floor(min(min(z))); uplim = ceil(max(max(z))); histint = 0.2;
    xedge = lowlim:histint:uplim;
    % calculate histogram and format output
    zed  = reshape(z,numel(z),1);
    zed  = zed(zed<am);
    zed  = zed(zed>=lowlim);
    zhist = histcounts(zed,xedge,'Normalization','probability');
    xcent = xedge(1:end-1)+histint/2;
    barh(xcent, zhist*100, 'histc');
    xlabel('Probability (%)')
    ylabel('Elevation (mOD)')
end