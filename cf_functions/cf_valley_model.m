function [xi,yi,zgrd,yz,Lv,Ls0] = cf_valley_model(obj,isfull)
%
%-------function help------------------------------------------------------
% NAME
%   cf_valley_model.m
% PURPOSE
%   construct idealised channel form using 3D CKFA model
% USAGE
%   [xi,yi,zgrd,Lv,Ls0] = cf_valley_model(obj,isfull)
% INPUTS
%   obj - CF_FormModel class instance
%   isfull - true returns full grid, false half-grid
% OUTPUTS
%   xi - x co-ordinate (m)
%   yi - y co-ordinate (m)
%   zgrd - bed elevation grid (m)
%   Lv - river valley convergence length for longitudinal elevation (m)
%   Ls0 - cross-valley convergene length between mtl and max elevation (m)
% NOTES
%   
% SEE ALSO
%   called by CF_ValleyModel as part of ChannelForm model
%
% Author: Ian Townend
% CoastalSEA (c) Jan 2022
%--------------------------------------------------------------------------
%
    if nargin<3
        isfull = true;
    end
    %get the required input parameter classes
    grdobj = obj.RunParam.GD_GridProps;
    hydobj = obj.RunParam.CF_HydroData;
    valobj = obj.RunParam.CF_ValleyData;
    sedobj = obj.RunParam.CF_SediData;
    
    %model run parameters
    Lt = diff(grdobj.XaxisLimits);   %length of model domain (m)

    %set-up co-ordinate system
    [xi,yi,delx] = getGridDimensions(grdobj);
    yi = yi(yi>=0);  %half the grid
    yix = repmat(yi',length(xi),1); 
    
    zMx = valobj.ValleyEle;          %Maximum elelvation (mOD)
    bv = valobj.ValleyWidth(1)/2;    %half-width of valley at mouth

    %sediment properties
    d50riv = sedobj.d50river;        %sediment grain size, D50 (m)
    tauriv = sedobj.tauriver;        %critical bed shear stress (Pa)

    %model constants
    cns = muiConstants.Evoke;
    rhow = cns.WaterDensity;      %density of water (default = 1025 kg/m^3)
    rhos = cns.SedimentDensity;   %density of sediment (default = 2650 kg/m^3) 
    
    %hydraulic parameters
%     zm0 = hydobj.zmt(end);
%     if isscalar(hydobj.zhw)
%         am = (hydobj.zhw-hydobj.zlw)/2; %tidal amplitude at mouth
%         zhw = zm0+am;
%     else  
%         zhw = hydobj.zhw(end);          %high water level at mouth
%     end
    zhw = hydobj.zhw(1);          %high water level at mouth
    %top of river valley
    coastlevel = 1;
    zmx = (zMx-(zhw+1))/Lt*xi'+zhw+coastlevel;
    zmx = repmat(zmx',1,length(yi));

    %valley properties
    z0 = valobj.ValleyDepth;       %depth of valley at mouth (mAD)
    xr = valobj.xTidalLimit;       %distance to tidal limit (m)
    ztl = valobj.zTidalLimit;      %water surface elevation at TL (mAD)
    xH = valobj.xValleyHead;       %distance to valley head (m)
    zH = valobj.zValleyHead;       %elevation at valley head (mAD)

    %river properties
    Qr = hydobj.Qr;                  %river discharge (m^3/s)
    %examined using river slope (too wide and shallow)
    [zm0,Lv] = CF_ValleyModel.findconvergencelength(xr,ztl-1,xH,zH,z0);
    zv = zm0*(exp(xi'/Lv)-1)+z0;
    gradV = gradient(zv)/delx;
    Sr = interp1(xi',gradV,xr);
    % Sr  = 2*am/xr;   %energy slope at tidal limit (-); **estimate**

    % River depth
    [hrv,Wrv,~] = river_regime(Qr,Sr,d50riv,tauriv,rhos,rhow);
    zr = ztl-hrv;                 %elevation of river bed at TL (mAD)

    %get convergence length for river valley longitudinal elevation
    [zm0,Lv] = CF_ValleyModel.findconvergencelength(xr,zr,xH,zH,z0);
    zv = (zm0*(exp(xi'/Lv)-1)+z0)';
    %get cross-valley convergence length between mtl and max elevation
    Ls0 = bv/log((zhw+zm0-z0)/zm0);

    %compute valley surface
    zi = (zm0*(exp(yix./Ls0)-1)+zv);
    %impose surrounding land surface
    idz = zi>zmx;        
    zi(idz) = zmx(idz);
            
    %initial attempt to add variable river as function of Sr (not
    %working)
    % Sr = interp1(xi,gradV,xr);
    % [hrv,Wrv,~] = river_regime(Qr,gradV,d50riv,tauriv,rhos,rhow);
    % for i=1:length(xi)
    %     [hrv(i),Wrv(i),~] = river_regime(Qr,gradV(i),d50riv,tauriv,rhos,rhow);
    % end
    % 
    % Lc = repmat(Wrv'/2,1,length(yi));
    % idy = yix<=Lc;
    % zbot = max(zi(idy),[],2);                 %elevation on edge of strip
    % zbot = repmat(zbot,1,length(yi));
    % zi(idy) = zbot(idy);                 %flatten bottom of channel

    %add river in bed
    Lc = Wrv/2;
    idy = yix<=Lc;
    zc = reshape(zi(idy),length(xi),[]); %strip width of river
    zbot = max(zc,[],2);                 %elevation on edge of strip
    zbot = repmat(zbot,1,length(yi));
    zi(idy) = zbot(idy);                 %flatten bottom of channel

    mu = 6*hrv/Lc;
    zi(idy) = zi(idy)-(mu*Lc/2*(1-(yix(idy)/Lc).^2)); %add river channel

    %plot half-sections along channel
    plotValley(xi,yi,zi)
    
    grid = struct('x',xi,'y',yi,'z',zi,'ishead',false);
    zwl = {hydobj.zhw(end),hydobj.zmt(end),hydobj.zlw(end)};
    yz = cf_plan_form(grid,zwl);
%     yz = num2cell(flipud(yz)',2);      %formatted to load into dstable
    %generate complete 3D channel form by mirroring half section
    if isfull                          %return full grid
        zgrd = cat(2,fliplr(zi(:,2:end)),zi);
        yi  = [-flipud(yi(2:end)); yi];
    else                               %return half grid
        zgrd = zi;
    end

    Lv = round(Lv); Ls0 = round(Ls0); %report to nearest metre
end
%%
function plotValley(xi,yi,zi)
    %contour plot of valley form
    hf = figure('Name','Valley Plot','Units','normalized','Tag','PlotFig');  
    hf.Position(1) = 0.7;
    plot(yi,zi(1,:))
    hold on
    nx = length(xi);
    for i=1:50:nx
        plot(yi,zi(i,:))
    end
    hold off
    xlabel('Width (m)')
    ylabel('Elevation (m)')
    title('Valley cross-section at intervals along channel';
end
% function [zm0,Lv] = findconvergencelength(xr,zr,xH,zH,z0)
%     %iterative solution for the convergence length of a valley with
%     %defined levels at the head, zH, and tidal limit, zr, and a 
%     %basal level, z0, at the mouth
%     % zm0 - elevation adjustment at mouth (offset to be applied to thalweg)
%     % Lv - along channel convergence length of valley bottom
%     %NB -same as findconvergencelength function in CF_ValleyModel
%     drp = zr-z0;       %depth between river invert at tidal limit and mouth
%     dHp = zH-z0;       %depth between valley head and mouth
%     lv_fun = @(lv) drp/dHp-(exp(xr/lv)-1)/(exp(xH/lv)-1);
%     Lv = fzero(lv_fun,xH/2);
%     zm0 = drp/(exp(xr/Lv)-1);
% end
%%
% function yz = planForm(hydobj,xi,yi,zi)
%     %compute controlling dimensions (y) for an exponential plan form
%     zwl(1) = hydobj.zhw(end); %constant so does not matter which end
%     zwl(2) = hydobj.zmt(end);
%     zwl(3) = hydobj.zlw(end);
%     
%     yz = zeros(length(xi),3);
%     for j=1:3
%         M = contour(yi,xi,zi,[zwl(j),zwl(j)]);
%         xj = M(2,2:end);
%         yj = M(1,2:end);
%         yz(:,j) = interp1(xj,yj,xi,'linear');
%     end
% end