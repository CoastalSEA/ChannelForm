function [xi,yi,zgrd,Wz,Rv] = cf_exp_models(obj,isfull)
%
%-------function help------------------------------------------------------
% NAME
%   cf_exp_models.m
% PURPOSE
%   construct idealised channel form using exponential functions in y to 
%   determine width and cross-section of various forms to determine z at each x interval
% USAGE
%   [xi,yi,zgrd,yz,Rv] = cf_exp_models(obj,isfull)
% INPUTS
%   obj - CF_FormModel class instance
%         obj.Selection.wlflag - indicates type of water surface to use
%            0=CSTmodel used to define water levels
%            1=constant HW tapering LW 
%            2=constant HW & LW
%   isfull - true returns full grid, false half-grid
% OUTPUT
%   xi - x co-ordinate (m)
%   yi - y co-ordinate (m)
%   zgrd - bed elevation grid (m)
%   Wz -  table of Whw, Wmt, Wlw for width at hw, mt, lw (m)
%   Rv - struct of river regime properties Hr, Wr, Ar
% NOTES
%   Cross-section comprises a channel using Cao&Knight section, or a 
%   rectangular prism and an intertidal using linear, stepped, uniform 
%   shear or muddy shore profiles and with an exponential plan form.
% SEE ALSO
%   used in CF_FormModel as part of ChannelForm model
%
% Author: Ian Townend
% CoastalSEA (c) Jan 2022
%--------------------------------------------------------------------------
%
    xi = []; yi = []; zgrd = []; Wz = []; Rv = [];
    if nargin<3
        isfull = true;
    end
    
    %channel properties
    if obj.Selection.wlflag==0 && isempty(obj.CSTparams)
        %provides initial guess of gross properties if cst_model called
        obj.Selection.wlflag = 1;       %fixed water level surface 
        obj= cf_set_hydroprops(obj);     
        obj.Selection.wlflag = 0;       %reset to cst surface
        obj = channel_properties(obj); 
    end
    
    %set the water level variations along the estuary
    [obj,ok] = cf_set_hydroprops(obj); 
    if ok<1
        fprintf('No water level surface found. Grid not defined\n')
        return; 
    end

    [xi,yi,zi,Wz,Rv] = channel_3D_form(obj);
    if isempty(xi),return; end    

    %model x-axis is from head. Reverse data for use in ChannelForm
    yzcell = num2cell(flipud(Wz)',2);  %formatted to load into table
    Wz = table(yzcell{:},'VariableNames',{'Whw','Wmt','Wlw'});
    %x is defined from head with origin at "shoulder". 
    %Change to origin at mouth with grid ordered from x=0
    xi = fliplr(max(xi)-xi);
    %generate complete 3D channel form by mirroring half section
    if isfull                          %return full grid
        zgrd = flipud(cat(2,fliplr(zi(:,2:end)),zi));
        yi  = [-flipud(yi(2:end)); yi];
    else                               %return half grid
        zgrd = zi;
    end
    
%     [X,Y] = meshgrid(xi,yi);
%     figure; surf(X,Y,zgrd);
end    
%%
function obj = channel_properties(obj)
    %compute the summary gross properties for the channel form
    expobj = obj.RunParam.CF_FormData;
    hydobj = obj.RunParam.CF_HydroData;
    
    Whw = expobj.HWmouthWidth;      %width of mouth at high water(m)
    Wlw = expobj.LWmouthWidth;      %width of mouth at low water level(m)
    nc = expobj.ChannelShapeParam;  %channel shape parameter (-)
    Lwu = expobj.HWwidthELength;    %width convergence length at high water (m)
    Lwl = expobj.LWwidthELength;    %width convergence length at low water (m)
    zm = obj.zMouthInvert;          %thalweg bed level at mouth to zero datum (m)
    Lt = hydobj.xTidalLimit;        %distance from mouth to tidal limit
    
    %aspect ratio and slope coefficient (mc) at mouth    
    dl = hydobj.zlw(1)-zm;          %depth of lower form at mouth to lw (m)
    ar = Wlw/dl;                    %aspect ratio of low water channel at mouth
    mc = 2*nc/ar;                   %submerged static coefficient of Coulomb 
                                    %friction (estimated from geometry)   
    %derive width, csa and convergence length at mean tide                                
    Wm = (Whw+Wlw)/2;               %estimate of mean tide width
    Lw = (Lwu+Lwl)/2;               %estimate of mean convergence length
    tr = (hydobj.zhw(1)-hydobj.zlw(1)); %tidal range at mouth 
    Am = Wm*tr+mc*Wlw^2/6;          %csa at the mouth (m^2)
    [~,~,Arv] = get_river_regime(obj,tr);    
    La = -Lt/log((Arv+1)/Am);       %assume mouth converges to river & min csa = 1
    %assign values
    obj.CSTparams.Wm = Wm;
    obj.CSTparams.Lw = Lw;
    obj.CSTparams.Am = Am;
    obj.CSTparams.La = La;  
end
%%
function [xi,yi,zi,Wz,Rv] = channel_3D_form(obj)
    %generate the 3D form
    xi = []; yi = []; zi = []; Wz = []; Rv = [];
    %get the required input parameter classes
    expobj = obj.RunParam.CF_FormData;
    grdobj = obj.RunParam.GD_GridProps;
    hydobj = obj.RunParam.CF_HydroData;
    sedobj = obj.RunParam.CF_SediData;
    
    %model run parameters
    Lm = diff(grdobj.XaxisLimits);   %length of model domain (m)
    offset = 2*grdobj.histint;       %offset from hw to supra-tidal form
    
    %channel form parameters    
    bu = expobj.HWmouthWidth/2;      %half-width of mouth at high water(m)
    bl = expobj.LWmouthWidth/2;      %half-width of mouth at low water level(m)
    nc = expobj.ChannelShapeParam;   %channel shape parameter (-)
    Lwu = expobj.HWwidthELength;     %width convergence length at high water (m)
    Lwl = expobj.LWwidthELength;     %width convergence length at low water (m)
    nu = expobj.HWwidthPower;        %width exponent at high water (-)
    nl = expobj.LWwidthPower;        %width exponent at low water (-)   
    zm = obj.zMouthInvert;           %thalweg bed level at mouth to zero datum (m)
    ki = expobj.FlatShapeParam;      %intertidal shape parameter[ki*100; range:0.01-0.5]
    Ll = hydobj.xTideRiver;          %distance from mouth to estuary/river switch
    Lt = hydobj.xTidalLimit;         %distance from mouth to tidal limit
    
    if strcmp(obj.Selection.planform,'Power')
        %negative exponents result in invalid complex forms
        if  any([nu,nl]<=0)
            warndlg('One or more exponents for power form not set or invalid')
            return; 
        end           
    elseif strcmp(obj.Selection.planform,'Exponential')
        %allow divergent landwards exponential at HW but not at LW
        if Lwl<=0
            warndlg('Low water exponents for exponential form not set or invalid')
            return; 
        end 
    end
    
    %sediment properties (if saltmarsh defined)
    dmax = 0;
    if ~isempty(sedobj.MaxMarshDepth)
        dm = sedobj.AvMarshDepth;    %average depth of marsh surface
        dmax = sedobj.MaxMarshDepth; %maximum depth of saltmarsh
    end

    fact = 2;
    switch obj.Selection.intertidalform
        case 'Uniform Shear'
            fact = 2.57;             %based on L=L*(pi/2+1) where pi/2+1=2.57
    end

    %set-up co-ordinate system
    [~,yi,delx] = getGridDimensions(grdobj);
    xi = -(Lm-Ll):delx:Ll;  %from head (-ve) to mouth (+ve)
    yi = yi(yi>=0);  %half the grid
  
    zi = zeros(length(xi),length(yi));
    Wz = zeros(length(xi),3);
    %water level properties based on amplitude+mtl or CST model (mAD)
    zHWxi = fliplr(hydobj.zhw);      %high water level(mAD)
    zLWxi = fliplr(hydobj.zlw);      %low water level(mAD)    
    amp0 = (hydobj.zhw(1)-hydobj.zlw(1))/2; %tidal amplitude at mouth 
    
    %river properties
    [Rv.Hr,Rv.Wr,Rv.Ar] = get_river_regime(obj,2*amp0);
    [hrv,br,mcr] = get_river_profile(obj,2*amp0,yi);   %#ok<ASGLU>

    %aspect ratio and slope coefficient (mc) at mouth    
    dl = hydobj.zlw(1)-zm;           %depth of lower form at mouth to lw (m)
    ar = 2*bl/dl;                    %aspect ratio of low water channel at mouth
    mct = 2*nc/ar;                   %submerged static coefficient of Coulomb 
                                     %friction (estimated from geometry)
    
    %calculate the transformation of the x co-ordinate for given x and y values
    %and then use this to obtain a vaule of z    
    for ix=1:length(xi)
        zhw = zHWxi(ix); 
        zlw = zLWxi(ix);
        ax = (zhw-zlw)/2;            %tidal amplitude(m)
        zo = zhw-ax;                 %mean tide level(m)  
        
        %------------------------------------------------------------------
        %Note - see how this relates to work of Uncles re along channel
        %slope variation. *************************************************
        % mc = mcr+(mct-mcr)*ax/amp0;%interpolate river and mouth slope 
        %                            %coefficients as function of tidal amplitude
        % mc = mct*amp0/ax;          %scale value at mouth based on tidal amplitude
        % if mc>mcr, mc = mcr; end   %river is limiting value
        %**************
        mc = mct;                    %force constant value
        %------------------------------------------------------------------
        switch obj.Selection.planform
            case 'Exponential'
                xexp = Ll-xi(ix);    %x defined from mouth for exponential          
                [yu,yo,yl] = expPlan(xexp,bu,bl,br,Lt,Lwu,Lwl,fact);
            case 'Power'
                [yu,yo,yl] = powPlan(xi(ix),Lt,Ll,bu,bl,br,nu,nl,fact);
        end
        ls = (yu-yl)/fact;           %lower intertidal width from lw to mtl 
        
        %supra-tidal form relative to mtl        
        if yu==0
            zds(1) = zhw-zo+offset;
        else
            zds = (zhw-zo+offset).*(yi>yu); %elevation of surrounding land
        end

        %define intertidal form
        switch obj.Selection.intertidalform            
            case 'Linear'
                %upper intertidal form
                zdu = (ax.*((yi-yl)./ls-1)).*(yi<=yu & yi>yo); 
                %lower intertidal form
                zdl = (ax.*((yi-yl)./ls-1)).*(yi<=yo & yi>=yl); 
            case 'Rectangular'
                %upper intertidal form
                zdu = 0.*(yi<=yu & yi>yo);                
                %lower intertidal form
                zdl = 0.*(yi<=yo & yi>=yl);
            case 'Stepped'
                %upper intertidal form
                zdu = ax/2.*(yi<=yu & yi>yo);                
                %lower intertidal form
                zdl = -ax/2.*(yi<=yo & yi>=yl);
            case 'Uniform Shear'
                if ls==0             %no intertidal present
                    zdu = yi*0;
                    zdl =yi*0;
                else
                    %upper intertidal form
                    zdu = (ax.*sin((yi-yl)./ls-1)).*(yi<=yu & yi>yo);
                    %lower intertidal form
                    zdl = (ax.*((yi-yl)./ls-1)).*(yi<=yo & yi>=yl);                    
                end
            case 'L&M muddy shore'
                if yu==yl             %no intertidal present
                    zdu = yi*0;
                    zdl =yi*0;
                else
                    yol = (yi-yl)/(yu-yl);
                    yol = yol.*(yol<=1);
                    %upper intertidal form
                    zdu = ax*(1 - 2*exp(4*ki.*yol).*(1-yol).^2).*(yi<=yu & yi>yo);                
                    %lower intertidal form
                    zdl = ax*(1- 2*exp(4*ki.*yol).*(1-yol).^2).*(yi<=yo & yi>=yl);
                end
        end

        %check for saltmarsh on upper intertidal form
        if dmax>0            
            idm = yi<=yu & yi>yo;    %width in which marsh can exist
            ismarsh = logical(zdu>(ax-dmax).*idm); %valid marsh depths
            zdu(ismarsh) = (ax-dm);  %assign mean marsh depth
        end
        
        %define low water channel form
        switch obj.Selection.channelform
            case 'Parabolic'
                if yl>0
                    zdco = mc*yl/nc.*(1-(yi/yl).^nc).*(yi<yl); %depth at yi
                    zdc = -ax.*(zdco>0) -zdco;
                    isriver = yi<br & hrv>-zdc;
                    if any(isriver)
                        zdcr = -hrv.*isriver;
                        zdc = zdc.*(~isriver)+zdcr;
                    end
                else
                    zdc = yi*0;
                end
            case 'Rectangular'
                if yl>0
                    %hd = mc.Wlw/2/(nc+1), or %hd = nc.hc/(nc+1),
                    hydep = mc*yl/(nc+1);        
                    zdco = ones(size(yi))*hydep.*(yi<yl); %depth at yl
                    zdc = -ax.*(zdco>0) -zdco;
                    isriver = yi<br & hrv>-zdc;
                    if any(isriver)
                        zdcr = -hrv.*isriver;
                        zdc = zdc.*(~isriver)+zdcr;
                    end
                else
                    zdc = yi*0;
                end
        end       
    
        %combine forms
        zix = zds+zdu+zdl+zdc;
        %adjust mtl datum to zero datum  
        zi(ix,:)  = zix + zo; 
        Wz(ix,:) = [yu, yo, yl]*2;   %full width       
    end
end
%%
function [yu,yo,yl] = expPlan(xexp,bu,bl,br,Lt,Lhw,Llw,fact)
    %compute controlling dimensions (y) for an exponential plan form
    Lmt = (Lhw+Llw)/2;               %width shape factor as average of upper and lower values
    Ls = (bu-bl)/fact;               %Lstar in F&A, lower intertial width (lw to mtl)(m) 
    bo = bl+Ls;                      %half-width at mtl(m)
    yu = (bu-br)*exp(-xexp/Lhw)+br;  %distance from centre line to hw
    yo = (bo-br)*exp(-xexp/Lmt)+br;  %distance from centre line to mtl
    yl = (bl-br)*exp(-xexp/Llw)+br;  %distance from centre line to lw 
    if xexp>Lt && br>0       
        if yu<br, yu = br; end        
        if yo<br, yo = br; end 
        if yl<br, yl = br; end 
    elseif xexp>Lt
        yu = 0; yo = 0; yl = 0;    
    end
end
%%
function [yu,yo,yl] = powPlan(xxi,Lt,Ll,bu,bl,br,mu,ml,fact)
    %compute controlling dimensions (y) for a power plan form
    mo = (mu+ml)/2; %width shape factor as average of upper and lower values
    %Lt = Ll+Lu;                     %total length of estuary (m)
    Lu = Lt-Ll;
    Lo = Ll+Lu/fact;                 %length of channel at mtl(m)
    Ls = (bu-bl)/fact;               %Lstar in F&A, lower intertial width (lw to mtl)(m) 
    bo = bl+Ls;                      %half-width at mtl(m)
    xu = xxi+Lu;                     %adjust to origin of high water plan form
    %use eps as an offset in x to avoid inf if x=0
    yu = ((bu-br)*((xu)/Lt)^mu)*(xu>0)+br;   %distance from centre line to hw
    xo = xxi+Lu/fact;                %adjust to origin of mtl plan form
    yo = ((bo-br)*((xo)/Lo)^mo)*(xo>0)+br;   %distance from centre line to mtl
    yl = ((bl-br)*((xxi)/Ll)^ml)*(xxi>0)+br; %distance from centre line to lw
    
    
    if xxi<Ll-Lt && br>0        
        if yu<br, yu = br; end        
        if yo<br, yo = br; end 
        if yl<br, yl = br; end 
    elseif xxi<Ll-Lt
        yu = 0; yo = 0; yl = 0;          
    end
end
