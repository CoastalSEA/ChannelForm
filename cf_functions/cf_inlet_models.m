function [xi,yi,zgrd,Wz,Rv] = cf_inlet_models(obj,isfull)
%
%-------function help------------------------------------------------------
% NAME
%   cf_inlet_models.m
% PURPOSE
%   construct idealised tidal inlet form using exponential function in y to 
%   determine LW channel width and cross-section of various forms to determine z at each x interval
% USAGE
%   [xi,yi,zgrd,yz,Rv] = cf_inlet_models(obj,isfull)
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
%   rectangular prism with an exponential plan for and an intertidal using 
%   linear, stepped, uniform shear or muddy shore profiles to idealised
%   basin plan forms - square, rectangle, circle, elipse
%   Only uses constant HW and constant or tapering????? LW levels
% SEE ALSO
%   used in CF_FormModel as part of ChannelForm model
%
% Author: Ian Townend
% CoastalSEA (c) June 2024
%--------------------------------------------------------------------------
%
    xi = []; yi = []; zgrd = []; Wz = [];
    if nargin<3
        isfull = true;
    end

    if obj.Selection.wlflag==0 
        fprintf('CSTmodel not currently supported for tidal inlets\n')
        return; 
    else
        %set the water level variations along the estuary
        [obj,ok] = cf_set_hydroprops(obj); 
        if ok<1
            fprintf('No water level surface found. Grid not defined\n')
            return; 
        end
    end

    [xi,yi,zi,Wz,Rv] = inlet_3D_form(obj);
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
%     figure; surf(Y,X,zgrd);
end
%%
function [xi,yi,zi,Wz,Rv] = inlet_3D_form(obj)
    %generate the 3D form
    
    %get the required input parameter classes
    inletobj = obj.RunParam.CF_FormData;
    grdobj = obj.RunParam.GD_GridProps;
    hydobj = obj.RunParam.CF_HydroData;
    sedobj = obj.RunParam.CF_SediData;

    %model run parameters
    Lm = diff(grdobj.XaxisLimits);   %length of model domain (m)
    offset = 2*grdobj.histint;       %offset from hw to supra-tidal form

    %channel form parameters     
    bu = inletobj.HWmouthWidth/2;      %half-width of mouth at high water(m)
    bl = inletobj.LWmouthWidth/2;      %half-width of mouth at low water level(m)
    bt = inletobj.HWbasinWidth/2;      %half-width of basin (m)
    Lwl = inletobj.LWwidthELength;     %width convergence length at low water (m)
    nl = inletobj.LWwidthPower;        %width exponent at low water (-)    
    Li = inletobj.xInletLength;        %length of mouth inlet (m)
    Ll = inletobj.xLWchannel;          %distance from mouth to end of LW channel (m)
    Lt = inletobj.xHWbasinLength;      %distance from mouth to tidal limit (m)
    nc = inletobj.ChannelShapeParam;   %channel shape parameter (-)
    ki = inletobj.FlatShapeParam;      %intertidal shape parameter[ki*100; range:0.01-0.5]
    zm = obj.zMouthInvert;             %thalweg bed level at mouth to zero datum (m)  
    
    if strcmp(obj.Selection.planform,'Power')
        %negative exponents result in invalid complex forms
        if nl<=0
            warndlg('LW exponent for power form not set or invalid')
            return; 
        end           
    elseif strcmp(obj.Selection.planform,'Exponential')
        %divergent landwards exponential at LW not possible
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

    %scaling for lower intertidal as proportion of full width
    fact = 2;
    switch obj.Selection.intertidalform
        case 'Uniform Shear'
            fact = 2.57;             %based on L=L*(pi/2+1) where pi/2+1=2.57
    end

    %set-up co-ordinate system
    [~,yi,delx,~] = getGridDimensions(grdobj);
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
    [hrv,br,mcr] = get_river_profile(obj,2*amp0,yi);  %#ok<ASGLU> %NB bh-half width

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
        axo = (zhw-zlw)/2;            %tidal amplitude(m)
        
        if xi(ix)<0 && strcmp(obj.Selection.lwform,'Yes') 
            %adjustment to get tidal flat to slope landwards as well as
            %laterally from limit of low water channel
            switch obj.Selection.intertidalform            
                case 'Linear'
                    zlw = zlw+(-xi(ix))*(zhw-zlw)/(Lt-Ll);
                    if zlw>zhw, zlw = zhw; end
                case 'Uniform Shear'
                    lsx = (Lt-Ll)/fact;
                    if xi(ix)>-lsx
                        zlw = axo.*(xi(ix)./lsx-1);    %lower intertidal form
                    else
                        zlw = axo.*sin(xi(ix)./lsx-1); %upper intertidal form
                    end
                case 'L&M muddy shore'
                    xol = xi(ix)/(Lt-Ll);
                    xol = xol*(xol<1);
                    zlw = axo*(1-2*exp(4*ki.*xol).*(1-xol).^2);
            end                    
        end
        
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
        basin = obj.Selection.basinform;
        switch obj.Selection.planform
            case 'Exponential'
                xexp = Ll-xi(ix);    %x defined from mouth for exponential          
                [yu,yo,yl] = expPlan(xexp,bu,bl,br,bt,Lt,Ll,Li,Lwl,fact,basin);
            case 'Power'
                [yu,yo,yl] = powPlan(xi(ix),bu,bl,br,bt,Lt,Ll,Li,nl,fact,basin);
        end
        
        if yu-yl==0
            ls = 1;                  %prevent divide by 0 when no flat
        else
            ls = (yu-yl)/fact;       %lower intertidal width from lw to mtl 
        end
        
        %supra-tidal form relative to mtl        
%         if yu==0
%             zds(1) = zhw-zo+offset;
%         else
%             zds = (zhw-zo+offset).*(yi>yu); %elevation of surrounding land
%         end
        zds = (zhw-zo+offset).*(yi>yu); %elevation of surrounding land
        if yu==0
            zds(1) = zhw-zo+offset;     %add in central value when no channel
        end
%         if xi(ix)<-(Lt-Ll)
%             %works for rectangle but not for ellipse!!!!****************
%             zds = zhw-zo+offset.*(yi>=yu);         %end of basin 
%         end

        %define intertidal form
        zdu = 0.*yi; zdl = 0.*yi;       %needed if yo=0 or yl=0
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
        
        if any(zdu>zhw-zo+offset)
            idu = zdu>zhw-zo+offset;
            zdu(idu) = 0;
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
function [yu,yo,yl] = expPlan(xexp,bu,bl,br,bt,Lt,Ll,Li,Llw,fact,basin)
    %compute controlling dimensions (y) for an exponential plan form at LW
    %xexp - distance from mouth; bu - mouth HW half-width; bl - mouth LW
    %half-width; bh - river half-width; bt - basin HW half-width;
    %Lt - tidal length from mouth; Ll - length of low water channel from mouth
    %Li - length of inlet; Lle - LW width convergence
    Ls = (bu-bl)/fact;               %Lstar in F&A, lower intertial width (lw to mtl)(m) 
    bo = bl+Ls;                      %half-width at mtl(m)    
    if xexp<=Li
        yl = bl;  yo = bo;  yu = bu;
    else
        yl = (bl-br)*exp(-(xexp-Li)/Llw)+br;  %distance from centre line to lw 
        if br==0 && xexp>Ll && yl<1
            yl = 0;                      %no river limit channel to low water length
        end
        %
        if strcmp(basin,'Rectangle')
            yu = bt;                     %distance from centre line to hw
        else
            xi = xexp-Li;
            yu = getHWForm(yl,xi,bu,bt,Lt,Ll,Li,basin);
        end
        %
        if xexp>Lt
            yu = yl;
        end
        yo = yl+(yu-yl)/fact;            %distance from centre line to mtl       
    end
    %
    if xexp>Lt && br>0   
        if yu<br, yu = br; end        
        if yo<br, yo = br; end 
        if yl<br, yl = br; end 
    elseif xexp>Lt
        yu = 0; yo = 0; yl = 0;
    end
end
%%
function [yu,yo,yl] = powPlan(xxi,bu,bt,Lt,Ll,Li,ml,basin)
    %compute controlling dimensions (y) for a power plan form
    %xxi - distance
    Ls = (bu-bl)/fact;               %Lstar in F&A, lower intertial width (lw to mtl)(m) 
    bo = bl+Ls;                      %half-width at mtl(m)    
    if xxi>Ll-Li                     %inlet reach 
        yl = bl;  yo = bo;  yu = bu;
    elseif xxi<Ll-Lt && br>0         %landward reach with river
        yl = br;  yu = br;  yo = br; 
    elseif xxi<Ll-Lt
        yu = 0;   yo = 0;    yl = 0; %landward reach no river  
    else                             %central channel reach
        yl = (((bl-br)*(xxi/(Ll-Li))^ml))*(xxi>0)+br; %distance from centre line to lw
        %
        if strcmp(basin,'Rectangle')
            yu = bt;                         %distance from centre line to hw
        else
            xi = Ll-xxi-Li;
            yu = getHWForm(yl,xi,bu,bl,br,bt,Lt,Ll,Li,fact,basin);
        end
        yo = yl+(yu-yl)/fact;            %distance from centre line to mtl           
    end
end
%%
function yu = getHWForm(yl,xi,bu,bt,Lt,Ll,Li,basin)
    %define high water width based in selected plan form
    %xi -distance upstream from landward end of inlet
    switch basin
        case 'Ellipse'
            %ellipse around tidal channel 
            if xi>=0 && 1-xi/(Lt-Li)>0                  
                yu = sqrt((4*bt^2)/(Lt-Li)*xi*(1-xi/(Lt-Li)));
            else 
                yu = yl;
            end           
        case 'Half-ellipse'
            %half-ellipse around tidal channel    
            if xi>=0 && (1-xi^2/(Lt-Li)^2)>0               
                yu = sqrt((bt^2)*(1-xi^2/(Lt-Li)^2));
            else 
                yu = yl;
             end    
        case 'Divergent-shore'
            %exponential divergence
            Lwu = (Ll-Li)/log(bt/bu);
            yu = bu*exp(xi/Lwu);
            if yu>bt, yu = bt; end
        case 'Divergent-bay'
            %exponential divergence
            Lwu = (Ll-Li)/log(bt/bu);
            yuexp = bu*exp(xi/Lwu);
            if yuexp>bt, yuexp = bt; end
            %parabolic head
            yupbl = sqrt(4*(Lt-Li-xi)*(Lt-Ll));
            %combined plan form
            yu = min(yuexp,yupbl);
        case 'Logistic-shore'
            %logistic variation around tidal limit
            x = xi/(Ll-Li);
            params = setGLstruct(bu,bt,Ll,Li);
            yu = general_logistic(x,params,false)*(bt-bu)+bu;
        case 'Logistic-bay'
            %logistic variation around tidal limit
            x = xi/(Ll-Li);
            params = setGLstruct(bu,bt,Ll,Li);
            yuexp = general_logistic(x,params,false)*(bt-bu)+bu;
            %parabolic head
            yupbl = sqrt(4*(Lt-Li-xi)*(Lt-Ll));
            %combined plan form
            yu = min(yuexp,yupbl);
    end
end
%%
function params = setGLstruct(bu,bt,Ll,Li)
    %   p - struct containing:
    A = 0.0;      %left horizontal asymptote
    B = 10;       %growth rate
    C = 1.0;      %upper asymptote of A+(K-A)/C^(1/nu) - typically 1
    K = 1.0;      %right horizontal asymptote when C=1
    nu = 1.0;     %>0 - affects proximity to which asymptote maximum growth occurs
    M = 0.5;      %start point (controls length of channel before expansion)
    Q = 1.0;      %related to the value y(0).

    params = struct('A',A,'B',B,'C',C,'K',K,'nu',nu,'M',M,'Q',Q);
end