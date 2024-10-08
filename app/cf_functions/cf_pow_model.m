function [xi,yi,zgrd,Wz,Rv] = cf_pow_model(obj,isfull)
%
%-------function help------------------------------------------------------
% NAME
%  cf_pow_model.m
% PURPOSE
%   function to compute 3D form of a creek or tidal channel
%   useing power laws to define width and hydraulic depth variations.
% USAGE
%   [xi,yi,zgrd,yz] = cf_pow_models(obj,isfull)
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
%   Wz -  table of Whw,Wmt,Wlw for width at hw,mt,lw (m)
%   Rv - struct of river regime properties Hr, Wr, Ar
% NOTES
%   Power law form is an adaptation to 3D of the form explored by Prandle &
%   Rahman, 1980, JPO.
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

    %channel properties
    if obj.Selection.wlflag==0
        %provides initial guess of gross properties if cst_model called
        obj.Selection.wlflag = 1;
        obj= cf_set_hydroprops(obj);     %fixed water level surface 
        obj.Selection.wlflag = 0;
        obj = pr_properties(obj); 
    end
    %set the water level variations along the estuary
    [obj,ok] = cf_set_hydroprops(obj);
    if ok<1, return; end
    
    [xi,yi,zi,Wz,Rv] = pr_3D_form(obj);
    if isempty(xi),return; end
    
    %model x-axis is from head. Reverse data for use in ChannelForm
    yzcell = num2cell(flipud(Wz)',2);  %formatted to load into dstable
    Wz = table(yzcell{:},'VariableNames',{'Whw','Wmt','Wlw'});
    %x is defined from head with origin at "shoulder". 
    %Change to origin at mouth
    xi = fliplr(max(xi)-xi);
    %generate complete 3D channel form by mirroring half section
    if isfull                        %return full grid
        zgrd = flipud(cat(2,fliplr(zi(:,2:end)),zi));
        yi  = [-flipud(yi(2:end)); yi];
    else                             %return half grid
        zgrd = zi;
    end
     
%     [X,Y] = meshgrid(xi,yi);
%     figure; surf(X,Y,zgrd);   
end
%%
function obj = pr_properties(obj)
    %compute the summary gross properties for the channel form
    [xi,yi,zi] = pr_3D_form(obj); 
    grid.x = fliplr(max(xi)-xi);
    grid.y  = [-flipud(yi(2:end)); yi];
    grid.z = flipud(cat(2,fliplr(zi(:,2:end)),zi));
    grid.t = years(0);
    grid.ishead = false; 
    grid.xM = 0;
    grid.Lt = max(xi);
    
    wl = obj.RunParam.CF_HydroData; 
    grdobj = obj.RunParam.GD_GridProps;
    [~,hypdst] = gd_basin_hypsometry(grid,wl,grdobj.histint,0);
    props = gd_section_properties(grid,wl,hypdst);
    gp = gd_gross_properties(grid,wl,props);
    %used in CSTmodel
    obj.CSTparams.Wm = gp.Wm;
    obj.CSTparams.Lw = gp.Lw;
    obj.CSTparams.Am = gp.Am;
    obj.CSTparams.La = gp.La;  
end
%%
function [xi,yi,zi,Wz,Rv] = pr_3D_form(obj)    
    %generate the 3D form
    xi = []; yi = []; zi = []; Wz = []; Rv = [];
    %get the required input parameter classes
    pwrobj = obj.RunParam.CF_FormData;
    grdobj = obj.RunParam.GD_GridProps;
    hydobj = obj.RunParam.CF_HydroData;
    sedobj = obj.RunParam.CF_SediData;

    %model run parameters
    Lt = diff(grdobj.XaxisLimits);   %length of model domain (m)
    offset = 2*grdobj.histint;       %offset from hw to supra-tidal form
    
    %channel form parameters
    bu = pwrobj.HWmouthWidth/2;    %half-width of mouth at high water(m)
    bl = pwrobj.LWmouthWidth/2;    %half-width of mouth at low water level(m)
    nu = pwrobj.HWwidthExponent;   %width exponent at high water (-)
    nl = pwrobj.LWwidthExponent;   %width exponent at low water (-)
    mu = pwrobj.HWdepthExponent;   %depth exponent at high water (-)
    ml = pwrobj.LWdepthExponent;   %depth exponent at low water (-)
    zm = obj.zMouthInvert;         %thalweg bed level at mouth to zero datum (m)    

    if any([nu,nl,mu,ml]<=0)
        %negative exponents result in invalid complex forms
        warndlg('One or more exponents for power form not set or invalid')
        return; 
    end           

    %sediment properties (if saltmarsh defined)
    dmax = 0;
    if ~isempty(sedobj.MaxMarshDepth)
        dm = sedobj.AvMarshDepth;    %average depth of marsh surface (m)
        dmax = sedobj.MaxMarshDepth; %maximum depth of saltmarsh (m)
    end
    
    Le = hydobj.xTidalLimit;       %length of tidal channel (m)
    Ll = hydobj.xTideRiver;        %distance from mouth to estuary/river switch (m)
    Lu = Le-Ll;                    %distance from estuary/river switch to head (m)
    
    %set-up co-ordinate system
    [~,yi,delx] = getGridDimensions(grdobj);
    xi = -(Lt-Ll):delx:Ll;
    yi = yi(yi>=0);    %half the grid
    nyi = length(yi);
    
    zi = zeros(length(xi),length(yi));
    Wz = zeros(length(xi),3);
    %water level properties based on amplitude+mtl or CST model (mAD)
    zHWxi = fliplr(hydobj.zhw);             %high water level(mAD)
    zLWxi = fliplr(hydobj.zlw);             %low water level(mAD)
    amp0 = (hydobj.zhw(1)-hydobj.zlw(1))/2; %tidal amplitude at mouth
    zso = hydobj.zmt(1);  %mean tide level at mouth used to define the level of the 'shoulder' in the power form (constant)
    du = hydobj.zhw(end)-zso;               %depth of upper form at head relative to zso (constant)
    dl = zm-zso;          %depth of lower form at mouth relative to zso (constant)

    %river properties
    [Rv.Hr,Rv.Wr,Rv.Ar] = get_river_regime(obj,2*amp0);
    [hrv,bh,~] = get_river_profile(obj,2*amp0,yi'); 
    
    %calculate the transformation of the x co-ordinate for given x and y values
    %and then use this to obtain a value of z
    nxi = length(xi);
    for ix=1:nxi
        zhw = zHWxi(ix); zlw = zLWxi(ix);
        am = (zhw-zlw)/2;     %tidal amplitude(m)
        zmt = zhw-am;         %mean tide level(m)
        
        % upper form
        %zu equation is for a coordinate system with x=0 at the 'shoulder'
        %where z=0 before the vertical offset for the 'shoulder', zso, 
        %is applied
        xLi = Lu+xi(ix);
        temp = Le/Lu.*((yi./bu).^(1/nu)-xi(ix)./Le);
        temp(temp<0) = 0;
        zu = du.*((temp.*(xLi>=0)).^mu); %equation only valid to tidal limit
        %elevation when above tidal limit set to be same as surrounding land
        if xLi<0, zu = ones(size(yi))*zhw+offset; end  
        yu = bu*(xLi/Le)^nu*(xLi>0);    

        % lower form
        %zl equation applies below and seaward of the 'shoulder' ie for x>0
        %and z<zso in the same coordinate system (see definitions sketch in
        %manual)
        temp = (xi(ix)./Ll-(yi./bl).^(1/nl));
        temp(temp<0) = 0;
        %apply along channel and cross-channel limits to validity of eqn
        limit = (xi(ix)>0) & (temp>=0);   
        zl = dl.*((temp.*limit).^ml); 
        yl = (bl*(xi(ix)/Ll)^nl)*(xi(ix)>0);
        yo = (yu+yl)/2;
        
        %check for saltmarsh on upper intertidal form        
        if dmax>0     
            idm = yi<=yu & yi>yo;    %width in which marsh can exist
            ismarsh = logical(zu>(am-dmax).*idm); %valid marsh depths
            zu(ismarsh) = (am-dm);   %assign mean marsh depth
        end
        zu = zu.*(yi>=yl);
        zl = zl.*(yi<yl);
        %elevations adjusted for variation in mtl. The model is strictly
        %adjusted by zso to define the domains relative to the 'shoulder'.
        %zmt is used here to generate an along channel vertical variation 
        %in the horizontal plane that defines the shoulder, so that system 
        %takes account of valley slope and link to the river system     
        zi(ix,:) = zu+zl+zmt;           
        isriver = yi'<bh & (du-hrv)<zi(ix,:);
        if any(isriver)
            if yu<bh, yu=bh; end
            if yl<bh, yl=bh; end
            zdcr = (zhw-hrv).*isriver;
            zi(ix,:) = zi(ix,:).*(~isriver)+zdcr;
        end
        Wz(ix,:) = [yu, yo, yl]*2;   %full width
    end
    
    %add mask to define surrounding land surface
    zhw = repmat(zHWxi',1,nyi);
    msk = (zhw+offset).*(zi>zhw);      %construct high water mask
    zi  = zi.*(zi<=zhw)+ msk;
end