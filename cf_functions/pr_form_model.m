function [xi,yi,zgrd,yz] = pr_form_model(~,mobj,isfull)
%
%-------function help------------------------------------------------------
% NAME
%   pr_form_model.m
% PURPOSE
%   function to compute 3D form of a creek or tidal channel
%   useing power laws to define width and hydraulic depth variations.
% USAGE
%   [xi,yi,zgrd,yz] = pr_form_model(obj,inputs,isfull)
% INPUTS
%   obj - CF_FormModel class instance (used to identify class method)
%   mobj - ChannelForm model UI instance
%   isfull - true returns full grid, false half-grid
% OUTPUT
%   xi - x co-ordinate (m)
%   yi - y co-ordinate (m)
%   zgrd - bed elevation grid (m)
%   yz - width at hw,mt,lw (m)
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
    if nargin<3
        isfull = true;
    end
    %get the required input parameter classes
    pwrobj = getClassObj(mobj,'Inputs','CF_PowerData');
    grdobj = getClassObj(mobj,'Inputs','CF_GridData');
    hydobj = getClassObj(mobj,'Inputs','CF_HydroData');
    sedobj = getClassObj(mobj,'Inputs','CF_SediData');

    %model run parameters
    Lt = diff(grdobj.XaxisLimits);   %length of model domain (m)
    nintx = grdobj.Xint;             %no of intervals in the x direction
    bt = diff(grdobj.YaxisLimits)/2; %half width of model domain (m)   
    ninty = grdobj.Yint/2;           %no of intervals in the y direction
    
    %channel form parameters
    Le = pwrobj.TotalLength;       %total length of channel (m)
    bu = pwrobj.HWmouthWidth/2;    %half-width of mouth at high water(m)
    bl = pwrobj.LWmouthWidth/2;    %half-width of mouth at low water level(m)
    nu = pwrobj.HWwidthExponent;   %width exponent at high water (-)
    nl = pwrobj.LWwidthExponent;   %width exponent at low water (-)
    mu = pwrobj.HWdepthExponent;   %depth exponent at high water (-)
    ml = pwrobj.LWdepthExponent;   %depth exponent at low water (-)
    zm = pwrobj.zMouthInvert;      %thalweg bed level at mouth to zero datum (m)    
    
    %sediment properties (if saltmarsh defined)
    dmax = 0;
    if ~isempty(sedobj.MaxMarshDepth)
        dm = sedobj.AvMarshDepth;    %average depth of marsh surface
        dmax = sedobj.MaxMarshDepth; %maximum depth of saltmarsh
    end
    
    Ll = hydobj.xTideRiver(1);  %distance from mouth to estuary/river switch
    Lu = Le-Ll;                 %distance from estuary/river switch to head
    
    %set-up co-ordinate system
    delx = Lt/nintx;
    dely = bt/ninty;
    Lu = Lu+delx;   
    xi = -(Lt-Ll):delx:Ll;
    yi = 0:dely:bt; yi(1) = 0.1;     %the offset ensures no duplicates
    
    z20 = zeros(length(xi),length(yi));
    yz = zeros(length(xi),3);
    %water level properties based on amplitude+mtl or CST model (mAD)
    if isscalar(hydobj.zhw)
        zHWxi = ones(length(xi),1)*hydobj.zhw;
        zLWxi = ones(length(xi),1)*hydobj.zlw;
    else
        %water level properties based on amplitude+mtl or CST model (mAD)
        zHWxi = hydobj.zhw;              %high water level(mAD)
        zLWxi = hydobj.zlw;              %low water level(mAD)
    end
    
    %river properties
    [hrv,bh,~] = get_river_profile(yi,mobj); 
    
    %calculate the transformation of the x co-ordinate for given x and y values
    %and then use this to obtain a value of z
    nxi = length(xi);
    for ix=1:nxi
        zhw = zHWxi(ix); zlw = zLWxi(ix);
        am = (zhw-zlw)/2;     %tidal amplitude(m)
        zo = zhw-am;          %mean tide level(m)
        du = zhw-zo;          %depth of upper form
        dl = zm-zo;           %depth of lower form
        %upper form
        zu = du.*(Lt/Lu.*((yi./bu).^(1/nu)-xi(ix)./Lt)).^mu;
        zu = zu.*(zu>=0);
        yu = (bu.*((Lu+xi(ix))./Lt).^nu).*(Lu+xi(ix)>0);
        %lower form
        temp = (dl.*(xi(ix)./Ll-(yi./bl).^(1/nl)).^ml);
        zl = temp.*(temp<=0);
        yl = (bl.*(xi(ix)./Ll).^nl).*(xi(ix)>0);
        yo = (yu+yl)/2;
        %check for saltmarsh on upper intertidal form        
        if dmax>0     
            idm = yi<=yu & yi>yo;    %width in which marsh can exist
            ismarsh = logical(zu>(am-dmax).*idm); %valid marsh depths
            zu(ismarsh) = (am-dm);   %assign mean marsh depth
        end
        zu = zu.*(yi>=yl);
        zl = zl.*(yi<yl);
        z20(ix,:) = zu+zl;           %elevations relative to zero datum
        yz(ix,:) = [yu, yo, yl];
        isriver = yi<bh & hrv>-z20(ix,:);
        if any(isriver)
            zdcr = -hrv.*isriver;
            z20(ix,:) = z20(ix,:).*(~isriver)+zdcr;
        end
    end
    
    yz = num2cell(yz',2);  %formatted to load into dstable 
    %combine upper and lower surface
    zi  = z20+zo;           %add vertical transformation to move origin  
    zhw = repmat(zHWxi,1,ninty+1);
    msk = (zhw+0.01).*(zi>zhw);      %construct high water mask
    zi  = zi.*(zi<=zhw)+ msk;

    %generate complete 3D channel form by mirroring half section
    if isfull                        %return full grid
        zgrd = cat(2,fliplr(zi),zi);
        yi  = [-fliplr(yi), yi];
    else                             %return half grid
        zgrd = zi;
    end
end