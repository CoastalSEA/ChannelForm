function formdims = ckfa_form_solver(inp,hme)
%
%-------function help------------------------------------------------------
% NAME
%   ckfa_form_solver.m 
% PURPOSE
%   find the hydraulic depth, e-folding length for width and the peak tidal
%   amplitude for defined tidal range and sediment properties
% USAGE
%   formdims = ckfa_form_solver(hme,params)
% INPUTS
%   hme = intial guess of hydraulic depth at mouth (m)
%   inp is a struct with fields (*)=not used in this function:
%       amp    = tidal amplitude (m)
%       omega =  angular frequency, 2pi/Tp (1/s)
%       Le   = estuary length (m)
%       Wrv  = width of regime river channel (m)
%       Arv  = CSA of regime river channel (m2)
%       Uw  = wind speed (m/s)
%       zw  = elevation of wind speed (m) - default is 10m
%       rhow = density of water (kg/m^3)
%     	rhoc = suspended sediment concentration (kg/m^3)
%       taucr= critical threshold bed shear stress (Pa)
%       d50  = median sediment grain size diameter (m)
%       ws   = sediment fall velocity (m/s)
%       me   = erosion rate coeficient (kg/N/s)
%       g    = acceleration due to gravity (m/s2)
%       gamma = Dronkers tidal asymetry coefficient (-)
% OUTPUTS
%   formdims is a struct containing:
%       La = e-folding length for area (m)
%       Lw = e-folding length for width (m)
%       hm = MTL hydraulic depth at mouth (m)
%       Wm = MTL Width at mouth (m)
%       Am = MTL CSA at mouth (m2)
%       Ucr = peak tidal amplitude (m/s)
%       Ph = tidal prism based on hydraulics (m3)
%       Ph = tidal prism based on geometry (m3)
%       r = hypsometry exponent (-)
% NOTES
%   v2 - See 'Tidal form solver.docx' for details of equations used.
%   v1 - replaces earlier version that used fzero to find root
% SEE ALSO
%   ckfa_wave_profile.m and ckfa_form_model.m
%
% Author: Ian Townend
% CoastalSEA (c)June 2018 updated 2020. converted for ChannelForm App 2022
%--------------------------------------------------------------------------
%
    formdims = [];
    fhL = @(hme) fun_hm(inp,hme);
    %unconstrained optimisation
    options = optimset('MaxIter',500,'TolFun',1e-9,'TolX',1e-3);
    [hm,~,ok] = fminsearch(fhL,hme,options);    %find minimum of fhL
    if ok<1, return; end    

    gamma = inp.gamma;
    % gamma = inp.Le/(sqrt(inp.g*hme)*12.4*3600);
    r = hypsometry_exponent(inp.amp,hm,gamma,true);
    Lw = widthConvergence(inp,hm,r);
    La = areaConvergence(inp,hm,Lw);
    Ucr= inp.amp*inp.omega*La/hm;
    [La,Lw] = waveConvergence(inp,hm,La,Lw,Ucr);    
    Wm = inp.Wrv*exp(inp.Le/Lw);
    Am = inp.Arv*exp(inp.Le/La);
    
    [Ph,Pf]= tidalPrism(inp,La,Lw,Ucr,hm,Am,Wm,r);
    formdims = table(La,Lw,hm,Wm,Am,Ucr,Ph,Pf,r); 
end
%%
function fy = fun_hm(inp,hme)
    % Find the hydraulic depth based on a balance of eorsion and deposition 
    % with Uc determined from the hydraulics and Lw from equating
    % hydraulic and geometric prism estimates
    
    %trap depths less that tidal amplitude
    if hme<inp.amp, fy = 1; return; end
    
    %current only shear stress
    Cdc = frictionCoeff(inp,hme);   

    gamma = inp.gamma;
    %gamma = inp.Le/(sqrt(inp.g*hme)*12.4*3600);
    r = hypsometry_exponent(inp.amp,hme,gamma,true);  %hypsometry exponent
    Lw = widthConvergence(inp,hme,r);            %width covergence length                  
    La = areaConvergence(inp,hme,Lw);            %area convergence length    
    %
    Uc  = inp.amp*inp.omega*La/hme;              %tidal velocity amplitude
    idu = Uc<0;
    Uc(idu) = 0;
    %
    tau = inp.rhow*Cdc*Uc^2;                     %bed shear stress

    if inp.Uw>0                                  %winds included
        Wm = inp.Wrv*exp(inp.Le/Lw);
        Fch = sqrt(2)*Wm; 
        [Hs, Tp, ~] = tma_spectrum(inp.Uw,inp.zw,Fch,hme,hme);
        Hrms = Hs/sqrt(2);
        [La,~] = waveConvergence(inp,hme,La,Lw,Uc);
        Uc  = inp.amp*inp.omega*La/hme;          %tidal velocity amplitude
        %shear stress under combined + aligned waves and current
        tauall = tau_bed(hme,inp.d50,inp.visc,inp.rhow,Uc,Hrms,Tp,0);
        tau = tauall.taur;
    end

    %estimate of erosion
    if tau>inp.taucr
        ero = (tau-2*inp.taucr)*(pi/2-asin(sqrt(inp.taucr/tau)));
        ero = inp.me/pi*(ero+sqrt(inp.taucr*(tau-inp.taucr)));
    else
        ero = 0;
    end
    %estimate of deposition
    sed = inp.rhoc*inp.ws;
    %
    fy = abs(ero-sed);                    %balance of erosion and depostion
end

%%
function Lw = widthConvergence(inp,hm,r)
    %approximate width or area convergence length for given depth 
    k  = inp.omega/sqrt(inp.g*hm);                  %wave number
    beta = ((r*hm-inp.amp)/(r*hm))^(r-1);           %stream width ratio
    eH = pi()*inp.amp*beta/4/hm;                    %amplitude-depth ratio

    Lw = 1/k*atan2(2*eH/(1+eH^2),(1-eH^2)/(1+eH^2));%width covergence length
    eL = 1-exp(-inp.Le/Lw);                         %length correction
    %inclusion of length adjustment
    if eH^2+eL^2-1>0      
        fact1 = sqrt(eH^4+eH^2*eL^2-eH^2);        
        fact2 = eH^2+eL^2;
        term1 = -(eL*(1+fact1)/fact2-1)/eH;
        term2 = (eL+fact1)/fact2;
        LwA = 1/k*atan2(term1,term2);

        term3 = -(eL*(1-fact1)/fact2-1)/eH;
        term4 = -(-eL+fact1)/fact2;
        LwB = 1/k*atan2(term3,term4);
        
        Lw = max(LwA,LwB);
    end   

    % Lw = 2/k*atan(eH);                              %width covergence length
    %for the cases examined this is the same as:
    % Lw = 1/k*atan2(2*eH/(1+eH^2),(1-eH^2)/(1+eH^2));
    %which simiplifies to Lw = 1/k*atan(2*eH/(1-eH^2)) for eH<1
    % if eH<1
    %     Lw = 1/k*atan(2*eH/(1-eH^2));      %width overgence length
    % else
    %     Lw = 2/k*atan(eH);                 %width covergence length     
    %inclusion of length adjustment ie only assume La=Lw
    % eL = 1-exp(-inp.Le/Lw);
    % if eH^2+eL^2-1>0       
    %     Lw = 2/k*atan((eH+sqrt(eH^2+eL^2-1))/(eL+1));
    % end

    Lw(Lw<0) = 0;
end
%%
function La = areaConvergence(inp,hm,Lw)
    %approximate area convergence length for given depth and width
    %convergence 
    k = inp.omega/sqrt(inp.g*hm);                   %wave number  
    fLa = @(LA) Lw*cos(mod(k*LA,pi/2))-LA;
    %unconstrained optimisation
    options = optimset('MaxIter',1000,'TolFun',1e-9,'TolX',1e-6);
    [La,~,ok] = fzero(fLa,Lw,options);    %find minimum of fhL
    if ok<1, La = Lw; end 
end

%%
function Cd = frictionCoeff(inp,h)
    %drag coefficient for water depth,h, using, d50, taucr, visc, rhow
    if h<=0, Cd=0; return; end
    a    = 0.0001615; b=6; c=-0.08;
    fact = inp.taucr/inp.rhow/a.*(h/inp.visc).^2;
    A    = b*c/2*(fact).^(c/2);
    LW   = lambertw(A);
    ucs  = sqrt(((inp.visc./h).^2).*exp(log(fact)-2*LW/c));
    cds  = a*exp(b*(ucs.*h/inp.visc).^c);
    %rough turbulent current
    zo  = inp.d50/12;
    cdr  = (0.4./(log(h/zo)-1)).^2;
    ucr = sqrt(inp.taucr/inp.rhow./cdr);
    if ucr<ucs, Cd = cdr; else, Cd = cds; end
end

%%
function [La,Lw] = waveConvergence(inp,hm,La,Lw,Uc)
    %adjust convergence lengths to account for the effect of waves
    Wm = inp.Wrv*exp(inp.Le/Lw);
    %use effective width We, to define angle and fetch
    %We = inp.Cwe*Wm;      %effective mean tide width (eg W @ Lw/2) 
    %phi = Wm/2/inp.Le*(1-exp(-inp.Le/Lw));
    phi = Wm/Lw;           %angle of shore to channel axis at mouth
    Lstar = Uc/inp.omega*tan(phi);      %forcLw one tidal flat
    Wlw = Wm-2*Lstar;
    if Wlw<0
        Wlw = 10;                       %minimum channel width
    end

    %get the increase in width and csa at mtl
    Fmt = sqrt(2)*Wm;                   %fetch at mean tide
    [~,ymt] = ckfa_wave_profile(inp,inp.Uw,Fmt,inp.rhoc,hm,inp.amp);

    %width and depth at low water
    Flw = sqrt(2)*Wlw;                  %fetch at low water
    if hm>inp.amp
        hlw = hm-inp.amp;               %hydraulic depth at low water
    else
        hlw = 0.5;                      %minimum value for drainage channel
    end
    conc = inp.rhoc*hm/hlw;             %modified concentration at low water
    [~,ylw,alw] = ckfa_wave_profile(inp,inp.Uw,Flw,conc,hlw,hlw);
    
    %update convergence lengths La and Lw
    Am = inp.Arv*exp(inp.Le/La);
    Amw = Am+alw+(ymt+ylw)*inp.amp;  %2 banks ->  2(dy1+dy2)/2*a*L
    La = -inp.Le/(log(inp.Arv/Amw));

    Wmw = Wm+2*ymt;
    Lw = -inp.Le/(log(inp.Arv/Wmw));
end

%%
function [Ph,Pf] = tidalPrism(inp,La,Lw,Uc,hm,Am,Wm,r)
    % Tidal prism
    w = inp.omega;
    amp = inp.amp;
    k = waveNumber(inp,hm,La,r);                    %wave number
    eA = mod(k*La,pi/2);                            %wave number ratio      
    % 
    eL = (1-exp(-inp.Le/Lw));                       %length correction
    beta = ((r*hm-amp)/(r*hm))^(r-1);               %stream width ratio
    eh = amp*beta/hm;                               %amplitude-depth ratio
    Pf = 2*amp*cos(eA)*Wm*Lw*eL;                    %form prism
    Ph = Uc*Am*(4-pi*eh*sin(eA))/2/w;               %hydraulic prism
    
    %the full integral using hypsometry gives very similar results
    % f_int = @(t) sin(w.*t-eA).*(1+amp.*cos(w.*t).*beta./hm).^r;
    % Ph = Uc*Am*integral(f_int,eA/w,(eA+pi)/w);   %eq.15
end

%%
function k = waveNumber(inp,hm,La,r)
    %compute the wave number based on Eq.23 in F&A'94
    Cr = (8/3/pi/inp.g)*(inp.amp*inp.omega^2);      %constants
    beta = ((r*hm-inp.amp)/(r*hm))^(r-1);           %stream width ratio
    Cd = frictionCoeff(inp,hm);                     %friction coefficient
    k = Cr*beta*Cd*La^2/hm^3;                       %wave number
end

%%
function [r,gma] = hypsometry_exponent(amp,Hm,gamma,isDronk)
    %get fucntion to set the hypsometry exponent and central depth
    %using hydraulic depth and tidal amplitude for reaches
    if isDronk  %uses Dronkers gamma
        func = @(r) abs(((r*Hm+amp)/(r*Hm-amp))^(3-r)-gamma);
    else        %uses length-wave length ratio
        func = @(r) abs(((r*Hm-amp)/(r*Hm+amp))^(r-1)-gamma);
    end
    options = optimset('TolFun',1e-9,'TolX',1e-9);
    r = fminbnd(func,1,3,options);
    gma = ((r*Hm+amp)/(r*Hm-amp))^(3-r);
end






