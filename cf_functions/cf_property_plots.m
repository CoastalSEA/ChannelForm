function cf_property_plots(obj,src)
%
%-------function help------------------------------------------------------
% NAME
%   cf_property_plots.m
% PURPOSE
%   plots displayed on Proprety tab in ChannelForm model
% USAGE
%   cf_model_tabs(obj,mobj,src)
% INPUTS
%   obj - instance of any form model that uses the GDinterface abstract class
%   src - handle to calling popup menu or create Figure button
% OUTPUT
%   generate plot or properties table on tab
% SEE ALSO
%   used to handle different types of plot accessed from the Properties tab
%
% Author: Ian Townend
% CoastalSEA (c) Jan 2022
%--------------------------------------------------------------------------
%
    tabsrc = src.Parent;           %handle to Form-Props tab
    zprop = obj.Data.Form;         %dstable of elevations and widths
    hprop = obj.Data.Hypsometry;   %dstable of hypsometry data
    sprop = obj.Data.SectionProps; %dstable of section properties
    cdesc = obj.Data.Form.Description;
    hydobj = obj.RunParam.CF_HydroData;

    if strcmp(src.Tag,'FigButton')
        hfig = figure('Tag','PlotFig');
        ax = axes('Parent',hfig,'Tag','PropertyPlot','Units','normalized');
        src = findobj(tabsrc,'Tag','PlotList'); 
    else
        ht = findobj(tabsrc,'Type','axes'); %remove any existing plot
        delete(ht); 
        ax = axes('Parent',tabsrc,'Tag','PropertyPlot',...
           'Units','normalized','Position', [0.08,0.11,0.86,0.58]);
    end
 
    %water levels at mouth
    wl_0 = zeros(1,3);
    wl_0(1) = hydobj.zhw(1);
    wl_0(2) = hydobj.zmt(1);    
    wl_0(3) = hydobj.zlw(1);   
    
    switch src.String{src.Value}
        case 'Hypsommetry'            
            zcentre = hprop.Dimensions.Z;
            zsurf = hprop.SurfaceArea;
            zvol = hprop.Volume;            
            hypsommetryplot(ax,zcentre,zsurf,zvol,wl_0,cdesc)
        case 'Cross-sections'
            xi = zprop.Dimensions.X;
            yi = zprop.Dimensions.Y;
            zi = squeeze(zprop.Z);    
            crossectionplot(ax,xi,yi,zi,wl_0,cdesc)            
        case 'Centre-line'
            ax.Position(3) = 0.82;
            xi = zprop.Dimensions.X;
            yz = reshape(zprop.DataTable{1,2:end},length(xi),3);
            zi = squeeze(zprop.Z);
            zi = zi(:,size(zi,2)/2+1:end);
            centrelineplot(ax,xi,yz,zi,cdesc)
%             sT = h_data.StepData{idr};
%             fprintf('Le0 = %.0f; Le = %.0f\n',sT.Length(1),sT.Length(end));
        case 'Form width'
            xj = sprop.Dimensions.X;  
            L = zeros(length(xj),3); W = L;
            for i=1:3
                W(:,i) = sprop.DataTable{:,i}; %selects Widths 
                L(i) = getconvergencelength(xj,W(:,i));
            end
            formwidthplot(ax,xj,W,L,cdesc)
        case 'Elevation-Area histogram'
            zcentre = hprop.Dimensions.Z;
            zhist = hprop.SAfreq; %SAfreq = histogram bin count x grid area
            elevationareaplot(ax,zcentre,zhist,wl_0,cdesc)
        case 'Hydraulic depth'
            zcentre = hprop.Dimensions.Z;
            zsurf = hprop.SurfaceArea;
            zvol = hprop.Volume;          
            hydraulicdepthplot(ax,zcentre,zsurf,zvol,wl_0,cdesc)
        case 'Area-Prism ratio'
            xj = sprop.Dimensions.X;
            AP = sprop.AoP;
            prismareaplot(ax,xj,AP,cdesc)
        case 'Cross-sectional area' 
            xj = sprop.Dimensions.X;   
            L = zeros(length(xj),3); W = L;
            for i=1:3
                A(:,i) = sprop.DataTable{:,3+i}; %selects CSAs 
                L(i) = getconvergencelength(xj,A(:,i));
            end
            csaplot(ax,xj,A,L,cdesc)
        case 'Prism'
            xj = sprop.Dimensions.X;  
            Pr = sprop.Prism;
            prismplot(ax,xj,Pr,cdesc)
        case 'a/h and Vs/Vc'
            xj = sprop.Dimensions.X';  
            Ahw = sprop.CSAhw;
            Amt = sprop.CSAmt;
            Alw = sprop.CSAlw;   
            Whw = sprop.Whw;
            Wmt = sprop.Wmt;
            Wlw = sprop.Wlw;
            amp = hydobj.zhw-hydobj.zlw;   %transient!!!!
            xi = zprop.Dimensions.X';
            amp = interp1(xi,amp,xj);
            
            ah(:,1) = amp./(Ahw./Whw);
            ah(:,2) = amp./(Amt./Wmt);
            ah(:,3) = amp./(Alw./Wlw);
            ah(ah>1.2) = NaN; ah(ah<0) = NaN;
            VsVc = (Ahw-Alw-2*amp.*Wlw)./(Alw+amp.*Wlw);
            VsVc(VsVc>3) = NaN;
            asymmratiosplot(ax,xj,ah,VsVc,cdesc);             
        case 'Hydraulics'
            xi = zprop.Dimensions.X;
            hydraulicsplot(ax,hydobj,xi,cdesc); 
        case 'Transgression'
            if ~isa(obj,'CF_TransModel')
                warndlg('Trangression model required for this plot');
                return;
            end
            steptable = obj.StepData;
            transgressionplot(ax,steptable)
        case 'xHydraulics'
%             cst = h_data.(aprop{4}){idr};
            xi = zprop.Dimensions.X;
            yi = zprop.Dimensions.Y;
            yz = reshape(zprop.DataTable{1,2:end},length(xi),3);
            zi = squeeze(zprop.Z);
%             zi = zi(:,size(zi,2)/2+1:end);
            projectionplot(ax,xi,hydobj,yi,yz,zi);
    end   
end
%%
function hypsommetryplot(ax,zcentre,zsurf,zvol,wl,casedesc)
    %plot volume and surface area as function of elevation 
    plot(ax,zvol,zcentre,'-r');
    hold on
    plot(ax,zsurf,zcentre,'-.b');
    
    %add water levels at mouth
    plot(xlim, wl(1)*[1 1],':','Color',[0.7,0.7,0.7]);
    plot(xlim, wl(2)*[1 1],'--','Color',[0.8,0.8,0.8]);
    plot(xlim, wl(3)*[1 1],':','Color',[0.7,0.7,0.7]);
    
    xlabel('Volume (m^3) and Area (m^2)'); 
    ylabel('Elevation (mAD)');
    legend('Volume','Surface area','Location','SouthEast');
    title(casedesc,'FontWeight','normal','FontSize',10);
    hold off
end
%%
function crossectionplot(ax,xi,yi,zgrd,wl,casedesc)
    %plot a series of cross-sections along length of channel
%     figure;  %over-ride plot to UI window and plot as stand alone figure
%     ax = axes;    
    noxi=length(xi);
    nx1=noxi;nx2=ceil(9*noxi/10);nx3=ceil(4*noxi/5);
    nx4=ceil(3*noxi/5);nx5=ceil(2*noxi/5);nx6=ceil(noxi/5);
    plot(ax,yi,zgrd(nx1,:),'-r');
    hold on
    plot(ax,yi,zgrd(nx2,:),'-.b');
    plot(ax,yi,zgrd(nx3,:),'--g');
    plot(ax,yi,zgrd(nx4,:),'-.c');
    plot(ax,yi,zgrd(nx5,:),'--y');
    plot(ax,yi,zgrd(nx6,:),'-.m');
    plot(ax,yi,zgrd(1,:),':k');
    
    %add water levels at mouth
    plot(xlim, wl(1)*[1 1],':','Color',[0.7,0.7,0.7]);
    plot(xlim, wl(2)*[1 1],'--','Color',[0.8,0.8,0.8]);
    plot(xlim, wl(3)*[1 1],':','Color',[0.7,0.7,0.7]);
    
    xlabel('Width (m)'); 
    ylabel('Elevation (mAD)');
    hL=legend('0','0.1L','0.2L','0.4L','0.6L','0.8L','R','Location','SouthEast');
    set(hL, 'Color', 'none');
    title(casedesc,'FontWeight','normal','FontSize',10);
    hold off
end
%%
function centrelineplot(ax,xi,yz,zi,casedesc)
    %plot centrel-line depth and width variation along channel
%     figure;  %over-ride plot to UI window and plot as stand alone figure
%     ax = axes;
    yyaxis left
    plot(ax,xi,zi(:,1),'Color','k');
    xlabel('<= head       Distance along channel (m)       mouth =>');  
    ylabel('Elevation (mAD)');
    yyaxis right
    plot(ax,xi,squeeze(yz(:,1)),'Color','g','LineStyle','--');
    hold on
    plot(ax,xi,squeeze(yz(:,2)),'Color','r','LineStyle',':');
    plot(ax,xi,squeeze(yz(:,3)),'Color','b','LineStyle','-.');
    hold off
    ylabel('Half-width (m)');
    hL=legend('Centre Line','High water', 'Mean tide level','Low water',...
                                                    'Location','West');
    set(hL, 'Color', 'none');
    title(casedesc,'FontWeight','normal','FontSize',10);
end
%%
function formwidthplot(ax,xi,W,L,casedesc)
    %plot centrel-line depth and width variation along channel
%     figure;  %over-ride plot to UI window and plot as stand alone figure
%     ax = axes;
    %order of W and L is HW, MT, LW
    ax.Position(1) = 0.1;
    plot(ax,xi,W(:,1),xi,W(:,2),xi,W(:,3));
    xlabel('<= head       Distance along channel (m)       mouth =>');
    ylabel('Width (m)');
    txt1 = sprintf('High water (L_W = %.0f)',L(1));
    txt2 = sprintf('Mean tide  (L_W = %.0f)',L(2));
    txt3 = sprintf('Low water  (L_W = %.0f)',L(3));
    hL=legend(txt1,txt2,txt3,'Location','North');                                                                                         
    set(hL, 'Color', 'none');
    title(casedesc,'FontWeight','normal','FontSize',10);
end
%%
function elevationareaplot(ax,zcentre,zhist,wl,casedesc)
    %histogram of amount of surface area at each elevation
    barh(ax,zcentre, zhist, 'histc'); 
    hold on   
    %add water levels at mouth
    plot(xlim, wl(1)*[1 1],':','Color',[0.7,0.7,0.7]);
    plot(xlim, wl(2)*[1 1],'--','Color',[0.8,0.8,0.8]);
    plot(xlim, wl(3)*[1 1],':','Color',[0.7,0.7,0.7]);
    hold off
    xlabel('Surface area (m^2)');
    ylabel('Elevation (mAD)');
    title(casedesc,'FontWeight','normal','FontSize',10);
end
%%
function hydraulicdepthplot(ax,zcentre,zsurf,zvol,wl,casedesc)
    %plot variation of hydraulic depth with elevation
    plot(ax,zvol./zsurf,zcentre);
    hold on   
    %add water levels at mouth
    plot(xlim, wl(1)*[1 1],':','Color',[0.7,0.7,0.7]);
    plot(xlim, wl(2)*[1 1],'--','Color',[0.8,0.8,0.8]);
    plot(xlim, wl(3)*[1 1],':','Color',[0.7,0.7,0.7]);
    hold off
    xlabel('Hydraulic depth (m)'); 
    ylabel('Elevation (mAD)');
    title(casedesc,'FontWeight','normal','FontSize',10);
end
%%
function prismareaplot(ax,xj,APratio,casedesc)
    %plot area-prism ratio along channel
    plot(ax,xj,APratio);
    xlabel('<= head       Distance along channel (m)       mouth =>'); 
    ylabel('Area-Prism ratio');
    title(casedesc,'FontWeight','normal','FontSize',10);
end
%%
function csaplot(ax,xj,A,L,casedesc)
    %plot cross-sectional area along channel  
    %order of A and L is HW, MT, LW
    plot(ax,xj,A(:,1),xj,A(:,2),xj,A(:,3));
    xlabel('<= head       Distance along channel (m)       mouth =>');
    ylabel('CSA');
    txt1 = sprintf('High water (L_A = %.0f)',L(1));
    txt2 = sprintf('Mean tide  (L_A = %.0f)',L(2));
    txt3 = sprintf('Low water  (L_A = %.0f)',L(3));
    hL=legend(txt1,txt2,txt3,'Location','North');
    set(hL, 'Color', 'none');   
    title(casedesc,'FontWeight','normal','FontSize',10);
end
%%
function prismplot(ax,xj,Pr,casedesc)
    %plot prism along channel
    plot(ax,xj,Pr);
    xlabel('<= head       Distance along channel (m)       mouth =>');
    ylabel('Prism');
    title(casedesc,'FontWeight','normal','FontSize',10);
end
%%
function asymmratiosplot(ax,xj,ah,VsVc,casedesc)
    %plot a/h and Vs/Vc along channel
%     figure;  %over-ride plot to UI window and plot as stand alone figure
%     ax = axes;
    ax.Position = [0.08,0.11,0.84,0.58];
    yyaxis left
    plot(ax,xj,ah(:,1),'Color','g','LineStyle','-');
    hold on
    plot(ax,xj,ah(:,2),'Color','r','LineStyle','-');
    plot(ax,xj,ah(:,3),'Color','b','LineStyle','-');
    hold off
    ylabel('a/h');
    yyaxis right
    plot(ax,xj,VsVc);
    ylabel('Vs/Vc');
    xlabel('<= head       Distance along channel (m)       mouth =>');
    hL=legend('a/h High water', 'a/h Mean tide level','a/h Low water',...
                                               'Vs/Vc','Location','best');
    set(hL, 'Color', 'none');
    title(casedesc,'FontWeight','normal','FontSize',10);
end
%%
function hydraulicsplot(ax,hydobj,x,casedesc)
    %plot the results from the CST hydraulic model or linear decay
    
    %over-ride plot to UI window and plot as stand alone figure
    % figure('Name','Tab plot','Tag','PlotFig');  
    % ax = axes; 
    green = mcolor('green');
    orange = mcolor('orange');
    
    zhw = hydobj.zhw;  %high water level
    zmt = hydobj.zmt;  %mean tide level
    zlw = hydobj.zlw;  %low water level
    
    incvelocity = false;  %used to plot water levels without velocities
    if isfield(hydobj.cstres,'U')
        U = hydobj.cstres.U;  %tidal velocity amplitude
        v = hydobj.cstres.v;  %river velocity 
        d = hydobj.cstres.d;  %hydraulic depth
        incvelocity = true;  %used to plot water levels without velocities
    end
    yyaxis(ax,'left')
    cla                                    %clear any existing plot lines
    plot(ax,x,zmt,'-r','LineWidth',1.0);   %plot time v elevation
    hold on
    plot(ax,x,zhw,'-.b','LineWidth',0.8)   %plot high water level
    plot(ax,x,zlw,'-.b','LineWidth',0.8)   %plot low water level    
    ax.YLim(2) = max(zhw)+0.1;
    ylabel('Elevation (mOD)'); 
    
    if incvelocity
        plot(ax,x,zmt-d,'-k','LineWidth',0.6); %hydraulic depth below mean tide level
        ax.YLim(1) = min(zmt-d);
        yyaxis(ax,'right')
        cla                                %clear any existing plot lines
        plot(ax,x,U,'--','Color',orange,'LineWidth',0.6)%plot tidal velocity
        plot(ax,x,v,'--','Color',green,'LineWidth',0.6) %plot river velocity
        ylabel('Velocity (m/s)'); 
            legend('MTL','HWL','LWL','Hydraulic depth',...
                'Tidal velocity','River velocity','Location','west');
    else
       legend('MTL','HWL','LWL','Location','west');
    end
    hold off
    ax.XLimMode = 'manual';
    ax.XLim = [min(x),max(x)];
    ax.Position(3) = ax.Position(3)*0.95;
    xlabel('<= head        Distance along channel (m)        mouth =>');
    title(casedesc,'FontWeight','normal','FontSize',10);
end
%%
function transgressionplot(ax,stepT)
    %plot transgression time step results
    t = stepT.Time;
    varnames = stepT.Properties.VariableNames(2:end);
    ok=1;
    hold on
    while ok>0
        [idx,ok] = listdlg('Name','Properties','PromptString','Select variable',...                         
                           'ListSize',[120,200],'SelectionMode','single',...
                           'ListString',varnames);
        if ok==0, continue; end
        var = varnames{idx};
        plot(ax,t,stepT.(var),'DisplayName',var);
    end
    hold off
    xlabel('Time (years)')
    legend
end
    