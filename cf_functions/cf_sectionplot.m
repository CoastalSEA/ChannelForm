function cf_sectionplot(obj)
%
%-------function help------------------------------------------------------
% NAME
%   cf_sectionplot.m
% PURPOSE
%   plot a set of along-channel cross-sections and the channel centre-line
%   (thalweg) for the initial and final grid
% USAGE
%   cf_sectionplot(obj)
% INPUTS
%   obj - instance of CF_TransModel class
% OUTPUT
%   summary plot
% SEE ALSO
%   CF_Transgression Model.sectionPlot and thalwegPlot, which generate
%   the same plot at the end of the model run using data held in transient
%   properties
%
% Author: Ian Townend
% CoastalSEA (c) Jan 2022
%--------------------------------------------------------------------------
%
    hf = figure('Name','Section Plot','Units','normalized',...
                                            'Tag','PlotFig');                                    
    hf.Position(1) = 0.1;
    hf.Position(3) = hf.Position(3)*2;
    
    dst = obj.Data.Transgression;
    slr = dst.SLR(end);
    p = uipanel('Parent',hf,'BorderType','none'); 
    p.Title = sprintf('Case: %s, slr=%0.3g m',obj.Data.Grid.Description,slr);
    p.TitlePosition = 'centertop'; 
    p.FontSize = 12;
    p.FontWeight = 'bold';

    ax1 = subplot(1,2,1,'Parent',p);
    crossectionPlot(obj,ax1,slr);

    ax2 = subplot(1,2,2,'Parent',p);
    thalwegPlot(obj,ax2,slr);

    %to replace the default background with a white background uncomment
    %used for figure plots for papers
    %hf.Color = [1,1,1];
    %ax1.Color = [1,1,1];
    %ax2.Color = [1,1,1];            
end       
%%
function crossectionPlot(obj,ax,slr)
    %plot cross-sections along length of channel at start and end of run
    % NB uses transient properties so only works in current session
    dst = obj.Data.Grid;
    
    grid0 = getGrid(obj,1);
    z0 = grid0.z; %intial grid

    grid1 = getGrid(obj,height(dst));
    xi = grid1.x; 
    yi = grid1.y;
    zi = grid1.z;%final grid

    gd_dir = gd_ax_dir(grid1);
    if gd_dir.x==1 || gd_dir.x==4                 
        %orientation of x-axis, x=0 is nearest the mouth if ishead=false
        zi = flipud(zi);
        z0 = flipud(z0);
    end  

    wl0 = obj.Data.WaterLevels.zmt(1);
    wli = obj.Data.WaterLevels.zmt(end);
    
    noxi=length(xi);
    ix0 = find(xi>=grid1.xM-eps,1,'first'); 
    noxi = noxi-ix0;
    nx1=ix0;  nx2=ceil(ix0+0.2*noxi);  
    nx3=ceil(ix0+0.4*noxi); nx4=ceil(ix0+0.6*noxi);

    green = mcolor('green');  orange = mcolor('orange');
    plot(ax,yi,z0(nx1,:),'-r','LineWidth',0.6);
    hold on
    plot(ax,yi,zi(nx1,:),'--r','LineWidth',0.6);

    plot(ax,yi,z0(nx2,:),'-b','LineWidth',0.58);
    plot(ax,yi,zi(nx2,:),'--b','LineWidth',0.58);

    plot(ax,yi,z0(nx3,:),'-','Color',green,'LineWidth',0.56);
    plot(ax,yi,zi(nx3,:),'--','Color',green,'LineWidth',0.56);

    plot(ax,yi,z0(nx4,:),'-','Color',orange,'LineWidth',0.56);
    plot(ax,yi,zi(nx4,:),'--','Color',orange,'LineWidth',0.56);

    %add water levels at mouth
    plot(xlim, wl0*[1 1],'-','Color',[0.7,0.7,0.7]);
    plot(xlim, wli*[1 1],'--','Color',[0.7,0.7,0.7]);

    hold off
    %reverse x-axis - sections are oriented looking from mouth up-estuary
    ax.XDir = 'reverse';
    %add meta-data
    xlabel('Width (m)'); 
    ylabel('Change in level (m)');
    hL=legend('0L-pre','0L-post','0.2L-pre','0.2L-post','0.4L-pre','0.4L-post','0.6L-pre','0.6L-post','Location','SouthEast');
    set(hL, 'Color', 'none');
    casedesc = sprintf('Cross-sections relative to mouth at end of run, slr=%0.2g m',slr);
    title(casedesc,'FontWeight','normal','FontSize',10);            
end
%%
function thalwegPlot(obj,ax,slr)
    %plot the thalweg at the start and end of run
    % NB uses transient properties so only works in current session
    dst = obj.Data.Grid;
    
    grid0 = getGrid(obj,1);
    z0 = grid0.z(:,ceil(size(grid0.z,2)/2)); %intial grid

    grid1 = getGrid(obj,height(dst));
    xi = grid1.x; 
    zi = grid1.z(:,ceil(size(grid1.z,2)/2));%final grid

    gd_dir = gd_ax_dir(grid1);
    if gd_dir.x==1 || gd_dir.x==4                 
        %orientation of x-axis, x=0 is nearest the mouth if ishead=false
        zi = flipud(zi);
        z0 = flipud(z0);
    end
            
    wl0 = obj.Data.WaterLevels.zmt(1);
    wli = obj.Data.WaterLevels.zmt(end);

    plot(ax,xi,z0(:,1),'-k');            
    hold on
    plot(ax,xi,zi(:,1),'--k');

    %add water levels at mouth
    xoffset = diff(ax.XLim)/100;    yoffset = diff(ax.YLim)/50; 
    plot(xlim, wl0*[1 1],'-','Color',[0.7,0.7,0.7]);
    plot(xlim, wli*[1 1],'--','Color',[0.7,0.7,0.7]);  
    text(xoffset,wl0+yoffset,'Zmt');
    %add position of mouth
    plot([1,1]*grid0.xM,ax.YLim,'-r');
    plot([1,1]*grid1.xM,ax.YLim,'--r');
    ytxt = ax.YLim(1)+yoffset;
    xtxt = grid0.xM+xoffset;
    text(xtxt,ytxt,'Mouth location');
    hold off
    xlabel('<= mouth       Distance along channel (m)       head =>'); 
    ylabel('Elevation (mAD)');
    hL=legend('Initial','Final','Location','SouthEast');
    set(hL, 'Color', 'none');
    casedesc = sprintf('Centre-line difference plot, slr=%0.2g m',slr);
    title(casedesc,'FontWeight','normal','FontSize',10); 
end