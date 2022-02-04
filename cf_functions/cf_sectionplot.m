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
    slr = dst.dSLR(end);
    p = uipanel('Parent',hf,'BorderType','none'); 
    p.Title = sprintf('Change plot, slr=%0.1g',slr);
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
    dst = obj.Data.Form;
    
    grid0 = getGrid(obj,1);
    z0 = grid0.z; %intial grid

    grid1 = getGrid(obj,height(dst));
    xi = grid1.x; 
    yi = grid1.y;
    zi = grid1.z;%final grid

    if grid1.ishead  %orientation of x-axis, x=0 is nearest the mouth if ishead=false
        zi = flipud(zi);                          
        z0 = flipud(z0);  
    end

    wl0 = obj.Data.WaterLevels.zmt(1);
    wli = obj.Data.WaterLevels.zmt(end);

%     figure('Name','XS Plot','Units','normalized','Tag','PlotFig');                                                    
%     ax = axes('Tag','PropertyPlot');
    noxi=length(xi);
    ix0 = find(xi>=grid1.xM-eps,1,'first'); 
    noxi = noxi-ix0;
    nx1=ix0;  nx2=ceil(ix0+0.1*noxi);  nx3=ceil(ix0+0.2*noxi);

    green = mcolor('green');
    plot(ax,yi,z0(nx1,:),'-r','LineWidth',0.6);
    hold on
    plot(ax,yi,zi(nx1,:),'--r','LineWidth',0.6);

    plot(ax,yi,z0(nx2,:),'-b','LineWidth',0.58);
    plot(ax,yi,zi(nx2,:),'--b','LineWidth',0.58);

    plot(ax,yi,z0(nx3,:),'-','Color',green,'LineWidth',0.56);
    plot(ax,yi,zi(nx3,:),'--','Color',green,'LineWidth',0.56);

    %add water levels at mouth
    plot(xlim, wl0*[1 1],'-','Color',[0.7,0.7,0.7]);
    plot(xlim, wli*[1 1],'--','Color',[0.7,0.7,0.7]);

    hold off
    xlabel('Width (m)'); 
    ylabel('Change in level (m)');
    hL=legend('0pre','0post','0.1pre','0.1post','0.2pre','0.2post','Location','SouthEast');
    set(hL, 'Color', 'none');
    casedesc = sprintf('Cross-sections difference plot, slr=%0.1g m',slr);
    title(casedesc,'FontWeight','normal','FontSize',10);            
end
%%
function thalwegPlot(obj,ax,slr)
    %plot the thalweg at the start and end of run
    % NB uses transient properties so only works in current session
    dst = obj.Data.Form;
    
    grid0 = getGrid(obj,1);
    z0 = grid0.z(:,ceil(size(grid0.z,2)/2)); %intial grid

    grid1 = getGrid(obj,height(dst));
    xi = grid1.x; 
    zi = grid1.z(:,ceil(size(grid1.z,2)/2));%final grid

    if grid1.ishead  %orientation of x-axis, x=0 is nearest the mouth if ishead=false
        zi = flipud(zi);                          
        z0 = flipud(z0);  
    end

    wl0 = obj.Data.WaterLevels.zmt(1);
    wli = obj.Data.WaterLevels.zmt(end);

    plot(ax,xi,z0(:,1),'-k');            
    hold on
    plot(ax,xi,zi(:,1),'--k');

    %add water levels at mouth
    plot(xlim, wl0*[1 1],'-','Color',[0.7,0.7,0.7]);
    plot(xlim, wli*[1 1],'--','Color',[0.7,0.7,0.7]);  

    hold off
    xlabel('<= mouth       Distance along channel (m)       head =>'); 
    ylabel('Elevation (mAD)');
    hL=legend('Initial','Final','Location','SouthEast');
    set(hL, 'Color', 'none');
    casedesc = sprintf('Centre-line difference plot, slr=%0.1g',slr);
    title(casedesc,'FontWeight','normal','FontSize',10); 
end