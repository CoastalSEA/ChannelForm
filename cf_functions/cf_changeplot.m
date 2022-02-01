function cf_changeplot(obj)
%
%-------function help------------------------------------------------------
% NAME
%   cf_sectionplot.m
% PURPOSE
%   plot volume difference over the model run and the variation of sediment
%   exchange with transgression distance
% USAGE
%   cf_sectionplot(obj)
% INPUTS
%   obj - instance of CF_TransModel class
% OUTPUT
%   summary plot
% SEE ALSO
%   CF_Transgression Model.changePlot, which generates the same plot
%   at the end of the model run using data held in transient properties
%
% Author: Ian Townend
% CoastalSEA (c) Jan 2022
%--------------------------------------------------------------------------
%
    hf = figure('Name','Change Plot','Units','normalized',...
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
    netchangePlot(obj,ax1,slr);

    ax2 = subplot(1,2,2,'Parent',p);
    rateofchangePlot(obj,ax2);

    %to replace the default background with a white background uncomment
    %used for figure plots for papers
    %hf.Color = [1,1,1];
    %ax1.Color = [1,1,1];
    %ax2.Color = [1,1,1];            
end
%%
function netchangePlot(obj,ax,slr)
    %plot the difference between the original and translated channel
    %origin of the plot is at the offset distance (amount translated) 
    % NB uses transient properties so only works in current session
    dst = obj.Data.Transgression;
    
    grid0 = getGrid(obj,1);
    z0 = grid0.z; %intial grid

    grid1 = getGrid(obj,height(dst));
    xs = grid1.x; 
    ys = grid1.y;
    zs = grid1.z;%final grid

    zdiff = zs-z0;

    wls = obj.Data.Plan;
    yzhw = wls.Whw/2;
    yzlw = wls.Wlw/2;
    zyz = ones(size(xs))*max(max(zdiff));
    
    
    adX = dst.estdX(end);  %cumulative estuary trangression
    cdX = dst.cstdX(end);  %cumulative coastal trangression

    %contourf(ax,xi,yi2,zgrd','LineStyle', 'none');
    surf(ax,xs,ys,zdiff','FaceColor','interp','EdgeColor', 'none');
    colormap(cmap_selection(20));
    caxis([-2*slr,2*slr])      %only show +/-2*slr change
    view(2);
    hold on
    plot3(xs,yzhw,zyz,'--k');  %high water line
    plot3(xs,-yzhw,zyz,'--k');
    plot3(xs,yzlw,zyz,'-.k');  %low water line
    plot3(xs,-yzlw,zyz,'-.k');
    hold off
    xlabel('Distance along channel (m)');
    ylabel('Width (m)');
    h_c = colorbar(ax);   
    h_c.Label.String = 'Change in elevation (m) for +/-2.slr';
    timetxt = cellstr(grid1.t);
    infotxt = sprintf('Dashed lines are Hw/Lw at T=%s',timetxt{1});
    tltxt = sprintf('dV=0 for: SLR = %0.1g m, Transgression = %d m, Coast erosion = %d m\n%s',...
                                slr,round(adX),round(cdX),infotxt);
    title(tltxt,'FontSize',10);
end  
%%
function rateofchangePlot(obj,ax)
    %plot the volume difference for different distanes of transgression
    % NB uses transient properties so only works in current session
    
    dst = obj.Data.Transgression;
    
    delX = dst.delX(end);
    dx = [0,delX/2,delX,delX+delX/2];  %x intervals used            
    dV = dst.vdiffx(end,:);      %volume change for dx
    trnobj = obj.RunParam.CF_TransData;
    sedvol = dst.sedVol(end);
    estdX = dst.estdX(end);
    
    plot(ax,dx,dV,'-k','LineWidth',1)
    ax.YGrid = 'on';
    xlabel('Landward transgression (m)');
    ylabel('Volume change (m^3)');
    %add explanatory text
    mxX = ax.XLim(2); mxY = ax.YLim(2);
    hold on
    plot([0,mxX],[0,0],'-.r');  %'Color',[0.96,0.96,0.96]
    plot([0,estdX],[sedvol,sedvol],'--b'); %transgression point
    plot([estdX,estdX],[0,sedvol],'--b'); %transgression point
    hold off
    txt1 = 'Channel only';
    if trnobj.inclFloodPlain
        txt1 = 'Include flood plain';
    end
    txt2 = 'No bed constraint';
    if trnobj.inclGeoConstraint
        txt2 = 'Including bed constraints';
    end
    txt3 = 'No HW constraint';
    if trnobj.inclHWConstraint
        txt3 = 'Constrained at HW';
    end

    title(sprintf('%s; %s; %s\n',txt1,txt2,txt3),'FontSize',10);
    text(mxX,mxY/10,'Sediment Import  ','HorizontalAlignment','right');
    text(mxX,-mxY/10,'Sediment Export  ','HorizontalAlignment','right');
end