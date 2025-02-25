function cf_summarygridplot(obj)   
%
%-------function help------------------------------------------------------
% NAME
%   cf_summarygridplot.m
% PURPOSE
%   plot initial grid, new grid and difference
% USAGE
%   cf_summarygridplot(obj)
% INPUTS
%   obj - instance of CF_TransModel class
% OUTPUT
%   summary plot
% SEE ALSO
%   CF_Transgression Model.summaryGridPlots, which generates the same plot
%   at the end of the model run using data held in transient properties
%
% Author: Ian Townend
% CoastalSEA (c) Jan 2022
%--------------------------------------------------------------------------
%
    dst = obj.Data.Transgression;
    slr = dst.SLR(end);
    
    grid0 = getGrid(obj,1);
    grid1 = getGrid(obj,height(dst));

    zdiff = grid1.z-grid0.z;
    startyr = obj.RunParam.RunProperties.StartYear;
    endyr = cellstr(grid1.t);
    xlimit = max(grid0.x);

    figure('Name','Grid Plot','Units','normalized','Tag','PlotFig');    
    s1 = subplot(3,1,1);
    surf(grid0.x,grid0.y,grid0.z','FaceColor','interp','EdgeColor', 'none');    
    view(2); 
    h1 = colorbar; 
    s1.XLim(2) = xlimit;
    hold on 
    plot(grid0.xM,0,'+r','DisplayName','xM')
    hold off
    h1.Label.String = 'Elevation (mAD)';
    title(sprintf('Intiial form, T = %s',num2str(startyr)))

    s2 = subplot(3,1,2);
    surf(grid1.x,grid1.y,grid1.z','FaceColor','interp','EdgeColor', 'none');
    view(2); 
    h2 = colorbar; 
    s2.XLim(2) = xlimit;
    hold on 
    plot(grid1.xM,0,'+r','DisplayName','xM')
    hold off
    h2.Label.String = 'Elevation (mAD)';
    title(sprintf('Form at time = %s',endyr{1}))

    s3 = subplot(3,1,3);
    surf(grid1.x,grid1.y,zdiff','FaceColor','interp','EdgeColor', 'none');            
    view(2); 
    colormap(s3,cmap_selection(20));
    clim([-2*slr,2*slr])      %only show +/-2*slr change
    h3 = colorbar; 
    s3.XLim(2) = xlimit;
    hold on 
    plot(grid1.xM,0,'+r','DisplayName','xM')
    hold off
    h3.Label.String = 'Change in elevation (m)';
    title('Difference plot (showing +/-2.slr)')
    casedesc = sprintf('%s, slr=%0.3g m',obj.Data.Grid.Description,slr);
    sgtitle(casedesc,'FontWeight','normal','FontSize',10); 
end