function cf_animation(obj)
%
%-------function help------------------------------------------------------
% NAME
%   cf_animation.m
% PURPOSE
%   animation of model run
% USAGE
%   cf_animation(obj)
% INPUTS
%   obj - instance of CF_TransModel class
% OUTPUT
%   animation
% SEE ALSO
%   example of direct call to muiPlots.newAnimation
% NOTES
%   obj is a data class instance and pobj is a plot class muiPlots instance
%
% Author: Ian Townend
% CoastalSEA (c) Jan 2022
%--------------------------------------------------------------------------
%
    gridData = obj.Data.Grid;
    hfig = figure('Name','Animation', ...
                    'Units','normalized', ...
                    'Resize','on','HandleVisibility','on', ...
                    'Visible','off','Tag','PlotFig');
    %create an instance of muiPlots and populate the properties that are
    %needed for the newAnimation method
    pobj = muiPlots.get_muiPlots();   %create new instance          
    pobj.Plot.CurrentFig = hfig;
    pobj.Plot.FigNum = hfig.Number;
    pobj.UIset.callTab = '3DT';
    pobj.UIset.Polar = false;
    pobj.UIset.Type.String = 'surf';
    pobj.UIset.iscmap = false;         %supress option to select colormap
    pobj.AxisLabels.X = 'Distance from mouth (m)';
    pobj.AxisLabels.Y = 'Width (m)';
    pobj.AxisLabels.Y = 'Elevation (mAD)';
    pobj.Legend = [];
    slr = obj.Data.Transgression.SLR;
    msl0 = obj.RunParam.WaterLevels.MSL0;
    pobj.Title = sprintf('Case: %s, slr=%0.3g m',gridData.Description,slr(end));
    %extract the plot data
    pobj.Data.X = gridData.Dimensions.X;
    pobj.Data.Y = gridData.Dimensions.Y;
    pobj.Data.Z = gridData.Z;
    pobj.Data.T = gridData.RowNames;
    xM = gridData.UserData.xM;
    qLt = round(gridData.UserData.Lt/4);
    nint = round(qLt/abs(pobj.Data.X(2)-pobj.Data.X(1)));
    
    ulim = max(pobj.Data.Z,[],'All');
    dlim = min(pobj.Data.Z,[],'All'); 
    msl = msl0+slr(end)-slr;     %timeseries of mean sea level

    %control range of z used in animation    
    prmptxt = {'Upper elevation limit','Lower elevation limit','View, plan=[0,90]'};
            dlgtitle = 'Animation';
            defaults = {num2str(ulim),num2str(dlim),'290 30'};
            answer = inputdlg(prmptxt,dlgtitle,1,defaults);
    if isempty(answer)
        return; 
    else
        uplimit = str2double(answer{1});    %z value of upper limit
        dnlimit = str2double(answer{2});    %z value of lower limit
        aview = str2num(answer{3}); %#ok<ST2NM> %azimuth and elevation for camera view
    end
    
    pobj.Data.Z(pobj.Data.Z>uplimit) = NaN;
    pobj.Data.Z(pobj.Data.Z<dnlimit) = NaN;
    for i=1:length(pobj.Data.T)
        pobj.Data.Z(i,:,:) = pobj.Data.Z(i,:,:)-msl(i);
        pobj.Data.Z(i,1,:) = dlim;
        [~,ixM(i,1)] = gd_basin_indices(getGrid(obj,i)); %nearest grid point 
    end
    %variables to store values that are modified by the function
    t = pobj.Data.T;  %pobj.Data.T is modified by call to convertTime in new3Dplot
    var = pobj.Data.Z;%pobj.Data.Z replaced by 3D plot calls
    
    [figax,hp] = setupAnimation(pobj,var);
    if ~isvalid(pobj.Plot.CurrentFig), return; end
    hold on
    line(figax,'XData',[xM(end),xM(end)],'YData',ylim,'ZData',[msl(end),msl(end)],...
                      'Color','k','LineStyle','-.','LineWidth',0.75);                    
    plot_sections(pobj,figax,var,ixM(1),nint(1),1)
    hold off
    view(aview(1),aview(2))
    
    %change color map (note this does not alter the saved animation)
    idsel = dlim:ulim;         %full range of data to preserve color intervals
    [cmap,~] = landsea(idsel,[dnlimit,uplimit]); %requires Mapping toolbox
    colormap(cmap)
    
    getAnimation(pobj,hp,hfig,t,var,ixM,nint);
    pobj.Data.T = t;     %restore datetime values
    pobj.Data.Z = var;   %restore Z values
    figax.UserData = pobj.Data;  %store data set in UserData to
                                %allow user to switch between plots
    %add replay and slider
    setControlPanel(pobj,hfig,length(t),string(t(1)));
end
%%
function [cmap,climits] = landsea(zi,range)
    %definition of land-sea cmap - requires Mapping toolbox
    if nargin<2, range = [min(zi,[],'All'),max(zi,[],'All')]; end
    cmapsea = [0,0,0.2;  0,0,1;  0,0.45,0.74;  0.30,0.75,0.93; 0.1,1,1];
    cmapland = [0.95,0.95,0.0;  0.1,0.7,0.2; 0,0.4,0.2; 0.8,0.9,0.7;  0.4,0.2,0];
    [cmap,climits] = demcmap(range,128,cmapsea,cmapland);
end
%%
function [figax,hp] = setupAnimation(pobj,var)
    %initialise 3Dplot and setup animation variables
    hfig = pobj.Plot.CurrentFig;
    vari = setTimeDependentVariable(pobj,var,1); 
    pobj.Data.Z = vari;  %first time step
    pobj.UIset.callTab = '3D';
    hfig.Visible = 'on';
    getAplot(pobj);
    pobj.UIset.callTab = '3DT';
    if ~isvalid(hfig), return; end
    %assign axes properties
    figax = gca;                
    figax.ZLimMode = 'manual'; %fix limits of z-axis
    figax.ZLim = minmax(var);   
    figax.NextPlot = 'replaceChildren';
    figax.Tag = 'PlotFigAxes';  
    %assign data source
    hp = figax.Children;
    hp.ZDataSource = 'vari';  
    %set limites of colorbar    
    hcb = findobj(hfig,'Type','colorbar');
    hcb.LimitsMode = 'manual'; %fix limits of contour bar
    hcb.Limits = figax.ZLim;
    %adjust tick labels and add title
    adjustAxisTicks(pobj,figax);  %adjust tick labels  if defined
    figax.Position = [0.16,0.16,0.65,0.75]; %make space for slider bar
    title(sprintf('%s \nTime = %s',pobj.Title,string(pobj.Data.T(1)))) 
end
%%
function getAnimation(pobj,hp,hfig,t,var,ixM,nint)
    %generate an animation for user selection.
    nrec = length(t);
    figax = gca;
    Mframes(nrec) = struct('cdata',[],'colormap',[]);
    Mframes(1) = getframe(gcf); %NB print function allows more control of 
    for i=2:nrec
        vari = setTimeDependentVariable(pobj,var,i); %#ok<NASGU>
        refreshdata(hp,'caller')
        title(sprintf('%s \nTime = %s',pobj.Title,string(t(i))))
        hline = findobj(figax,'Tag','msection');
        delete(hline)
        hold on
        plot_sections(pobj,figax,var,ixM(i),nint(i),i)
        hold off
        drawnow;                 
        Mframes(i) = getframe(gcf); 
        %NB print function allows more control of resolution 
    end
    idm = size(pobj.ModelMovie,1);            
    pobj.ModelMovie{idm+1,1} = hfig.Number;
    pobj.ModelMovie{idm+1,2} = Mframes;   %save movie to class property
end    
%%
function plot_sections(pobj,figax,var,ixM,nint,i)
    %plot sections at intervals along channel
    line(figax,'XData',repmat(pobj.Data.X(ixM),1,length(pobj.Data.Y)),...
              'YData',pobj.Data.Y,'ZData',var(i,ixM,:),...
              'Color','k','LineStyle',':','LineWidth',0.6,'Tag','msection');
    line(figax,'XData',repmat(pobj.Data.X(ixM+nint),1,length(pobj.Data.Y)),...
              'YData',pobj.Data.Y,'ZData',var(i,ixM+nint,:),...
              'Color','k','LineStyle',':','LineWidth',0.6,'Tag','msection'); 
    line(figax,'XData',repmat(pobj.Data.X(ixM+2*nint),1,length(pobj.Data.Y)),...
              'YData',pobj.Data.Y,'ZData',var(i,ixM+2*nint,:),...
              'Color','k','LineStyle',':','LineWidth',0.6,'Tag','msection');     
          
end
    