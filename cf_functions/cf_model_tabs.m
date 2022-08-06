function cf_model_tabs(obj,src)
%
%-------function help------------------------------------------------------
% NAME
%   cf_model_tabs.m
% PURPOSE
%   generate tabPlot and tabProperties for any GDinterface model
% USAGE
%   cf_model_tabs(obj,mobj,src)
% INPUTS
%   obj - instance of any form model that uses the GDinterface abstract class
%   src - handle to calling tab
% OUTPUT
%   generate plot or properties table on tab
% SEE ALSO
%   used to handle tabPlot calls for the different models which all have
%   the same output data stucture defined by GDinterface
%
% Author: Ian Townend
% CoastalSEA (c) Jan 2022
%--------------------------------------------------------------------------
%
    dst = obj.Data.Form;
    
    if height(dst)>1
        %propmpt user to select timestep
        list = dst.DataTable.Properties.RowNames;
        irec = listdlg('PromptString','Select timestep:',...
                       'Name','Tab plot','SelectionMode','single',...
                       'ListSize',[250,100],'ListString',list);
        if isempty(irec), return; end
    else
        irec = 1;
    end

    switch src.Tag
        case {'Plot','FigButton'}
            tabPlot(obj,src,irec,false);
        case 'FormProps'
            tabProperties(obj,src,irec);
    end
end
%%
function tabPlot(obj,src,irec,iswidth)
    %generate plot for display on Plot tab. To produce stand alone plots 
    %set the flags to true as required. Can also use button on tab
    if nargin<4
        %to use the high and low water widths rather than interpolated values
        iswidth = false; 
    end

    grid = getGrid(obj,irec);
    if isempty(grid.z), return; end
    
    if grid.x(1)<0
        xi = max(grid.x)-grid.x;
    else
        xi = (grid.x);
    end
    yi = grid.y;
    zi = grid.z';
    %can use zlevels to plot contours or recover the widths
    widths = reshape(obj.Data.Plan.DataTable{1,:},length(xi),3);
    hydobj = obj.Data.WaterLevels;
    zo = hydobj.zmt(1);
    zhw = hydobj.zhw(1);
    zlw = hydobj.zlw(1);
    %set button for rotate and standalone figure
    tabcb = @(src,evdat)cf_model_tabs(obj,src);
    ax = tabfigureplot(obj,src,tabcb,true);
    %force ax to be the current axes (Otherwise if standalone figure 
    %generated in one tab and user changes tab, it can plot on the wrong 
    %figure if switching between FormProps and Q-Plot    
    axes(ax); 
    
    %plot form as a surface
    surfc(ax,xi,yi,zi, 'FaceColor','interp','EdgeColor', 'none');

    % caxis(ax,[-15,10]);  %constrain the range of the colormap used by surfc
    % zlim(ax,[-10,10]);   %constrain the z-axis limits
    zlimits = [min(zi,[],'All'),10];   %constrain the z-color range limits
    gd_colormap(zlimits);
    % view(180,90); %plan view
	view(315,30);
    hold on
    if iswidth
        line(ax,'XData',xi,'YData',widths(:,1),'ZData',repmat(zhw,size(xi)),'Color','k','LineStyle','-.','LineWidth',0.75);
        line(ax,'XData',xi,'YData',-widths(:,1),'ZData',repmat(zhw,size(xi)),'Color','k','LineStyle','-.');
        line(ax,'XData',xi,'YData',widths(:,2),'ZData',repmat(zo,size(xi)),'Color','k','LineStyle','--');
        line(ax,'XData',xi,'YData',-widths(:,2),'ZData',repmat(zo,size(xi)),'Color','k','LineStyle','--');
        line(ax,'XData',xi,'YData',widths(:,3),'ZData',repmat(zlw,size(xi)),'Color','k','LineStyle','-.');
        line(ax,'XData',xi,'YData',-widths(:,3),'ZData',repmat(zlw,size(xi)),'Color','k','LineStyle','-.');
    else
        contour3(ax,xi,yi,zi,'LineColor',[0.8,0.8,0.8]);
        c_vals = [zlw,zo,zhw];
        [C,h] = contour3(ax,xi,yi,zi,c_vals,'-.k');
        clabel(C,h); 
        contour3(ax,xi,yi,zi,[zo zo],'--c');
    end
    xlabel('Distance from mouth (m)'); 
    ylabel('Width (m)');  
    zlabel('Elevation (mAD)');
    ttltxt = sprintf('%s (%s)',grid.desc,char(grid.t));
    title(ttltxt,'FontWeight','normal','FontSize',10);
    hold off
    cb = colorbar;
    %c.Limits = [-10,10]; %constrain the range of the colorbar
    cb.Label.String = 'Elevation (mAD)';
    ax.Color = [0.96,0.96,0.96];  %needs to be set after plot
end
%%
function tabProperties(obj,tabsrc,irec)
    %generate table and plot for display on Properties tab
    if ~isfield(obj.Data,'GrossProps') || isempty(obj.Data.GrossProps)
        getdialog('No Form Properties available for selected grid')
        return;
    end
    T = getDSTable(obj.Data.GrossProps,irec,[]);
    %generate table of gross properties
    uitable('Parent',tabsrc,'Data',T.DataTable{:,:},...
            'ColumnName',T.VariableNames,...
            'ColumnWidth',{80},...
            'RowName',T.DataTable.Properties.RowNames,...
            'Units','Normalized','Position',[0.05,0.84,0.9,0.14],...
            'Tag','grossprops');
    %user popup to select a type of plot 
    popup = findobj(tabsrc,'Style','popup','Tag','PlotList');
    if isempty(popup)
        plotlist = {'Hypsommetry','Cross-sections','Thalweg + Plan width',...
            'Along-channel properties','Hydraulic depth',...
            'Prism-Area ratio','Elevation-Area histogram',...
            'a/h and Vs/Vc','Hydraulics','Transgression'};    
        uicontrol('Parent',tabsrc,'Style','text',...
           'Units','Normalized','Position', [0.05 0.79 0.1 0.04],...
           'String', 'Select plot:');  
        popup = uicontrol('Parent',tabsrc,'Style','popup',...
           'String',plotlist,'Tag','PlotList',...
           'Units','Normalized','Position', [0.05 0.74 0.9 0.05],...
           'Callback', @(src,evdat)gd_property_plots(obj,irec,src)); %#ok<NASGU>

        %Create push button to copy data to clipboard
        uicontrol('Parent',tabsrc,'Style','pushbutton',...                    
            'String','>Table','UserData',T,'Tag','CopyButton',...
            'TooltipString','Copy grossproperties table to clipboard',...
            'Units','normalized','Position',[0.75 0.793 0.10 0.044],...                    
            'Callback',@copydata2clip);  
        
        %create push button to create tab plot as a stand alone figure
        uicontrol('Parent',tabsrc,'Style','pushbutton',...                    
            'String','>Figure','Tag','FigButton',...
            'TooltipString','Create plot as stand alone figure',...
            'Units','normalized','Position',[0.86 0.793 0.10 0.044],...                    
            'Callback',@(src,evdat)gd_property_plots(obj,irec,src));  
        
    else
        %update obj in Callbacks and table in UserData
        popup.Callback = @(src,evdat)gd_property_plots(obj,irec,src);
        hb = findobj(tabsrc,'Tag','CopyButton');           
        hb.UserData = T;
        hf = findobj(tabsrc,'Tag','FigButton');
        hf.Callback = @(src,evdat)gd_property_plots(obj,irec,src);
        gd_property_plots(obj,irec,popup); %set plot to current popup selection
    end           
end 