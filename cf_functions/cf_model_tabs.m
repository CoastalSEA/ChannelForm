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
    switch src.Tag
        case 'Plot'
            tabPlot(obj,src);
        case 'FormProps'
            tabProperties(obj,src);
    end
end
%%
function tabPlot(obj,src,iswidth,isfig)
    %generate plot for display on Plot tab. To produce stand alone plots 
    %set the flags to true as required. Can also use button on tab
    if nargin<4
        %to create stand alone figure set isfig=true   
        isfig = false; 
        %to use the high and low water widths rather than interpolated values
        iswidth = false; 
    elseif nargin<5
        isfig = false; 
    end

    selecteddata = obj.Data.Form;
    xi = selecteddata.Dimensions.X;
    if (xi(1)<0)
        xi = max(xi)-xi;
    else
        xi = flipud(xi);
    end
    yi = selecteddata.Dimensions.Y;
    zi = squeeze((selecteddata.Z))';
    %can use zlevels to plot contours or recover the widths
    widths = reshape(selecteddata.DataTable{1,2:end},length(xi),3);
    hydobj = obj.RunParam.CF_HydroData;
    zo = hydobj.zmt(1);
    zhw = hydobj.zhw(1);
    zlw = hydobj.zlw(1);

    %surface plot of 3D form
    ht = findobj(src,'Type','axes');
    delete(ht);
    %
    if isfig
        hf = figure;
        ax = axes('Parent',hf,'Tag','Surface');
    else
        ax = axes('Parent',src,'Tag','Surface');
    end
    
    surfc(ax,xi,yi,zi, 'FaceColor','interp','EdgeColor', 'none');

    % caxis(ax,[-15,10]);  %constrain the range of the colormap used by surfc
    % zlim(ax,[-10,10]);   %constrain the z-axis limits
    if license('test','MAP_Toolbox')
        cmapsea = [0,0,0.2;  0,0,1;  0,0.45,0.74;  0.30,0.75,0.93; 0.1,1,1];
        cmapland = [0.95,0.95,0.0;  0.1,0.7,0.2; 0,0.4,0.2; 0.8,0.9,0.7;  0.4,0.2,0];
        demcmap([min(zi,[],'All'),10],128,cmapsea,cmapland) %mapping toolbox
    else
        cmap = cmap_selection;
        if isempty(cmap), cmap = 'parula'; end
        colormap(cmap)        
    end

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
    title(selecteddata.Description,'FontWeight','normal','FontSize',10);
    hold off
    cb = colorbar;
    %c.Limits = [-10,10]; %constrain the range of the colorbar
    cb.Label.String = 'Elevation (mAD)';
    ax.Color = [0.96,0.96,0.96];  %needs to be set after plot
    
    hb = findobj(src,'Style','pushbutton');
    delete(hb) %delete button so that new axes is assigned to callback
    uicontrol('Parent',src,...  %callback button
        'Style','pushbutton',...
        'String', 'Rotate off',...
        'Units','normalized', ...
        'Position', [0.02,0.92,0.1,0.05],...
        'TooltipString','Turn OFF when finished, otherwise tabs do not work',...
        'Callback',@(src,evtdat)rotatebutton(ax,src,evtdat));  
end
%%
function tabProperties(obj,tabsrc)
    %generate table and plot for display on Properties tab
    T = obj.Data.GrossProps;
    %generate table of gross properties
    uitable('Parent',tabsrc,'Data',T.DataTable{:,:},...
            'ColumnName',T.VariableNames,...
            'ColumnWidth',{80},...
            'RowName',T.DataTable.Properties.RowNames,...
            'Units','Normalized','Position',[0.05,0.84,0.9,0.14],...
            'Tag','grossprops');
    %user popup to select a type of plot 
    popup = findobj(tabsrc,'Style','popup');
    if isempty(popup)
        plotlist = {'Hypsommetry','Cross-sections','Centre-line',...
                'Form width','Elevation-Area histogram','Hydraulic depth',...
                'Area-Prism ratio','Cross-sectional area','Prism',...
                'a/h and Vs/Vc','Hydraulics','Transgression'};    
        uicontrol('Parent',tabsrc,'Style','text',...
           'Units','Normalized','Position', [0.05 0.79 0.1 0.04],...
           'String', 'Select plot:');  
        popup = uicontrol('Parent',tabsrc,'Style','popup',...
           'String',plotlist,'Tag','PlotList',...
           'Units','Normalized','Position', [0.05 0.74 0.9 0.05],...
           'Callback', @(src,evdat)cf_property_plots(obj,src));

        %Create push button to copy data to clipboard
        uicontrol('Parent',tabsrc,'Style','pushbutton',...                    
            'String','>Table','UserData',T,'Tag','CopyButton',...
            'TooltipString','Copy grossproperties table to clipboard',...
            'Units','normalized','Position',[0.75 0.793 0.10 0.044],...                    
            'Callback',@copydata2clip);  
        
        %create push button to create tab plot as a stand alone figure
        uicontrol('Parent',tabsrc,'Style','pushbutton',...                    
            'String','>Figure','UserData',popup.Value,'Tag','FigButton',...
            'TooltipString','Create plot as stand alone figure',...
            'Units','normalized','Position',[0.86 0.793 0.10 0.044],...                    
            'Callback',@(src,evdat)cf_property_plots(obj,src));  
        
    else
        %update obj in Callbacks and table in UserData
        popup.Callback = @(src,evdat)cf_property_plots(obj,src);
        hb = findobj(tabsrc,'Tag','CopyButton');           
        hb.UserData = T;
        hf = findobj(tabsrc,'Tag','FigButton');
        hf.Callback = @(src,evdat)cf_property_plots(obj,src);
        cf_property_plots(obj,popup); %set plot to current popup selection
    end           
end 