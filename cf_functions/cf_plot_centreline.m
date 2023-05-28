function cf_plot_centreline()    
%
%-------function help------------------------------------------------------
% NAME
%   cf_plot_centreline.m
% PURPOSE
%   plot 2 or3D line e.g., for valley/channel thalweg
% USAGE
%   cf_plot_centreline()
% INPUTS
%   user prompted to select a text file with x,y,(z) coordinates
% OUTPUT
%   plot of 2 or 3D line
% SEE ALSO
%   used in ChannelForm
%
% Author: Ian Townend
% CoastalSEA (c) June 2023
%--------------------------------------------------------------------------
%
        %get data from file
        userprompt = 'Select data file(s)>';
        [fname, path]=uigetfile('*.txt',userprompt,'MultiSelect','off');
        if isequal(fname,0), return; end
        cline = load([path,fname]);
        %initialise figure
        hf = figure('Name','CentreLine','Tag','PlotFig');
        ax = axes(hf);
        %plots for 2 or 3 dimensions
        if size(cline,2)==2
            plot(ax,cline(:,1),cline(:,2),'DisplayName','Centreline')            
            %axis(ax,'equal')
        else
            plot3(ax,cline(:,1),cline(:,2),cline(:,3))
            zlabel('Elevation (m)')
        end
        xlabel('X-distance (m)')
        ylabel('Y-distance (m)')
        legend

        %code to generate migrating meander figure
        % hold on
        % ylim = ax.YLim;
        % answer = inputdlg('Migration distance','Meander',1,{'10000'});
        % xdist = str2double(answer{1});
        % plot(ax,[xdist,xdist],ylim,'-g','DisplayName','Eroded shore')
        % xx = cline(:,1)+xdist;
        % plot(ax,xx,cline(:,2),'--b','DisplayName','Migrating meander')   
        % idx = cline(:,1)>=xdist;
        % xx = cline(idx,1); yy =cline(idx,2);
        % plot(ax,xx,yy,'--r','DisplayName','Stationary meander')   
        % hold off
        % legend
end

