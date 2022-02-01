function cf_summarytransplot(obj)
%
%-------function help------------------------------------------------------
% NAME
%   cf_summarytransplot.m
% PURPOSE
%   composite set of plots to show results form Transgression table
% USAGE
%   cf_summarytransplot(obj)
% INPUTS
%   obj - instance of CF_TransModel class
% OUTPUT
%   summary plot
% SEE ALSO
%   cf_property_plots - Transgression option
%
% Author: Ian Townend
% CoastalSEA (c) Jan 2022
%--------------------------------------------------------------------------
%
    %composite set of plots to show results form Transgression table
    dst = obj.Data.Transgression;
    slr = dst.dSLR(end);
    t = dst.RowNames;
    vardesc = dst.VariableDescriptions;
    %varnames: 'delX','estdX','cstdX','dSLR','Lt','FPA','sedVol','vdiffx'
    varnames =  dst.VariableNames;
    figure('Name','Trans Plot','Units','normalized','Tag','PlotFig'); 
    %
    subplot(2,2,1)
    %plot transgression distances idx=1-3
    plot(t,dst.(varnames{1}),'DisplayName',vardesc{1});
    hold on 
    for i=2:3
        plot(t,dst.(varnames{i}),'DisplayName',vardesc{i});
    end
    hold off
    plot_labels(dst,3)
    %
    subplot(2,2,2)
    % 
    yyaxis left
    plot(t,dst.(varnames{4}),'DisplayName',vardesc{4});
    ylabel(dst.VariableLabels{4})
    yyaxis right
    plot(t,dst.(varnames{7}),'DisplayName',vardesc{7});
    plot_labels(dst,7)
    %
    subplot(2,2,3)
    %
    yyaxis left
    plot(t,dst.(varnames{5}),'DisplayName',vardesc{5});
    ylabel(dst.VariableLabels{5})
    yyaxis right
    plot(t,dst.(varnames{6}),'DisplayName',vardesc{6});
    plot_labels(dst,6)
    
    %
    subplot(2,2,4)
    %plot Volume change for [0,dx/2,dx,3dx/2]
    vardescs = {[varnames{8},'-0'],[varnames{8},'-dx/2'],...
                [varnames{8},'-dx'],[varnames{8},'-3dx/2']};
    hp = plot(t,dst.(varnames{8}));
    for i=1:4
        hp(i).DisplayName = vardescs{i};
    end
    plot_labels(dst,8)
    
    sgtitle(sprintf('%s, slr=%0.1g',dst.Description,slr),'FontSize',12);
end
%%
function plot_labels(dst,idx)
    xlabel('Time (years)')
    ylabel(dst.VariableLabels{idx})
    ytickformat('%3.2g')
    legend
end