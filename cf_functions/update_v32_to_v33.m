function update_v32_to_v33(obj)
%
%-------header-------------------------------------------------------------
% NAME
%   update_v32_to_v33.m 
% PURPOSE
%   update saved models from v3.2.
% USAGE
%   update_v32_to_v33(sm)
% INPUTS
%   obj - instance of model
% RESULTS
%   saved model updated from v3.2.
% NOTES
%   Called in muiModelUI.loadModel when old and new version numbers do not
%   match.
% To use from command line, open ASMITA using:
% >>sm=ChannelForm;     and then call
% >>update_v32_to_v33(sm)
%
% Author: Ian Townend
% CoastalSEA (c) Mar 2024
%--------------------------------------------------------------------------
%

    %add run time plot option to RunProperties and sadjust tab position
    if isfield(obj.Inputs,'RunProperties') && ...
                            ~isfield(obj.Inputs.RunProperties,'isRunPlot')
        obj.Inputs.RunProperties.PropertyLabels = {' Time Step (years)',...
                          ' Number of Time Steps',...
                          ' Output Interval (No. of time steps)', ...
                          ' Start Year',...
                          ' Run time plot (0/1)'};

        obj.Inputs.RunProperties.isRunPlot = false;      
    end

    getdialog('Project updated from v3.2 to v3.3')
end
