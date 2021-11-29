%% SedTools
% Toolbox for sediment related analyses

%% Licence
% The code is provided as Open Source code (issued under a GNU General 
% Public License).

%% Requirements
% ModelUI is written in Matlab(TM) and requires v2016b, or later. In addition, 
% ModelUI requires both the <matlab:doc('dstoolbox') dstoolbox> and the 
% <matlab:doc('muitoolbox') muitoolbox>

%% Background

%% SedTools classes
% *ModelUI* - defines the behaviour of the main UI.

%% SedTools functions
%  Currently provides:
%      1) tool to analyse settling column data
% 	
%  Quick Guide to setting up a new analysis:
% 
%  File>New 
%      To create a new project space.
% 	
%  Setup>Setup>Settling Parameters
%      Define input parameters for settling column analysis. To see the parameters that are currently
%      set for the model use the Setup tab.
% 	
%  File>Save as
%      Save project to a *.mat file.
% 	
%  Run>Run Settling Data 
%      Runs the model and prompts user for a description of the scenario. The results can be viewed
%      on the Plot tab.
% 	
%  Completed scenarios are listed based on the user descriptions on the Scenarios tab.
% 
%  Plot>Plot Menu
%      Select data and generate plots of results
%      Use Plot tab to generate customised analysis plot and see statistical output
%% Manual
% The <matlab:open_manual manual> provides further details of setup and 
% configuration of the model. The files for the example use case can be found in
% the example folder <matlab:example_folder here>. 

%% See Also
% <matlab:doc('muitoolbox') muitoolbox>, <matlab:doc('dstoolbox') dstoolbox>.
	