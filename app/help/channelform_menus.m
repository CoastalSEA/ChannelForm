%% Menu Options
% Summary of the options available for each drop down menu.

%% File
% * *New*: clears any existing model (prompting to save if not already saved) and a popup dialog box prompts for Project name and Date (default is current date). 
% * *Open*: existing Asmita models are saved as *.mat files. User selects a model from dialog box.
% * *Save*: save a file that has already been saved.
% * *Save as*: save a file with a new or different name.
% * *Exit*: exit the program. The close window button has the same effect.

%% Tools
% * *Refresh*: updates Cases tab.
% * *Clear all > Project*: deletes the current project, including all Setup data and all Cases.
% * *Clear all > Figures*: deletes all results plot figures (useful if a large number of plots have been produced).
% * *Clear all > Cases*: deletes all Cases listed on the Cases tab but does not affect the model setup.

%% Project
% * *Project Info*: edit the Project name and Date
% * *Cases > Edit Description*: user selects a Case to edit the Case description.
% * *Cases > Edit DS properties*: initialises the  UI for editing Data Set properties (DSproperties).
% * *Cases > Edit Data Set*: initialises the Edit Data UI for editing data sets.
% * *Cases > Save*: user selects a data set to be saved from a list box of Cases and the is then prompted to name the file. The data are written to an Excel spreadsheet. 
% * *Cases > Delete*: user selects Case(s) to be deleted from a list box of Cases and results are then deleted (model setup is not changed).
% * *Cases > Reload*: user selects a Case to reload as the current parameter settings.
% * *Cases > View settings*: user selects a Case to display a table listing the parameters used for the selected Case. 
% * *Export/Import > Export*: user selects a Case class instance to export as a mat file.
% * *Export/Import > Import*: user selects an exported Case class instance (mat file) to be loaded.
%%
% <html>
% <table border=1><tr><td><u>Note</u>: to export the data from a Case for use in another application 
% (eg text file, Excel, etc), use the <b>Project>Cases>Edit Data Set</b> option 
% to make a selection and then use the <i>Copy to Clipboard</i> button to paste 
% the selection to the clipboard.
% </td></tr></table>
% </html>

%% Setup > Form Parameters
% * *Exp Form Parameters*: UI to input the parameters needed for the
% exponential form model. 
% * *Power Form Parameters*: UI to input the parameters needed for the
% power form model. 
% * *Valley Parameters*: UI to input the parameters needed for the
% valley form model. 
% * *Shore Parameters*: UI to input parameters needed for the shore profile,
% when adding a shoreline strip to Channel or Valley models.
%%
% <html>
% <table border=1><tr><td><u>Note</u>: the CKFA model input parameters are 
% derived from the water level, hydrodynamic and sediment parameter settings.
% </td></tr></table>
% </html>

%% Setup > System Parameters
% * *Hydraulic Parameters > Tidal Forcing*: UI to input water level parameters including 
% sea level rise and any cyclic variations in tidal range. 
% * *Hydraulic Parameters > Hydraulic Model*: UI to input parameters needed
% for the CST hydraulic model.
% * *Sediment Parameters*: UI to input parameters needed to determine
% external sediment exchange.
% * *Transgression Parameters*: UI to set inclusion/exclusion of various
% constraints, computation parameters and open coast erosion settings.
% * *Morphological Modifications*: UI to input the x,y,z definitions of
% the, one or more, boxout forms that can be imposed on model forms. 

%% Setup > Run Parameters
% * *Grid Parameters*: define dimensions of default grid
% * *Run Time Parameters*: settings for time step and duration of model run

%% Setup > Import Grid
% Load x,y,z data from a file. To create a new instance (e.g. for a
% different location or data source) use Load. To add data to an
% existing data set, use Add.
%%
% * *Load*: prompts user for file to be loaded. Once files have been read, user is prompted for a description (working title) for the data set. 
% * *Add*: prompts user for file to be loaded (only one file at a time can be added). Only files with the same dimensions as the existing data set can be used to Add data to a data record.
% * *Delete*: delete a grid from an existing Case table (ie a row).

%% Setup > Grid Tools
% * *Grid Tools > Translate Grid*: interactively translate grid x-y
% coordinates;
% * *Grid Tools > Rotate Grid*: interactively flip or rotate grid;   
% * *Grid Tools > Re-Grid*: regrid a gridded dataset to match another grid or to user
% specified dimensions;
% * *Grid Tools > Sub-Grid*: interactively define a subgrid and save grid as a new Case;               
% * *Grid Tools > Combine Grids*: superimpose one grid on another based on maximum
% or minimum set of values;
% * *Grid Tools > Add Surface*: add horizontal surface to an extisting
% grid;
% * *Grid Tools > To curvilinear*: map grid from cartesian to curvilinear coordinates; 
% * *Grid Tools > From curvilinear*: map grid from curvilinear to cartesian
% coordinates;
% * *Grid Tools > Display Dimensions*: display a table with the dimensions
% of a selected grid;
% * *Grid Tools > Difference Plot*: generate a plot of the difference
% between two grids;
% * *Grid Tools > Plot Sections*: interactively define sections and plot
% them on a figure;
% * *Grid Tools > Digitise Line*: interactively digitise a line (with
% option to add z values) using selected grid as base map;
% * *Grid Tools > Export xyz Grid*: select a Case and export grid as xyz
% tuples;

%% Setup (other)
% * *Add Properties*: add form properties to a gridded data
% set.
% * *Delete Properties*: delete ALL property tables associated with a 
% selected gridded data set.
% * *Edit Inlet Definition*: provides access to the definition of the 
% position of the head of the inlet, distance from the x-axis origin to 
% the mouth and any definition of the channel centre-line (if used).
% * *Model Constants*: various constants are defined for use in models, such as the acceleration due to gravity, viscosity and density of sea water, and density of sediment. 

%% Utilities 
% * *Hydraulic Model*: runs the CSTmodel (if CSTmodel App is installed).
% * *Add Form to Valley*: combines channel and valley grids.
% * *Utilities>Add Meander*: adds a user defined meander to an idealised 
% channel form model. After selecting a model, the user is prompted to 
% load a text file containing the x-y co-ordinates of the channel 
% centreline. These are used to map the straight channel into curvilinear 
% co-ordinates that follow the centreline and then interpolate the results 
% onto the cartesian model grid.
% * *Utilities>Add Shoreline>Model Shore*: add a shoreline strip to a 
% channel or valley model. This extends the grid seawards by the distance 
% to the closure depth specified in Setup>Form Parameters>Shore Parameters. 
% This comprises a beach profile using Dean’s equilibrium model on the open 
% coast and the defined offshore bed slope across the mouth of the channel.
% * *Utilities>Add Shoreline>Extrapolate Shore*: add a shoreline strip to a 
% grid (imported or model). The user is prompted to define the width of the 
% extension, the seaward or offshore depth and an exponent. An exponent 
% of 1, gives linear extrapolation between the grid shore values and the 
% defined offshore value, a value greater than one gives a shore that is 
% concave down and a value greater than one gives a shore that is convex up.
% * *Utilities>Add Modifications*: adds defined modifications to a grid.
% * *Utilities>Add Thalweg to Valley*: uses an xyz definition of a channel 
% thalweg, to create a valley base in an existing model or imported grid. 
% This is used to modify an existing valley grid, to extend the valley 
% form down to the pre-Holocene surface (i.e. below the existing channel 
% bed) using the defined the thalweg data. User is prompted a thalweg file. 
% This has a header of %f %f %f, followed by rows of x, y, z values that 
% define the pre-Holocene valley bed. The upper and lower cut-offs must 
% then be defined.  These are used to mask out the channel in the valley 
% grid in order to introduce the new channel base (thalweg). Typically 
% values just above high water and below the deepest depth in the valley
% grid channel should suffice.
% * *River Dimensions*: calculates the river dimensions based on the
% current input parameters.
% * *Valley Thalweg*: plots the valley centre-line along the x-axis.
% * *Area of Flood Plain*: calculates the area of the flood plain.
% * *CKFA Channel Dimensions*: calculates the governing dimensions that
% would be used in the CKFA model based on the current input parameters.
% * *Morphological Timescale*: calculate the morphological timescale using
% the ASMITA model parameterisation and the current input parameters.
% * *Channel-Valley Sub-Plot*: create a plot showing the channel, valley
% and combined form as 3 sub-plots.

%% Run
% * *Exponential form model*: run the exponential model. Prompts to select
% type of water surface and additional form selection properties.
% * *Power form model*:  run the power model. Prompts to select
% type of water surface
% * *CKFA form model*:  run the CKFA model. Prompts to select
% type of water surface.
% * *Valley form model*: run the Valley model
% * *Transgression model*: run the transgression model. Prompts to select
% source channel and valley models to use.
% * *Derive Output*: initialises the Derive Output UI to select and define manipulations of the data or call external functions and load the result as new data set.

%% Analysis
% * *Plots*: initialises the Plot UI to select variables and produce various types of plot. The user selects the Case, Dataset and Variable to used, along with the Plot type and any Scaling to be applied from a series of drop down lists, 
% * *Statistics*: initialiss the Statistics UI to select data and run a range of standard statistical methods.
% * *Trangression Plots*: summary plots to present the results of the
% transgression model or export the GrossProps and Transgression output
% tables to a file (see <matlab:doc('channelform_output') Output help> for details).
% * *Transgression Plots>Grid Plot*: plot initial grid, new grid and
% difference.
% * *Transgression Plots>Change Plot*: plot volume difference over the 
% model run and the variation of sediment exchange with transgression
% distance.
% * *Transgression Plots>Section Plot*: plot a set of along-channel cross-sections and the channel centre-line
% (thalweg) for the initial and final grid.
% * *Transgression Plots>Transgression Plot*: composite set of plots to show 
% results form Transgression table.
% * *Transgression Plots>Export Table*: write GrossProps and Transgression 
% output tables from the Transgression model to a file.

%% Help
% * *Help*: access the online documentation for CoastalTools.

%% See Also
% The <matlab:cfm_open_manual manual> provides further details of setup and 
% configuration of the model.