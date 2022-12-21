%% ChannelForm classes
% The ChannelForm App is built using the <matlab:doc('muitoolbox') muitoolbox>
% and a number of model specific classes.

%% Model Classes
% The ChannelForm App uses the following classes:
%%
% * *ChannelForm* - main UI for the ChannelForm interface, which implements the 
%   <matlab:doc('muimodelui') muiModelUI> abstract class to define main menus.

%%
% Data input classes:
%%
% * *CF_ExpData* - handles data input for exponential form parameters.
% * *CF_PowerData* - handles data input for power form parameters. 
% * *CF_ValleyData* - handles data input for valley form parameters.
% * *CF_HydroData* - handles data input of the hydraulic parameters used by CSTmodel App.
% * *CF_SediData* - handles data input for sediment flux parameters.
% * *CF_ShoreData* – handles data input for shoreline form parameters.
% * *CF_ModsData* - handles data input for imposed modifications to the form. 
% * *CF_TransData* - handles data input for channel transgression
% parameter settings.
% * *RunProperties* - define the model run time parameters.
% * *WaterLevels* - setup and access to water level definitions, tidal 
% parameters and changes over time.
%%
% Model classes:
%%
% * *CF_FormModel* - methods needed to create and manipulate exponential,
% power, or CKFA form models.
% * *CF_ValleyModel* - methods needed to create and manipulate an idealised
% valley form.
% * *CF_TransModel* - transgression model to compute the movement of a
% channel within a valley in response to sea level rise using a simple
% kinetic model.


%% Grid Classes
% Classes used to manipulate cartesian grids can be found in the
% _muiAppGridFcns_ folder and include the following:
%%
% * *GD_GridProps*: class inherits <matlab:doc('muipropertyui') muiPropertyUI> 
% abstract class, providing an interface to define the extent and intervals
% of a cartesian grid. 
% * *GDinterface*: an abstract class to support classes that need additional
% functionality to handle grids. The class inherits <matlab:doc('muidataset') muiDataSet> 
% and together they provide an extensive set of methods to handle datasets
% of various types (eg from models or imported files). 
% * *GD_ImportData*: class inherits <matlab:doc('gdinterface') GDinterface> abstract class (see above)
% to load xyz data from a file.
%%
% Further details can be found in <matlab:doc('grid_class_fcns') Grid classes and functions>
% 

%% See Also 
% <matlab:doc('channelform_functions') ChannelForm functions> used in 
% ChannelForm App and the <matlab:cfm_open_manual manual>, which provides further details 
% of setup and configuration of the model.
