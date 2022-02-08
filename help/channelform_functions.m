%% ChannelForm functions
% Summary of functions available in the _cf_functions_ folder. Use the Matlab(TM)
% help function in the command window to get further details of each
% function.

%% Model functions
% Functions that implement various form models

%%
% * *cf_changeplot* - plot volume difference over the model run and the
% variation of sediment exchange with transgression distance.
% * *cf_cst2grid* - interpolate the CST water levels onto the model grid.
% * *cf_exp_models* - construct idealised channel form using exponential 
% functions in y to determine width and CKFA cross-section to determine 
% z at each x interval.
% * *cf_model_tabs* - generate tabPlot and tabProperties for any
% GDinterface model.
% * *cf_offset_wls* - translate water levels landwards when mouth is some distance from x=0 
% and pad the water level vectors with values at the mouth.
% * *cf_plan_form* - compute planform variation along the x-axis at
% specified planar levels.
% * *cf_pow_model* - function to compute 3D form of a creek or tidal channel
%   using power laws to define width and hydraulic depth variations.
% * *cf_property_plots* - plots displayed on Proprety tab in ChannelForm
% model.
% * *cf_sectionplot* - plot a set of along-channel cross-sections and the channel centre-line
% (thalweg) for the initial and final grid.
% * *cf_set_hydroprops* - set water levels for form model using either the surface defined by the
% cst_model or high water at the mouth and a reducing tidal amplitude.
% * *cf_summarygridplot* - plot initial grid, new grid and difference.
% * *cf_summarytransplot* - composite set of plots to show 
% results form Transgression table.
% * *cf_valley_model* -  construct idealised channel form using 3D
% exponential form model.
% * *cf_writetable* - write GrossProps and Transgression output tables to a file.
% * *ckfa_3D_form* - construct 3D form of an idealised channel based on the CKFA model.
% * *ckfa_dimensions* - utility to get the CKFA channel form properties:
% hydraulic depth, hm, width convergence length, Lw, critical tidal velocity amplitude, Ucr.
% * *ckfa_form_model* - construct idealised channel form using 3D CKFA model
% * *ckfa_form_solver* - find the hydraulic depth, e-folding length for width and the peak tidal
%   amplitude for defined tidal range and sediment properties.
% * *ckfa_wave_form* - obtain depth and width of wave formed profile at high and low water
% * *ckfa_wave_profile* - depth and width of wave formed profile
% d=dw*(1-y/yw)^2/3.
% * *get_river_profile* - get river profile using regime theoretical
% profile (Cao & Knight,1996).
% * *get_river_regime* - get river regime properties using regime
% theoretical profile (Cao & Knight,1996).
% * *get_sed_flux* - single element ASMITA model to compute amount of
% sedimentation.


%% Grid functions
% Functions used to manipulate cartesian grids can be found in the
% _muiAppGridFcns_ folder and includes the following.
%%
% * *gd_channel_hypsometry*
% - compute area and volume hypsometry from gridded elevation data.
% * *gd_colormap*
% - check if Mapping toolbox is installed to use land/sea colormap, or call
% _cmap_selection_ (see <matlab:doc('psfunctions') Plotting and statistical functions> 
% in the <matlab:doc('muitoolbox') muitoolbox>) if not available.
% * *gd_gross_properties*
% - compute the gross properties of a gridded bathymetry.
% * *gd_section_properties*
% - compute the width, cross-sectional area and prism along channel.
% * *getconvergencelength* -  least squares fit using fminsearch to
% find the convergence length of a channel from a distance-width xy data set 
% * *InterX* - intersection of two curves. MATLAB Central File Exchange, 
% Author: NS, 2010, https://www.mathworks.com/matlabcentral/fileexchange/22441-curve-intersections.
% * *xy2sn* - Juernjakob Dugge, 2015, jdugge/xy2sn, https://github.com/jdugge/xy2sn 
% * *sn2xy* - Juernjakob Dugge, 2015, jdugge/xy2sn, https://github.com/jdugge/xy2sn 
% * *arclength* - John D'Errico, 2012, https://www.mathworks.com/matlabcentral/fileexchange/34871-arclength 
% * *distance2curve* - John D'Errico, 2013,
% http://www.mathworks.de/matlabcentral/fileexchange/34869-distance2curve
% * *interparc* - John D'Errico, 2012, https://www.mathworks.com/matlabcentral/fileexchange/34874-interparc
%%
% Further details can be found in <matlab:doc('grid_class_fcns') Grid classes and functions>, 
% or by using the >>help _function name_ command in the Command Window.

%% Generic functions
% Other functions used in ChannelForm, include:
%%
% * *river_regime* 
% - function to compute the width and hydraulic depth of a river section given the
%   discharge and energy slope.
% * *sealevelrise.m*
% - function to compute sea level rise using linear or exponential rate of
% change.
% * *settling_velocity* - calculate the settling velocity using Soulsby
% equation.
% * *ucrit* -  compute the critical flow velocity for a given critical shear stress and
% wave conditions in the wave-current case.


%% See Also
% The <matlab:doc('asmita_classes') additional classes> used in Asmita and
% the <matlab:open_manual manual>, which provides further details of setup and 
% configuration of the model.
