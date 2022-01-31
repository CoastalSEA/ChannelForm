%% Output
% The ChannelForm model generates gridded bathymetry/topography terrain
% models as snapshots or timeseries datasets. Similar externally generated
% data sets can also be imported. The grids and various derived properties
% are held in a set of <matlab:doc('dstable') dstables>.

%% Output Assignment
% The ChannelForm model model and data import classes use the 
% <matlab:doc('gdinterface') GDinterface> abstract class to manage the data. 
% The data are stored in a set of <matlab:doc('dstable') dstables>, which
% are assigned as a <matlab:doc('struct') struct> to the <matlab:doc('muicatalogue') muiCatalogue>
% _Data_ property (typically as part of a UI based on the <matlab:doc('muitoolbox') muitoolbox>). 
% The metatdata for the dstable content are defined in the _*setDSproperties*_ 
% method in <matlab:doc('gdinterface') GDinterface>. 

%% Output Tables
% Imported and model grids are both held in the _'Form'_ <matlab:doc('dstable') dstable>.
% In addition, data related to the properties of the channel or inlet and
% any data on water levels are held in a number of addtional tables, as
% detailed below.
%%
% A number of additional tables support applications that represent
% inlet, estuary, river and valley formations. These include _'Hypsometry'_ for the
% variation of surface area and volume as a function of elevation,
% _'SectionProps'_ holds the variation of the width and cross-sectional area
% along the x-axis derived from the gridded data, _'Plan'_ holds the widths along 
% the x-axis specified in the model, or based on an exponential fit to the
% gridded data, _'WaterLevels'_ holds the variation in water levels at high
% water, mean tide and low water along the x-axis, and _'GrossProps'_ holds
% various summary properties of the form related to length, area, volume,
% rate of convergence, etc. 

%%% Form
% Elevations of grid along with the x and y dimensions, with the 
% following variables:
%%
% * Z - Elevation (mAD)
%%
% _Dimensions_ of time (table rows), X and Y.
%
% _UserData.ishead_ - orientation of x-axis relative to mouth (true if
% minimum x is at the head of the estuary/inlet, default is false).
%
% _UserData.xM_ - distance from x=0 to mouth of estuary/inlet (default is
% 0).

%%% Hypsometry
% Variation of surface area and volume as a function of
% elevation, with the following variables:
%%
% * Volume - Volume (m^3)
% * SurfaceArea - Surface area (m^2)
% * SAfreq - Surface area frequency
%%
% _Dimensions_ of time (table rows), and Z.

%%% SectionProps
% Variation of the width and cross-sectional area
% along the x-axis derived from the gridded data, with the following variables:
%%
% * Whw - High Water Width (m)
% * Wmt - Mean Tide Width (m)
% * Wlw - Low Water Width (m)
% * CSAhw - High Water Cross-sectional Area (m^2)
% * CSAmt - Mean Tide Cross-sectional Area (m^2)
% * CSAlw - Low Water Cross-sectional Area (m^2)
% * Prism - Tidal Prism (m^3)
% * AoP - Prism to CSA ratio (m)
%%
% _Dimensions_ of time (table rows), and X.

%%% Plan
% Widths along the x-axis, specified in the model, or based on an exponential 
% fit to the gridded data, with the following variables:
%%
% * Whw - High Water Width (m)
% * Wmt - Mean Tide Width (m)
% * Wlw - Low Water Width (m)
%%
% _Dimensions_ of time (table rows), and X.

%%% GrossProps 
% Various summary properties of the form related to length, area, volume,
% rate of convergence, etc, with the following variables:
%%
% * Shw - High Water Surface Area (m^2)
% * Slw - Low Water Surface Area (m^2)
% * Vhw - High Water Volume (m^3)
% * Vlw - Low Water Volume (m^3)
% * Pr - Tidal Prism (m^3)
% * Gamma - Dronkers' Gamma (-)
% * Vs - Storage Volume (m^3)
% * Vc - Channel Volume (m^3) 
% * VsoVc - Storage to Channel Volume ratio (-)
% * amp - Tidal amplitude at mouth (m)
% * hyd - Hydraulic depth of channel (m)
% * aoh - Amplitude to Depth ratio (-)
% * Wm - Mean Tide Width  at mouth (m)
% * Am - Mean Tide Cross-sectional Area at mouth (m^2)
% * PoA - Prism to Cross-sectional Area ratio at mouth (m)
% * Lw - Width convergence length
% * La - Cross-sectional Area convergence length
%%
% _Dimensions_ of time (table rows).

%%% WaterLevels
% Variation in water levels at high water, mean tide and low water along
% the x-axis, with the following variables:
%%
% * zhw - High Water level (mAD), where mAD = metres above Datum
% * zmt - Mean Tide level (mAD)
% * zlw - Low Water level (mAD)
%%
% _Dimensions_ of time (table rows), and X.

%% See Also
% The <matlab:open_manual manual> provides further details of setup and 
% configuration of the model.
%%
% Townend I H, Zhou Z, Guo L and Coco G, 2021, A morphological investigation 
% of marine transgression in estuaries. Earth Surf. Process. Landf., 46, 
% 626–641, https://doi.org/10.1002/esp.5050.

