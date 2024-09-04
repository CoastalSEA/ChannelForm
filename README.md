# ChannelForm
ChannelForm is a Matlab(TM) App to model the interaction of estuary channels and the surrounding flood plain and valley landscape.

## Licence
The code is provided as Open Source code (issued under a BSD 3-clause License).

## Requirements
ChannelForm is written in Matlab(TM) and requires v2016b, or later. In addition, ChannelForm requires the _dstoolbox_ and the _muitoolbox_. To include hydrodynamic water surfaces, the CSTmodel App also needs to be installed.

## Background
ChannelForm is a Matlab(TM) App that provides a framework for studying estuary and inlet morphological landforms in the context of the surrounding landscape (valley and flood plain). A set of tools allow digital terrain models (DTMs) to be created from models or imported. The models/data can define the channel form and the surrounding or antecedent landscape. A set of grid manipulation tools allow the grids to be manipulated and combined to generate composite DTMs of the channel in the landscape. Several idealised models are provided that use some combination of the plan form and cross-section to generate a 3D form, or generate the form using parameters that are not a function of the estuary itself (i.e. they are exogenous). In both cases, the water level surface can be defined using plane surfaces, or derived from a quasi-analytical hydrodynamic model that accounts for tides and river flows. Using channel and valley forms, marine transgression as a function of sea level rise can also be examined. This identifies the horizontal transgression of the initial form to achieve a net mass balance when any exchange with the external environment is accounted for. The output includes the gross properties of the channel, various well-known parameter relationships, and a time series of DTMs

## Bibliography
Townend I H and Pethick J, 2002, Estuarine flooding and managed retreat. Phil.Trans.R.Soc.Lond.A, 360 (1796), 1477-1495.

Townend I H, 2010, An exploration of equilibrium in Venice Lagoon using an idealised form model. Continental Shelf Research, 30 (8), 984-999.

Townend I H, 2012, The estimation of estuary dimensions using a simplified form model and the exogenous controls. Earth Surface Processes and Landforms, 37, 1573-1583.

Townend I H, Zhou Z, Guo L and Coco G, 2021, A morphological investigation of marine transgression in estuaries. Earth Surf. Process. 
Landf., 46, 626–641, https://doi.org/10.1002/esp.5050.

## Acknowledgements
The ChannelForm App makes use of functions from the Matlab(TM) File Exchange Form, including:
##
* InterX - NS (2010). https://www.mathworks.com/matlabcentral/fileexchange/22441-curve-intersections 
* xy2sn and sn2xy - Juernjakob Dugge (2015). jdugge/xy2sn, https://github.com/jdugge/xy2sn  

## Manual
The ChannelForm manual in the app/doc folder provides further details of setup and configuration of the model. The files for the example use case can be found in the app/example folder. 

## See Also
The repositories for _dstoolbox_, _muitoolbox_ and _muiAppLIb_.
