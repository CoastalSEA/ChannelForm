Cases are those presented in 
Townend, I. H., Z. Zhou, L. Guo, and G. Coco (2020), A morphological investigation of marine transgression in estuaries, Earth Surf. Process. Landf., 46, 626–641, doi:https://doi.org/10.1002/esp.5050.
Case 1 is a short steep valley (e.g. the Dyffi)
Case 2 is a long steep valley (e.g. the Humber)
Case 3 is a long flat valley (e.g. the Thames)

The parameter settings in the *.mat file provided are as detailed in the paper, with the exception of the equilbrium concentration. Changes in the definition of the hydraulic depth and equilibrium coefficients that are inputs the sediment flux function mean that the calibration of the model has changed. The following table illusrates the changes between the original and current version of the model. 


					cE	Case 1	cE		Case 2	cE	Case 3
original kinematic	-	1119	-		4331	-	4540
constant kinematic	-	1122	-		4384	-	4845
cst kinematic		-	1361	-		3812	-	4044

original dynamic	0.5	992		2		3363	1		3324
constant dynamic	0.5	734		0.2		1635	0.1		18
								0.05	2884	0.01	554
					0	924		0		3545	0		3930
cst dynamic			0.5	884		0.05	1249	0.05	-2217
								0.02	2016	0.01	1006
								0.02	2290	0.005	1741
					0	1005	0		2592	0		3152
