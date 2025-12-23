# The health and economic repercussions of declining MMR coverage in the United States

Chad R. Wells <sup>1</sup>, Abhishek Pandey <sup>1</sup>, Yang Ye <sup>1</sup>, Carolyn Bawden <sup>1</sup>, Rebecca Giglio <sup>2</sup>,  Charlene Wong <sup>2</sup>, Velda Wang <sup>2</sup>, Chelsea Cipriano <sup>2</sup>, Lamia Ayaz <sup>1</sup>, Seyed M. Moghadas <sup>3</sup>, Meagan C. Fitzpatrick <sup>4</sup>, Dave A. Chokshi <sup>2</sup>, Burton H. Singer <sup>5</sup>, Alison P. Galvani <sup>1</sup>

<sup>1</sup> Center for Infectious Disease Modeling and Analysis (CIDMA), Yale School of Public Health, New Haven, Connecticut, 06510, USA<br /> 
<sup>2</sup> The Common Health Coalition, New York, NY, 10016<br /> 
<sup>3</sup> Agent-Based Modelling Laboratory, York University, Toronto, Ontario, M3J 1P3, Canada<br /> 
<sup>4</sup> Center for Vaccine Development and Global Health, University of Maryland School of Medicine, Baltimore, Maryland, 21201, USA<br /> 
<sup>5</sup> Emerging Pathogens Institute, University of Florida, Gainesville, Florida, 32610, USA<br /> 

Copyright 2025, Chad Wells et al. All rights reserved. Released under the GNU GENERAL PUBLIC LICENSE v3.

The MATLAB code provided here will run the fitting and analysis for the modelling portion of the manuscript.

## Vaccination scripts and functions
Age_Adjustment_Factor_25_plus - Computes the adjustment factors for the age groups 25 years and older
Age_Adjustment_Factor - Computes the adjustment factors for the age groups 0-24 years of age
Compute_Age_Group_Reduction - Computes the level of reduction in the 0-4; 5-9; 10-14 age groups for a specified annual reduction
Decline_Adjustment_Factor - computes the adjustment factors required to get the specified level of reduction
Estimated_NE_Health_District_Vaccine_Uptake - computes the vaccine coverage among the Nebraska Health Districts
Estimated_State_Vaccine_Uptake - computes the vaccine coverage among the states
Fit_Model - Runs the optimization for the specified training dataset for the vaccination model
Load_Data_Adjustment - loads the vaccination data required for computing the necessary adjustments to get the county and state level coverage
Load_Data - Loads the data for the vaccination model and the data used in the training of the model
Objective_Adjustment_Coverage_County - the objective function used for paramterizing the adjustments factor at the county level
Objective_Adjustment_Coverage_State - the objective function used for paramterizing the adjustments factor at the state level
Objective_Function_Coverage - The obejective fucntion used in the training of the vaccination model
Parameters - returns the parameters for the specified vaccination model
Refine_Fit_Model - Refines the fitting of the vaccination model; This is run after "Fit_Model"
Vaccination_Coverage_Adjusted - Computes the adjusted vaccination coverage
Vaccination_Coverage - Computes the vaccination coverage
F1 - Runs fitting for a spatial stratification of data
F2 - Runs fitting for a spatial stratification of data
F3 - Runs fitting for a spatial stratification of data
F4 - Runs fitting for a spatial stratification of data
Fit_Final_Model - Fits the final vaccination model using all the avaialble data
