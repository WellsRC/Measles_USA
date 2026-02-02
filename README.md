# The health and economic repercussions of declining MMR coverage in the United States

Chad R. Wells <sup>1</sup>, Abhishek Pandey <sup>1</sup>, Yang Ye <sup>1</sup>, Carolyn Bawden <sup>1</sup>, Rebecca Giglio <sup>2</sup>,  Charlene Wong <sup>2</sup>, Velda Wang <sup>2</sup>, Chelsea Cipriano <sup>2</sup>, Lamia Ayaz <sup>1</sup>, Gergely Röst <sup>3</sup>, Seyed M. Moghadas <sup>4</sup>, Meagan C. Fitzpatrick <sup>5</sup>, Burton H. Singer <sup>6</sup>, Alison P. Galvani <sup>1</sup>

<sup>1</sup> Center for Infectious Disease Modeling and Analysis (CIDMA), Yale School of Public Health, New Haven, Connecticut, 06510, USA<br /> 
<sup>2</sup> The Common Health Coalition, New York, NY, 10016<br /> 
<sup>3</sup> National Laboratory for Health Security, University of Szeged, Szeged, Hungary<br /> 
<sup>4</sup> Agent-Based Modelling Laboratory, York University, Toronto, Ontario, M3J 1P3, Canada<br /> 
<sup>5</sup> Center for Vaccine Development and Global Health, University of Maryland School of Medicine, Baltimore, Maryland, 21201, USA<br /> 
<sup>6</sup> Emerging Pathogens Institute, University of Florida, Gainesville, Florida, 32610, USA<br /> 

Copyright 2025, Chad Wells et al. All rights reserved. Released under the GNU GENERAL PUBLIC LICENSE v3.

The MATLAB code provided here will run the fitting and analysis for the modelling portion of the manuscript.

## Scripts and functions
Age_Adjustment_Factor_25_plus - Computes the adjustment factors for the age groups 25 years and older <br /> 
Age_Adjustment_Factor - Computes the adjustment factors for the age groups 0-24 years of age <br /> 
Case_Importation_Sample Constructs the sample of imported measles cases into the US <br /> 
Chain_Size_Distribution_CDF - Cumulative distribution function for the stuttering chain <br />
Chain_Size_Distribution - Probability distribution function for the stuttering chain  <br />
Compute_Age_Group_Reduction - Computes the level of reduction in the 0-4; 5-9; 10-14 age groups for a specified annual reduction <br /> 
County_Level_Costs - Computes the total costs at the county level <br />
Compute_Direct_Medical_Costs_County - Computes the direct medical costs at the county level <br />
Compute_Direct_Medical_Costs - Computes the direct medical costs at the national level  <br />
Compute_Productivity_Losses_County - Computes productivty losses at the county level  <br />
Compute_Productivity_Losses - Computes productivty losses at the national level  <br />
County_Transmission - Computes the transmission rate for the county <br />
Death_Probability - Returns the probability of death for specified age groups <br />
Decline_Adjustment_Factor - computes the adjustment factors required to get the specified level of reduction <br /> 
Determine_Model_County_Case_Count - Interpolation function used to approximate the final size for a specified transmission rate during the optimization process <br />
Estimate_HDI - Estimates the mode and the HDI for the provided psoterior kernel distribtuon <br />
Estimated_NE_Health_District_Vaccine_Uptake - computes the vaccine coverage among the Nebraska Health Districts <br /> 
Estimated_State_Vaccine_Uptake - computes the vaccine coverage among the states <br /> 
Fit_Model - Runs the optimization for the specified training dataset for the vaccination model <br /> 
Hospitalization_Probability - Returns the probability of hospitalization for specified age groups <br />
Hurdle_Parameters - Returns the parameters for the Hurdle Model <br />
Load_Data_Adjustment - loads the vaccination data required for computing the necessary adjustments to get the county and state level coverage <br /> 
Load_Transmission_Covariates - Loads the covariates used in determining the county-level transmission <br /> 
Load_Data - Loads the data for the vaccination model and the data used in the training of the model <br /> 
Objective_Adjustment_Coverage_County - the objective function used for paramterizing the adjustments factor at the county level <br /> 
Objective_Adjustment_Coverage_State - the objective function used for paramterizing the adjustments factor at the state level <br /> 
Objective_Function_Coverage - The obejective fucntion used in the training of the vaccination model <br /> 
Parameters - returns the parameters for the specified vaccination model <br /> 
Refine_Fit_Model - Refines the fitting of the vaccination model; This is run after "Fit_Model" <br /> 
Vaccination_Coverage_Adjusted - Computes the adjusted vaccination coverage <br /> 
Vaccination_Coverage - Computes the vaccination coverage <br /> 
F1 - Runs fitting for a spatial stratification of data <br /> 
F2 - Runs fitting for a spatial stratification of data <br /> 
F3 - Runs fitting for a spatial stratification of data <br /> 
F4 - Runs fitting for a spatial stratification of data <br /> 
Fit_Final_Model - Fits the final vaccination model using all the avaialble data <br /> 
