function J = Objective_Adjustment_Coverage_County(x,beta_x,beta_insurance,County_Data_model)

dZ_County = x;

[Estimated_Vaccination_Coverage] = Vaccination_Coverage_Adjusted([County_Data_model.X County_Data_model.XI County_Data_model.X2],beta_x,beta_insurance,County_Data_model,dZ_County);

v_county=Estimated_Vaccination_Coverage.Overall;

Z_county=10^4.*(County_Data_model.Vaccine_Uptake-v_county);
J = [Z_county];
end