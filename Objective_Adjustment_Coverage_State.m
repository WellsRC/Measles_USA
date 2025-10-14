function J = Objective_Adjustment_Coverage_State(x,beta_x,beta_insurance,County_Data_model,State_Data)

dZ_County = x.*ones(size(County_Data_model,1),1);

[Estimated_Vaccination_Coverage] = Vaccination_Coverage_Adjusted([County_Data_model.X County_Data_model.XI County_Data_model.X2],beta_x,beta_insurance,County_Data_model,dZ_County);

v_county=Estimated_Vaccination_Coverage.Overall;
v_state = Estimated_State_Vaccine_Uptake(v_county,County_Data_model.Weight,County_Data_model.State_FIP,State_Data.State_FIP);
   
VAC_STAT=State_Data.Vaccine_Uptake;

Z_State=10^4.*(VAC_STAT-v_state);
J = [Z_State];
end