function J = Objective_Redution_National_Coverage(x,beta_x,beta_insurance,County_Data_model_0_to_4,County_Data_model_5_to_9,dZ_County,National_Reduction)

dZ_Reduction = x;

[Estimated_Vaccination_Coverage] = Vaccination_Coverage_Adjusted([County_Data_model_0_to_4.X County_Data_model_0_to_4.XI County_Data_model_0_to_4.X2],beta_x,beta_insurance,County_Data_model_0_to_4,dZ_County(:,1));
v_county_0_to_4=Estimated_Vaccination_Coverage.Overall;

[Estimated_Vaccination_Coverage] = Vaccination_Coverage_Adjusted([County_Data_model_5_to_9.X County_Data_model_5_to_9.XI County_Data_model_5_to_9.X2],beta_x,beta_insurance,County_Data_model_5_to_9,dZ_County(:,2));
v_county_5_to_9=Estimated_Vaccination_Coverage.Overall;


[Estimated_Vaccination_Coverage] = Vaccination_Coverage_Adjusted([County_Data_model_0_to_4.X County_Data_model_0_to_4.XI County_Data_model_0_to_4.X2],beta_x,beta_insurance,County_Data_model_0_to_4,dZ_County(:,1)+dZ_Reduction);
Reduced_v_county_0_to_4=Estimated_Vaccination_Coverage.Overall;

[Estimated_Vaccination_Coverage] = Vaccination_Coverage_Adjusted([County_Data_model_5_to_9.X County_Data_model_5_to_9.XI County_Data_model_5_to_9.X2],beta_x,beta_insurance,County_Data_model_5_to_9,dZ_County(:,2)+dZ_Reduction);
Reduced_v_county_5_to_9=Estimated_Vaccination_Coverage.Overall;


v_baseline=sum(v_county_0_to_4(:).*County_Data_model_0_to_4.Weight(:)+v_county_5_to_9(:).*County_Data_model_5_to_9.Weight(:))./sum(County_Data_model_0_to_4.Weight(:)+County_Data_model_5_to_9.Weight(:));
v_reduction=sum(Reduced_v_county_0_to_4(:).*County_Data_model_0_to_4.Weight(:)+Reduced_v_county_5_to_9(:).*County_Data_model_5_to_9.Weight(:))./sum(County_Data_model_0_to_4.Weight(:)+County_Data_model_5_to_9.Weight(:));

J=10^4.*(National_Reduction-(v_baseline-v_reduction));
end