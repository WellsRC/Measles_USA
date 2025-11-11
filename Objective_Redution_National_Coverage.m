function J = Objective_Redution_National_Coverage(x,beta_x,beta_insurance,County_Data_model_Age_Group,dZ_County,VC_Reduction)

dZ_Reduction = x;

[Estimated_Vaccination_Coverage] = Vaccination_Coverage_Adjusted([County_Data_model_Age_Group.X County_Data_model_Age_Group.XI County_Data_model_Age_Group.X2],beta_x,beta_insurance,County_Data_model_Age_Group,dZ_County(:,1));
v_county_Age_Group=Estimated_Vaccination_Coverage.Overall(:).*County_Data_model_Age_Group.Weight(:);
Pop_Age_Group=County_Data_model_Age_Group.Weight(:);

[Estimated_Vaccination_Coverage] = Vaccination_Coverage_Adjusted([County_Data_model_Age_Group.X County_Data_model_Age_Group.XI County_Data_model_Age_Group.X2],beta_x,beta_insurance,County_Data_model_Age_Group,dZ_County(:,1)+dZ_Reduction);
Reduced_v_county_Age_Group=Estimated_Vaccination_Coverage.Overall(:).*County_Data_model_Age_Group.Weight(:);


v_baseline=sum(v_county_Age_Group(:))./sum(Pop_Age_Group(:));
v_reduction=sum(Reduced_v_county_Age_Group(:))./sum(Pop_Age_Group(:));

J=10^4.*(VC_Reduction-(v_baseline-v_reduction)).^2;
end