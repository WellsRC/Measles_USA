function J = Objective_Redution_National_Coverage(x,beta_x,beta_insurance,County_Data_model_0_to_4,County_Data_model_5_to_9,County_Data_model_10_to_14,County_Data_model_15_to_19,County_Data_model_20_to_24,dZ_County,National_Reduction,Age_0_6)

dZ_Reduction = x;

[Estimated_Vaccination_Coverage] = Vaccination_Coverage_Adjusted([County_Data_model_0_to_4.X County_Data_model_0_to_4.XI County_Data_model_0_to_4.X2],beta_x,beta_insurance,County_Data_model_0_to_4,dZ_County(:,1));
v_county_0_to_4=Estimated_Vaccination_Coverage.Overall(:).*County_Data_model_0_to_4.Weight(:);
Pop_0_to_4=County_Data_model_0_to_4.Weight(:);

[Estimated_Vaccination_Coverage] = Vaccination_Coverage_Adjusted([County_Data_model_0_to_4.X County_Data_model_0_to_4.XI County_Data_model_0_to_4.X2],beta_x,beta_insurance,County_Data_model_0_to_4,dZ_County(:,1)+dZ_Reduction);
Reduced_v_county_0_to_4=Estimated_Vaccination_Coverage.Overall(:).*County_Data_model_0_to_4.Weight(:);

if(~isempty(County_Data_model_5_to_9))
    [Estimated_Vaccination_Coverage] = Vaccination_Coverage_Adjusted([County_Data_model_5_to_9.X County_Data_model_5_to_9.XI County_Data_model_5_to_9.X2],beta_x,beta_insurance,County_Data_model_5_to_9,dZ_County(:,2));
    if(~Age_0_6)
        v_county_5_to_9=Estimated_Vaccination_Coverage.Overall(:).*County_Data_model_5_to_9.Weight(:);
        Pop_5_to_9=County_Data_model_5_to_9.Weight(:);
    else
        v_county_5_to_9=(2/5).*Estimated_Vaccination_Coverage.Overall(:).*County_Data_model_5_to_9.Weight(:);
        Pop_5_to_9=(2/5).*County_Data_model_5_to_9.Weight(:);
    end

    [Estimated_Vaccination_Coverage] = Vaccination_Coverage_Adjusted([County_Data_model_5_to_9.X County_Data_model_5_to_9.XI County_Data_model_5_to_9.X2],beta_x,beta_insurance,County_Data_model_5_to_9,dZ_County(:,2)+dZ_Reduction);
    if(~Age_0_6)
        Reduced_v_county_5_to_9=Estimated_Vaccination_Coverage.Overall(:).*County_Data_model_5_to_9.Weight(:);
    else
        Reduced_v_county_5_to_9=(2/5).*Estimated_Vaccination_Coverage.Overall(:).*County_Data_model_5_to_9.Weight(:);
    end
else
    v_county_5_to_9=0;
    Reduced_v_county_5_to_9=0;
    Pop_5_to_9=0;
end

if(~isempty(County_Data_model_10_to_14))
    [Estimated_Vaccination_Coverage] = Vaccination_Coverage_Adjusted([County_Data_model_10_to_14.X County_Data_model_10_to_14.XI County_Data_model_10_to_14.X2],beta_x,beta_insurance,County_Data_model_10_to_14,dZ_County(:,3));
    v_county_10_to_14=Estimated_Vaccination_Coverage.Overall(:).*County_Data_model_10_to_14.Weight(:);
    Pop_10_to_14=County_Data_model_10_to_14.Weight(:);

    [Estimated_Vaccination_Coverage] = Vaccination_Coverage_Adjusted([County_Data_model_10_to_14.X County_Data_model_10_to_14.XI County_Data_model_10_to_14.X2],beta_x,beta_insurance,County_Data_model_10_to_14,dZ_County(:,3)+dZ_Reduction);
    Reduced_v_county_10_to_14=Estimated_Vaccination_Coverage.Overall(:).*County_Data_model_10_to_14.Weight(:);
else
    v_county_10_to_14=0;
    Reduced_v_county_10_to_14=0;
    Pop_10_to_14=0;
end

if(~isempty(County_Data_model_15_to_19))
    [Estimated_Vaccination_Coverage] = Vaccination_Coverage_Adjusted([County_Data_model_15_to_19.X County_Data_model_15_to_19.XI County_Data_model_15_to_19.X2],beta_x,beta_insurance,County_Data_model_15_to_19,dZ_County(:,4));
    v_county_15_to_19=Estimated_Vaccination_Coverage.Overall(:).*County_Data_model_15_to_19.Weight(:);
    Pop_15_to_19=County_Data_model_15_to_19.Weight(:);

    [Estimated_Vaccination_Coverage] = Vaccination_Coverage_Adjusted([County_Data_model_15_to_19.X County_Data_model_15_to_19.XI County_Data_model_15_to_19.X2],beta_x,beta_insurance,County_Data_model_15_to_19,dZ_County(:,4)+dZ_Reduction);
    Reduced_v_county_15_to_19=Estimated_Vaccination_Coverage.Overall(:).*County_Data_model_15_to_19.Weight(:);
else
    v_county_15_to_19=0;
    Reduced_v_county_15_to_19=0;
    Pop_15_to_19=0;
end

if(~isempty(County_Data_model_20_to_24))
    [Estimated_Vaccination_Coverage] = Vaccination_Coverage_Adjusted([County_Data_model_20_to_24.X County_Data_model_20_to_24.XI County_Data_model_20_to_24.X2],beta_x,beta_insurance,County_Data_model_20_to_24,dZ_County(:,5));
    v_county_20_to_24=Estimated_Vaccination_Coverage.Overall(:).*County_Data_model_20_to_24.Weight(:);
    Pop_20_to_24=County_Data_model_20_to_24.Weight(:);

    [Estimated_Vaccination_Coverage] = Vaccination_Coverage_Adjusted([County_Data_model_20_to_24.X County_Data_model_20_to_24.XI County_Data_model_20_to_24.X2],beta_x,beta_insurance,County_Data_model_20_to_24,dZ_County(:,5)+dZ_Reduction);
    Reduced_v_county_20_to_24=Estimated_Vaccination_Coverage.Overall(:).*County_Data_model_20_to_24.Weight(:);
else
    v_county_20_to_24=0;
    Reduced_v_county_20_to_24=0;
    Pop_20_to_24=0;
end

v_baseline=sum(v_county_0_to_4(:)+v_county_5_to_9(:)+v_county_10_to_14(:)+v_county_15_to_19(:)+v_county_20_to_24(:))./sum(Pop_0_to_4(:)+Pop_5_to_9(:)+Pop_10_to_14(:)+Pop_15_to_19(:)+Pop_20_to_24(:));
v_reduction=sum(Reduced_v_county_0_to_4(:)+Reduced_v_county_5_to_9(:)+Reduced_v_county_10_to_14(:)+Reduced_v_county_15_to_19(:)+Reduced_v_county_20_to_24(:))./sum(Pop_0_to_4(:)+Pop_5_to_9(:)+Pop_10_to_14(:)+Pop_15_to_19(:)+Pop_20_to_24(:));

J=10^4.*(National_Reduction-(v_baseline-v_reduction));
end