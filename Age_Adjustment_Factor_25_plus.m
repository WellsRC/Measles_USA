function [County_Info,dZ_County]=Age_Adjustment_Factor_25_plus(Vaccine,Model_Num,Age_Group)

% Load the data
[County_Data,~] = Load_Data_Adjustment(Vaccine,Age_Group);

% Construct the covariates for the selected mode
Model_Var=false(2^7,43);


Model_Var(:,[1 2 3])=true;

temp_var=false(7,43);
temp_var(1,4:8)=true;
temp_var(2,9:13)=true;
temp_var(3,14)=true;
temp_var(4,15)=true;
temp_var(5,16:22)=true;
temp_var(6,23:34)=true;
temp_var(7,35:43)=true;

bin_model=dec2bin([0:2^7-1]',7)=='1';
[~,srt_indx]=sort(sum(bin_model,2),'ascend');
bin_model=bin_model(srt_indx,:);

for ii=1:2^7
    temp_overall=false(1,43);
    for jj=1:7
        if(bin_model(ii,jj))
            temp_overall=temp_overall|temp_var(jj,:);
        end
    end
    Model_Var(ii,:)=Model_Var(ii,:) | temp_overall;
end

Model_Var=Model_Var(Model_Num,:);

% Contruct the interaction terms and covariates for the selected moel
County_Data_model=County_Data;
County_Data_model.X=County_Data_model.X(:,Model_Var);
County_Data_model.X2=County_Data_model.X;
if(Model_Var(end))
    County_Data_model.X2=County_Data_model.X2(:,3:end-9);
else
    County_Data_model.X2=County_Data_model.X2(:,3:end);
end
XT=[];
for ii=1:size(County_Data_model.X2,2)
    for jj=(ii+1):size(County_Data_model.X2,2)
        XT=[XT County_Data_model.X2(:,ii).*County_Data_model.X2(:,jj)];
    end
end

County_Data_model.X2=County_Data_model.X2.^2;
County_Data_model.XI=XT;   

% Load the trained parameters for the selected model
load([Vaccine '_Model_' num2str(Model_Num) '.mat'],'final_model_par');
[beta_x,beta_insurance,~,~] = Parameters(final_model_par);


dZ_County = zeros(size(County_Data.Weight));

opts_lsq=optimoptions('lsqnonlin','ConstraintTolerance',10^(-12),'FunctionTolerance',10^(-12),'MaxFunctionEvaluations',10^6,'MaxIterations',10^4,'OptimalityTolerance',10^(-12),'StepTolerance',10^(-12));

f_val=zeros(size(County_Data.Weight));

for cc=1:length(f_val)
    trim_County_Data_model.X=County_Data_model.X(cc,:);
    trim_County_Data_model.XI=County_Data_model.XI(cc,:);
    trim_County_Data_model.X2=County_Data_model.X2(cc,:);

    trim_County_Data_model.Uninsured=County_Data_model.Uninsured(cc);
    trim_County_Data_model.Private=County_Data_model.Private(cc);
    trim_County_Data_model.Public=County_Data_model.Public(cc);
    trim_County_Data_model.Immunity=County_Data_model.Immunity(cc);

    [Estimated_Vaccination_Coverage] = Vaccination_Coverage_Adjusted([trim_County_Data_model.X trim_County_Data_model.XI trim_County_Data_model.X2],beta_x,beta_insurance,trim_County_Data_model,0);
        
    v_county=Estimated_Vaccination_Coverage.Overall;
    if(0.97.*v_county>trim_County_Data_model.Immunity)
        trim_County_Data_model.Vaccine_Uptake=trim_County_Data_model.Immunity./0.97;
        X0=log(trim_County_Data_model.Vaccine_Uptake./(1-trim_County_Data_model.Vaccine_Uptake))-log(v_county./(1-v_county));
        [dZ_County(cc),f_val(cc)]=lsqnonlin(@(x)Objective_Adjustment_Coverage_County(x,beta_x,beta_insurance,trim_County_Data_model),X0,-40,40,[],[],[],[],[],opts_lsq);
    end
end

[Estimated_Vaccination_Coverage] = Vaccination_Coverage_Adjusted([County_Data_model.X County_Data_model.XI County_Data_model.X2],beta_x,beta_insurance,County_Data_model,dZ_County);

County_Info.Vaccine_Uptake_Uninsured=Estimated_Vaccination_Coverage.Uninsured;
County_Info.Vaccine_Uptake_Public=Estimated_Vaccination_Coverage.Public;
County_Info.Vaccine_Uptake_Private=Estimated_Vaccination_Coverage.Private;
County_Info.Vaccine_Uptake=Estimated_Vaccination_Coverage.Overall;
County_Info.County=County_Data.County;
County_Info.State=County_Data.State;
County_Info.GEOID=County_Data.GEOID;
County_Info.Population=County_Data.Population;
County_Info.Uninsured=County_Data.Uninsured;
County_Info.Public=County_Data.Public;
County_Info.Private=County_Data.Private;
County_Info.Total_Population=County_Data.Total_Population;

end

