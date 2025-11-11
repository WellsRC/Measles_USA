function [County_Info,dZ_Reduction]=Decline_Adjustment_Factor(Vaccine,Model_Num,dZ_County,Age_Group_Reduction)

% Load the data
[County_Data_0_to_4,~] = Load_Data_Adjustment(Vaccine,'Age_0_to_4');

County_Data_0_to_4.Vaccine_Uptake(County_Data_0_to_4.Vaccine_Uptake==1)=1-10^(-4);
County_Data_0_to_4.Vaccine_Uptake(County_Data_0_to_4.Vaccine_Uptake==0)=10^(-4);

    [County_Data_5_to_9,~] = Load_Data_Adjustment(Vaccine,'Age_5_to_9');
    
    County_Data_5_to_9.Vaccine_Uptake(County_Data_5_to_9.Vaccine_Uptake==1)=1-10^(-4);
    County_Data_5_to_9.Vaccine_Uptake(County_Data_5_to_9.Vaccine_Uptake==0)=10^(-4);

    [County_Data_10_to_14,~] = Load_Data_Adjustment(Vaccine,'Age_10_to_14');

    County_Data_10_to_14.Vaccine_Uptake(County_Data_10_to_14.Vaccine_Uptake==1)=1-10^(-4);
    County_Data_10_to_14.Vaccine_Uptake(County_Data_10_to_14.Vaccine_Uptake==0)=10^(-4);


    % [County_Data_15_to_19,~] = Load_Data_Adjustment(Vaccine,'Age_15_to_19');
    % 
    % County_Data_15_to_19.Vaccine_Uptake(County_Data_15_to_19.Vaccine_Uptake==1)=1-10^(-4);
    % County_Data_15_to_19.Vaccine_Uptake(County_Data_15_to_19.Vaccine_Uptake==0)=10^(-4);
    % 
    % [County_Data_20_to_24,~] = Load_Data_Adjustment(Vaccine,'Age_20_to_24');
    % 
    % County_Data_20_to_24.Vaccine_Uptake(County_Data_20_to_24.Vaccine_Uptake==1)=1-10^(-4);
    % County_Data_20_to_24.Vaccine_Uptake(County_Data_20_to_24.Vaccine_Uptake==0)=10^(-4);


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
County_Data_model_0_to_4=County_Data_0_to_4;
County_Data_model_0_to_4.X=County_Data_model_0_to_4.X(:,Model_Var);
County_Data_model_0_to_4.X2=County_Data_model_0_to_4.X;
if(Model_Var(end))
    County_Data_model_0_to_4.X2=County_Data_model_0_to_4.X2(:,3:end-9);
else
    County_Data_model_0_to_4.X2=County_Data_model_0_to_4.X2(:,3:end);
end
XT=[];
for ii=1:size(County_Data_model_0_to_4.X2,2)
    for jj=(ii+1):size(County_Data_model_0_to_4.X2,2)
        XT=[XT County_Data_model_0_to_4.X2(:,ii).*County_Data_model_0_to_4.X2(:,jj)];
    end
end

County_Data_model_0_to_4.X2=County_Data_model_0_to_4.X2.^2;
County_Data_model_0_to_4.XI=XT;   

% 5 to 9
    County_Data_model_5_to_9=County_Data_5_to_9;
    County_Data_model_5_to_9.X=County_Data_model_5_to_9.X(:,Model_Var);
    County_Data_model_5_to_9.X2=County_Data_model_5_to_9.X;
    if(Model_Var(end))
        County_Data_model_5_to_9.X2=County_Data_model_5_to_9.X2(:,3:end-9);
    else
        County_Data_model_5_to_9.X2=County_Data_model_5_to_9.X2(:,3:end);
    end
    XT=[];
    for ii=1:size(County_Data_model_5_to_9.X2,2)
        for jj=(ii+1):size(County_Data_model_5_to_9.X2,2)
            XT=[XT County_Data_model_5_to_9.X2(:,ii).*County_Data_model_5_to_9.X2(:,jj)];
        end
    end
    
    County_Data_model_5_to_9.X2=County_Data_model_5_to_9.X2.^2;
    County_Data_model_5_to_9.XI=XT;   

% 10 to 14
    County_Data_model_10_to_14=County_Data_10_to_14;
    County_Data_model_10_to_14.X=County_Data_model_10_to_14.X(:,Model_Var);
    County_Data_model_10_to_14.X2=County_Data_model_10_to_14.X;
    if(Model_Var(end))
        County_Data_model_10_to_14.X2=County_Data_model_10_to_14.X2(:,3:end-9);
    else
        County_Data_model_10_to_14.X2=County_Data_model_10_to_14.X2(:,3:end);
    end
    XT=[];
    for ii=1:size(County_Data_model_10_to_14.X2,2)
        for jj=(ii+1):size(County_Data_model_10_to_14.X2,2)
            XT=[XT County_Data_model_10_to_14.X2(:,ii).*County_Data_model_10_to_14.X2(:,jj)];
        end
    end
    
    County_Data_model_10_to_14.X2=County_Data_model_10_to_14.X2.^2;
    County_Data_model_10_to_14.XI=XT;   


%     County_Data_model_15_to_19=County_Data_15_to_19;
%     County_Data_model_15_to_19.X=County_Data_model_15_to_19.X(:,Model_Var);
%     County_Data_model_15_to_19.X2=County_Data_model_15_to_19.X;
%     if(Model_Var(end))
%         County_Data_model_15_to_19.X2=County_Data_model_15_to_19.X2(:,3:end-9);
%     else
%         County_Data_model_15_to_19.X2=County_Data_model_15_to_19.X2(:,3:end);
%     end
%     XT=[];
%     for ii=1:size(County_Data_model_15_to_19.X2,2)
%         for jj=(ii+1):size(County_Data_model_15_to_19.X2,2)
%             XT=[XT County_Data_model_15_to_19.X2(:,ii).*County_Data_model_15_to_19.X2(:,jj)];
%         end
%     end
% 
%     County_Data_model_15_to_19.X2=County_Data_model_15_to_19.X2.^2;
%     County_Data_model_15_to_19.XI=XT;   

%     County_Data_model_20_to_24=County_Data_20_to_24;
%     County_Data_model_20_to_24.X=County_Data_model_20_to_24.X(:,Model_Var);
%     County_Data_model_20_to_24.X2=County_Data_model_20_to_24.X;
%     if(Model_Var(end))
%         County_Data_model_20_to_24.X2=County_Data_model_20_to_24.X2(:,3:end-9);
%     else
%         County_Data_model_20_to_24.X2=County_Data_model_20_to_24.X2(:,3:end);
%     end
%     XT=[];
%     for ii=1:size(County_Data_model_20_to_24.X2,2)
%         for jj=(ii+1):size(County_Data_model_20_to_24.X2,2)
%             XT=[XT County_Data_model_20_to_24.X2(:,ii).*County_Data_model_20_to_24.X2(:,jj)];
%         end
%     end
% 
%     County_Data_model_20_to_24.X2=County_Data_model_20_to_24.X2.^2;
%     County_Data_model_20_to_24.XI=XT;   


% Load the trained parameters for the selected model
load([Vaccine '_Model_' num2str(Model_Num) '.mat'],'final_model_par');
[beta_x,beta_insurance,~,~] = Parameters(final_model_par);


opts_lsq=optimoptions('fmincon','FunctionTolerance',10^(-12),'MaxFunctionEvaluations',10^6,'MaxIterations',10^4,'OptimalityTolerance',10^(-12),'StepTolerance',10^(-12));

[dZ_Reduction_0_4,f_val_0_4]=fmincon(@(x)Objective_Redution_National_Coverage(x,beta_x,beta_insurance,County_Data_model_0_to_4,dZ_County(:,1),Age_Group_Reduction(1)),0,[],[],[],[],-40,40,[],opts_lsq);    

[Estimated_Vaccination_Coverage] = Vaccination_Coverage_Adjusted([County_Data_model_0_to_4.X County_Data_model_0_to_4.XI County_Data_model_0_to_4.X2],beta_x,beta_insurance,County_Data_model_0_to_4,dZ_County(:,1)+dZ_Reduction_0_4);
County_Info.Vaccine_Uptake_0_to_4=Estimated_Vaccination_Coverage.Overall;

[dZ_Reduction_5_9,f_val_5_9]=fmincon(@(x)Objective_Redution_National_Coverage(x,beta_x,beta_insurance,County_Data_model_5_to_9,dZ_County(:,2),Age_Group_Reduction(2)),0,[],[],[],[],-40,40,[],opts_lsq);    

[Estimated_Vaccination_Coverage] = Vaccination_Coverage_Adjusted([County_Data_model_5_to_9.X County_Data_model_5_to_9.XI County_Data_model_5_to_9.X2],beta_x,beta_insurance,County_Data_model_5_to_9,dZ_County(:,2)+dZ_Reduction_5_9);
County_Info.Vaccine_Uptake_5_to_9=Estimated_Vaccination_Coverage.Overall;

[dZ_Reduction_10_14,f_val_10_14]=fmincon(@(x)Objective_Redution_National_Coverage(x,beta_x,beta_insurance,County_Data_model_10_to_14,dZ_County(:,3),Age_Group_Reduction(3)),0,[],[],[],[],-40,40,[],opts_lsq);    

[Estimated_Vaccination_Coverage] = Vaccination_Coverage_Adjusted([County_Data_model_10_to_14.X County_Data_model_10_to_14.XI County_Data_model_10_to_14.X2],beta_x,beta_insurance,County_Data_model_10_to_14,dZ_County(:,3)+dZ_Reduction_10_14);
County_Info.Vaccine_Uptake_10_to_14=Estimated_Vaccination_Coverage.Overall;

County_Info.County=County_Data_0_to_4.County;
County_Info.State=County_Data_0_to_4.State;
County_Info.GEOID=County_Data_0_to_4.GEOID;
County_Info.Population=County_Data_0_to_4.Population;
County_Info.Total_Population=County_Data_0_to_4.Total_Population;

dZ_Reduction=[dZ_Reduction_0_4 dZ_Reduction_5_9 dZ_Reduction_10_14];

end

