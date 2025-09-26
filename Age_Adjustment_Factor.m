function Age_Adjustment_Factor(Vaccine,Model_Num,Age_Group)
Year=2023;
% Load the data
[County_Data,State_Data] = Load_Data_Adjustment(Vaccine,Age_Group);

County_Data.Vaccine_Uptake(County_Data.Vaccine_Uptake==1)=1-10^(-4);
County_Data.Vaccine_Uptake(County_Data.Vaccine_Uptake==0)=10^(-4);

State_Data.Vaccine_Uptake(State_Data.Vaccine_Uptake==1)=1-10^(-4);
State_Data.Vaccine_Uptake(State_Data.Vaccine_Uptake==0)=10^(-4);

% Construct the covariates for the selected mode
Model_Var=false(2^7,42);

Model_Var(:,[1 2])=true;

temp_var=false(7,42);
temp_var(1,3:7)=true;
temp_var(2,8:12)=true;
temp_var(3,13)=true;
temp_var(4,14)=true;
temp_var(5,15:21)=true;
temp_var(6,22:33)=true;
temp_var(7,34:42)=true;


bin_model=dec2bin([0:2^7-1]',7)=='1';
[~,srt_indx]=sort(sum(bin_model,2),'ascend');
bin_model=bin_model(srt_indx,:);

for ii=1:2^7
    temp_overall=false(1,42);
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
County_Data_model.X2=County_Data_model.X.^2;
if(Model_Var(end))
    County_Data_model.X2=County_Data_model.X2(:,3:end-9);
else
    County_Data_model.X2=County_Data_model.X2(:,3:end);
end
XT=[];
for ii=3:size(County_Data_model.X2,2)
    for jj=(ii+1):size(County_Data_model.X2,2)
        XT=[XT County_Data_model.X(:,ii).*County_Data_model.X(:,jj)];
    end
end
County_Data_model.XI=XT;  

% Load the trained parameters for the selected model
load([Vaccine '_Model_' num2str(Model_Num) '.mat'],'model_par');
[beta_x,beta_insurance,~,~] = Parameters(model_par);

% Find the coutnies where we have known values of vaccine uptale
known_county=~isnan(County_Data.Vaccine_Uptake);
% Fidn states where we have known info about state-level vaccine uptake
state_fill=ismember(County_Data.State_FIP,State_Data.State_FIP) & isnan(County_Data.Vaccine_Uptake);
u_state_fill=unique(County_Data.State_FIP(state_fill));
% Specificy the bounds for the adjustment
lb=-30.*ones(sum(known_county)+length(u_state_fill),1);
ub=30.*ones(sum(known_county)+length(u_state_fill),1);

opts=optimoptions('lsqnonlin','ConstraintTolerance',10^(-12),'FunctionTolerance',10^(-12),'MaxFunctionEvaluations',10^6,'MaxIterations',10^4,'OptimalityTolerance',10^(-12),'StepTolerance',10^(-12));
[dZ_County_est,test_val]=lsqnonlin(@(x)Objective_Adjustment_Coverage(x,beta_x,beta_insurance,County_Data_model,known_county,state_fill,u_state_fill,State_Data),0.*lb,lb,ub,[],[],[],[],[],opts);


dZ_County = zeros(size(known_county));

dZ_County(known_county)=dZ_County_est(1:sum(known_county));

for ss=1:length(u_state_fill)
    t_f=state_fill & County_Data_model.State_FIP==u_state_fill(ss);
    dZ_County(t_f)=dZ_County_est(sum(known_county)+ss);
end

save([Vaccine '_Adjustment_' Age_Group '_Year=' num2str(Year) '_Model=' num2str(Model_Num) '.mat'],'dZ_County');

end

