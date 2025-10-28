function [County_Info,dZ_County]=Age_Adjustment_Factor(Vaccine,Model_Num,Age_Group)
Year=2023;
% Load the data
[County_Data,State_Data] = Load_Data_Adjustment(Vaccine,Age_Group);

County_Data.Vaccine_Uptake(County_Data.Vaccine_Uptake==1)=1-10^(-4);
County_Data.Vaccine_Uptake(County_Data.Vaccine_Uptake==0)=10^(-4);

State_Data.Vaccine_Uptake(State_Data.Vaccine_Uptake==1)=1-10^(-4);
State_Data.Vaccine_Uptake(State_Data.Vaccine_Uptake==0)=10^(-4);

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

% Find the coutnies where we have known values of vaccine uptale
known_county=~isnan(County_Data.Vaccine_Uptake);
u_state_known=unique(County_Data.State_FIP(known_county));
% Fidn states where we have known info about state-level vaccine uptake
state_fill=ismember(County_Data.State_FIP,State_Data.State_FIP) & isnan(County_Data.Vaccine_Uptake);
u_state_fill=unique(County_Data.State_FIP(state_fill));
u_state_fill=u_state_fill(~ismember(u_state_fill,u_state_known));  
% Specificy the bounds for the adjustment


dZ_County = NaN.*zeros(size(known_county));

opts_lsq=optimoptions('lsqnonlin','ConstraintTolerance',10^(-12),'FunctionTolerance',10^(-12),'MaxFunctionEvaluations',10^6,'MaxIterations',10^4,'OptimalityTolerance',10^(-12),'StepTolerance',10^(-12));

dZ_County_est=zeros(sum(known_county)+length(u_state_fill),1);
f_val=zeros(sum(known_county)+length(u_state_fill),1);

indx_known_county=find(known_county);

for cc=1:length(indx_known_county)
    t_f=indx_known_county(cc);
    trim_County_Data_model.X=County_Data_model.X(t_f,:);
    trim_County_Data_model.XI=County_Data_model.XI(t_f,:);
    trim_County_Data_model.X2=County_Data_model.X2(t_f,:);

    trim_County_Data_model.Uninsured=County_Data_model.Uninsured(t_f);
    trim_County_Data_model.Private=County_Data_model.Private(t_f);
    trim_County_Data_model.Public=County_Data_model.Public(t_f);

    trim_County_Data_model.Vaccine_Uptake=County_Data_model.Vaccine_Uptake(t_f);

    [Estimated_Vaccination_Coverage] = Vaccination_Coverage_Adjusted([trim_County_Data_model.X trim_County_Data_model.XI trim_County_Data_model.X2],beta_x,beta_insurance,trim_County_Data_model,0);
        
    v_county=Estimated_Vaccination_Coverage.Overall;

    X0=log(trim_County_Data_model.Vaccine_Uptake./(1-trim_County_Data_model.Vaccine_Uptake))-log(v_county./(1-v_county));

    [dZ_County_est(cc),f_val(cc)]=lsqnonlin(@(x)Objective_Adjustment_Coverage_County(x,beta_x,beta_insurance,trim_County_Data_model),X0,-40,40,[],[],[],[],[],opts_lsq);
  test=0;  
end

for ss=1:length(u_state_fill)
    t_f=County_Data.State_FIP==u_state_fill(ss);
    trim_County_Data_model.X=County_Data_model.X(t_f,:);
    trim_County_Data_model.XI=County_Data_model.XI(t_f,:);
    trim_County_Data_model.X2=County_Data_model.X2(t_f,:);

    trim_County_Data_model.Uninsured=County_Data_model.Uninsured(t_f);
    trim_County_Data_model.Private=County_Data_model.Private(t_f);
    trim_County_Data_model.Public=County_Data_model.Public(t_f);

    trim_County_Data_model.Vaccine_Uptake=County_Data_model.Vaccine_Uptake(t_f);

    trim_County_Data_model.Weight=County_Data_model.Weight(t_f);
    trim_County_Data_model.State_FIP=County_Data_model.State_FIP(t_f);

    t_f=State_Data.State_FIP==u_state_fill(ss);
    trim_State_Data.Vaccine_Uptake=State_Data.Vaccine_Uptake(t_f);
    trim_State_Data.State_FIP=State_Data.State_FIP(t_f);

    [Estimated_Vaccination_Coverage] = Vaccination_Coverage_Adjusted([trim_County_Data_model.X trim_County_Data_model.XI trim_County_Data_model.X2],beta_x,beta_insurance,trim_County_Data_model,0);
    
    v_county=Estimated_Vaccination_Coverage.Overall;
    v_state = Estimated_State_Vaccine_Uptake(v_county,trim_County_Data_model.Weight,trim_County_Data_model.State_FIP,trim_State_Data.State_FIP);

    X0=log(trim_State_Data.Vaccine_Uptake./(1-trim_State_Data.Vaccine_Uptake))-log(v_state./(1-v_state));

    [dZ_County_est(sum(known_county)+ss),f_val(sum(known_county)+ss)]=lsqnonlin(@(x)Objective_Adjustment_Coverage_State(x,beta_x,beta_insurance,trim_County_Data_model,trim_State_Data),X0,-40,40,[],[],[],[],[],opts_lsq);
    test=0;
end



dZ_County(known_county)=dZ_County_est(1:sum(known_county));
for ss=1:length(u_state_fill)
    t_f=state_fill & County_Data_model.State_FIP==u_state_fill(ss);
    dZ_County(t_f)=dZ_County_est(sum(known_county)+ss);
end

temp_dZ=dZ_County;
f_nan=find(isnan(dZ_County));

for nn=1:length(f_nan)
    t_samp=County_Data_model.State_FIP==County_Data_model.State_FIP(f_nan(nn)) & ~isnan(temp_dZ);

    if(sum(t_samp)==0)
        dZ_County(f_nan(nn))=0;
    else
        dZ_County(f_nan(nn))=median(temp_dZ(t_samp));
    end
end

[Estimated_Vaccination_Coverage] = Vaccination_Coverage_Adjusted([County_Data_model.X County_Data_model.XI County_Data_model.X2],beta_x,beta_insurance,County_Data_model,dZ_County);

County_Info.Vaccine_Uptake=Estimated_Vaccination_Coverage.Overall;
County_Info.County=County_Data.County;
County_Info.State=County_Data.State;
County_Info.GEOID=County_Data.GEOID;
County_Info.Population=County_Data.Population;
County_Info.Total_Population=County_Data.Total_Population;

end

