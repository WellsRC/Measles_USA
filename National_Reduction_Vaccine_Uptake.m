clear;

Vaccine='MMR';
L_F=zeros(2^7,4);
L_S=zeros(2^7,4);
for Spatial_Valiation=1:4
    load(['Refine_' Vaccine '_Set_Spatial_Validation=' num2str(Spatial_Valiation) '.mat']);
    L_F(:,Spatial_Valiation)=L_fit;
    L_S(:,Spatial_Valiation)=L_spatial_val;
end

L_T=sum(L_F,2)+sum(L_S,2);
    
[~,Model_Num]=max(L_T);

clearvars -except Model_Num Vaccine

Age_Class={'Age_0_to_4','Age_5_to_9'};

dZ_County=zeros(3143,2);
for aa=1:length(Age_Class)
    [~,dZ_County(:,aa)]=Age_Adjustment_Factor(Vaccine,Model_Num,Age_Class{aa});
end

National_Reduction=0.10;
[County_Info,dZ_Reduction]=Decline_Adjustment_Factor(Vaccine,Model_Num,dZ_County,National_Reduction);


load([Vaccine '_Immunity.mat'],'County_Data')

County_Data_Vaccine_Reduction=County_Data;

clear County_Data

epsv1=0.93;
epsv2=0.97;

Immunity_0_to_4=epsv1.*County_Info.Vaccine_Uptake_0_to_4;
Immunity_5_to_9=epsv2.*County_Info.Vaccine_Uptake_5_to_9;

County_Data_Vaccine_Reduction.Immunity(:,1)=array2table(Immunity_0_to_4);
County_Data_Vaccine_Reduction.Immunity(:,2)=array2table(Immunity_5_to_9);

County_Data_Vaccine_Reduction.Total_Immunity=sum(table2array(County_Data_Vaccine_Reduction.Population).*table2array(County_Data_Vaccine_Reduction.Immunity),2);

load('Baseline_Estimate_Measles_Incidence.mat','beta_j','beta_seed');

Final_Size_Est=zeros(length(County_Data_Vaccine_Reduction.State),1);
R_eff=zeros(length(County_Data_Vaccine_Reduction.State),1);
R_eff_Seed=zeros(length(County_Data_Vaccine_Reduction.State),1);

opts=optimoptions('lsqnonlin','OptimalityTolerance',10^(-12),'StepTolerance',10^(-12),'FunctionTolerance',10^(-16),'MaxFunctionEvaluations',10^6,'MaxIterations',10^(6));

for s_indx=1:length(County_Data_Vaccine_Reduction.State)
    State_temp=County_Data_Vaccine_Reduction.State{s_indx};
    State_temp(State_temp==' ')='_';
    M=readtable([pwd '/State_Contact_Matrix/United_States_subnational_' State_temp '_M_overall_contact_matrix_18']);
    M=table2array(M);
    
    Tot_Pop=County_Data_Vaccine_Reduction.Total_Population(s_indx);
    Pop=table2array(County_Data_Vaccine_Reduction.Population(s_indx,:));
    Immunity=table2array(County_Data_Vaccine_Reduction.Immunity(s_indx,:));
    Pop=Tot_Pop.*Pop;
    S_Pop=(1-Immunity).*Pop; % Susceptible population
    
    % Equation 13
    % (https://pmc.ncbi.nlm.nih.gov/articles/PMC7088810/pdf/11538_2010_Article_9623.pdf)
    
    A_eff=beta_j.*repmat(S_Pop,18,1).*M./repmat(Pop,18,1);
    A_eff(repmat(Pop,18,1)==0)=0;
    R_eff(s_indx)=max(eig(A_eff));
    if(max(eig(A_eff))>1)            
        z=lsqnonlin(@(x) A_eff*x(:)+log(1-x(:)), 0.999.*ones(1,18),zeros(1,18),ones(1,18),[],[],[],[],[],opts);
        Final_Size_Est(s_indx)=S_Pop*z(:);
    end
    
    A_eff=beta_seed.*repmat(S_Pop,18,1).*M./repmat(Pop,18,1);
    A_eff(repmat(Pop,18,1)==0)=0;
    R_eff_Seed(s_indx)=max(eig(A_eff));
end


County_Data_Vaccine_Reduction.R_eff_Seed=R_eff_Seed;
County_Data_Vaccine_Reduction.R_eff=R_eff;
County_Data_Vaccine_Reduction.Final_Size_Est=Final_Size_Est;
save(['National_Reduction_' num2str(100*National_Reduction) '.mat'],'County_Data_Vaccine_Reduction')

