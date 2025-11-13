function REMOVE_Reduction_Vaccine_Uptake()

load(['MMR_Immunity.mat'],'County_Data')

County_Data_Vaccine_Reduction=County_Data;

clear County_Data

County_Data_Vaccine_Reduction.Immunity=0.65.*County_Data_Vaccine_Reduction.Immunity;


County_Data_Vaccine_Reduction.Total_Immunity=sum(table2array(County_Data_Vaccine_Reduction.Population).*table2array(County_Data_Vaccine_Reduction.Immunity),2);

load('Baseline_Estimate_Measles_Incidence.mat','beta_j','beta_seed');

Final_Size_Est=zeros(length(County_Data_Vaccine_Reduction.State),1);
Proportion_Size_Age_Unvaccinated=zeros(length(County_Data_Vaccine_Reduction.State),18);
Proportion_Size_Age_Vaccinated=zeros(length(County_Data_Vaccine_Reduction.State),18);
R_eff=zeros(length(County_Data_Vaccine_Reduction.State),1);
R_eff_Seed=zeros(length(County_Data_Vaccine_Reduction.State),1);

R_0=zeros(length(County_Data_Vaccine_Reduction.State),1);

All_Contacts=zeros(length(County_Data_Vaccine_Reduction.State),18);
Unvaccinated_Contacts=zeros(length(County_Data_Vaccine_Reduction.State),18);

opts=optimoptions('lsqnonlin','OptimalityTolerance',10^(-12),'StepTolerance',10^(-12),'FunctionTolerance',10^(-16),'MaxFunctionEvaluations',10^6,'MaxIterations',10^(6));

for s_indx=1:length(County_Data_Vaccine_Reduction.State)
    State_temp=County_Data_Vaccine_Reduction.State{s_indx};
    State_temp(State_temp==' ')='_';
    M=readtable([pwd '/State_Contact_Matrix/United_States_subnational_' State_temp '_M_overall_contact_matrix_18']);
    M=table2array(M);

    All_Contacts(s_indx,:)=sum(M,2);
    
    Tot_Pop=County_Data_Vaccine_Reduction.Total_Population(s_indx);
    Pop=table2array(County_Data_Vaccine_Reduction.Population(s_indx,:));
    Immunity=table2array(County_Data_Vaccine_Reduction.Immunity(s_indx,:));
    Pop=Tot_Pop.*Pop;
    S_Pop=(1-Immunity).*Pop; % Susceptible population
    
    % Equation 13
    % (https://pmc.ncbi.nlm.nih.gov/articles/PMC7088810/pdf/11538_2010_Article_9623.pdf)
    A0=beta_j.*repmat(Pop,18,1).*M./repmat(Pop,18,1);
    A0(repmat(Pop,18,1)==0)=0;
    R_0(s_indx)=max(abs(eig(A0)));

    A_eff=beta_j.*repmat(S_Pop,18,1).*M./repmat(Pop,18,1);
    A_eff(repmat(Pop,18,1)==0)=0;
    R_eff(s_indx)=max(abs(eig(A_eff)));
    if(max(abs(eig(A_eff)))>1)            
        z=lsqnonlin(@(x) A_eff*x(:)+log(1-x(:)), 0.999.*ones(1,18),zeros(1,18),ones(1,18),[],[],[],[],[],opts);
        Final_Size_Est(s_indx)=S_Pop*z(:);
    else
        n_s=S_Pop./sum(S_Pop);
    end

    A_eff=beta_seed.*repmat(S_Pop,18,1).*M./repmat(Pop,18,1);
    A_eff(repmat(Pop,18,1)==0)=0;
    R_eff_Seed(s_indx)=max(abs(eig(A_eff)));
end

County_Data_Vaccine_Reduction.R_eff_Seed=R_eff_Seed;
County_Data_Vaccine_Reduction.R_eff=R_eff;
County_Data_Vaccine_Reduction.R_0=R_0;
County_Data_Vaccine_Reduction.Final_Size_Est=Final_Size_Est;

save(['TEST_MODEL_IMM=60.mat'],'County_Data_Vaccine_Reduction')

end