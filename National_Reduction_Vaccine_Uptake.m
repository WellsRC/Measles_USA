function National_Reduction_Vaccine_Uptake(National_Reduction,Age_Reduction)

for ii=2:length(Age_Reduction)
    Age_Reduction(ii)=Age_Reduction(ii-1) & Age_Reduction(ii);
end
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

clearvars -except Model_Num Vaccine National_Reduction Age_Reduction

Age_Class={'Age_0_to_4','Age_5_to_9','Age_10_to_14','Age_15_to_19','Age_20_to_24'};
Age_Class=Age_Class(Age_Reduction);

FN_Age_Class={'Age_0_to_4','_Age_5_to_9','_Age_10_to_14','_Age_15_to_19','_Age_20_to_24'};
FN_Age_Class=FN_Age_Class(Age_Reduction);

dZ_County=zeros(3143,sum(Age_Reduction));
for aa=1:length(Age_Class)
    [~,dZ_County(:,aa)]=Age_Adjustment_Factor(Vaccine,Model_Num,Age_Class{aa});
end

[County_Info,~]=Decline_Adjustment_Factor(Vaccine,Model_Num,dZ_County,National_Reduction,Age_Reduction);


load([Vaccine '_Immunity.mat'],'County_Data')

County_Data_Vaccine_Reduction=County_Data;

clear County_Data

epsv1=0.93;
epsv2=0.97;

Immunity_0_to_4=epsv1.*County_Info.Vaccine_Uptake_0_to_4;
Vaccine_Uptake_0_to_4=County_Info.Vaccine_Uptake_0_to_4;

if(Age_Reduction(2))
    Immunity_5_to_9=epsv2.*County_Info.Vaccine_Uptake_5_to_9;
    Vaccine_Uptake_5_to_9=County_Info.Vaccine_Uptake_5_to_9;
else
    Immunity_5_to_9=table2array(County_Data_Vaccine_Reduction.Immunity(:,2));
    Vaccine_Uptake_5_to_9=table2array(County_Data_Vaccine_Reduction.Vaccine_Uptake(:,2));
end

if(Age_Reduction(3))
    Immunity_10_to_14=epsv2.*County_Info.Vaccine_Uptake_10_to_14;
    Vaccine_Uptake_10_to_14=County_Info.Vaccine_Uptake_10_to_14;
else
    Immunity_10_to_14=table2array(County_Data_Vaccine_Reduction.Immunity(:,3));
    Vaccine_Uptake_10_to_14=table2array(County_Data_Vaccine_Reduction.Vaccine_Uptake(:,3));
end

if(Age_Reduction(4))
    Immunity_15_to_19=epsv2.*County_Info.Vaccine_Uptake_15_to_19;
    Vaccine_Uptake_15_to_19=County_Info.Vaccine_Uptake_15_to_19;
else
    Immunity_15_to_19=table2array(County_Data_Vaccine_Reduction.Immunity(:,4));
    Vaccine_Uptake_15_to_19=table2array(County_Data_Vaccine_Reduction.Vaccine_Uptake(:,4));
end

if(Age_Reduction(5))
    Immunity_20_to_24=epsv2.*County_Info.Vaccine_Uptake_20_to_24;
    Vaccine_Uptake_20_to_24=County_Info.Vaccine_Uptake_20_to_24;
else
    Immunity_20_to_24=table2array(County_Data_Vaccine_Reduction.Immunity(:,5));
    Vaccine_Uptake_20_to_24=table2array(County_Data_Vaccine_Reduction.Vaccine_Uptake(:,5));
end

Vac_Susceptible=zeros(size(County_Data_Vaccine_Reduction.Immunity));
Unvac_Susceptible=1-table2array(County_Data_Vaccine_Reduction.Immunity);

Vac_Susceptible(:,1)=(1-epsv1).*Vaccine_Uptake_0_to_4;
Vac_Susceptible(:,2)=(1-epsv2).*Vaccine_Uptake_5_to_9;
Vac_Susceptible(:,3)=(1-epsv2).*Vaccine_Uptake_10_to_14;
Vac_Susceptible(:,4)=(1-epsv2).*Vaccine_Uptake_15_to_19;
Vac_Susceptible(:,5)=(1-epsv2).*Vaccine_Uptake_20_to_24;

for ss=6:18
    tf=table2array(County_Data_Vaccine_Reduction.Immunity(:,ss))<=epsv2;
    Vac_Susceptible(tf,ss)=table2array(County_Data_Vaccine_Reduction.Immunity(tf,ss)).*((1-epsv2)./epsv2);
end

Unvac_Susceptible(:,1)=(1-Vaccine_Uptake_0_to_4);
Unvac_Susceptible(:,2)=(1-Vaccine_Uptake_5_to_9);
Unvac_Susceptible(:,3)=(1-Vaccine_Uptake_10_to_14);
Unvac_Susceptible(:,4)=(1-Vaccine_Uptake_15_to_19);
Unvac_Susceptible(:,5)=(1-Vaccine_Uptake_20_to_24);

for ss=6:18
    tf=table2array(County_Data_Vaccine_Reduction.Immunity(:,ss))<=epsv2;
    Unvac_Susceptible(tf,ss)=1-table2array(County_Data_Vaccine_Reduction.Immunity(tf,ss))./epsv2;
end

County_Data_Vaccine_Reduction.Immunity(:,1)=array2table(Immunity_0_to_4);
County_Data_Vaccine_Reduction.Immunity(:,2)=array2table(Immunity_5_to_9);
County_Data_Vaccine_Reduction.Immunity(:,3)=array2table(Immunity_10_to_14);
County_Data_Vaccine_Reduction.Immunity(:,4)=array2table(Immunity_15_to_19);
County_Data_Vaccine_Reduction.Immunity(:,5)=array2table(Immunity_20_to_24);

County_Data_Vaccine_Reduction.Vaccine_Uptake(:,1)=array2table(Vaccine_Uptake_0_to_4);
County_Data_Vaccine_Reduction.Vaccine_Uptake(:,2)=array2table(Vaccine_Uptake_5_to_9);
County_Data_Vaccine_Reduction.Vaccine_Uptake(:,3)=array2table(Vaccine_Uptake_10_to_14);
County_Data_Vaccine_Reduction.Vaccine_Uptake(:,4)=array2table(Vaccine_Uptake_15_to_19);
County_Data_Vaccine_Reduction.Vaccine_Uptake(:,5)=array2table(Vaccine_Uptake_20_to_24);


County_Data_Vaccine_Reduction.Total_Immunity=sum(table2array(County_Data_Vaccine_Reduction.Population).*table2array(County_Data_Vaccine_Reduction.Immunity),2);

load('Baseline_Estimate_Measles_Incidence.mat','beta_j','beta_seed');

Final_Size_Est=zeros(length(County_Data_Vaccine_Reduction.State),1);
Proportion_Size_Age_Unvaccinated=zeros(length(County_Data_Vaccine_Reduction.State),18);
Proportion_Size_Age_Vaccinated=zeros(length(County_Data_Vaccine_Reduction.State),18);
R_eff=zeros(length(County_Data_Vaccine_Reduction.State),1);
R_eff_Seed=zeros(length(County_Data_Vaccine_Reduction.State),1);

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
    
    v_temp=Unvac_Susceptible(s_indx,:);
    Unvaccinated_Contacts(s_indx,:)=M*v_temp(:);
    % Equation 13
    % (https://pmc.ncbi.nlm.nih.gov/articles/PMC7088810/pdf/11538_2010_Article_9623.pdf)
    
    A_eff=beta_j.*repmat(S_Pop,18,1).*M./repmat(Pop,18,1);
    A_eff(repmat(Pop,18,1)==0)=0;
    R_eff(s_indx)=max(eig(A_eff));
    if(max(eig(A_eff))>1)            
        z=lsqnonlin(@(x) A_eff*x(:)+log(1-x(:)), 0.999.*ones(1,18),zeros(1,18),ones(1,18),[],[],[],[],[],opts);
        Final_Size_Est(s_indx)=S_Pop*z(:);
        w_vac=Vac_Susceptible(s_indx,:)./(Vac_Susceptible(s_indx,:)+Unvac_Susceptible(s_indx,:));
        Total_Cases=S_Pop*z(:);
        Proportion_Size_Age_Unvaccinated(s_indx,:)=S_Pop(:).*(z(:).*(1-w_vac(:)))./Total_Cases;
        Proportion_Size_Age_Vaccinated(s_indx,:)=S_Pop(:).*(z(:).*(w_vac(:)))./Total_Cases;
    else
        n_s=S_Pop./sum(S_Pop);
        w_vac=Vac_Susceptible(s_indx,:)./(Vac_Susceptible(s_indx,:)+Unvac_Susceptible(s_indx,:));
        Proportion_Size_Age_Unvaccinated(s_indx,:)=n_s(:).*(1-w_vac(:));
        Proportion_Size_Age_Vaccinated(s_indx,:)=n_s(:).*(w_vac(:));
    end
    
    A_eff=beta_seed.*repmat(S_Pop,18,1).*M./repmat(Pop,18,1);
    A_eff(repmat(Pop,18,1)==0)=0;
    R_eff_Seed(s_indx)=max(eig(A_eff));
end


County_Data_Vaccine_Reduction.R_eff_Seed=R_eff_Seed;
County_Data_Vaccine_Reduction.R_eff=R_eff;
County_Data_Vaccine_Reduction.Final_Size_Est=Final_Size_Est;
County_Data_Vaccine_Reduction.All_Contacts=All_Contacts;
County_Data_Vaccine_Reduction.Unvaccinated_Contacts=Unvaccinated_Contacts;
save(['National_Reduction=' num2str(100*National_Reduction) '_' FN_Age_Class{:} '.mat'],'County_Data_Vaccine_Reduction','Proportion_Size_Age_Unvaccinated','Proportion_Size_Age_Vaccinated')

end