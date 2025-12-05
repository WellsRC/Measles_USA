function National_Reduction_Vaccine_Uptake(National_Annual_Reduction,Year_Reduced)

if(National_Annual_Reduction>0)
    Age_Group_Reduction=Compute_Age_Group_Reduction(National_Annual_Reduction,Year_Reduced);
elseif(National_Annual_Reduction==0)
    Age_Group_Reduction=zeros(1,3);
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

clearvars -except Model_Num Vaccine National_Annual_Reduction Year_Reduced Age_Group_Reduction

Age_Class={'Age_0_to_4','Age_5_to_9','Age_10_to_14'}; 


if(National_Annual_Reduction~=0)
    dZ_County=zeros(3143,length(Age_Class));
    for aa=1:length(Age_Class)
        [~,dZ_County(:,aa)]=Age_Adjustment_Factor(Vaccine,Model_Num,Age_Class{aa});
    end
    
    [County_Info,~]=Decline_Adjustment_Factor(Vaccine,Model_Num,dZ_County,Age_Group_Reduction);
end

load([Vaccine '_Immunity.mat'],'County_Data')

County_Data_Vaccine_Reduction=County_Data;

clear County_Data

epsv1=0.93;
epsv2=0.97;

if(National_Annual_Reduction~=0)
    Immunity_0_to_4=epsv1.*County_Info.Vaccine_Uptake_0_to_4.Overall;
    Vaccine_Uptake_0_to_4=County_Info.Vaccine_Uptake_0_to_4.Overall;
    Vaccine_Uptake_0_to_4_Uninsured=County_Info.Vaccine_Uptake_0_to_4.Uninsured;
    Vaccine_Uptake_0_to_4_Public=County_Info.Vaccine_Uptake_0_to_4.Public;
    Vaccine_Uptake_0_to_4_Private=County_Info.Vaccine_Uptake_0_to_4.Private;
    
    Immunity_5_to_9=epsv2.*County_Info.Vaccine_Uptake_5_to_9.Overall;
    Vaccine_Uptake_5_to_9=County_Info.Vaccine_Uptake_5_to_9.Overall;
    Vaccine_Uptake_5_to_9_Uninsured=County_Info.Vaccine_Uptake_5_to_9.Uninsured;
    Vaccine_Uptake_5_to_9_Public=County_Info.Vaccine_Uptake_5_to_9.Public;
    Vaccine_Uptake_5_to_9_Private=County_Info.Vaccine_Uptake_5_to_9.Private;
    
    Immunity_10_to_14=epsv2.*County_Info.Vaccine_Uptake_10_to_14.Overall;
    Vaccine_Uptake_10_to_14=County_Info.Vaccine_Uptake_10_to_14.Overall;
    Vaccine_Uptake_10_to_14_Uninsured=County_Info.Vaccine_Uptake_10_to_14.Uninsured;
    Vaccine_Uptake_10_to_14_Public=County_Info.Vaccine_Uptake_10_to_14.Public;
    Vaccine_Uptake_10_to_14_Private=County_Info.Vaccine_Uptake_10_to_14.Private;
    
else
    Immunity_0_to_4=table2array(County_Data_Vaccine_Reduction.Immunity(:,1));
    Vaccine_Uptake_0_to_4=table2array(County_Data_Vaccine_Reduction.Vaccine_Uptake(:,1));
    Vaccine_Uptake_0_to_4_Uninsured=(County_Data_Vaccine_Reduction.Vaccine_Uptake_Uninsured(:,1));
    Vaccine_Uptake_0_to_4_Public=(County_Data_Vaccine_Reduction.Vaccine_Uptake_Public(:,1));
    Vaccine_Uptake_0_to_4_Private=(County_Data_Vaccine_Reduction.Vaccine_Uptake_Private(:,1));
    
    Immunity_5_to_9=table2array(County_Data_Vaccine_Reduction.Immunity(:,2));
    Vaccine_Uptake_5_to_9=table2array(County_Data_Vaccine_Reduction.Vaccine_Uptake(:,2));
    Vaccine_Uptake_5_to_9_Uninsured=(County_Data_Vaccine_Reduction.Vaccine_Uptake_Uninsured(:,2));
    Vaccine_Uptake_5_to_9_Public=(County_Data_Vaccine_Reduction.Vaccine_Uptake_Public(:,2));
    Vaccine_Uptake_5_to_9_Private=(County_Data_Vaccine_Reduction.Vaccine_Uptake_Private(:,2));

    Immunity_10_to_14=table2array(County_Data_Vaccine_Reduction.Immunity(:,3));
    Vaccine_Uptake_10_to_14=table2array(County_Data_Vaccine_Reduction.Vaccine_Uptake(:,3));
    Vaccine_Uptake_10_to_14_Uninsured=(County_Data_Vaccine_Reduction.Vaccine_Uptake_Uninsured(:,3));
    Vaccine_Uptake_10_to_14_Public=(County_Data_Vaccine_Reduction.Vaccine_Uptake_Public(:,3));
    Vaccine_Uptake_10_to_14_Private=(County_Data_Vaccine_Reduction.Vaccine_Uptake_Private(:,3));
end


w_u=County_Data_Vaccine_Reduction.Uninsured./(County_Data_Vaccine_Reduction.Uninsured+County_Data_Vaccine_Reduction.Public+County_Data_Vaccine_Reduction.Private);
w_pu=County_Data_Vaccine_Reduction.Public./(County_Data_Vaccine_Reduction.Uninsured+County_Data_Vaccine_Reduction.Public+County_Data_Vaccine_Reduction.Private);
w_pr=County_Data_Vaccine_Reduction.Private./(County_Data_Vaccine_Reduction.Uninsured+County_Data_Vaccine_Reduction.Public+County_Data_Vaccine_Reduction.Private);

Vac_Susceptible=zeros(size(County_Data_Vaccine_Reduction.Immunity));
Unvac_Susceptible=1-table2array(County_Data_Vaccine_Reduction.Immunity);

Vac_Susceptible_Uninsured=zeros(size(County_Data_Vaccine_Reduction.Immunity));
Unvac_Susceptible_Uninsured=1-table2array(County_Data_Vaccine_Reduction.Immunity);

Vac_Susceptible_Public=zeros(size(County_Data_Vaccine_Reduction.Immunity));
Unvac_Susceptible_Public=1-table2array(County_Data_Vaccine_Reduction.Immunity);

Vac_Susceptible_Private=zeros(size(County_Data_Vaccine_Reduction.Immunity));
Unvac_Susceptible_Private=1-table2array(County_Data_Vaccine_Reduction.Immunity);

Vac_Susceptible(:,1)=(1-epsv1).*Vaccine_Uptake_0_to_4;
Vac_Susceptible(:,2)=(1-epsv2).*Vaccine_Uptake_5_to_9;
Vac_Susceptible(:,3)=(1-epsv2).*Vaccine_Uptake_10_to_14;

Vac_Susceptible_Uninsured(:,1)=w_u(:,1).*(1-epsv1).*Vaccine_Uptake_0_to_4_Uninsured;
Vac_Susceptible_Uninsured(:,2)=w_u(:,2).*(1-epsv2).*Vaccine_Uptake_5_to_9_Uninsured;
Vac_Susceptible_Uninsured(:,3)=w_u(:,3).*(1-epsv2).*Vaccine_Uptake_10_to_14_Uninsured;

Vac_Susceptible_Public(:,1)=w_pu(:,1).*(1-epsv1).*Vaccine_Uptake_0_to_4_Public;
Vac_Susceptible_Public(:,2)=w_pu(:,2).*(1-epsv2).*Vaccine_Uptake_5_to_9_Public;
Vac_Susceptible_Public(:,3)=w_pu(:,3).*(1-epsv2).*Vaccine_Uptake_10_to_14_Public;

Vac_Susceptible_Private(:,1)=w_pr(:,1).*(1-epsv1).*Vaccine_Uptake_0_to_4_Private;
Vac_Susceptible_Private(:,2)=w_pr(:,2).*(1-epsv2).*Vaccine_Uptake_5_to_9_Private;
Vac_Susceptible_Private(:,3)=w_pr(:,3).*(1-epsv2).*Vaccine_Uptake_10_to_14_Private;

Vac_Susceptible(:,4)=table2array(County_Data_Vaccine_Reduction.Vaccine_Uptake(:,4)).*((1-epsv2));

Vac_Susceptible_Uninsured(:,4)=w_u(:,4).*(County_Data_Vaccine_Reduction.Vaccine_Uptake_Uninsured(:,4)).*((1-epsv2));
Vac_Susceptible_Private(:,4)=w_pu(:,4).*(County_Data_Vaccine_Reduction.Vaccine_Uptake_Private(:,4)).*((1-epsv2));
Vac_Susceptible_Public(:,4)=w_pr(:,4).*(County_Data_Vaccine_Reduction.Vaccine_Uptake_Public(:,4)).*((1-epsv2));

Vac_Susceptible(:,5)=table2array(County_Data_Vaccine_Reduction.Vaccine_Uptake(:,5)).*((1-epsv2));

Vac_Susceptible_Uninsured(:,5)=w_u(:,5).*(County_Data_Vaccine_Reduction.Vaccine_Uptake_Uninsured(:,5)).*((1-epsv2));
Vac_Susceptible_Private(:,5)=w_pu(:,5).*(County_Data_Vaccine_Reduction.Vaccine_Uptake_Private(:,5)).*((1-epsv2));
Vac_Susceptible_Public(:,5)=w_pr(:,5).*(County_Data_Vaccine_Reduction.Vaccine_Uptake_Public(:,5)).*((1-epsv2));


Unvac_Susceptible(:,1)=(1-Vaccine_Uptake_0_to_4);
Unvac_Susceptible(:,2)=(1-Vaccine_Uptake_5_to_9);
Unvac_Susceptible(:,3)=(1-Vaccine_Uptake_10_to_14);
Unvac_Susceptible(:,4)=(1-table2array(County_Data_Vaccine_Reduction.Vaccine_Uptake(:,4)));
Unvac_Susceptible(:,5)=(1-table2array(County_Data_Vaccine_Reduction.Vaccine_Uptake(:,5)));

Unvac_Susceptible_Uninsured(:,1)=w_u(:,1).*(1-Vaccine_Uptake_0_to_4_Uninsured);
Unvac_Susceptible_Uninsured(:,2)=w_u(:,2).*(1-Vaccine_Uptake_5_to_9_Uninsured);
Unvac_Susceptible_Uninsured(:,3)=w_u(:,3).*(1-Vaccine_Uptake_10_to_14_Uninsured);
Unvac_Susceptible_Uninsured(:,4)=w_u(:,4).*(1-(County_Data_Vaccine_Reduction.Vaccine_Uptake_Uninsured(:,4)));
Unvac_Susceptible_Uninsured(:,5)=w_u(:,5).*(1-(County_Data_Vaccine_Reduction.Vaccine_Uptake_Uninsured(:,5)));

Unvac_Susceptible_Public(:,1)=w_pu(:,1).*(1-Vaccine_Uptake_0_to_4_Public);
Unvac_Susceptible_Public(:,2)=w_pu(:,2).*(1-Vaccine_Uptake_5_to_9_Public);
Unvac_Susceptible_Public(:,3)=w_pu(:,3).*(1-Vaccine_Uptake_10_to_14_Public);
Unvac_Susceptible_Public(:,4)=w_pu(:,4).*(1-(County_Data_Vaccine_Reduction.Vaccine_Uptake_Public(:,4)));
Unvac_Susceptible_Public(:,5)=w_pu(:,5).*(1-(County_Data_Vaccine_Reduction.Vaccine_Uptake_Public(:,5)));

Unvac_Susceptible_Private(:,1)=w_pr(:,1).*(1-Vaccine_Uptake_0_to_4_Private);
Unvac_Susceptible_Private(:,2)=w_pr(:,2).*(1-Vaccine_Uptake_5_to_9_Private);
Unvac_Susceptible_Private(:,3)=w_pr(:,3).*(1-Vaccine_Uptake_10_to_14_Private);
Unvac_Susceptible_Private(:,4)=w_pr(:,4).*(1-(County_Data_Vaccine_Reduction.Vaccine_Uptake_Private(:,4)));
Unvac_Susceptible_Private(:,5)=w_pr(:,5).*(1-(County_Data_Vaccine_Reduction.Vaccine_Uptake_Private(:,5)));

for ss=6:18
    
    s_temp=max(((table2array(County_Data_Vaccine_Reduction.Immunity(:,ss)))-epsv2.*table2array(County_Data_Vaccine_Reduction.Vaccine_Uptake(:,ss)))./(1-table2array(County_Data_Vaccine_Reduction.Vaccine_Uptake(:,ss))+(1-epsv2).*table2array(County_Data_Vaccine_Reduction.Vaccine_Uptake(:,ss))),0);
    Unvac_Susceptible(:,ss)=(1-s_temp).*(1-table2array(County_Data_Vaccine_Reduction.Vaccine_Uptake(:,ss)));

    Unvac_Susceptible_Uninsured(:,ss)=w_u(:,ss).*(1-s_temp).*(1-(County_Data_Vaccine_Reduction.Vaccine_Uptake_Uninsured(:,ss)));
    Unvac_Susceptible_Public(:,ss)=w_pu(:,ss).*(1-s_temp).*(1-(County_Data_Vaccine_Reduction.Vaccine_Uptake_Public(:,ss)));
    Unvac_Susceptible_Private(:,ss)=w_pr(:,ss).*(1-s_temp).*(1-(County_Data_Vaccine_Reduction.Vaccine_Uptake_Private(:,ss)));


    Vac_Susceptible(:,ss)=(1-s_temp).*table2array(County_Data_Vaccine_Reduction.Vaccine_Uptake(:,ss)).*((1-epsv2));

    Vac_Susceptible_Uninsured(:,ss)=w_u(:,ss).*(1-s_temp).*(County_Data_Vaccine_Reduction.Vaccine_Uptake_Uninsured(:,ss)).*((1-epsv2));
    Vac_Susceptible_Private(:,ss)=w_pu(:,ss).*(1-s_temp).*(County_Data_Vaccine_Reduction.Vaccine_Uptake_Private(:,ss)).*((1-epsv2));
    Vac_Susceptible_Public(:,ss)=w_pr(:,ss).*(1-s_temp).*(County_Data_Vaccine_Reduction.Vaccine_Uptake_Public(:,ss)).*((1-epsv2));
end

County_Data_Vaccine_Reduction.Immunity(:,1)=array2table(Immunity_0_to_4);
County_Data_Vaccine_Reduction.Immunity(:,2)=array2table(Immunity_5_to_9);
County_Data_Vaccine_Reduction.Immunity(:,3)=array2table(Immunity_10_to_14);

County_Data_Vaccine_Reduction.Vaccine_Uptake(:,1)=array2table(Vaccine_Uptake_0_to_4);
County_Data_Vaccine_Reduction.Vaccine_Uptake(:,2)=array2table(Vaccine_Uptake_5_to_9);
County_Data_Vaccine_Reduction.Vaccine_Uptake(:,3)=array2table(Vaccine_Uptake_10_to_14);


County_Data_Vaccine_Reduction.Vaccine_Uptake_Uninsured(:,1)=(Vaccine_Uptake_0_to_4_Uninsured);
County_Data_Vaccine_Reduction.Vaccine_Uptake_Uninsured(:,2)=(Vaccine_Uptake_5_to_9_Uninsured);
County_Data_Vaccine_Reduction.Vaccine_Uptake_Uninsured(:,3)=(Vaccine_Uptake_10_to_14_Uninsured);

County_Data_Vaccine_Reduction.Vaccine_Uptake_Public(:,1)=(Vaccine_Uptake_0_to_4_Public);
County_Data_Vaccine_Reduction.Vaccine_Uptake_Public(:,2)=(Vaccine_Uptake_5_to_9_Public);
County_Data_Vaccine_Reduction.Vaccine_Uptake_Public(:,3)=(Vaccine_Uptake_10_to_14_Public);

County_Data_Vaccine_Reduction.Vaccine_Uptake_Private(:,1)=(Vaccine_Uptake_0_to_4_Private);
County_Data_Vaccine_Reduction.Vaccine_Uptake_Private(:,2)=(Vaccine_Uptake_5_to_9_Private);
County_Data_Vaccine_Reduction.Vaccine_Uptake_Private(:,3)=(Vaccine_Uptake_10_to_14_Private);


County_Data_Vaccine_Reduction.Total_Immunity=sum(table2array(County_Data_Vaccine_Reduction.Population).*table2array(County_Data_Vaccine_Reduction.Immunity),2);

load('Baseline_Estimate_Measles_Incidence.mat','beta_j');

Final_Size_Est=zeros(length(County_Data_Vaccine_Reduction.State),1);
Proportion_Size_Age_Unvaccinated=zeros(length(County_Data_Vaccine_Reduction.State),18);
Proportion_Size_Age_Vaccinated=zeros(length(County_Data_Vaccine_Reduction.State),18);

Proportion_Age_Unvaccinated_Uninsured=zeros(length(County_Data_Vaccine_Reduction.State),18);
Proportion_Age_Vaccinated_Uninsured=zeros(length(County_Data_Vaccine_Reduction.State),18);

Proportion_Age_Unvaccinated_Public=zeros(length(County_Data_Vaccine_Reduction.State),18);
Proportion_Age_Vaccinated_Public=zeros(length(County_Data_Vaccine_Reduction.State),18);

Proportion_Age_Unvaccinated_Private=zeros(length(County_Data_Vaccine_Reduction.State),18);
Proportion_Age_Vaccinated_Private=zeros(length(County_Data_Vaccine_Reduction.State),18);

R_eff=zeros(length(County_Data_Vaccine_Reduction.State),1);

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
    
    A_eff=beta_j(s_indx).*repmat(S_Pop,18,1).*M./repmat(Pop,18,1);
    A_eff(repmat(Pop,18,1)==0)=0;
    R_eff(s_indx)=max(abs(eig(A_eff)));
    if(max(abs(eig(A_eff)))>1)            
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

    w_ins=Unvac_Susceptible_Uninsured(s_indx,:)./(Unvac_Susceptible_Uninsured(s_indx,:)+Unvac_Susceptible_Private(s_indx,:)+Unvac_Susceptible_Public(s_indx,:));
    Proportion_Age_Unvaccinated_Uninsured(s_indx,:)=w_ins;
    w_ins=Unvac_Susceptible_Public(s_indx,:)./(Unvac_Susceptible_Uninsured(s_indx,:)+Unvac_Susceptible_Private(s_indx,:)+Unvac_Susceptible_Public(s_indx,:));
    Proportion_Age_Unvaccinated_Public(s_indx,:)=w_ins;
    Proportion_Age_Unvaccinated_Private(s_indx,:)=1-Proportion_Age_Unvaccinated_Uninsured(s_indx,:)-Proportion_Age_Unvaccinated_Public(s_indx,:);

    w_ins=Vac_Susceptible_Uninsured(s_indx,:)./(Vac_Susceptible_Uninsured(s_indx,:)+Vac_Susceptible_Private(s_indx,:)+Vac_Susceptible_Public(s_indx,:));
    Proportion_Age_Vaccinated_Uninsured(s_indx,:)=w_ins;
    w_ins=Vac_Susceptible_Public(s_indx,:)./(Vac_Susceptible_Uninsured(s_indx,:)+Vac_Susceptible_Private(s_indx,:)+Vac_Susceptible_Public(s_indx,:));
    Proportion_Age_Vaccinated_Public(s_indx,:)=w_ins;
    Proportion_Age_Vaccinated_Private(s_indx,:)=1-Proportion_Age_Vaccinated_Uninsured(s_indx,:)-Proportion_Age_Vaccinated_Public(s_indx,:);
end

County_Data_Vaccine_Reduction.R_eff=R_eff;
County_Data_Vaccine_Reduction.Final_Size_Est=Final_Size_Est;
County_Data_Vaccine_Reduction.All_Contacts=All_Contacts;
County_Data_Vaccine_Reduction.Unvaccinated_Contacts=Unvaccinated_Contacts;

save(['National_Reduction=' num2str(100*National_Annual_Reduction) '_Year=' num2str(Year_Reduced) '.mat'],'County_Data_Vaccine_Reduction','Proportion_Size_Age_Unvaccinated','Proportion_Size_Age_Vaccinated','Proportion_Age_Unvaccinated_Uninsured','Proportion_Age_Unvaccinated_Public','Proportion_Age_Unvaccinated_Private','Proportion_Age_Vaccinated_Uninsured','Proportion_Age_Vaccinated_Public','Proportion_Age_Vaccinated_Private')

end