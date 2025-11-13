clear;
clc;
parpool(48)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Immunity
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

Age_Class={'Age_0_to_4','Age_5_to_9','Age_10_to_14','Age_15_to_19','Age_20_to_24','Age_25_to_29','Age_30_to_34','Age_35_to_39','Age_40_to_44','Age_45_to_49','Age_50_to_54','Age_55_to_59','Age_60_to_64','Age_65_to_69','Age_70_to_74','Age_75_to_79','Age_80_to_84','Age_85_plus'};
V_Age_Class={'Age_0_to_4','Age_5_to_9','Age_10_to_14','Age_15_to_19','Age_20_to_24'};
T=readtable('Immunity_Age.xlsx');
epsv1=0.93;
epsv2=0.97;

for aa=1:length(Age_Class)
    if(ismember(Age_Class{aa},V_Age_Class))
        [temp_county,~]=Age_Adjustment_Factor(Vaccine,Model_Num,Age_Class{aa});
        if aa==1
            County_Data.State=temp_county.State;
            County_Data.County=temp_county.County;
            County_Data.GEOID=temp_county.GEOID;
            County_Data.Population=temp_county.Population;
            County_Data.Total_Population=temp_county.Total_Population;
            Immunity=epsv1.*temp_county.Vaccine_Uptake;
            Vaccine_Uptake=temp_county.Vaccine_Uptake;
        else
            Immunity=[Immunity epsv2.*temp_county.Vaccine_Uptake];
            Vaccine_Uptake=[Vaccine_Uptake temp_county.Vaccine_Uptake];
        end
    else
        temp_immunity=zeros(length(County_Data.County),1);
        T_temp=T(strcmp(T.Age_Group,Age_Class{aa}),:);
        for ss=1:height(T_temp)
            tf=strcmp(County_Data.State,T_temp.State{ss});
            temp_immunity(tf)=T_temp.Immunity(ss);
        end
        Immunity=[Immunity temp_immunity];
    end
end

Vaccine_Uptake=array2table(Vaccine_Uptake);
Vaccine_Uptake.Properties.VariableNames=V_Age_Class;
Immunity=array2table(Immunity);
Immunity.Properties.VariableNames=Age_Class;
County_Data.Vaccine_Uptake=Vaccine_Uptake;
County_Data.Immunity=Immunity;
County_Data.Total_Immunity=sum(table2array(County_Data.Population).*table2array(County_Data.Immunity),2);

beta_j=10.^(linspace(-1.1,log10(1.5),10^3));
Final_Size_Est=zeros(length(County_Data.State),length(beta_j));
opts=optimoptions('lsqnonlin','OptimalityTolerance',10^(-12),'StepTolerance',10^(-12),'FunctionTolerance',10^(-16),'MaxFunctionEvaluations',10^6,'MaxIterations',10^(6));
    
[b_indx,s_indx]=meshgrid([1:10^3],1:length(County_Data.State));
old_size=size(s_indx);
s_indx=s_indx(:);
b_indx=b_indx(:);
Final_Size_Est=Final_Size_Est(:);
R_eff=Final_Size_Est(:);
R_0=Final_Size_Est(:);

parfor ii=1:length(Final_Size_Est)
    State_temp=County_Data.State{s_indx(ii)};
    State_temp(State_temp==' ')='_';
    M=readtable([pwd '/State_Contact_Matrix/United_States_subnational_' State_temp '_M_overall_contact_matrix_18']);
    M=table2array(M);
    
    Tot_Pop=County_Data.Total_Population(s_indx(ii));
    Pop=table2array(County_Data.Population(s_indx(ii),:));
    Immunity=table2array(County_Data.Immunity(s_indx(ii),:));
    Pop=Tot_Pop.*Pop;
    S_Pop=(1-Immunity).*Pop; % Susceptible population

    % Equation 13
    % (https://pmc.ncbi.nlm.nih.gov/articles/PMC7088810/pdf/11538_2010_Article_9623.pdf)
    A_0=beta_j(b_indx(ii)).*repmat(Pop,18,1).*M./repmat(Pop,18,1);
    A_0(repmat(Pop,18,1)==0)=0;
    R_0(ii)=max(abs(eig(A_0)));
    
    A_eff=beta_j(b_indx(ii)).*repmat(S_Pop,18,1).*M./repmat(Pop,18,1);
    A_eff(repmat(Pop,18,1)==0)=0;
    R_eff(ii)=max(abs(eig(A_eff)));
    if(max(abs(eig(A_eff)))>1)            
        z=lsqnonlin(@(x) A_eff*x(:)+log(1-x(:)), 0.999.*ones(1,18),zeros(1,18),ones(1,18),[],[],[],[],[],opts);
        Final_Size_Est(ii)=S_Pop*z(:);
    end
end

Final_Size_Est=reshape(Final_Size_Est,old_size);
County_Data.Final_Size_Est=Final_Size_Est;
County_Data.beta_j=beta_j;
County_Data.R_eff=reshape(R_eff,old_size);
County_Data.R_0=reshape(R_0,old_size);

save([Vaccine '_Immunity.mat'],'County_Data')



