clear;
clc;
parpool(10);
T=readtable('National_Measles_Cases_Weekly_2023_2025.xlsx');
Week_Nat_Case_Count_2025=T.total_cases(T.Year==2025);
Vaccine='MMR';
load([Vaccine '_Immunity.mat'],'County_Data')

Max_Outbreak=County_Data.Total_Population(:).*(1-County_Data.Total_Immunity(:));
load('County_Matrix_Gravity_Covariates.mat',"Distance_Matrix_ij",'Population_j','Population_i','County_GEOID')

indx_G=zeros(length(County_Data.County),1);
for cc=1:length(County_Data.County)
   indx_G(cc)=find(strcmp(County_GEOID,County_Data.GEOID{cc}));
end

Distance_Matrix_ij=Distance_Matrix_ij(indx_G,:);
Distance_Matrix_ij=Distance_Matrix_ij(:,indx_G);

Population_j=Population_j(indx_G,:);
Population_j=Population_j(:,indx_G);

Population_i=Population_i(indx_G,:);
Population_i=Population_i(:,indx_G);

Measles_Cases=readtable('County_Level_Measles_Cases_Adjusted.csv');

Imported_Case=zeros(length(County_Data.County),1);

Known_Ind_Cases=zeros(length(County_Data.County),1);
Unknown_Ind_Cases=NaN.*zeros(length(County_Data.County),2);
Unknown_Ind_Cases_Weight=NaN.*zeros(length(County_Data.County),2);

% Known imported cases
for cc=1:length(Imported_Case)
    t_f=str2double(County_Data.GEOID{cc})==Measles_Cases.GEOID & strcmp(Measles_Cases.type,'imported') & ~isnan(Measles_Cases.case_count);
    if(sum(t_f)>0)
        Imported_Case(cc)=Measles_Cases.case_count(t_f);
    end

    t_f=str2double(County_Data.GEOID{cc})==Measles_Cases.GEOID & strcmp(Measles_Cases.type,'local') & ~isnan(Measles_Cases.case_count);
    if(sum(t_f)>0)
        Known_Ind_Cases(cc)=Measles_Cases.case_count(t_f);
    end
end

for indx=1:max(Measles_Cases.ID_Unknown)
    t_f=Measles_Cases.ID_Unknown==indx;

    t_county=ismember(str2double(County_Data.GEOID),Measles_Cases.GEOID(t_f));

    w_pop=County_Data.Total_Population(t_county)-County_Data.Total_Immunity(t_county);
    w_pop=w_pop./sum(w_pop);

    imp_pop=County_Data.Total_Population(t_county)./sum(County_Data.Total_Population(t_county));
    if(sum(~isnan(Unknown_Ind_Cases(t_county,1)))>1)
        sc_unkown=2;
    else
        sc_unkown=1;
    end
    if(ismember('local',Measles_Cases.type(t_f)))
        Unknown_Ind_Cases(t_county,sc_unkown)=Measles_Cases.unkown_case_count(t_f);
        Unknown_Ind_Cases_Weight(t_county,sc_unkown)=w_pop;
    else
        Imported_Case(t_county)=Imported_Case(t_county)+imp_pop.*Measles_Cases.unkown_case_count(t_f);
    end
end



X_Samp=[];
L_Samp=0;
load('Explore_Beta_Relation_Parameters.mat')

[X_Samp,indx_s]=unique(X_Samp,"rows");
L_Samp=L_Samp(indx_s);



A=[];
X0=[];


% Bounds for parameters for Gravity model https://link.springer.com/article/10.1007/s42001-025-00414-7

opts=optimoptions('surrogateopt','PlotFcn',[],'MaxFunctionEvaluations',10^3,'UseParallel',false,'InitialPoints',X0);
lb=[ -6 -9 -9 -9 log10(0.05) -7 1 0 0 0];
ub=[3 0  0  log10(3)  log10(0.5) -3 length(L_Samp) 30 5 5];

rng(20251009)
r_samp_pc_2025=rand(length(Known_Ind_Cases),100);
r_samp_outbreak_2025=rand(length(Known_Ind_Cases),100);

[par_0,fval_0]=surrogateopt(@(x)Objective_Estimate_R0(x,County_Data,Imported_Case,Known_Ind_Cases,Unknown_Ind_Cases,Unknown_Ind_Cases_Weight,Population_i,Population_j,Distance_Matrix_ij,Week_Nat_Case_Count_2025,r_samp_pc_2025,r_samp_outbreak_2025,Max_Outbreak,X_Samp,L_Samp),lb,ub,[6 7 8],[],[],[],[],opts);

% Need to adjust since pattern search does not do integer constraints
lb(end-3:end)=lb(end-3:end)-0.499;
ub(end-3:end)=ub(end-3:end)+0.499;

opts_ps=optimoptions('patternsearch','UseParallel',false,'FunctionTolerance',10^(-9),'MaxIterations',10^3,'MaxFunctionEvaluations',10^4,'PlotFcn',[],'UseCompleteSearch',true,'UseCompletePoll',true,'Cache','on');
[par,fval]=patternsearch(@(x)Objective_Estimate_R0(x,County_Data,Imported_Case,Known_Ind_Cases,Unknown_Ind_Cases,Unknown_Ind_Cases_Weight,Population_i,Population_j,Distance_Matrix_ij,Week_Nat_Case_Count_2025,r_samp_pc_2025,r_samp_outbreak_2025,Max_Outbreak,X_Samp,L_Samp),par_0,[],[],[],[],lb,ub,[],opts_ps);

if(fval_0<fval)
    par=par_0;
end
lambda_0=-10.^par(1);
lambda_i=10.^par(2);
lambda_j=10.^par(3);
lambda_d=10.^par(4);
k_mealses=10.^par(5); 
lambda_out=10.^par(6); 


indx_beta=round(par(7));
beta_j=Transmission_Relation(1-County_Data.Total_Immunity,X_Samp(indx_beta,1:4));

R_NHG=round(par(8));
Import_Gaines=round(par(9));
Import_Kansas=round(par(10));

save('Baseline_Estimate_Measles_Incidence.mat','lambda_out',"R_NHG","lambda_0","lambda_i","lambda_j","lambda_d",'k_mealses','Import_Gaines','Import_Kansas','beta_j','indx_beta');