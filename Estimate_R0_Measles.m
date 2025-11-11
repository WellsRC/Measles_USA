clear;
clc;

load('Turncated_Negative_Binomial_Parameter.mat');
F_NB = scatteredInterpolant(kv(:),avg_fs(:),log(pv(:)./(1-pv(:))));

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

% [Importation_Cases_County_2024] = Case_Importation_Sample('Sample_2024',500);
% [Importation_Cases_County_2023] = Case_Importation_Sample('Sample_2023',500);

A=[];
X0=[-0.900585602899950	2.61176764317632	-5.75974820979949	-4.38964539471613	-0.416995137899064	-3.68399151986741	-0.0456643792169190; 
    -0.900585602899950	2.61176764317632	-5.75974820979949	-2.38964539471613	-0.416995137899064	-1.68399151986741	-0.0456643792169190;
    -0.919105557593447	2.44932597042403	-5.39855770747366	-2.27514071705288	-0.452726356477612	-3.69668885902970	-0.0438257808097951;
    -0.919105557593447	2.44932597042403	-5.39855770747366	-4.27514071705288	-0.452726356477612	-3.69668885902970	-0.0438257808097951;
    -1.03591091251131	-2.73390885648704	-3.60928486970385	-3.31907054831215	-0.817456549781539	-1.71044810642087	-0.569053623227815];

opts=optimoptions('surrogateopt','PlotFcn','surrogateoptplot','MaxFunctionEvaluations',10^3,'UseParallel',false,'InitialPoints',X0);
lb=[-1.1 -6 -6 -6 -3 -5 -1.1];
ub=[log10(1.5) 3 1  1  1 1 log10(1.5)];

rng(20251009)
r_samp_pc_2025=rand(length(Known_Ind_Cases),100);
r_samp_outbreak_2025=rand(length(Known_Ind_Cases),100);

[par_0,fval_0]=surrogateopt(@(x)Objective_Estimate_R0(x,County_Data,Imported_Case,Known_Ind_Cases,Unknown_Ind_Cases,Unknown_Ind_Cases_Weight,Population_i,Population_j,Distance_Matrix_ij,Week_Nat_Case_Count_2025,r_samp_pc_2025,r_samp_outbreak_2025,F_NB,Max_Outbreak),lb,ub,[],[],[],[],[],opts);

opts_ps=optimoptions('patternsearch','UseParallel',false,'FunctionTolerance',10^(-9),'MaxIterations',10^3,'MaxFunctionEvaluations',10^4,'PlotFcn','psplotbestf','UseCompleteSearch',true,'UseCompletePoll',true,'Cache','on');
[par,fval]=patternsearch(@(x)Objective_Estimate_R0(x,County_Data,Imported_Case,Known_Ind_Cases,Unknown_Ind_Cases,Unknown_Ind_Cases_Weight,Population_i,Population_j,Distance_Matrix_ij,Week_Nat_Case_Count_2025,r_samp_pc_2025,r_samp_outbreak_2025,F_NB,Max_Outbreak),par_0,[],[],[],[],lb,ub,[],opts_ps);

if(fval_0<fval)
    par=par_0;
end

beta_seed=10.^par(1);
lambda_0=par(2);
lambda_i=10.^par(3);
lambda_j=10.^par(4);
lambda_d=10.^par(5);
k_nbin=10.^par(6);
beta_j=10.^par(7);
k_mealses=0.23; 
save('Baseline_Estimate_Measles_Incidence.mat',"beta_seed","k_nbin","beta_j","lambda_0","lambda_i","lambda_j","lambda_d",'k_mealses');