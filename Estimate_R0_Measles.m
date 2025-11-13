clear;
clc;
parpool(10);

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

A=[];
X0=[0.0990970910821299	-4.35381196369646	0.0143402710041440	-1.81826367097497	-3	-2.47434793088613	-5.76929090582909 log10(0.85) 1 0;
    0.0990970910821299	-4.35381196369646	0.0143402710041440	-1.81826367097497	-3	-2.47434793088613	-5.76929090582909 log10(0.85) 0 0;
    0.0990970910821299	-4.35381196369646	0.0143402710041440	-1.81826367097497	-3	-2.47434793088613	-5.76929090582909 log10(0.85) 0 1;
    0.0990970910821299	-4.35381196369646	0.0143402710041440	-1.81826367097497	-3	-2.47434793088613	-5.76929090582909 log10(0.85) 2 0;
    0.0990970910821299	-4.35381196369646	0.0143402710041440	-1.81826367097497	-3	-2.47434793088613	-5.76929090582909 log10(0.85) 0 2;
    0.0990970910821299	-4.35381196369646	0.0143402710041440	-1.81826367097497	-3	-2.47434793088613	-5.76929090582909 log10(0.85) 1 1;
    0.0990970910821299	-4.35381196369646	0.0143402710041440	-1.81826367097497	-3	-2.47434793088613	-5.76929090582909 log10(0.85) 1 2;
    0.0990970910821299	-4.35381196369646	0.0143402710041440	-1.81826367097497	-3	-2.47434793088613	-5.76929090582909 log10(0.85) 2 1;
    0.0990970910821299	-4.35381196369646	0.0143402710041440	-1.81826367097497	-3	-2.47434793088613	-5.76929090582909 log10(0.85) 2 2;];

opts=optimoptions('surrogateopt','PlotFcn',[],'MaxFunctionEvaluations',10^3,'UseParallel',false,'InitialPoints',X0);
lb=[-1.1      -6 -6 -6 -6 -5 -8 log10(0.65) 0 0];
ub=[log10(1.5) 3 2  2  2 2   0   log10(1.35) 2 2];

rng(20251009)
r_samp_pc_2025=rand(length(Known_Ind_Cases),100);
r_samp_outbreak_2025=rand(length(Known_Ind_Cases),100);

[par_0,fval_0]=surrogateopt(@(x)Objective_Estimate_R0(x,County_Data,Imported_Case,Known_Ind_Cases,Unknown_Ind_Cases,Unknown_Ind_Cases_Weight,Population_i,Population_j,Distance_Matrix_ij,Week_Nat_Case_Count_2025,r_samp_pc_2025,r_samp_outbreak_2025,F_NB,Max_Outbreak),lb,ub,[9 10],[],[],[],[],opts);

% Need to adjust since pattern search does not do integer constraints
lb(end-1:end)=lb(end-1:end)-0.5;
ub(end-1:end)=ub(end-1:end)+0.5;

opts_ps=optimoptions('patternsearch','UseParallel',false,'FunctionTolerance',10^(-9),'MaxIterations',10^3,'MaxFunctionEvaluations',10^4,'PlotFcn',[],'UseCompleteSearch',true,'UseCompletePoll',true,'Cache','on');
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
k_mealses=10.^par(7); 
beta_j=10.^par(8);
Import_Gaines=round(par(9));
Import_Kansas=round(par(10));
save('Baseline_Estimate_Measles_Incidence.mat',"beta_seed","k_nbin","beta_j","lambda_0","lambda_i","lambda_j","lambda_d",'k_mealses','Import_Gaines','Import_Kansas');