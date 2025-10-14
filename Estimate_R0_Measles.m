clear;
clc;
% Estimate additional cases for 2025
X=[2	1	1	11	19	27	65	85	75	90	84	105	89	115	74	86	48	49 35	25	55	41	24	20	29	35	18	15	12	17	11	11	20	37	32	18 32	15	19	12	4];
opts=optimoptions('ga','UseParallel',true,'MaxGenerations',10^4,'FunctionTolerance',10^(-9));
[par_nat,fval]=ga(@(x) -sum(log(nbinpdf(X,x(5).*ones(size(1:length(X))),x(5)./(x(5)+x(3).*gampdf([1:length(X)]-x(4),x(1),x(2)))))),5,[],[],[],[],[],[],[],[],opts);
[par_nat,fval]=fmincon(@(x) -sum(log(nbinpdf(X,x(5).*ones(size(1:length(X))),x(5)./(x(5)+x(3).*gampdf([1:length(X)]-x(4),x(1),x(2)))))),par_nat);
Ex_2025=round(par_nat(3).*(gamcdf(52-par_nat(4),par_nat(1),par_nat(2))-gamcdf(41-par_nat(4),par_nat(1),par_nat(2))));
Nat_Case_Count=[49
121
59
285
1563+Ex_2025]; % Estimated an additional 56 cases for 2025
Vaccine='MMR';
load([Vaccine '_Immunity.mat'],'County_Data')

load('County_Matrix_Gravity_Covariate.mat',"Gravity_Model",'County_GEOID');

indx_G=zeros(length(County_Data.County),1);
for cc=1:length(County_Data.County)
   indx_G(cc)=find(strcmp(County_GEOID,County_Data.GEOID{cc}));
end

G=Gravity_Model(indx_G,:);
G=G(:,indx_G);

load('Meales_Epi_Weight_2025.mat')
epi_weight=ones(length(County_Data.County),1);
for ss=1:length(Unique_State)
    tf=strcmp(Unique_State{ss},County_Data.State);
    epi_weight(tf)=weights_epi_Week(ss);
end

epi_weight(epi_weight==0)=1;

Measles_Cases=readtable('County_Level_Measles_Cases_Adjusted.csv');

Imported_Case=zeros(length(County_Data.County),1);

Known_Ind_Cases=zeros(length(County_Data.County),1);
Unknown_Ind_Cases=NaN.*zeros(length(County_Data.County),1);
Unknown_Ind_Cases_Weight=NaN.*zeros(length(County_Data.County),1);

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

    w_pop=County_Data.Total_Immunity(t_county);
    w_pop=w_pop./sum(w_pop);
    if(ismember('local',Measles_Cases.type(t_f)))
        Unknown_Ind_Cases(t_county)=Measles_Cases.unkown_case_count(t_f);
        Unknown_Ind_Cases_Weight(t_county)=w_pop;
    else
        Imported_Case(t_county)=Imported_Case(t_county)+w_pop.*Measles_Cases.unkown_case_count(t_f);
    end
end

% X=[-0.771984187426246	-2.79619300864409	-1.35704097440345	-0.868815640301847;
%    -0.780492869626304	-2.86929076365012	-1.17423174257773	-0.827384134109794
%    -0.827441614212315	-2.83704207090879	-1.18325755894387	-0.829283651370020];
A=[-1 0 0 1];
opts=optimoptions('surrogateopt','PlotFcn','surrogateoptplot','MaxFunctionEvaluations',5.*10^2,'UseParallel',false,'InitialPoints',[]);
[par_0,fval_0]=surrogateopt(@(x)Objective_Estimate_R0(x,County_Data,Imported_Case,Known_Ind_Cases,Unknown_Ind_Cases,Unknown_Ind_Cases_Weight,G,Nat_Case_Count,epi_weight),[-1.1 -3.5 -3.5 -1.1],[log10(1.5) -0.5 -0.5 log10(1.5)],[],A,0,[],[],opts);

opts_ps=optimoptions('fmincon','FunctionTolerance',10^(-6),'MaxFunctionEvaluations',10^3,'MaxIterations',10^3,'PlotFcn','optimplotfval','UseParallel',true);
[par,fval]=fmincon(@(x)Objective_Estimate_R0(x,County_Data,Imported_Case,Known_Ind_Cases,Unknown_Ind_Cases,Unknown_Ind_Cases_Weight,G,Nat_Case_Count,epi_weight),par_0,A,0,[],[],[-1.1 -3.5 -3.5 -1.1],[log10(1.5) -0.5 -0.5 log10(1.5)],[],opts_ps);


% lb_BM=[-0.95 -2.95 -1.45 -0.95];
% ub_BM=[-0.7 -2.70 -1.05 -0.7];
% r=lhsdesign(10^6,4);
% x_samp=repmat(lb_BM(1:3),10^6,1)+repmat(ub_BM(1:3)-lb_BM(1:3),10^6,1).*r(:,1:3);
% x_samp=[x_samp lb_BM(4)+(x_samp(:,1)-lb_BM(4)).*r(:,4)];
% 
% L=zeros(10^6,1);
% for jj=1:10^6
%     L(jj)=-Objective_Estimate_R0(x_samp(jj,:),County_Data,Imported_Case,Known_Ind_Cases,Unknown_Ind_Cases,Unknown_Ind_Cases_Weight,G,Nat_Case_Count,epi_weight);
% end

psi_c=County_Data.Total_Immunity;

beta_seed= 10.^par(1);
lambda_g= -10.^par(2);
k_nbin=10.^par(3);
beta_j=10.^par(4);
k_mealses=0.23; 
save('Baseline_Estimate_Measles_Incidence.mat',"beta_seed","k_nbin","beta_j","lambda_g",'k_mealses');