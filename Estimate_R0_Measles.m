clear;
clc;
% parpool(24);
Vaccine='MMR';
load([Vaccine '_Immunity.mat'],'County_Data')

Max_Outbreak=County_Data.Total_Population(:).*(1-County_Data.Total_Immunity(:));
load('County_Matrix_Gravity_Covariates.mat',"Distance_Matrix_ij",'Population_i','Population_j')

load('Prior_log_Regression_Transmission.mat','prior_mean_transmission','transmission_bnds')
% Data as of December 30
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

Nat_Case_Count_2025=sum(Imported_Case)+sum(Known_Ind_Cases);

for indx=1:max(Measles_Cases.ID_Unknown)
    t_f=Measles_Cases.ID_Unknown==indx;

    Nat_Case_Count_2025=Nat_Case_Count_2025+unique(Measles_Cases.unkown_case_count(t_f));
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


[County_Transmission_X,RUCC_j] = Load_Transmission_Covariates(County_Data.GEOID);



[beta_j] = County_Transmission(prior_mean_transmission',County_Transmission_X.X);

[Case_Count,Reff]=Determine_Model_County_Case_Count(County_Data,beta_j);

t_f=strcmp(County_Data.County,'Gaines') & strcmp(County_Data.State,'Texas');

q_0=integral(@(x)nbinpdf(0,0.23,0.23./(0.23+Reff(t_f).*x)),0,1);

Import_Gaines=round(log(1-0.5)/log(q_0));



A=[];
load('Test_X0.mat','XN');
X0=[-0.682650734356208	-0.634984685061989	1.000000000000000	2.292915774596826	-3.885716365515917	1.349212745863400	1.689635478876136	-1.179180475574071	4.911991738864590	1.558763055599865	1.723791136646484	9.169578307464704	-2.406373602562711	56.513257866276874	9.019585705315873	-16.091994759357902	-1.674559337355726	23.993468450650575	8.671364425223091	-0.632640158570092	-20.523870128039437	7.792293135294987	-8.218379805515285	-11.513240556970029	-11.151701503897414	-10.735080698329510	-11.320892140350258	-11.305759531359817	-10.910912412765109	-11.125292001115980	-11.548730396773925	-10.908308911007506	-8.996586156804241;
    -0.745150734356208	-0.642125798343239	1.000000000000000	2.292915774596826	-3.885716365515917	1.349212745863400	1.689635478876136	-1.179241510730321	4.911991738864590	1.558396844662365	0.973302855396484	10.403953307464704	-2.406373602562711	56.512281303776874	9.019585705315873	-15.338088509357902	-1.674559337355726	25.993468450650575	8.854958175223091	-0.632640158570092	-20.523870128039437	7.792293135294987	-8.218379805515285	-11.385555010095029	-11.162931972647414	-10.735080698329510	-11.371185109100258	-11.225803476672317	-10.910912412765109	-11.029588876115980	-11.048730396773925	-10.900740551632506	-8.996586156804241];
X0=[X0;XN];
lb=[-5 log10(0.05) 1  transmission_bnds(1,:) -16];
ub=[1 log10(0.5)  30 transmission_bnds(2,:)  -6];


XS=[X0];

NS=10;
for xx=1:size(X0,1)
    Xt=X0(xx,:).*(1+0.01.*(0.5-rand(NS,length(lb))));
    Xt(:,3)=round(Xt(:,3)+(randi(3,size(Xt(:,3)))-2));
    Xt(Xt(:,3)<lb(3),3)=lb(3);
    Xt(Xt(:,3)>ub(3),3)=ub(3);
    XS=[XS;Xt];
end

% Bounds for parameters for Gravity model https://link.springer.com/article/10.1007/s42001-025-00414-7

rng(20251009)
r_samp_pc_2025=rand(length(Known_Ind_Cases),2500);
r_samp_outbreak_2025=rand(length(Known_Ind_Cases),2500);

Lt=zeros(size(XS,1),1);

parfor ii=1:size(XS,1)
    Lt(ii) = Objective_Estimate_R0(XS(ii,:), County_Data, Imported_Case, Known_Ind_Cases, Unknown_Ind_Cases, Unknown_Ind_Cases_Weight, Population_i,Population_j, Distance_Matrix_ij, Nat_Case_Count_2025, r_samp_pc_2025, r_samp_outbreak_2025, Max_Outbreak, County_Transmission_X.X,RUCC_j,Import_Gaines);
end
XS=XS(~isnan(Lt) & ~isinf(Lt),:);
opts=optimoptions('surrogateopt','PlotFcn','surrogateoptplot','MaxFunctionEvaluations',5.*10^3,'UseParallel',true,'InitialPoints',XS);


[par_0,fval_0,exitflag,output,trials]=surrogateopt(@(x)Objective_Estimate_R0(x,County_Data,Imported_Case,Known_Ind_Cases,Unknown_Ind_Cases,Unknown_Ind_Cases_Weight,Population_i,Population_j,Distance_Matrix_ij,Nat_Case_Count_2025,r_samp_pc_2025,r_samp_outbreak_2025,Max_Outbreak,County_Transmission_X.X,RUCC_j,Import_Gaines),lb,ub,[3],[],[],[],[],opts);

% Need to adjust since pattern search does not do integer constraints
lb(3)=lb(3)-0.499;
ub(3)=ub(3)+0.499;

opts_ps=optimoptions('patternsearch','UseParallel',true,'FunctionTolerance',10^(-9),'MaxIterations',10^3,'MaxFunctionEvaluations',10^4,'PlotFcn','psplotbestf','UseCompleteSearch',true,'UseCompletePoll',true,'Cache','on');
[par,fval]=patternsearch(@(x)Objective_Estimate_R0(x,County_Data,Imported_Case,Known_Ind_Cases,Unknown_Ind_Cases,Unknown_Ind_Cases_Weight,Population_i,Population_j,Distance_Matrix_ij,Nat_Case_Count_2025,r_samp_pc_2025,r_samp_outbreak_2025,Max_Outbreak,County_Transmission_X.X,RUCC_j,Import_Gaines),par_0,[],[],[],[],lb,ub,[],opts_ps);

if(fval_0<fval)
    par=par_0;
end

lambda_d=10.^par(1);
k_mealses=10.^par(2); 

R_NHG=round(par(3));

beta_transmission=par(4:32);

lambda_RUCC=10.^par(33).*ones(1,9); %par(33:41);
lambda_RUCC=lambda_RUCC(:);

[beta_j] = County_Transmission(beta_transmission,County_Transmission_X.X);

save('Baseline_Estimate_Measles_Incidence.mat',"R_NHG","lambda_d",'k_mealses','Import_Gaines','beta_j','beta_transmission','lambda_RUCC');