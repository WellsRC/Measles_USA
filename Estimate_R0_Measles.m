clear;
clc;
% parpool(24);
Vaccine='MMR';
load([Vaccine '_Immunity.mat'],'County_Data')

Max_Outbreak=County_Data.Total_Population(:).*(1-County_Data.Total_Immunity(:));
load('County_Matrix_Gravity_Covariates.mat',"Distance_Matrix_ij",'Population_j','Population_i','County_GEOID')

load('Prior_log_Regression_Transmission.mat','prior_mean_transmission','prior_covariance_transmission','transmission_bnds')
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
r1=log10(0.1+2.9.*rand(50,1));
X0=[ -7.509357645172464	-0.671900491142074	1.000000000000000	2.382228572256517	-3.992938004280498	1.038313605948621	1.447321814202811	-1.356604596770767	4.627256585738012	1.384426569213439	1.849183193965266	8.970985805067208	-1.062541489401979	55.956262530004757	9.504328740938426	-16.289161950564722	-3.664324358869183	20.767182276344538	9.308744986245594	-0.826455742962981	-20.159275604606936	7.601866462353867	-8.439388100948729	-11.523649761995596	-11.420230005687026	-10.875662739990188	-11.369685109942699	-11.107570842565682	-11.010864624844476	-11.321766299037595	-11.263185198214151	-11.400073942579734	-0.079906114293943	-0.118784110739932	-1.699355073087321	-0.663507810648139	-0.313236195282526	-1.987015296897200	-0.055476699328921	-2.639471288884554	-1.520089299947547
    -7.009357645172464	-0.601587991142074	1.000000000000000	2.382228572256517	-3.992938004280498	1.038313605948621	1.447321814202811	-1.356604596770767	4.627256585738012	1.384426569213439	1.880433193965266	8.970985805067208	-1.156291489401979	55.956262530004757	9.504328740938426	-16.289161950564722	-3.789324358869183	20.767182276344538	9.308744986245594	-0.826455742962981	-20.159275604606936	7.601866462353867	-8.439388100948729	-11.523649761995596	-11.420230005687026	-10.875662739990188	-11.369685109942699	-10.857570842565682	-11.010864624844476	-11.196766299037595	-11.138185198214151	-11.400073942579734	-3.079906114293943	-2.243784110739932	-3.199355073087321	-2.163507810648139	-3.063236195282526	-2.987015296897200	-2.805476699328921	-2.889471288884554	-2.645089299947547
    -8.973430561225456	-0.658750861907905	1.000000000000000	2.371085282501771	-4.004700277203306	1.012052722417525	1.441720009506884	-1.357960886090466	4.607314705045557	1.245435912334789	1.425685759219667	9.341979122722583	-1.065429324831885	55.959785573285778	9.733907098205830	-15.984189985691822	-3.652761002607996	27.057716952277527	9.821538369366044	-0.816035262922418	-20.108402750655383	7.461632261655628	-8.444223810337862	-11.498017296843088	-11.219323464731801	-10.806413349550006	-11.357814883707455	-11.397659732005700	-10.995018942722920	-11.410246852529314	-11.034425430951124	-11.176142741719314	-0.003191391623338	-0.055057638888049	-1.141804093148454	-0.337209440690931	-0.007926563747852	-0.869763574099851	-0.107084797348463	-1.794361546837253	-0.007402372390477];

lb=[-9 log10(0.05) 1  transmission_bnds(1,:)  -6.*ones(1,9)];
ub=[-5 log10(0.5)  30 transmission_bnds(2,:)  0.*ones(1,9)];

NS=1000;
XS=X0;
for xx=1:size(X0,1)
    Xt=X0(xx,:).*(1+0.05.*(0.5-rand(NS,length(lb))));
    Xt(:,3)=round(Xt(:,3)+(randi(3,size(Xt(:,3)))-2));
    Xt(Xt(:,3)<lb(3),3)=lb(3);
    Xt(Xt(:,3)>ub(3),3)=ub(3);
    XS=[XS;Xt];
end

% Bounds for parameters for Gravity model https://link.springer.com/article/10.1007/s42001-025-00414-7

rng(20251009)
r_samp_pc_2025=rand(length(Known_Ind_Cases),100);
r_samp_outbreak_2025=rand(length(Known_Ind_Cases),100);

Lt=zeros(size(XS,1),1);
parfor ii=1:size(XS,1)
    Lt(ii) = Objective_Estimate_R0(XS(ii,:), County_Data, Imported_Case, Known_Ind_Cases, Unknown_Ind_Cases, Unknown_Ind_Cases_Weight, Population_i, Distance_Matrix_ij, Nat_Case_Count_2025, r_samp_pc_2025, r_samp_outbreak_2025, Max_Outbreak, County_Transmission_X.X,RUCC_j,Import_Gaines);
end
XS=XS(~isnan(Lt) & ~isinf(Lt),:);
opts=optimoptions('surrogateopt','PlotFcn','surrogateoptplot','MaxFunctionEvaluations',5.*10^3,'UseParallel',true,'InitialPoints',XS);


[par_0,fval_0]=surrogateopt(@(x)Objective_Estimate_R0(x,County_Data,Imported_Case,Known_Ind_Cases,Unknown_Ind_Cases,Unknown_Ind_Cases_Weight,Population_i,Distance_Matrix_ij,Nat_Case_Count_2025,r_samp_pc_2025,r_samp_outbreak_2025,Max_Outbreak,County_Transmission_X.X,RUCC_j,Import_Gaines),lb,ub,[3],[],[],[],[],opts);

% Need to adjust since pattern search does not do integer constraints
lb(3)=lb(3)-0.499;
ub(3)=ub(3)+0.499;

opts_ps=optimoptions('patternsearch','UseParallel',true,'FunctionTolerance',10^(-9),'MaxIterations',10^3,'MaxFunctionEvaluations',10^4,'PlotFcn','psplotbestf','UseCompleteSearch',true,'UseCompletePoll',true,'Cache','on');
[par,fval]=patternsearch(@(x)Objective_Estimate_R0(x,County_Data,Imported_Case,Known_Ind_Cases,Unknown_Ind_Cases,Unknown_Ind_Cases_Weight,Population_i,Distance_Matrix_ij,Nat_Case_Count_2025,r_samp_pc_2025,r_samp_outbreak_2025,Max_Outbreak,County_Transmission_X.X,RUCC_j,Import_Gaines),par_0,[],[],[],[],lb,ub,[],opts_ps);

if(fval_0<fval)
    par=par_0;
end

lambda_d=10.^par(1);
k_mealses=10.^par(2); 

R_NHG=round(par(3));

beta_transmission=par(4:32);

lambda_RUCC=par(33:41);
lambda_RUCC=lambda_RUCC(:);

[beta_j] = County_Transmission(beta_transmission,County_Transmission_X.X);

save('Baseline_Estimate_Measles_Incidence.mat',"R_NHG","lambda_d",'k_mealses','Import_Gaines','beta_j','beta_transmission','lambda_RUCC');