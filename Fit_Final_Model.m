clear;
clc;
Vaccine='MMR';
L_F=zeros(2^7,4);
L_S=zeros(2^7,4);
for Spatial_Valiation=1:4
    load([Vaccine '_Set_Spatial_Validation=' num2str(Spatial_Valiation) '.mat']);
    L_F(:,Spatial_Valiation)=L_fit;
    L_S(:,Spatial_Valiation)=L_spatial_val;
end

L_T=sum(L_F,2)+sum(L_S,2);

[~,Model_Num]=max(L_T);

[County_Data,State_Data] = Load_Data(Vaccine);

filter_county_train= ismember(County_Data.Spatial_Identifier,[1:4]);
filter_state_train= ismember(State_Data.Spatial_Identifier,[1:4]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
% Model index for fitting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
Model_Var=false(2^7,43);

Var_Names={'Age','Race','Income','Gini_Index','Education','Poverty_Income_Ratio','Rural_Urban_Code'};

Model_Var(:,[1 2 3])=true;

temp_var=false(7,43);
temp_var(1,4:8)=true;
temp_var(2,9:13)=true;
temp_var(3,14)=true;
temp_var(4,15)=true;
temp_var(5,16:22)=true;
temp_var(6,23:34)=true;
temp_var(7,35:43)=true;

bin_model=dec2bin([0:2^7-1]',7)=='1';
[~,srt_indx]=sort(sum(bin_model,2),'ascend');
bin_model=bin_model(srt_indx,:);

for ii=1:2^7
    temp_overall=false(1,43);
    for jj=1:7
        if(bin_model(ii,jj))
            temp_overall=temp_overall|temp_var(jj,:);
        end
    end
    Model_Var(ii,:)=Model_Var(ii,:) | temp_overall;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Run fitting and validation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
County_Data.Vaccine_Uptake(County_Data.Vaccine_Uptake==1)=1-10^(-4);
County_Data.Vaccine_Uptake(County_Data.Vaccine_Uptake==0)=10^(-4);

        
County_Data_model=County_Data;
County_Data_model.X=County_Data_model.X(:,Model_Var(Model_Num,:));
County_Data_model.X2=County_Data_model.X.^2;
if(Model_Var(Model_Num,end))
    County_Data_model.X2=County_Data_model.X2(:,3:end-9);
else
    County_Data_model.X2=County_Data_model.X2(:,3:end);
end
XT=[];
for ii=3:size(County_Data_model.X2,2)
    for jj=(ii+1):size(County_Data_model.X2,2)
        XT=[XT County_Data_model.X(:,ii).*County_Data_model.X(:,jj)];
    end
end

County_Data_model.XI=XT;   

num_par=size(County_Data_model.X,2)+size(County_Data_model.X2,2)+size(County_Data_model.XI,2);

t_fit=~isnan(County_Data_model.Vaccine_Uptake) & County_Data_model.Vaccine_Uptake>0 & County_Data_model.Vaccine_Uptake<1;
v_temp=County_Data_model.Vaccine_Uptake;
z_fit=log(v_temp(t_fit)./(1-v_temp(t_fit)));    


if(~isempty(County_Data_model.XI))

    lX=length(County_Data_model.X(1,:));
    mdl_temp=fitlm([County_Data_model.X(t_fit,:) County_Data_model.X2(t_fit,:)],z_fit);
    par_0=[mdl_temp.Coefficients.Estimate(1).*ones(10,3) repmat(mdl_temp.Coefficients.Estimate(2:(2+lX-1))',10,1) zeros(10,length(County_Data_model.XI(1,:))) repmat(mdl_temp.Coefficients.Estimate((2+lX):end)',10,1) flip(linspace(-1,3,10))' flip(linspace(-1,3,10))' (linspace(1.01,3,10))' (linspace(1.01,3,10))'];

    mdl_temp=fitlm([County_Data_model.X(t_fit,:) County_Data_model.XI(t_fit,:)],z_fit);
    par_0=[par_0; mdl_temp.Coefficients.Estimate(1).*ones(10,2) repmat(mdl_temp.Coefficients.Estimate',10,1) zeros(10,length(County_Data_model.X2(1,:))) flip(linspace(-1,3,10))' flip(linspace(-1,3,10))' (linspace(1.01,3,10))' (linspace(1.01,3,10))'];

    mdl_temp=fitlm([County_Data_model.X(t_fit,:)],z_fit);
    par_0=[par_0; mdl_temp.Coefficients.Estimate(1).*ones(10,2) repmat(mdl_temp.Coefficients.Estimate',10,1) zeros(10,length(County_Data_model.XI(1,:))) zeros(10,length(County_Data_model.X2(1,:))) flip(linspace(-1,3,10))' flip(linspace(-1,3,10))' (linspace(1.01,3,10))' (linspace(1.01,3,10))'];

    mdl_temp=fitlm([County_Data_model.X(t_fit,:) County_Data_model.XI(t_fit,:) County_Data_model.X2(t_fit,:)],z_fit);
    par_0=[par_0;mdl_temp.Coefficients.Estimate(1).*ones(10,2) repmat(mdl_temp.Coefficients.Estimate',10,1) flip(linspace(-1,3,10))' flip(linspace(-1,3,10))' (linspace(1.01,3,10))' (linspace(1.01,3,10))'];

elseif(~isempty(County_Data_model.X2))
    mdl_temp=fitlm([County_Data_model.X(t_fit,:) County_Data_model.X2(t_fit,:)],z_fit);
    par_0=[mdl_temp.Coefficients.Estimate(1).*ones(10,2) repmat(mdl_temp.Coefficients.Estimate',10,1) flip(linspace(-1,3,10))' flip(linspace(-1,3,10))' (linspace(1.01,3,10))' (linspace(1.01,3,10))'];
else
    mdl_temp=fitlm([County_Data_model.X(t_fit,:)],z_fit);
    par_0=[mdl_temp.Coefficients.Estimate(1).*ones(10,2) repmat(mdl_temp.Coefficients.Estimate',10,1) flip(linspace(-1,3,10))' flip(linspace(-1,3,10))' (linspace(1.01,3,10))' (linspace(1.01,3,10))'];
end


test=coefCI(mdl_temp)';

test=[test(:,1) test(:,1) test];

par_samp=[repmat(test(1,:),100,1)+repmat((test(2,:)-test(1,:)),100,1).*rand(100,length(test(1,:))) -1+4.*rand(100,2) 1.01+2.99.*rand(100,2)];

par_0=[par_0; par_samp];


lb=-5000.*ones(1,num_par+7); % 3 interept and 4 hyper paramters
lb(end-3:end-2)=-1;
lb(end-1:end)=1;
ub=5000.*ones(1,num_par+7);  % 3 interept and 4 hyper paramters
ub(end-3:end-2)=3;
ub(end-1:end)=3;

opts_ga=optimoptions("ga","PlotFcn",'gaplotbestf','UseParallel',true,"MaxGenerations",100,"FunctionTolerance",10^(-9),'CrossoverFcn','crossoverheuristic','MigrationInterval',25,'SelectionFcn',{@selectiontournament,8},'PopulationSize',250,'InitialPopulationMatrix',par_0);
[par_est]=ga(@(x)Objective_Function_Coverage(x,County_Data_model,State_Data,filter_county_train,filter_state_train),length(lb),[],[],[],[],lb,ub,[],[],opts_ga);

opts_pats=optimoptions('patternsearch','UseParallel',true,"PlotFcn",'psplotbestf','FunctionTolerance',10^(-12),'MaxIterations',1000,'StepTolerance',10^(-9),'MaxFunctionEvaluations',10^4,'Cache','on');
[par_final,f_final]=patternsearch(@(x)Objective_Function_Coverage(x,County_Data_model,State_Data,filter_county_train,filter_state_train),par_est,[],[],[],[],lb,ub,[],opts_pats);

model_par=par_final;
L_fit=-f_final;

save([Vaccine '_Model_' num2str(Model_Num) '.mat'],'model_par','L_fit','Var_Names');


