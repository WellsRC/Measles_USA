function  Refine_Fit_Model(Vaccine,Spatial_Validation)

Spatial_Set=[1:4];
Spatial_Training=Spatial_Set(~ismember(Spatial_Set,Spatial_Validation)); 

[County_Data,State_Data,NE_Data] = Load_Data(Vaccine);

filter_county_train= ismember(County_Data.Spatial_Identifier,Spatial_Training);
filter_state_train= ismember(State_Data.Spatial_Identifier,Spatial_Training);
filter_NE_train= ismember(NE_Data.Spatial_Identifier,Spatial_Training);

filter_county_spatial_validation= ismember(County_Data.Spatial_Identifier,Spatial_Validation);
filter_state_spatial_validation= ismember(State_Data.Spatial_Identifier,Spatial_Validation);
filter_NE_spatial_validation= ismember(NE_Data.Spatial_Identifier,Spatial_Validation);

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
    load([ Vaccine '_Set_Spatial_Validation=' num2str(Spatial_Validation) '.mat'],'model_par');
    filter_paramters=false(2^7,571);
    
    parfor model_num=1:2^7
        x_temp=Model_Var(model_num,:);
        x2_temp=Model_Var(model_num,3:end-9);
        xi_temp=[];
        for ii=1:length(x2_temp)
            for jj=(ii+1):length(x2_temp)
                xi_temp=[xi_temp x2_temp(ii).*x2_temp(jj)];
            end
        end
        filter_paramters(model_num,:)=[x_temp xi_temp x2_temp];
    end
    County_Data.Vaccine_Uptake(County_Data.Vaccine_Uptake==1)=1-10^(-4);
    County_Data.Vaccine_Uptake(County_Data.Vaccine_Uptake==0)=10^(-4);


    refine_model_par=cell(2^7,1);
    L_fit=zeros(2^7,1);
    L_spatial_val=zeros(2^7,1);
    parfor model_num=1:2^7
            
        County_Data_model=County_Data;
        County_Data_model.X=County_Data_model.X(:,Model_Var(model_num,:));
        County_Data_model.X2=County_Data_model.X;
        if(Model_Var(model_num,end))
            County_Data_model.X2=County_Data_model.X2(:,3:end-9);
        else
            County_Data_model.X2=County_Data_model.X2(:,3:end);
        end
        XT=[];
        for ii=1:size(County_Data_model.X2,2)
            for jj=(ii+1):size(County_Data_model.X2,2)
                XT=[XT County_Data_model.X2(:,ii).*County_Data_model.X2(:,jj)];
            end
        end
    
        County_Data_model.X2=County_Data_model.X2.^2;
        County_Data_model.XI=XT;   
        
        num_par=size(County_Data_model.X,2)+size(County_Data_model.X2,2)+size(County_Data_model.XI,2);
    
        par_0=[];
        for kk=1:model_num
            par_temp=zeros(1,571);
            dx=filter_paramters(model_num,:)-filter_paramters(kk,:);
            if(min(dx)>=0)
                par_temp([true true true filter_paramters(kk,:) true true true true true])=model_par{kk};
                par_0=[par_0; par_temp];
            end
        end
        par_0=par_0(:,[true true true filter_paramters(model_num,:) true true true true true]);
        
        lb=-5*10^4.*ones(1,num_par+8); % 3 intercept and 5 hyper paramters
        lb(end-4:end-2)=-2;
        lb(end-1:end)=0;
        lb(1:5)=-50;
        ub=5*10^4.*ones(1,num_par+8);  % 3 intercept and 5 hyper paramters
        ub(end-4:end-2)=4;
        ub(end-1:end)=4;
        ub(1:5)=50;
        
        opts_ga=optimoptions("ga","PlotFcn",[],'UseParallel',false,"MaxGenerations",100,"FunctionTolerance",10^(-9),'CrossoverFcn','crossoverheuristic','MigrationInterval',25,'SelectionFcn',{@selectiontournament,8},'PopulationSize',250,'InitialPopulationMatrix',par_0);
        [par_est]=ga(@(x)Objective_Function_Coverage(x,County_Data_model,State_Data,NE_Data,filter_county_train,filter_state_train,filter_NE_train),length(lb),[],[],[],[],lb,ub,[],[],opts_ga);
    
        opts_pats=optimoptions('patternsearch','UseParallel',false,"PlotFcn",[],'FunctionTolerance',10^(-12),'MaxIterations',1000,'StepTolerance',10^(-9),'MaxFunctionEvaluations',10^4,'Cache','on');
        [par_final,f_final]=patternsearch(@(x)Objective_Function_Coverage(x,County_Data_model,State_Data,NE_Data,filter_county_train,filter_state_train,filter_NE_train),par_est,[],[],[],[],lb,ub,[],opts_pats);

        refine_model_par{model_num}=par_final;
        L_fit(model_num)=-f_final;
        L_spatial_val(model_num)=-Objective_Function_Coverage(par_final,County_Data_model,State_Data,NE_Data,filter_county_spatial_validation,filter_state_spatial_validation,filter_NE_spatial_validation);
    end
    save(['Refine_' Vaccine '_Set_Spatial_Validation=' num2str(Spatial_Validation) '.mat'],'refine_model_par','L_fit','L_spatial_val','Var_Names');
end

