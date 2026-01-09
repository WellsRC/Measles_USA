clear;
clc;

Vaccine='MMR';
load([Vaccine '_Immunity.mat'],'County_Data')


% Data as of December 30
Measles_Cases=readtable('County_Level_Measles_Cases_Adjusted.csv');

Imported_Case=zeros(length(County_Data.County),1);

Known_Ind_Cases=zeros(length(County_Data.County),1);
Unknown_Ind_Cases=NaN.*zeros(length(County_Data.County),2);
Unknown_Ind_Cases_Weight=NaN.*zeros(length(County_Data.County),2);

% Known imported cases
for cc=1:length(Known_Ind_Cases)
    t_f=str2double(County_Data.GEOID{cc})==Measles_Cases.GEOID & strcmp(Measles_Cases.type,'local') & ~isnan(Measles_Cases.case_count);
    if(sum(t_f)>0)
        Known_Ind_Cases(cc)=Measles_Cases.case_count(t_f);
    end
end

[County_Transmission_X] = Load_Transmission_Covariates(County_Data.GEOID);

X=County_Transmission_X.X;
Y=NaN.*zeros(size(X,1),1);
for ii=1:length(Y)
    if(Known_Ind_Cases(ii)>=3)
        tz=County_Data.Final_Size_Est(ii,:)>0;
        Y(ii)=pchip(County_Data.Final_Size_Est(ii,tz),County_Data.beta_j(tz),Known_Ind_Cases(ii));
    end
end

I=Known_Ind_Cases(~isnan(Y));

X=X(~isnan(Y),:);
Y=log(Y(~isnan(Y)));

w=I./sum(I);
mdl=fitlm(X,Y,Weights=w,Intercept=false);

ypred=exp(predict(mdl,X));

scatter(exp(Y),ypred);
hold on;
plot(linspace(min(exp(Y)),max(exp(Y)),1001),linspace(min(exp(Y)),max(exp(Y)),1001))
ylabel('predicted')

prior_mean_transmission=mdl.Coefficients.Estimate';
prior_covariance_transmission=mdl.CoefficientCovariance;
transmission_bnds = coefCI(mdl)';

save('Prior_log_Regression_Transmission.mat','prior_mean_transmission','transmission_bnds');