function [County_Transmission_X] = Load_Transmission_Covariates(GEOID)

Year=2023.*ones(length(GEOID),1);
State_FIP=cell(length(GEOID),1);
State_Name=cell(length(GEOID),1);
County_Name=cell(length(GEOID),1);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%555
% State Data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%555

% Var_Names={'Age','Race','Income','Gini_Index','Education','Poverty_Income_Ratio','Rural_Urban_Code'};
yy=2023;
T=readtable([pwd '/County_Data.xlsx'],'Sheet',['Year_' num2str(yy)]);
T.Physicians_per_capita(isnan(T.Physicians_per_capita))=0;
GEO_temp=T.GEOID;
X_temp=[log10(table2array(T(:,41))) table2array(T(:,42)) table2array(T(:,53:58)) table2array(T(:,59:70))]; % Education less than grade nine was removed such that the covariance matrix would be positive definite
temp_RUCC=T.Rural_Urban_Continum_Code;
RUCC=zeros(height(T),9);
for ss=1:9
    RUCC(temp_RUCC==ss,ss)=1;
end

RUCC(sum(RUCC,2)==0,:)=NaN.*RUCC(sum(RUCC,2)==0,:);
X_temp=table2array([table(X_temp) array2table(RUCC)]);

X=zeros(length(GEOID),size(X_temp,2));

for ii=1:length(GEOID)
    tf=strcmp(GEO_temp,GEOID{ii});
    X(ii,:)=X_temp(tf,:);
    State_FIP{ii}=T.State_FP{tf};
    State_Name{ii}=T.State{tf};
    County_Name{ii}=T.County{tf};
end

County_Transmission_X.Year=Year;
County_Transmission_X.State=State_Name;
County_Transmission_X.County=County_Name;
County_Transmission_X.GEOID=GEOID;

County_Transmission_X.State_FIP=State_FIP;
County_Transmission_X.X=X;
end
