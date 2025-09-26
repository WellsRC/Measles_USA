function [County_Data,State_Data] = Load_Data_Adjustment(Vaccine,Age_Group)
year_data=2023;
year_data_A_15_19=year_data-2;
% Var_Names={'Age','Race','Income','Gini_Index','Education','Poverty_Income_Ratio','Rural_Urban_Code'};

    T=readtable([pwd '/County_Data.xlsx'],'Sheet',['Year_' num2str(year_data)]);
    
    X=[10.^4.*T.Physicians_per_capita table2array(T(:,22:26)) table2array(T(:,36:40)) log10(table2array(T(:,41))) table2array(T(:,42)) table2array(T(:,52:58)) table2array(T(:,59:70))];
    temp_RUCC=T.Rural_Urban_Continum_Code;
    RUCC=zeros(height(T),9);
    for ss=1:9
        RUCC(temp_RUCC==ss)=1;
    end
    
    Spatial_Identifier=[T.Spatial_Identifier];
    if(strcmp(Age_Group,'Age_0_to_4'))
        Weight=[T.Population_0_to_4.*T.Total_Population];

        Uninsured=[T.Uninsured_Under_6];
        Private=[T.Private_insured_Under_6];
        Public=[T.Public_insured_Under_6];
    elseif(strcmp(Age_Group,'Age_5_to_9'))
        Weight=[T.Population_5_to_9.*T.Total_Population];
        
        Uninsured=[T.Uninsured_Under_6];
        Private=[T.Private_insured_Under_6];
        Public=[T.Public_insured_Under_6];
    elseif(strcmp(Age_Group,'Age_10_to_14'))
        Weight=[T.Population_10_to_14.*T.Total_Population];
        
        Uninsured=[T.Uninsured_6_to_18];
        Private=[T.Private_insured_6_to_18];
        Public=[T.Public_insured_6_to_18];
    elseif(strcmp(Age_Group,'Age_15_to_19'))
        Weight=[T.Population_15_to_19.*T.Total_Population];

        Uninsured=[T.Uninsured_6_to_18];
        Private=[T.Private_insured_6_to_18];
        Public=[T.Public_insured_6_to_18];

    elseif(strcmp(Age_Group,'Age_20_to_24'))
        Weight=[T.Population_20_to_24.*T.Total_Population];

        Uninsured=[T.Uninsured_19_to_25];
        Private=[T.Private_insured_19_to_25];
        Public=[T.Public_insured_19_to_25];
    end
    State_Name=[T.State];
    County_Name=[T.County];
    GEOID=[T.GEOID];

    State_FIP=[str2double(T.State_FP)];
    Year=[year_data.*ones(height(T),1)];
    Vaccine_Uptake=[T.(Vaccine)];
    if(strcmp(Vaccine,'MMR'))
        Religious_Exemption=[T.MMR_Religious_Exemption];
        Philosophial_Exemption=[T.MMR_Philosophical_Exemption];
    else
        Religious_Exemption=[T.Other_Religious_Exemption];
        Philosophial_Exemption=[T.Other_Philosophical_Exemption];
    end


RUCC(sum(RUCC,2)==0,:)=NaN.*RUCC(sum(RUCC,2)==0,:);
X=[table(Religious_Exemption,Philosophial_Exemption) table(X) array2table(RUCC)];
t_not_nan=~isnan(sum(table2array(X),2));

X=X(t_not_nan,:);
County_Data.Year=Year(t_not_nan);
County_Data.State=State_Name(t_not_nan);
County_Data.County=County_Name(t_not_nan);
County_Data.GEOID=GEOID(t_not_nan);
County_Data.Spatial_Identifier=Spatial_Identifier(t_not_nan);
County_Data.Uninsured=Uninsured(t_not_nan);
County_Data.Private=Private(t_not_nan);
County_Data.Public=Public(t_not_nan);

County_Data.State_FIP=State_FIP(t_not_nan);
County_Data.Weight=Weight(t_not_nan);
if(strcmp(Age_Group,'Age_5_to_9'))
    County_Data.Vaccine_Uptake=Vaccine_Uptake(t_not_nan);
else
    County_Data.Vaccine_Uptake=NaN.*Vaccine_Uptake(t_not_nan);
end
County_Data.X=table2array(X);




Spatial_Stratification=cell(4,1);
Spatial_Stratification{1}={'ME','VT','NH','MA','RI','CT','NY','PA','NJ'};

Spatial_Stratification{2}={'DE','MD','WV','VA','KY','NC','TN','AR','OK','SC','GA','AL','MS','LA','TX','FL','DC'};

Spatial_Stratification{3}={'ND','MN','WI','MI','SD','IA','IL','IN','OH','NE','KS','MO'};

Spatial_Stratification{4}={'WA','MT','OR','ID','WY','CA','NV','UT','CO','AZ','NM'};

if(strcmp(Age_Group,'Age_0_to_4'))
    T=readtable(['Estimated_State_' Vaccine '_Uptake_0_to_4.xlsx']);


    Spatial_Identifier_temp=zeros(height(T),1);
    for ss=1:4
        tf=ismember(T.STUSPS,Spatial_Stratification{ss});
        Spatial_Identifier_temp(tf)=ss;
    end

    Spatial_Identifier=[Spatial_Identifier_temp];
    State_FIP=[T.State_FIPs];
    Year=[year_data.*ones(height(T),1)];
    
    Vaccine_Uptake=T.Vaccine_Uptake;

elseif(strcmp(Age_Group,'Age_5_to_9'))
    T=readtable([pwd '/State_Vaccination_Data.xlsx'],'Sheet',Vaccine);
   
    Spatial_Identifier_temp=zeros(height(T),1);
    for ss=1:4
        tf=ismember(T.STUSPS,Spatial_Stratification{ss});
        Spatial_Identifier_temp(tf)=ss;
    end

    Spatial_Identifier=[Spatial_Identifier_temp];
    State_FIP=[T.State_FIPs];
    Year=[year_data.*ones(height(T),1)];
    
    V=table2array(T(:,4+(year_data-2017)));

    Vaccine_Uptake=V;
elseif(strcmp(Age_Group,'Age_10_to_14'))
    T=readtable([pwd '/Teen_State_Vaccination_Data.xlsx'],'Sheet',Vaccine);
   
    Spatial_Identifier_temp=zeros(height(T),1);
    for ss=1:4
        tf=ismember(T.STUSPS,Spatial_Stratification{ss});
        Spatial_Identifier_temp(tf)=ss;
    end

    Spatial_Identifier=[Spatial_Identifier_temp];
    State_FIP=[T.State_FIPs];
    Year=[year_data.*ones(height(T),1)];
    
    V=table2array(T(:,4+(year_data-2016)));

    Vaccine_Uptake=V;

elseif(strcmp(Age_Group,'Age_15_to_19'))
    T=readtable([pwd '/Teen_State_Vaccination_Data.xlsx'],'Sheet',Vaccine);
   
    Spatial_Identifier_temp=zeros(height(T),1);
    for ss=1:4
        tf=ismember(T.STUSPS,Spatial_Stratification{ss});
        Spatial_Identifier_temp(tf)=ss;
    end

    Spatial_Identifier=[Spatial_Identifier_temp];
    State_FIP=[T.State_FIPs];
    Year=[year_data.*ones(height(T),1)];
    
    V=table2array(T(:,4+(year_data_A_15_19-2016)));

    Vaccine_Uptake=V;
elseif(strcmp(Age_Group,'Age_20_to_24'))
    T=readtable([pwd '/Teen_State_Vaccination_Data.xlsx'],'Sheet',Vaccine);
   
    Spatial_Identifier_temp=zeros(height(T),1);
    for ss=1:4
        tf=ismember(T.STUSPS,Spatial_Stratification{ss});
        Spatial_Identifier_temp(tf)=ss;
    end

    Spatial_Identifier=[Spatial_Identifier_temp];
    State_FIP=[T.State_FIPs];
    Year=[year_data.*ones(height(T),1)];
    
    V=table2array(T(:,4)); % Takes year 2016

    Vaccine_Uptake=V;
end


State_Data=table(Year,State_FIP,Spatial_Identifier,Vaccine_Uptake);
State_Data=State_Data(~isnan(Vaccine_Uptake),:);
end
