function [County_Data,State_Data] = Load_Data(Vaccine)

Year=[];
Vaccine_Uptake=[];
Religious_Exemption=[];
Philosophial_Exemption=[];
State_FIP=[];
Age_5_to_9=[];
State_Name=[];
County_Name=[];
GEOID=[];
Spatial_Identifier=[];
Uninsured=[];
Private=[];
Public=[];

% Var_Names={'Age','Race','Income','Gini_Index','Education','Poverty_Income_Ratio','Rural_Urban_Code'};
for yy=2017:2023
    T=readtable([pwd '/County_Data.xlsx'],'Sheet',['Year_' num2str(yy)]);
    T.Physicians_per_capita(isnan(T.Physicians_per_capita))=0;
    if(yy==2017)
        X=[10.^4.*T.Physicians_per_capita table2array(T(:,22:26)) table2array(T(:,36:40)) log10(table2array(T(:,41))) table2array(T(:,42)) table2array(T(:,52:58)) table2array(T(:,59:70))];
        temp_RUCC=T.Rural_Urban_Continum_Code;
        RUCC=zeros(height(T),9);
        for ss=1:9
            RUCC(temp_RUCC==ss)=1;
        end
    else
        X=[X;10.^4.*T.Physicians_per_capita table2array(T(:,22:26)) table2array(T(:,36:40)) log10(table2array(T(:,41))) table2array(T(:,42)) table2array(T(:,52:58)) table2array(T(:,59:70))];
        temp_RUCC=T.Rural_Urban_Continum_Code;
        RUCC_t=zeros(height(T),9);
        for ss=1:9
            RUCC_t(temp_RUCC==ss,ss)=1;
        end
        RUCC=[RUCC;RUCC_t];
    end
    Uninsured=[Uninsured; T.Uninsured_Under_6];
    Private=[Private;T.Private_insured_Under_6];
    Public=[Public; T.Public_insured_Under_6];
    Spatial_Identifier=[Spatial_Identifier;T.Spatial_Identifier];
    Age_5_to_9=[Age_5_to_9;T.Population_5_to_9.*T.Total_Population];
    State_Name=[State_Name;T.State];
    County_Name=[County_Name;T.County];
    GEOID=[GEOID;T.GEOID];

    State_FIP=[State_FIP;str2double(T.State_FP)];
    Year=[Year; yy.*ones(height(T),1)];
    Vaccine_Uptake=[Vaccine_Uptake;T.(Vaccine)];
    if(strcmp(Vaccine,'MMR'))
        Religious_Exemption=[Religious_Exemption;T.MMR_Religious_Exemption];
        Philosophial_Exemption=[Philosophial_Exemption;T.MMR_Philosophical_Exemption];
    else
        Religious_Exemption=[Religious_Exemption;T.Other_Religious_Exemption];
        Philosophial_Exemption=[Philosophial_Exemption;T.Other_Philosophical_Exemption];
    end
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
County_Data.Age_5_to_9=Age_5_to_9(t_not_nan);
County_Data.Vaccine_Uptake=Vaccine_Uptake(t_not_nan);
County_Data.X=table2array(X);




Spatial_Stratification=cell(4,1);
Spatial_Stratification{1}={'ME','VT','NH','MA','RI','CT','NY','PA','NJ'};

Spatial_Stratification{2}={'DE','MD','WV','VA','KY','NC','TN','AR','OK','SC','GA','AL','MS','LA','TX','FL','DC'};

Spatial_Stratification{3}={'ND','MN','WI','MI','SD','IA','IL','IN','OH','NE','KS','MO'};

Spatial_Stratification{4}={'WA','MT','OR','ID','WY','CA','NV','UT','CO','AZ','NM'};

Year=[];
a_Beta_Parameters_Vaccine_Uptake=[];
b_Beta_Parameters_Vaccine_Uptake=[];
Spatial_Identifier=[];
State_FIP=[];
OR_Uninsured_Private=[];
OR_Public_Private=[];

for yy=2017:2023
    T=readtable([pwd '/State_Vaccination_Data.xlsx'],'Sheet',Vaccine);
    Surveyed=readtable([pwd '/State_Vaccination_Data.xlsx'],'Sheet','Surveyed_Population');
    Kindergarten=readtable([pwd '/State_Vaccination_Data.xlsx'],'Sheet','Kindergarten_Population');
    Vac_Coverage_Type=readtable([pwd '/Percent of children age 35 months who received recommended vaccines by Coverage Type.csv']);
    Spatial_Identifier_temp=zeros(height(T),1);
    for ss=1:4
        tf=ismember(T.STUSPS,Spatial_Stratification{ss});
        Spatial_Identifier_temp(tf)=ss;
    end

    Spatial_Identifier=[Spatial_Identifier;Spatial_Identifier_temp];
    State_FIP=[State_FIP;T.State_FIPs];
    Year=[Year; yy.*ones(height(T),1)];
    N=table2array(Kindergarten(:,3+(yy-2017)));
    V=table2array(T(:,4+(yy-2017)));
    P=table2array(Surveyed(:,3+(yy-2017)));

    a_Beta_Parameters_Vaccine_Uptake=[a_Beta_Parameters_Vaccine_Uptake; N.*V.*P];
    b_Beta_Parameters_Vaccine_Uptake=[b_Beta_Parameters_Vaccine_Uptake; N.*P-N.*V.*P];

    Vac_Uninsured=NaN.*ones(height(T.State_FIPs),1);
    Vac_Public=NaN.*ones(height(T.State_FIPs),1);
    Vac_Private=NaN.*ones(height(T.State_FIPs),1);
    for ss=1:height(T.State_FIPs)
        tf=T.State_FIPs(ss)==Vac_Coverage_Type.Fips & Vac_Coverage_Type.Year==yy & strcmp(Vac_Coverage_Type.CoverageType,'Uninsured');
        if(sum(tf)>0)
            Vac_Uninsured(ss)=Vac_Coverage_Type.Data(tf);
        end

        tf=T.State_FIPs(ss)==Vac_Coverage_Type.Fips & Vac_Coverage_Type.Year==yy & strcmp(Vac_Coverage_Type.CoverageType,'Medicaid');
        if(sum(tf)>0)
            Vac_Public(ss)=Vac_Coverage_Type.Data(tf);
        end

        tf=T.State_FIPs(ss)==Vac_Coverage_Type.Fips & Vac_Coverage_Type.Year==yy & strcmp(Vac_Coverage_Type.CoverageType,'Private');
        if(sum(tf)>0)
            Vac_Private(ss)=Vac_Coverage_Type.Data(tf);
        end
    end
    OR_Uninsured_Private=[OR_Uninsured_Private; Vac_Uninsured.*(1-Vac_Private)./(Vac_Private.*(1-Vac_Uninsured)) ];
    OR_Public_Private=[OR_Public_Private; Vac_Public.*(1-Vac_Private)./(Vac_Private.*(1-Vac_Public)) ];
end

State_Data=table(Year,State_FIP,Spatial_Identifier,a_Beta_Parameters_Vaccine_Uptake,b_Beta_Parameters_Vaccine_Uptake,OR_Uninsured_Private,OR_Public_Private);
State_Data=State_Data(~isnan(b_Beta_Parameters_Vaccine_Uptake),:);
end
