% clear;
close all;
Age_0_to_6=true;
National_Reduction=0.01;

if(~Age_0_to_6)
    Age_Reduction=[true(1,1) false(1,4)];
else
    Age_Reduction=[true(1,2) false(1,3)];
end

S=shaperead([pwd '\Shapefile\cb_2023_us_county_20m.shp'],"UseGeoCoords",true);

FN_Age_Class={'Age_0_to_4','_Age_5_to_9','_Age_10_to_14','_Age_15_to_19','_Age_20_to_24'};
FN_Age_Class=FN_Age_Class(Age_Reduction);

if(Age_0_to_6)
    load(['National_Reduction=' num2str(100*0) '_Ages_0_to_6.mat'],'County_Data_Vaccine_Reduction')
else
    load(['National_Reduction=' num2str(100*0) '_' FN_Age_Class{:} '.mat'],'County_Data_Vaccine_Reduction')
end

Uptake_Baseline=(table2array(County_Data_Vaccine_Reduction.Vaccine_Uptake(:,1)).*table2array(County_Data_Vaccine_Reduction.Population(:,1))+(2/5).*table2array(County_Data_Vaccine_Reduction.Vaccine_Uptake(:,2)).*table2array(County_Data_Vaccine_Reduction.Population(:,2)))./(table2array(County_Data_Vaccine_Reduction.Population(:,1))+(2/5).*table2array(County_Data_Vaccine_Reduction.Population(:,2)));
V_Baseline=NaN.*zeros(length(S),1);
V_Reduction=NaN.*zeros(length(S),1);


for cc=1:length(S)
    tf=strcmp(County_Data_Vaccine_Reduction.GEOID,S(cc).GEOID);
    if(sum(tf)>0)
        if(~Age_0_to_6)
            V_Baseline(cc)=table2array(County_Data_Vaccine_Reduction.Vaccine_Uptake(tf,Age_Reduction));
        else
            V_Baseline(cc)=(table2array(County_Data_Vaccine_Reduction.Vaccine_Uptake(tf,1)).*table2array(County_Data_Vaccine_Reduction.Population(tf,1))+(2/5).*table2array(County_Data_Vaccine_Reduction.Vaccine_Uptake(tf,2)).*table2array(County_Data_Vaccine_Reduction.Population(tf,2)))./(table2array(County_Data_Vaccine_Reduction.Population(tf,1))+(2/5).*table2array(County_Data_Vaccine_Reduction.Population(tf,2)));
        end
    end
end

if(Age_0_to_6)
    load(['National_Reduction=' num2str(100*National_Reduction) '_Ages_0_to_6.mat'],'County_Data_Vaccine_Reduction')
else
    load(['National_Reduction=' num2str(100*National_Reduction) '_' FN_Age_Class{:} '.mat'],'County_Data_Vaccine_Reduction')
end

Uptake_Reduction=(table2array(County_Data_Vaccine_Reduction.Vaccine_Uptake(:,1)).*table2array(County_Data_Vaccine_Reduction.Population(:,1))+(2/5).*table2array(County_Data_Vaccine_Reduction.Vaccine_Uptake(:,2)).*table2array(County_Data_Vaccine_Reduction.Population(:,2)))./(table2array(County_Data_Vaccine_Reduction.Population(:,1))+(2/5).*table2array(County_Data_Vaccine_Reduction.Population(:,2)));

Pop_0_to_6=table2array(County_Data_Vaccine_Reduction.Population(:,1).*County_Data_Vaccine_Reduction.Total_Population(:))+(2/5).*table2array(County_Data_Vaccine_Reduction.Population(:,2).*County_Data_Vaccine_Reduction.Total_Population(:));
State_c=County_Data_Vaccine_Reduction.State;
UState=unique(County_Data_Vaccine_Reduction.State);

for cc=1:length(S)
    tf=strcmp(County_Data_Vaccine_Reduction.GEOID,S(cc).GEOID);
     if(sum(tf)>0)
        if(~Age_0_to_6)
            V_Reduction(cc)=table2array(County_Data_Vaccine_Reduction.Vaccine_Uptake(tf,Age_Reduction));
        else
            V_Reduction(cc)=(table2array(County_Data_Vaccine_Reduction.Vaccine_Uptake(tf,1)).*table2array(County_Data_Vaccine_Reduction.Population(tf,1))+(2/5).*table2array(County_Data_Vaccine_Reduction.Vaccine_Uptake(tf,2)).*table2array(County_Data_Vaccine_Reduction.Population(tf,2)))./(table2array(County_Data_Vaccine_Reduction.Population(tf,1))+(2/5).*table2array(County_Data_Vaccine_Reduction.Population(tf,2)));
        end
    end
end


V_State_Baseline=zeros(length(UState),1);
V_State_Reducton=zeros(length(UState),1);

for ss=1:length(UState)
    tf=strcmp(UState{ss},State_c) & ~isnan(Uptake_Baseline);
    V_State_Baseline(ss)=sum(Pop_0_to_6(tf).*Uptake_Baseline(tf))./sum(Pop_0_to_6(tf));
    V_State_Reducton(ss)=sum(Pop_0_to_6(tf).*Uptake_Reduction(tf))./sum(Pop_0_to_6(tf));
end

D_Reduction_State=V_State_Baseline-V_State_Reducton;

D_Reduction_County=V_Baseline-V_Reduction;

tf=~isnan(Uptake_Baseline);
US_Baseline=sum(Pop_0_to_6(tf).*Uptake_Baseline(tf))./sum(Pop_0_to_6(tf));
US_Reduction=sum(Pop_0_to_6(tf).*Uptake_Reduction(tf))./sum(Pop_0_to_6(tf));

[~,SD]=sort(D_Reduction_State);
L_Ustate=UState(SD);


