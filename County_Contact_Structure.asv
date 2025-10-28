clear;


clc;
clear;
clear;
close all;
Vaccine='MMR';
L_F=zeros(2^7,4);
L_S=zeros(2^7,4);
for Spatial_Valiation=1:4
    load(['Refine_' Vaccine '_Set_Spatial_Validation=' num2str(Spatial_Valiation) '.mat']);
    L_F(:,Spatial_Valiation)=L_fit;
    L_S(:,Spatial_Valiation)=L_spatial_val;
end

L_T=sum(L_F,2)+sum(L_S,2);
    
[~,Model_Num]=max(L_T);

clearvars -except Model_Num Vaccine Age_Class


[temp_county]=Age_Adjustment_Factor(Vaccine,Model_Num,'Age_0_to_4');

County_Contacts.County=temp_county.County;
County_Contacts.State=temp_county.State;
County_Contacts.GEOID=temp_county.GEOID;

clear temp_county

Age_Class={'Age_0_to_4','Age_5_to_9','Age_10_to_14','Age_15_to_19','Age_20_to_24','Age_25_to_29','Age_30_to_34','Age_35_to_39','Age_40_to_44','Age_45_to_49','Age_50_to_54','Age_55_to_59','Age_60_to_64','Age_65_to_69','Age_70_to_74','Age_75_to_79','Age_80_to_84','Age_85_plus'};

County_Contacts.Age_Class=Age_Class;

Contacts=zeros(length(County_Contacts.County),length(Age_Class));

for cc=1:length(County_Contacts.County)
    State_temp=County_Contacts.State{cc};
    State_temp(State_temp==' ')='_';
    M=readtable([pwd '/State_Contact_Matrix/United_States_subnational_' State_temp '_M_overall_contact_matrix_18']);
    M=table2array(M);
    Contacts(cc,:)=sum(M,2);
end

County_Contacts.Age_Class=Age_Class;

County_Contacts=Contacts;

save('County_Contacts_Age.mat','County_Contacts');
