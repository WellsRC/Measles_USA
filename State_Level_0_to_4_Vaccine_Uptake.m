clear;

A_0_to_4=readtable('Child population by single age.xlsx');

A_0_to_4=A_0_to_4(A_0_to_4.TimeFrame==2023 & A_0_to_4.SingleAge<5,:);
% First dose of MMR does not happen until 13 months of age so those in the
% first year will have 0 uptake

w0=A_0_to_4(A_0_to_4.SingleAge==0,:);
w1=A_0_to_4(A_0_to_4.SingleAge==1,:);
w2=A_0_to_4(A_0_to_4.SingleAge==2,:);
w3=A_0_to_4(A_0_to_4.SingleAge==3,:);
w4=A_0_to_4(A_0_to_4.SingleAge==4,:);

V_13=readtable('13months_State_Vaccination_Data.xlsx','Sheet','MMR');
V_24=readtable('24months_State_Vaccination_Data.xlsx','Sheet','MMR');
V_35=readtable('35months_State_Vaccination_Data.xlsx','Sheet','MMR');

State=V_13.State;
STUSPS=V_13.STUSPS;
State_FIPs=V_13.State_FIPs;

Vaccine_Uptake=zeros(length(State),1);

for ss=1:length(State)

    w_NB=w0.Data(strcmp(w0.Location,State{ss}));
    temp_13=V_13.Birth_Year_2021(strcmp(V_13.State,State{ss}));
    w_13=w1.Data(strcmp(w1.Location,State{ss}));

    temp_24=V_24.Birth_Year_2021(strcmp(V_24.State,State{ss}));
    w_24=w2.Data(strcmp(w2.Location,State{ss}));

    temp_35=V_35.Birth_Year_2020(strcmp(V_35.State,State{ss}));
    w_35=w3.Data(strcmp(w3.Location,State{ss}));

    temp_48=V_35.Birth_Year_2019(strcmp(V_35.State,State{ss}));
    w_48=w4.Data(strcmp(w4.Location,State{ss}));

    Vaccine_Uptake(ss)=(temp_13.*w_13+temp_24.*w_24+temp_35.*w_35+temp_48.*w_48)./(w_13+w_24+w_35+w_NB+w_48);
end

T=table(State,STUSPS,State_FIPs,Vaccine_Uptake);

writetable(T,'Estimated_State_MMR_Uptake_0_to_4.xlsx');




