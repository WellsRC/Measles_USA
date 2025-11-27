clear;
clc;
close all;

Vaccine='MMR';
load([Vaccine '_Immunity.mat'],'County_Data')

Known_Ind_Cases=zeros(length(County_Data.County),1);

Measles_Cases=readtable('County_Level_Measles_Cases_Adjusted.csv');

Imported_Case=zeros(length(County_Data.County),1);

% Known imported cases

beta_Cap=zeros(length(County_Data.GEOID),1);
for cc=1:length(County_Data.GEOID)
    t_f=str2double(County_Data.GEOID{cc})==Measles_Cases.GEOID & strcmp(Measles_Cases.type,'imported') & ~isnan(Measles_Cases.case_count);
    if(sum(t_f)>0)
        Imported_Case(cc)=Measles_Cases.case_count(t_f);
    end

    t_f=str2double(County_Data.GEOID{cc})==Measles_Cases.GEOID & strcmp(Measles_Cases.type,'local') & ~isnan(Measles_Cases.case_count);
    if(sum(t_f)>0)
        Known_Ind_Cases(cc)=Measles_Cases.case_count(t_f);
    end
    beta_Cap(cc)=pchip(County_Data.R_0(cc,:),County_Data.beta_j,18.1);
end


f_out=find(Known_Ind_Cases>=3);
b_temp=zeros(length(f_out),1);
R0=zeros(length(b_temp),1);


for jj=1:length(b_temp)
    xx=County_Data.Final_Size_Est(f_out(jj),:);
    bb=County_Data.beta_j;
    bb=bb(xx>0);
    xx=xx(xx>0);
    b_temp(jj)=pchip(xx,bb,Known_Ind_Cases(f_out(jj)));
    R0(jj)=interp1(County_Data.beta_j',County_Data.R_0(f_out(jj),:)',b_temp(jj))';
end
M_beta=median(beta_Cap);
N_temp=1-County_Data.Total_Immunity(f_out);


load('Baseline_Estimate_Measles_Incidence.mat','indx_beta');
load('Explore_Beta_Relation_Parameters.mat')

[X_Samp,indx_s]=unique(X_Samp,"rows");
L_Samp=L_Samp(indx_s);

X_Samp=X_Samp(indx_beta,:);


figure(2)
ii=linspace(0,0.5,1001);
plot(ii,Transmission_Relation(ii,X_Samp(1:4)),'k','LineWidth',2); hold on
scatter((N_temp),b_temp,30,'r','filled');
set(gca,'tickdir','out','FontSize',16,'XTick',[0:0.1:0.5],'LineWidth',2)
xlabel('Susceptible','FontSize',18)
ylabel('\beta_c','FontSize',18)
box off;

print(gcf,['Figure_S3.png'],'-dpng','-r300');
