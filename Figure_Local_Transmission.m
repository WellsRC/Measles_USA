clear;
clc;
close all;

load('National_Reduction=1_Year=1.mat','County_Data_Vaccine_Reduction');
load('Baseline_Estimate_Measles_Incidence_OLD_KEEP.mat','k_mealses');

I=[0:10];

y1=zeros(length(County_Data_Vaccine_Reduction.R_eff),length(I));

for jj=1:length(y1)
    y1(jj,:)=1-integral(@(x)nbinpdf(0,k_mealses,k_mealses./(k_mealses+County_Data_Vaccine_Reduction.R_eff(jj).*x)),0,1).^I;
end

load('National_Reduction=1_Year=5.mat','County_Data_Vaccine_Reduction');


y5=zeros(length(County_Data_Vaccine_Reduction.R_eff),length(I));

for jj=1:length(y1)
    y5(jj,:)=1-integral(@(x)nbinpdf(0,k_mealses,k_mealses./(k_mealses+County_Data_Vaccine_Reduction.R_eff(jj).*x)),0,1).^I;
end