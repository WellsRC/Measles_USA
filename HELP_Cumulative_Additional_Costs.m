clear;
close all;
clc;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Baseline Calculations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Age_0_to_6=true;
National_Reduction=0.01;
[~,~,pd_1_cost]=National_Outcome_Distribution(National_Reduction,'Sample_2025',Age_0_to_6);
National_Reduction=0.02;
[~,~,pd_2_cost]=National_Outcome_Distribution(National_Reduction,'Sample_2025',Age_0_to_6);
National_Reduction=0.03;
[~,~,pd_3_cost]=National_Outcome_Distribution(National_Reduction,'Sample_2025',Age_0_to_6);
National_Reduction=0.04;
[~,~,pd_4_cost]=National_Outcome_Distribution(National_Reduction,'Sample_2025',Age_0_to_6);
National_Reduction=0.05;
[~,~,pd_5_cost]=National_Outcome_Distribution(National_Reduction,'Sample_2025',Age_0_to_6);

Cost=cell(1,5);
% Cost
x0=icdf(pd_1_cost,0.5);
lb=icdf(pd_1_cost,0.05);
ub=icdf(pd_1_cost,0.95);
temp_X=10.^fmincon(@(z)-log(pdf(pd_1_cost,10.^z)),log10(x0),[],[],[],[],log10(lb),log10(ub));
Cost{1} = ['$' num2str(temp_X,'%5.1f') '(95% CrI: $' num2str(icdf(pd_1_cost,0.025),'%5.1f') char(8211) '$' num2str(icdf(pd_1_cost,0.975),'%5.1f') ')'];


Y2=random(pd_1_cost,10^4,1)+random(pd_2_cost,10^4,1);
pd=fitdist(Y2,'Kernel','Support','positive');
x0=icdf(pd,0.5);
lb=icdf(pd,0.05);
ub=icdf(pd,0.95);
temp_X=10.^fmincon(@(z)-log(pdf(pd,10.^z)),log10(x0),[],[],[],[],log10(lb),log10(ub));
Cost{2} =['$' num2str(temp_X,'%5.1f') '(95% CrI: $' num2str(icdf(pd,0.025),'%5.1f') char(8211) '$' num2str(icdf(pd,0.975),'%5.1f') ')'];

Y3=Y2+random(pd_3_cost,10^4,1);
pd=fitdist(Y3,'Kernel','Support','positive');
x0=icdf(pd,0.5);
lb=icdf(pd,0.05);
ub=icdf(pd,0.95);
temp_X=10.^fmincon(@(z)-log(pdf(pd,10.^z)),log10(x0),[],[],[],[],log10(lb),log10(ub));
Cost{3} = ['$' num2str(temp_X,'%5.1f') '(95% CrI: $' num2str(icdf(pd,0.025),'%5.1f') char(8211) '$' num2str(icdf(pd,0.975),'%5.1f') ')'];


Y4=Y3+random(pd_3_cost,10^4,1);
pd=fitdist(Y4,'Kernel','Support','positive');
x0=icdf(pd,0.5);
lb=icdf(pd,0.05);
ub=icdf(pd,0.95);
temp_X=10.^fmincon(@(z)-log(pdf(pd,10.^z)),log10(x0),[],[],[],[],log10(lb),log10(ub));
Cost{4} = ['$' num2str(temp_X,'%5.1f') '(95% CrI: $' num2str(icdf(pd,0.025),'%5.1f') char(8211) '$' num2str(icdf(pd,0.975),'%5.1f') ')'];


Y5=Y4+random(pd_4_cost,10^4,1);
pd=fitdist(Y5,'Kernel','Support','positive');
x0=icdf(pd,0.5);
lb=icdf(pd,0.05);
ub=icdf(pd,0.95);
temp_X=10.^fmincon(@(z)-log(pdf(pd,10.^z)),log10(x0),[],[],[],[],log10(lb),log10(ub));
Cost{5} = ['$' num2str(temp_X,'%5.1f') '(95% CrI: $' num2str(icdf(pd,0.025),'%5.1f') char(8211) '$' num2str(icdf(pd,0.975),'%5.1f') ')'];

[~,~,pd_0_cost]=National_Outcome_Distribution(0,'Sample_2025',Age_0_to_6);

Cost_Baseline=cell(1,5);
% Cost
x0=icdf(pd_0_cost,0.5);
lb=icdf(pd_0_cost,0.05);
ub=icdf(pd_0_cost,0.95);
temp_X=10.^fmincon(@(z)-log(pdf(pd_0_cost,10.^z)),log10(x0),[],[],[],[],log10(lb),log10(ub));
Cost_Baseline{1} = ['$' num2str(temp_X,'%5.1f') '(95% CrI: $' num2str(icdf(pd_0_cost,0.025),'%5.1f') char(8211) '$' num2str(icdf(pd_0_cost,0.975),'%5.1f') ')'];


Y2=random(pd_0_cost,10^4,1)+random(pd_0_cost,10^4,1);
pd=fitdist(Y2,'Kernel','Support','positive');
x0=icdf(pd,0.5);
lb=icdf(pd,0.05);
ub=icdf(pd,0.95);
temp_X=10.^fmincon(@(z)-log(pdf(pd,10.^z)),log10(x0),[],[],[],[],log10(lb),log10(ub));
Cost_Baseline{2} =['$' num2str(temp_X,'%5.1f') '(95% CrI: $' num2str(icdf(pd,0.025),'%5.1f') char(8211) '$' num2str(icdf(pd,0.975),'%5.1f') ')'];

Y3=Y2+random(pd_0_cost,10^4,1);
pd=fitdist(Y3,'Kernel','Support','positive');
x0=icdf(pd,0.5);
lb=icdf(pd,0.05);
ub=icdf(pd,0.95);
temp_X=10.^fmincon(@(z)-log(pdf(pd,10.^z)),log10(x0),[],[],[],[],log10(lb),log10(ub));
Cost_Baseline{3} = ['$' num2str(temp_X,'%5.1f') '(95% CrI: $' num2str(icdf(pd,0.025),'%5.1f') char(8211) '$' num2str(icdf(pd,0.975),'%5.1f') ')'];


Y4=Y3+random(pd_0_cost,10^4,1);
pd=fitdist(Y4,'Kernel','Support','positive');
x0=icdf(pd,0.5);
lb=icdf(pd,0.05);
ub=icdf(pd,0.95);
temp_X=10.^fmincon(@(z)-log(pdf(pd,10.^z)),log10(x0),[],[],[],[],log10(lb),log10(ub));
Cost_Baseline{4} = ['$' num2str(temp_X,'%5.1f') '(95% CrI: $' num2str(icdf(pd,0.025),'%5.1f') char(8211) '$' num2str(icdf(pd,0.975),'%5.1f') ')'];


Y5=Y4+random(pd_0_cost,10^4,1);
pd=fitdist(Y5,'Kernel','Support','positive');
x0=icdf(pd,0.5);
lb=icdf(pd,0.05);
ub=icdf(pd,0.95);
temp_X=10.^fmincon(@(z)-log(pdf(pd,10.^z)),log10(x0),[],[],[],[],log10(lb),log10(ub));
Cost_Baseline{5} = ['$' num2str(temp_X,'%5.1f') '(95% CrI: $' num2str(icdf(pd,0.025),'%5.1f') char(8211) '$' num2str(icdf(pd,0.975),'%5.1f') ')'];
